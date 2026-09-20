#include "n2s/fit/layout.hpp"
#include "n2s/fit/locate.hpp"
#include "n2s/fit/refine.hpp"
#include "n2s/io/case_json.hpp"
#include "n2s/tolerances.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <optional>
#include <vector>

using Catch::Approx;
using n2s::fit::DomainLayout;
using n2s::fit::LayoutLocator;

namespace {

/// The bilinear map, written out here rather than called from the library, so
/// that the round-trip tests below have an oracle the implementation cannot
/// quietly agree with by sharing a bug.
Eigen::Vector2d blend(const Eigen::Vector2d& a,
                      const Eigen::Vector2d& b,
                      const Eigen::Vector2d& c,
                      const Eigen::Vector2d& d,
                      double u,
                      double v) {
    return (1.0 - u) * (1.0 - v) * a + u * (1.0 - v) * b + u * v * c + (1.0 - u) * v * d;
}

DomainLayout thesis_layout() {
    const std::filesystem::path path = std::filesystem::path{N2S_DATA_DIR} / "thesis_doublevb.json";
    const n2s::io::Case loaded = n2s::io::read_case(path);
    REQUIRE(loaded.layout.has_value());
    return *loaded.layout;
}

} // namespace

// ---------------------------------------------------------------------------
// The inversion itself, against closed-form answers.
// ---------------------------------------------------------------------------

TEST_CASE("an axis-aligned rectangle inverts to the obvious parameters", "[fit][locate]") {
    // The one case where the answer can be written down with no algebra: on a
    // rectangle the bilinear map is separable and (u, v) is just the point's
    // position along each side.
    const Eigen::Vector2d a{2.0, 5.0};
    const Eigen::Vector2d b{6.0, 5.0};
    const Eigen::Vector2d c{6.0, 11.0};
    const Eigen::Vector2d d{2.0, 11.0};

    for (const double u : {0.0, 0.25, 0.5, 1.0}) {
        for (const double v : {0.0, 1.0 / 3.0, 0.75, 1.0}) {
            const Eigen::Vector2d p{2.0 + 4.0 * u, 5.0 + 6.0 * v};

            const std::optional<Eigen::Vector2d> uv = n2s::fit::invert_bilinear(a, b, c, d, p);

            INFO("u " << u << " v " << v);
            REQUIRE(uv.has_value());
            CHECK(uv->x() == Approx(u).margin(1e-14));
            CHECK(uv->y() == Approx(v).margin(1e-14));
        }
    }
}

TEST_CASE("a parallelogram inverts through the linear branch", "[fit][locate]") {
    // A sheared quad still has a vanishing quadratic term, so this exercises
    // the branch that solves a linear system. The oracle is the inverse of the
    // 2x2 matrix [B D], which is exact here because the map really is affine.
    const Eigen::Vector2d a{0.0, 0.0};
    const Eigen::Vector2d b{2.0, 1.0};
    const Eigen::Vector2d c{3.0, 4.0};
    const Eigen::Vector2d d{1.0, 3.0};

    Eigen::Matrix2d basis;
    basis.col(0) = b - a;
    basis.col(1) = d - a;

    for (const double u : {0.1, 0.5, 0.9}) {
        for (const double v : {0.2, 0.6, 1.0}) {
            const Eigen::Vector2d p = a + basis * Eigen::Vector2d{u, v};

            const std::optional<Eigen::Vector2d> uv = n2s::fit::invert_bilinear(a, b, c, d, p);

            INFO("u " << u << " v " << v);
            REQUIRE(uv.has_value());
            CHECK(uv->x() == Approx(u).margin(1e-13));
            CHECK(uv->y() == Approx(v).margin(1e-13));
        }
    }
}

TEST_CASE("a trapezoid inverts through the quadratic branch", "[fit][locate]") {
    // Now the quadratic term is genuinely non-zero, so the quadratic has two
    // roots and the right one has to be picked. The oracle is the forward map
    // written out in this file.
    const Eigen::Vector2d a{0.0, 0.0};
    const Eigen::Vector2d b{4.0, 0.0};
    const Eigen::Vector2d c{3.0, 2.0};
    const Eigen::Vector2d d{1.0, 2.0};

    for (const double u : {0.0, 0.2, 0.5, 0.8, 1.0}) {
        for (const double v : {0.0, 0.3, 0.7, 1.0}) {
            const Eigen::Vector2d p = blend(a, b, c, d, u, v);

            const std::optional<Eigen::Vector2d> uv = n2s::fit::invert_bilinear(a, b, c, d, p);

            INFO("u " << u << " v " << v);
            REQUIRE(uv.has_value());
            CHECK(uv->x() == Approx(u).margin(1e-12));
            CHECK(uv->y() == Approx(v).margin(1e-12));
        }
    }
}

TEST_CASE("the inversion picks the root inside the quad", "[fit][locate]") {
    // A strongly tapered quad: the second root of the quadratic is real and
    // finite, and a solver that simply took the first one would return a
    // parameter pair outside [0,1] for an interior point.
    const Eigen::Vector2d a{0.0, 0.0};
    const Eigen::Vector2d b{8.0, 0.0};
    const Eigen::Vector2d c{4.2, 3.0};
    const Eigen::Vector2d d{3.8, 3.0};

    for (const double u : {0.05, 0.35, 0.65, 0.95}) {
        for (const double v : {0.05, 0.4, 0.9}) {
            const Eigen::Vector2d p = blend(a, b, c, d, u, v);

            const std::optional<Eigen::Vector2d> uv = n2s::fit::invert_bilinear(a, b, c, d, p);

            INFO("u " << u << " v " << v);
            REQUIRE(uv.has_value());
            CHECK(uv->x() == Approx(u).margin(1e-10));
            CHECK(uv->y() == Approx(v).margin(1e-10));
        }
    }
}

// ---------------------------------------------------------------------------
// Locating a point in a whole layout.
// ---------------------------------------------------------------------------

TEST_CASE("locating inverts the layout's own forward map", "[fit][locate]") {
    // The property that actually matters, stated as a round trip: whatever
    // location comes back, mapping it forward has to land on the point asked
    // for. Stated this way it stays true on a shared edge, where two different
    // locations are both correct answers.
    const DomainLayout layout = n2s::fit::grid_domain_layout(5, 4);
    const LayoutLocator locator{layout};

    for (std::size_t face = 0; face < layout.num_quads(); ++face) {
        for (const double u : {0.0, 0.25, 0.5, 0.75, 1.0}) {
            for (const double v : {0.0, 1.0 / 3.0, 0.5, 1.0}) {
                const Eigen::Vector2d p =
                    n2s::fit::bilinear_point(layout, static_cast<int>(face), u, v);

                const std::optional<n2s::LimitLocation> found = locator.locate(p);

                INFO("face " << face << " u " << u << " v " << v);
                REQUIRE(found.has_value());
                CHECK(found->u >= -n2s::tol::kLayoutContainment);
                CHECK(found->u <= 1.0 + n2s::tol::kLayoutContainment);
                CHECK(found->v >= -n2s::tol::kLayoutContainment);
                CHECK(found->v <= 1.0 + n2s::tol::kLayoutContainment);

                const Eigen::Vector2d back =
                    n2s::fit::bilinear_point(layout, found->face, found->u, found->v);
                CHECK((back - p).norm() < 1e-12);
            }
        }
    }
}

TEST_CASE("locating round-trips on the thesis layout", "[fit][locate]") {
    // The grid layout above is every quad a rectangle. The thesis layout is
    // hand-authored, with quads that are neither rectangular nor convex-looking
    // and four extraordinary vertices, which is what the inversion will
    // actually meet.
    const DomainLayout layout = n2s::fit::refine_quads(thesis_layout(), 3);
    const LayoutLocator locator{layout};

    std::size_t checked = 0;
    for (std::size_t face = 0; face < layout.num_quads(); ++face) {
        for (const double u : {0.1, 0.5, 0.9}) {
            for (const double v : {0.15, 0.45, 0.85}) {
                const Eigen::Vector2d p =
                    n2s::fit::bilinear_point(layout, static_cast<int>(face), u, v);

                const std::optional<n2s::LimitLocation> found = locator.locate(p);

                INFO("face " << face << " u " << u << " v " << v);
                REQUIRE(found.has_value());

                const Eigen::Vector2d back =
                    n2s::fit::bilinear_point(layout, found->face, found->u, found->v);
                CHECK((back - p).norm() < 1e-10);
                ++checked;
            }
        }
    }
    REQUIRE(checked > 0);
}

TEST_CASE("an interior sample is located in exactly one quad", "[fit][locate]") {
    // Interior means strictly inside one quad, so the answer is unambiguous and
    // the face is pinned, not merely the point it maps back to.
    const DomainLayout layout = n2s::fit::grid_domain_layout(4, 4);
    const LayoutLocator locator{layout};

    for (std::size_t face = 0; face < layout.num_quads(); ++face) {
        const Eigen::Vector2d p =
            n2s::fit::bilinear_point(layout, static_cast<int>(face), 0.5, 0.5);

        const std::optional<n2s::LimitLocation> found = locator.locate(p);
        INFO("face " << face);
        REQUIRE(found.has_value());
        CHECK(found->face == static_cast<int>(face));
        CHECK(found->u == Approx(0.5).margin(1e-12));
        CHECK(found->v == Approx(0.5).margin(1e-12));
    }
}

TEST_CASE("a point outside the layout is reported as outside", "[fit][locate]") {
    // Reported, not snapped to the nearest quad. A sample outside the layout is
    // a sample the fit has nothing to say about, and silently pulling it to the
    // boundary would add a row claiming the surface should pass somewhere it
    // was never asked to.
    const LayoutLocator locator{n2s::fit::grid_domain_layout(4, 4)};

    for (const Eigen::Vector2d& outside : {Eigen::Vector2d{-0.5, 0.5},
                                           Eigen::Vector2d{1.5, 0.5},
                                           Eigen::Vector2d{0.5, -0.25},
                                           Eigen::Vector2d{0.5, 2.0},
                                           Eigen::Vector2d{-3.0, -3.0}}) {
        INFO("point " << outside.x() << ", " << outside.y());
        CHECK_FALSE(locator.locate(outside).has_value());
    }
}

TEST_CASE("a point on a shared edge is located deterministically", "[fit][locate]") {
    // Both neighbouring quads are correct answers. Which one comes back must
    // not depend on anything but the layout and the point, or a fit built twice
    // from the same samples would produce two different systems.
    const DomainLayout layout = n2s::fit::grid_domain_layout(4, 4);
    const LayoutLocator locator{layout};

    const Eigen::Vector2d on_edge =
        n2s::fit::bilinear_point(layout, 0, 1.0, 0.5); // the edge face 0 shares

    const std::optional<n2s::LimitLocation> first = locator.locate(on_edge);
    REQUIRE(first.has_value());

    for (int repeat = 0; repeat < 4; ++repeat) {
        const LayoutLocator rebuilt{layout};
        const std::optional<n2s::LimitLocation> again = rebuilt.locate(on_edge);
        REQUIRE(again.has_value());
        CHECK(again->face == first->face);
        CHECK(again->u == Approx(first->u));
        CHECK(again->v == Approx(first->v));
    }

    const Eigen::Vector2d back = n2s::fit::bilinear_point(layout, first->face, first->u, first->v);
    CHECK((back - on_edge).norm() < 1e-12);
}

TEST_CASE("the located parameters agree with the metrics correspondence", "[fit][locate]") {
    // The locator and metrics::DomainMap have to be inverses of each other, not
    // merely similar. If they drift, the fit optimises against one
    // correspondence while the error is measured against another, and the
    // resulting parametric error is partly a disagreement between conventions.
    const DomainLayout layout = n2s::fit::grid_domain_layout(5, 5);
    const LayoutLocator locator{layout};
    const n2s::metrics::DomainMap forward = n2s::fit::bilinear_domain_map(layout);

    for (std::size_t face = 0; face < layout.num_quads(); ++face) {
        for (const double u : {0.2, 0.7}) {
            for (const double v : {0.3, 0.8}) {
                const n2s::LimitLocation asked{static_cast<int>(face), u, v};
                const Eigen::Vector2d p = forward(asked);

                const std::optional<n2s::LimitLocation> found = locator.locate(p);
                INFO("face " << face << " u " << u << " v " << v);
                REQUIRE(found.has_value());
                CHECK((forward(*found) - p).norm() < 1e-12);
            }
        }
    }
}

TEST_CASE("a layout with no quads is rejected", "[fit][locate]") {
    CHECK_THROWS_AS(LayoutLocator{DomainLayout{}}, std::invalid_argument);
}

TEST_CASE("the spatial index finds everything a full scan would", "[fit][locate]") {
    // The bucket grid is an optimisation, and the failure mode of a wrong one
    // is silent: a sample it fails to find is a row the least-squares system
    // never gets, so the fit is built from fewer points than it reports and
    // nothing anywhere says so. Compared here against scanning every quad.
    const DomainLayout layout = n2s::fit::refine_quads(thesis_layout(), 2);
    const LayoutLocator locator{layout};

    Eigen::Vector2d low = layout.vertices.front();
    Eigen::Vector2d high = layout.vertices.front();
    for (const Eigen::Vector2d& v : layout.vertices) {
        low = low.cwiseMin(v);
        high = high.cwiseMax(v);
    }

    // A deterministic lattice rather than random points, so a failure is
    // reproducible without recording a seed.
    const int side = 60;
    std::size_t inside = 0;
    for (int i = 0; i <= side; ++i) {
        for (int j = 0; j <= side; ++j) {
            const Eigen::Vector2d p{low.x() + (high.x() - low.x()) * static_cast<double>(i) / side,
                                    low.y() + (high.y() - low.y()) * static_cast<double>(j) / side};

            bool scan_found = false;
            for (std::size_t f = 0; f < layout.num_quads(); ++f) {
                const std::array<int, 4>& quad = layout.quads[f];
                const std::optional<Eigen::Vector2d> uv =
                    n2s::fit::invert_bilinear(layout.vertices[static_cast<std::size_t>(quad[0])],
                                              layout.vertices[static_cast<std::size_t>(quad[1])],
                                              layout.vertices[static_cast<std::size_t>(quad[2])],
                                              layout.vertices[static_cast<std::size_t>(quad[3])],
                                              p);
                if (uv.has_value() && uv->x() >= -n2s::tol::kLayoutContainment &&
                    uv->x() <= 1.0 + n2s::tol::kLayoutContainment &&
                    uv->y() >= -n2s::tol::kLayoutContainment &&
                    uv->y() <= 1.0 + n2s::tol::kLayoutContainment) {
                    scan_found = true;
                    break;
                }
            }

            const std::optional<n2s::LimitLocation> found = locator.locate(p);
            INFO("point " << p.x() << ", " << p.y());
            REQUIRE(found.has_value() == scan_found);

            if (found.has_value()) {
                ++inside;
                const Eigen::Vector2d back =
                    n2s::fit::bilinear_point(layout, found->face, found->u, found->v);
                CHECK((back - p).norm() < 1e-9);
            }
        }
    }

    // The layout is not a rectangle, so a lattice over its bounding box must
    // land both inside and outside it. If everything landed one way the test
    // above would be vacuous.
    INFO("inside " << inside << " of " << (side + 1) * (side + 1));
    CHECK(inside > 0);
    CHECK(inside < static_cast<std::size_t>((side + 1) * (side + 1)));
}
