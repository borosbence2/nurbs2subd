#include "n2s/fit/layout.hpp"
#include "n2s/trim/cases.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <stdexcept>
#include <vector>

using Catch::Approx;
using n2s::LimitLocation;
using n2s::fit::DomainLayout;

TEST_CASE("the grid layout agrees with the grid domain map", "[fit][layout]") {
    // Two independent routes to the same correspondence. They have to agree,
    // or a metric computed through one and a layout built through the other
    // would quietly disagree about which point corresponds to which.
    const int rows = 5;
    const int columns = 4;

    const DomainLayout layout = n2s::fit::grid_domain_layout(rows, columns);
    const n2s::metrics::DomainMap from_layout = n2s::fit::bilinear_domain_map(layout);
    const n2s::metrics::DomainMap from_grid = n2s::metrics::grid_domain_map(rows, columns);

    REQUIRE(layout.num_quads() == static_cast<std::size_t>((rows - 1) * (columns - 1)));

    for (int face = 0; face < static_cast<int>(layout.num_quads()); ++face) {
        for (const double u : {0.0, 0.25, 0.5, 1.0}) {
            for (const double v : {0.0, 0.5, 0.75, 1.0}) {
                const LimitLocation location{face, u, v};
                INFO("face " << face << " at (" << u << ", " << v << ")");
                REQUIRE((from_layout(location) - from_grid(location)).norm() < 1e-15);
            }
        }
    }
}

TEST_CASE("the bilinear map sends face corners to layout vertices", "[fit][layout]") {
    DomainLayout layout;
    layout.vertices = {{0.1, 0.2}, {0.8, 0.1}, {0.9, 0.7}, {0.2, 0.9}};
    layout.quads = {{0, 1, 2, 3}};

    const n2s::metrics::DomainMap map = n2s::fit::bilinear_domain_map(layout);

    // (0,0) is the first vertex, u runs to the second, v to the fourth. If this
    // ordering is wrong the whole correspondence is transposed, and a
    // transposed correspondence produces a plausible-looking parametric error
    // that measures nothing.
    CHECK((map({0, 0.0, 0.0}) - layout.vertices[0]).norm() < 1e-15);
    CHECK((map({0, 1.0, 0.0}) - layout.vertices[1]).norm() < 1e-15);
    CHECK((map({0, 1.0, 1.0}) - layout.vertices[2]).norm() < 1e-15);
    CHECK((map({0, 0.0, 1.0}) - layout.vertices[3]).norm() < 1e-15);

    // The centre is the average of the four corners.
    const Eigen::Vector2d centre =
        (layout.vertices[0] + layout.vertices[1] + layout.vertices[2] + layout.vertices[3]) / 4.0;
    CHECK((map({0, 0.5, 0.5}) - centre).norm() < 1e-15);
}

TEST_CASE("the correspondence is continuous across a shared edge", "[fit][layout]") {
    // Two quads meeting along an edge must agree along it, or the layout's
    // correspondence has a seam and every metric computed through it jumps
    // there.
    DomainLayout layout;
    layout.vertices = {{0.0, 0.0}, {0.5, 0.05}, {0.5, 0.95}, {0.0, 1.0}, {1.0, 0.1}, {1.0, 0.9}};
    layout.quads = {{0, 1, 2, 3}, {1, 4, 5, 2}};

    const n2s::metrics::DomainMap map = n2s::fit::bilinear_domain_map(layout);

    for (int i = 0; i <= 8; ++i) {
        const double t = static_cast<double>(i) / 8.0;
        // The shared edge is quad 0's u = 1 side and quad 1's u = 0 side, both
        // walked by v.
        INFO("t = " << t);
        CHECK((map({0, 1.0, t}) - map({1, 0.0, t})).norm() < 1e-15);
    }
}

TEST_CASE("lifting places control points on the surface", "[fit][layout]") {
    const n2s::NurbsSurface surface = n2s::cases::saddle(1.0, 0.6);
    const DomainLayout layout = n2s::fit::grid_domain_layout(4, 4);

    const n2s::ControlMesh mesh = n2s::fit::lift(layout, surface);

    REQUIRE(mesh.num_vertices() == layout.num_vertices());
    REQUIRE(mesh.num_quads() == layout.num_quads());

    for (std::size_t i = 0; i < layout.num_vertices(); ++i) {
        const Eigen::Vector2d& p = layout.vertices[i];
        INFO("vertex " << i);
        CHECK((mesh.vertices()[i] - surface.evaluate(p.x(), p.y())).norm() < 1e-15);
    }

    // grid_layout is exactly lift(grid_domain_layout(...)).
    const n2s::ControlMesh direct = n2s::fit::grid_layout(surface, 4, 4);
    for (std::size_t i = 0; i < mesh.num_vertices(); ++i) {
        CHECK((mesh.vertices()[i] - direct.vertices()[i]).norm() < 1e-15);
    }
}

TEST_CASE("an irregular layout still yields a usable correspondence", "[fit][layout]") {
    // The layouts that matter are not grids. This one is a three-quad fan with
    // an interior vertex, the shape a hand-authored layout actually takes.
    DomainLayout layout;
    layout.vertices = {{0.5, 0.5},
                       {1.0, 0.5},
                       {1.0, 1.0},
                       {0.5, 1.0},
                       {0.0, 1.0},
                       {0.0, 0.5},
                       {0.0, 0.0},
                       {0.5, 0.0}};
    layout.quads = {{0, 1, 2, 3}, {0, 3, 4, 5}, {0, 5, 6, 7}};

    const n2s::metrics::DomainMap map = n2s::fit::bilinear_domain_map(layout);
    const n2s::NurbsSurface surface = n2s::cases::plane();
    const n2s::ControlMesh mesh = n2s::fit::lift(layout, surface);

    CHECK(mesh.num_quads() == 3);
    // The shared centre vertex is (0,0) of every quad.
    for (int face = 0; face < 3; ++face) {
        CHECK((map({face, 0.0, 0.0}) - layout.vertices[0]).norm() < 1e-15);
    }
}

TEST_CASE("invalid layouts are rejected", "[fit][layout]") {
    SECTION("empty") {
        CHECK_THROWS_AS(n2s::fit::bilinear_domain_map(DomainLayout{}), std::invalid_argument);
    }

    SECTION("an out-of-range vertex index") {
        DomainLayout layout;
        layout.vertices = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
        layout.quads = {{0, 1, 2, 9}};
        CHECK_THROWS_AS(n2s::fit::bilinear_domain_map(layout), std::invalid_argument);
    }

    SECTION("a face index the layout does not have") {
        const n2s::metrics::DomainMap map =
            n2s::fit::bilinear_domain_map(n2s::fit::grid_domain_layout(2, 2));
        CHECK_THROWS_AS(map({5, 0.5, 0.5}), std::invalid_argument);
    }
}
