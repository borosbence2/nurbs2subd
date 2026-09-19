#include "n2s/nurbs/differential.hpp"
#include "n2s/nurbs/surface.hpp"

#include "support/analytic_nurbs.hpp"

#include <Eigen/Geometry>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <optional>
#include <vector>

using Catch::Approx;
using n2s::KnotVector;
using n2s::NurbsSurface;
using n2s::SurfaceDerivatives;

TEST_CASE("a plane is flat everywhere", "[nurbs][curvature]") {
    const KnotVector knots{1, {0, 0, 1, 1}};
    const std::vector<Eigen::Vector3d> net{{0, 0, 1}, {0, 3, 1}, {2, 0, 1}, {2, 3, 1}};
    const NurbsSurface plane = NurbsSurface::bspline(knots, knots, net);

    for (int i = 0; i <= 5; ++i) {
        for (int j = 0; j <= 5; ++j) {
            const double u = static_cast<double>(i) / 5.0;
            const double v = static_cast<double>(j) / 5.0;
            const SurfaceDerivatives d = plane.derivatives(u, v, 2);

            const auto curvature = n2s::surface_curvature(d);
            REQUIRE(curvature.has_value());

            INFO("u = " << u << ", v = " << v);
            CHECK(curvature->mean == Approx(0.0).margin(1e-14));
            CHECK(curvature->gaussian == Approx(0.0).margin(1e-14));

            const auto normal = n2s::surface_normal(d);
            REQUIRE(normal.has_value());
            CHECK(std::abs(normal->z()) == Approx(1.0).margin(1e-14));
        }
    }
}

TEST_CASE("a sphere has curvature 1/r everywhere", "[nurbs][curvature]") {
    // The analytic oracle from the plan. Signs follow the orientation of
    // (Su x Sv), which for this surface of revolution points inward, so the
    // magnitudes are what is asserted.
    const double radius = 1.4;
    const NurbsSurface ball = n2s::testing::sphere(radius);

    for (int i = 0; i <= 24; ++i) {
        for (int j = 1; j < 24; ++j) { // j = 0 and j = 24 are the poles
            const double u = static_cast<double>(i) / 24.0;
            const double v = static_cast<double>(j) / 24.0;
            const SurfaceDerivatives d = ball.derivatives(u, v, 2);

            const auto curvature = n2s::surface_curvature(d);
            REQUIRE(curvature.has_value());

            INFO("u = " << u << ", v = " << v);
            CHECK(std::abs(curvature->mean) == Approx(1.0 / radius).epsilon(1e-9));
            CHECK(curvature->gaussian == Approx(1.0 / (radius * radius)).epsilon(1e-9));

            // Both principal curvatures are equal on a sphere: every point is
            // umbilic. That is also the worst case for their accuracy --
            // k = H +- sqrt(H^2 - K) with H^2 - K rounding to about machine
            // epsilon, so the square root returns roughly sqrt(eps) * |H| and
            // only half the digits of H and K survive. A looser tolerance here
            // is the honest one, not a fudge.
            CHECK(std::abs(curvature->k1) == Approx(1.0 / radius).epsilon(1e-7));
            CHECK(std::abs(curvature->k2) == Approx(1.0 / radius).epsilon(1e-7));

            // The unit normal is radial.
            const auto normal = n2s::surface_normal(d);
            REQUIRE(normal.has_value());
            const Eigen::Vector3d radial = ball.evaluate(u, v).normalized();
            CHECK(std::abs(normal->dot(radial)) == Approx(1.0).epsilon(1e-9));
        }
    }
}

TEST_CASE("a cylinder has mean curvature 1/(2r) and zero Gaussian curvature",
          "[nurbs][curvature]") {
    const double radius = 2.0;
    const NurbsSurface cyl = n2s::testing::cylinder(radius, 3.0);

    for (int i = 0; i <= 24; ++i) {
        for (int j = 0; j <= 6; ++j) {
            const double u = static_cast<double>(i) / 24.0;
            const double v = static_cast<double>(j) / 6.0;
            const SurfaceDerivatives d = cyl.derivatives(u, v, 2);

            const auto curvature = n2s::surface_curvature(d);
            REQUIRE(curvature.has_value());

            INFO("u = " << u << ", v = " << v);
            CHECK(std::abs(curvature->mean) == Approx(1.0 / (2.0 * radius)).epsilon(1e-9));
            CHECK(curvature->gaussian == Approx(0.0).margin(1e-10));

            // A cylinder is developable: one principal curvature is exactly zero.
            const double smaller = std::min(std::abs(curvature->k1), std::abs(curvature->k2));
            const double larger = std::max(std::abs(curvature->k1), std::abs(curvature->k2));
            CHECK(smaller == Approx(0.0).margin(1e-10));
            CHECK(larger == Approx(1.0 / radius).epsilon(1e-9));
        }
    }
}

TEST_CASE("the first fundamental form measures lengths on the surface", "[nurbs][curvature]") {
    const NurbsSurface ball = n2s::testing::sphere(1.0);
    const SurfaceDerivatives d = ball.derivatives(0.3, 0.45, 2);

    const n2s::FirstFundamentalForm form = n2s::first_fundamental_form(d);

    CHECK(form.e == Approx(d.du().squaredNorm()));
    CHECK(form.f == Approx(d.du().dot(d.dv())));
    CHECK(form.g == Approx(d.dv().squaredNorm()));

    // EG - F^2 is the squared area of the parallelogram spanned by the partials.
    CHECK(form.determinant() == Approx(d.du().cross(d.dv()).squaredNorm()));
}

TEST_CASE("curvature is undefined at a degenerate point", "[nurbs][curvature]") {
    // v = 0 is the north pole of the sphere of revolution: the whole u-row of
    // the control net collapses to one point, Su vanishes, and there is no
    // tangent plane to speak of.
    const NurbsSurface ball = n2s::testing::sphere(1.0);
    const SurfaceDerivatives d = ball.derivatives(0.3, 0.0, 2);

    REQUIRE(d.du().norm() < 1e-12);
    CHECK_FALSE(n2s::surface_normal(d).has_value());
    CHECK_FALSE(n2s::surface_curvature(d).has_value());
}

TEST_CASE("principal curvatures reproduce mean and Gaussian curvature", "[nurbs][curvature]") {
    // On an arbitrary patch the identities H = (k1 + k2)/2 and K = k1 * k2 must
    // hold whatever the values are.
    const NurbsSurface patch = n2s::testing::bicubic_bezier(n2s::testing::wavy_bicubic_net());

    for (int i = 1; i < 8; ++i) {
        for (int j = 1; j < 8; ++j) {
            const double u = static_cast<double>(i) / 8.0;
            const double v = static_cast<double>(j) / 8.0;
            const SurfaceDerivatives d = patch.derivatives(u, v, 2);

            const auto c = n2s::surface_curvature(d);
            REQUIRE(c.has_value());

            INFO("u = " << u << ", v = " << v);
            CHECK((c->k1 + c->k2) / 2.0 == Approx(c->mean).epsilon(1e-12));
            CHECK(c->k1 * c->k2 == Approx(c->gaussian).margin(1e-12));
            CHECK(c->k1 >= c->k2);
        }
    }
}
