#include "n2s/nurbs/surface.hpp"

#include "support/analytic_nurbs.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

using Catch::Approx;
using n2s::KnotVector;
using n2s::NurbsSurface;
using n2s::SurfaceDerivatives;

TEST_CASE("a bilinear patch reproduces bilinear interpolation", "[nurbs][surface]") {
    const KnotVector knots{1, {0, 0, 1, 1}};
    const std::vector<Eigen::Vector3d> net{
        {0, 0, 0}, {0, 1, 1}, {1, 0, 2}, {1, 1, 5}}; // (0,0) (0,1) (1,0) (1,1)
    const NurbsSurface patch = NurbsSurface::bspline(knots, knots, net);

    for (int i = 0; i <= 10; ++i) {
        for (int j = 0; j <= 10; ++j) {
            const double u = static_cast<double>(i) / 10.0;
            const double v = static_cast<double>(j) / 10.0;

            const Eigen::Vector3d expected = (1 - u) * (1 - v) * net[0] + (1 - u) * v * net[1] +
                                             u * (1 - v) * net[2] + u * v * net[3];

            INFO("u = " << u << ", v = " << v);
            CHECK((patch.evaluate(u, v) - expected).norm() < 1e-14);
        }
    }
}

TEST_CASE("a bicubic Bezier patch matches closed-form Bernstein evaluation", "[nurbs][surface]") {
    const std::vector<Eigen::Vector3d> net = n2s::testing::wavy_bicubic_net();
    const NurbsSurface patch = n2s::testing::bicubic_bezier(net);

    REQUIRE_FALSE(patch.is_rational());
    REQUIRE(patch.num_u() == 4);
    REQUIRE(patch.num_v() == 4);

    for (int i = 0; i <= 20; ++i) {
        for (int j = 0; j <= 20; ++j) {
            const double u = static_cast<double>(i) / 20.0;
            const double v = static_cast<double>(j) / 20.0;

            INFO("u = " << u << ", v = " << v);
            CHECK((patch.evaluate(u, v) - n2s::testing::bernstein_patch(net, u, v)).norm() < 1e-14);
        }
    }
}

TEST_CASE("a rational cylinder is exact", "[nurbs][surface]") {
    const double radius = 3.0;
    const double height = 2.0;
    const NurbsSurface cyl = n2s::testing::cylinder(radius, height);

    double worst_radius = 0.0;
    double worst_height = 0.0;

    for (int i = 0; i <= 120; ++i) {
        for (int j = 0; j <= 10; ++j) {
            const double u = static_cast<double>(i) / 120.0;
            const double v = static_cast<double>(j) / 10.0;
            const Eigen::Vector3d p = cyl.evaluate(u, v);

            worst_radius = std::max(worst_radius, std::abs(std::hypot(p.x(), p.y()) - radius));
            worst_height = std::max(worst_height, std::abs(p.z() - v * height));
        }
    }

    INFO("worst radius error " << worst_radius << ", worst height error " << worst_height);
    CHECK(worst_radius < 1e-12);
    CHECK(worst_height < 1e-12);
}

TEST_CASE("a rational sphere of revolution is exact", "[nurbs][surface]") {
    const double radius = 1.75;
    const NurbsSurface ball = n2s::testing::sphere(radius);

    double worst = 0.0;
    for (int i = 0; i <= 60; ++i) {
        for (int j = 0; j <= 60; ++j) {
            const double u = static_cast<double>(i) / 60.0;
            const double v = static_cast<double>(j) / 60.0;
            worst = std::max(worst, std::abs(ball.evaluate(u, v).norm() - radius));
        }
    }

    INFO("worst radial error " << worst);
    CHECK(worst < 1e-12);
}

TEST_CASE("surface derivatives match central finite differences", "[nurbs][surface]") {
    // A rational, non-symmetric patch: uniform weights would hide sign errors
    // in the quotient rule, and symmetry would hide a u/v transposition.
    const KnotVector knots_u{2, {0, 0, 0, 0.5, 1, 1, 1}};
    const KnotVector knots_v{3, {0, 0, 0, 0, 0.4, 1, 1, 1, 1}};

    std::vector<Eigen::Vector3d> net;
    std::vector<double> weights;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 5; ++j) {
            const double x = 0.7 * static_cast<double>(i);
            const double y = 0.5 * static_cast<double>(j);
            net.emplace_back(x, y, 0.3 * std::sin(1.1 * x + 0.6) + 0.2 * std::cos(0.9 * y));
            weights.push_back(1.0 + 0.35 * static_cast<double>((i * 5 + j) % 4));
        }
    }
    const NurbsSurface patch{knots_u, knots_v, net, weights};
    REQUIRE(patch.is_rational());

    const double h = 1e-5;
    for (int i = 1; i < 12; ++i) {
        for (int j = 1; j < 12; ++j) {
            // Offset off the grid, and skip knots outright: across a simple
            // interior knot of a degree p surface only the first p-1
            // derivatives are continuous, so a central difference straddling
            // one measures the discontinuity, not the analytic derivative.
            const double u = (static_cast<double>(i) + 0.41) / 12.0;
            const double v = (static_cast<double>(j) + 0.29) / 12.0;
            if (knots_u.multiplicity(u) > 0 || knots_v.multiplicity(v) > 0) {
                continue;
            }

            const SurfaceDerivatives d = patch.derivatives(u, v, 2);

            const Eigen::Vector3d fd_u =
                (patch.evaluate(u + h, v) - patch.evaluate(u - h, v)) / (2.0 * h);
            const Eigen::Vector3d fd_v =
                (patch.evaluate(u, v + h) - patch.evaluate(u, v - h)) / (2.0 * h);
            const Eigen::Vector3d fd_uu =
                (patch.evaluate(u + h, v) - 2.0 * patch.evaluate(u, v) + patch.evaluate(u - h, v)) /
                (h * h);
            const Eigen::Vector3d fd_vv =
                (patch.evaluate(u, v + h) - 2.0 * patch.evaluate(u, v) + patch.evaluate(u, v - h)) /
                (h * h);
            const Eigen::Vector3d fd_uv =
                (patch.evaluate(u + h, v + h) - patch.evaluate(u + h, v - h) -
                 patch.evaluate(u - h, v + h) + patch.evaluate(u - h, v - h)) /
                (4.0 * h * h);

            INFO("u = " << u << ", v = " << v);
            CHECK((d.position() - patch.evaluate(u, v)).norm() < 1e-14);
            CHECK((d.du() - fd_u).norm() < 1e-6);
            CHECK((d.dv() - fd_v).norm() < 1e-6);
            CHECK((d.duu() - fd_uu).norm() < 1e-4);
            CHECK((d.dvv() - fd_vv).norm() < 1e-4);
            CHECK((d.duv() - fd_uv).norm() < 1e-4);
        }
    }
}

TEST_CASE("a clamped patch interpolates its corner control points", "[nurbs][surface]") {
    const std::vector<Eigen::Vector3d> net = n2s::testing::wavy_bicubic_net();
    const NurbsSurface patch = n2s::testing::bicubic_bezier(net);

    CHECK((patch.evaluate(0.0, 0.0) - patch.control_point(0, 0)).norm() < 1e-14);
    CHECK((patch.evaluate(0.0, 1.0) - patch.control_point(0, 3)).norm() < 1e-14);
    CHECK((patch.evaluate(1.0, 0.0) - patch.control_point(3, 0)).norm() < 1e-14);
    CHECK((patch.evaluate(1.0, 1.0) - patch.control_point(3, 3)).norm() < 1e-14);
}

TEST_CASE("invalid surface inputs are rejected", "[nurbs][surface]") {
    const KnotVector knots = KnotVector::uniform_clamped(2, 4);
    const std::vector<Eigen::Vector3d> net(16, Eigen::Vector3d::Zero());

    SECTION("net size must match the knot vectors") {
        const std::vector<Eigen::Vector3d> wrong(15, Eigen::Vector3d::Zero());
        CHECK_THROWS_AS(NurbsSurface::bspline(knots, knots, wrong), std::invalid_argument);
    }

    SECTION("weights must be positive") {
        std::vector<double> weights(16, 1.0);
        weights[7] = -0.5;
        CHECK_THROWS_AS((NurbsSurface{knots, knots, net, weights}), std::invalid_argument);
    }
}
