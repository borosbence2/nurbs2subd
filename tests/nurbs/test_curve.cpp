#include "n2s/nurbs/curve.hpp"

#include "support/analytic_nurbs.hpp"

#include <Eigen/Geometry> // Vector3d::cross
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

using Catch::Approx;
using n2s::KnotVector;
using n2s::NurbsCurve;

TEST_CASE("a rational quadratic reproduces a circle exactly", "[nurbs][curve]") {
    // The headline oracle from the plan: the conic is exact under the rational
    // representation, so the radius error must be at machine precision, not at
    // approximation precision.
    for (const double radius : {1.0, 0.25, 137.5}) {
        const NurbsCurve circle = n2s::testing::full_circle(radius);
        REQUIRE(circle.is_rational());

        double worst = 0.0;
        for (int i = 0; i <= 500; ++i) {
            const double u = static_cast<double>(i) / 500.0;
            const double r = circle.evaluate(u).norm();
            worst = std::max(worst, std::abs(r - radius) / radius);
        }
        INFO("radius " << radius << ", worst relative error " << worst);
        CHECK(worst < 1e-12);
    }
}

TEST_CASE("the circle parameterisation hits the expected quadrant corners", "[nurbs][curve]") {
    const NurbsCurve circle = n2s::testing::full_circle(1.0);

    // Each quarter of the parameter range is one 90 degree arc.
    const Eigen::Vector3d at_0 = circle.evaluate(0.0);
    const Eigen::Vector3d at_quarter = circle.evaluate(0.25);
    const Eigen::Vector3d at_half = circle.evaluate(0.5);

    CHECK(at_0.x() == Approx(1.0).margin(1e-14));
    CHECK(at_0.y() == Approx(0.0).margin(1e-14));
    CHECK(at_quarter.x() == Approx(0.0).margin(1e-14));
    CHECK(at_quarter.y() == Approx(1.0).margin(1e-14));
    CHECK(at_half.x() == Approx(-1.0).margin(1e-14));
    CHECK(at_half.y() == Approx(0.0).margin(1e-14));
}

TEST_CASE("a quarter circle has constant speed times curvature", "[nurbs][curve]") {
    // For any regular parameterisation of a circle of radius r the curvature
    // |C' x C''| / |C'|^3 equals 1/r everywhere, even though a rational
    // quadratic is emphatically not arc-length parameterised.
    const double radius = 2.5;
    const NurbsCurve arc = n2s::testing::quarter_circle(radius);

    for (int i = 0; i <= 50; ++i) {
        const double u = static_cast<double>(i) / 50.0;
        const std::vector<Eigen::Vector3d> d = arc.derivatives(u, 2);

        const double speed = d[1].norm();
        REQUIRE(speed > 1e-9);
        const double curvature = d[1].cross(d[2]).norm() / (speed * speed * speed);

        INFO("u = " << u);
        CHECK(curvature == Approx(1.0 / radius).epsilon(1e-10));
    }
}

TEST_CASE("curve derivatives match central finite differences", "[nurbs][curve]") {
    const KnotVector knots{3, {0, 0, 0, 0, 0.3, 0.6, 1, 1, 1, 1}};
    const std::vector<Eigen::Vector3d> points{
        {0, 0, 0}, {1, 2, 0}, {2, -1, 1}, {3, 0.5, 2}, {4, 1, -1}, {5, 0, 0}};
    const std::vector<double> weights{1.0, 2.0, 0.5, 1.5, 1.0, 3.0};
    const NurbsCurve curve{knots, points, weights};

    const double h = 1e-5;
    for (int i = 1; i < 40; ++i) {
        // Offset off the grid, and skip knots outright: a degree p curve is only
        // C^(p-1) at a simple interior knot, so a central difference straddling
        // one compares two different polynomials and reports an error that has
        // nothing to do with the correctness of the analytic derivative.
        const double u = (static_cast<double>(i) + 0.37) / 40.0;
        if (knots.multiplicity(u) > 0) {
            continue;
        }

        const std::vector<Eigen::Vector3d> d = curve.derivatives(u, 2);

        const Eigen::Vector3d fd1 = (curve.evaluate(u + h) - curve.evaluate(u - h)) / (2.0 * h);
        const Eigen::Vector3d fd2 =
            (curve.evaluate(u + h) - 2.0 * curve.evaluate(u) + curve.evaluate(u - h)) / (h * h);

        INFO("u = " << u);
        CHECK((d[0] - curve.evaluate(u)).norm() < 1e-14);
        CHECK((d[1] - fd1).norm() < 1e-6);
        CHECK((d[2] - fd2).norm() < 1e-4);
    }
}

TEST_CASE("a clamped curve interpolates its end control points", "[nurbs][curve]") {
    const KnotVector knots = KnotVector::uniform_clamped(3, 6);
    const std::vector<Eigen::Vector3d> points{
        {0, 0, 0}, {1, 2, 0}, {2, -1, 1}, {3, 0.5, 2}, {4, 1, -1}, {5, 0, 0}};
    const NurbsCurve curve = NurbsCurve::bspline(knots, points);

    CHECK((curve.evaluate(0.0) - points.front()).norm() < 1e-14);
    CHECK((curve.evaluate(1.0) - points.back()).norm() < 1e-14);
    CHECK_FALSE(curve.is_rational());
}

TEST_CASE("invalid curve inputs are rejected", "[nurbs][curve]") {
    const KnotVector knots = KnotVector::uniform_clamped(3, 5);
    const std::vector<Eigen::Vector3d> five(5, Eigen::Vector3d::Zero());

    SECTION("control point count must match the knot vector") {
        const std::vector<Eigen::Vector3d> four(4, Eigen::Vector3d::Zero());
        CHECK_THROWS_AS(NurbsCurve::bspline(knots, four), std::invalid_argument);
    }

    SECTION("weight count must match the control point count") {
        CHECK_THROWS_AS((NurbsCurve{knots, five, std::vector<double>(4, 1.0)}),
                        std::invalid_argument);
    }

    SECTION("weights must be positive") {
        std::vector<double> weights(5, 1.0);
        weights[2] = 0.0;
        CHECK_THROWS_AS((NurbsCurve{knots, five, weights}), std::invalid_argument);

        weights[2] = -1.0;
        CHECK_THROWS_AS((NurbsCurve{knots, five, weights}), std::invalid_argument);
    }
}
