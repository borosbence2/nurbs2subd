#include "n2s/nurbs/curve.hpp"
#include "n2s/nurbs/knot_insertion.hpp"
#include "n2s/nurbs/surface.hpp"

#include "support/analytic_nurbs.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

using Catch::Approx;
using n2s::Direction;
using n2s::KnotVector;
using n2s::NurbsCurve;
using n2s::NurbsSurface;

namespace {

/// The defining property of knot insertion: the geometry does not move.
void check_same_curve(const NurbsCurve& before, const NurbsCurve& after, double tolerance) {
    for (int i = 0; i <= 200; ++i) {
        const double t = static_cast<double>(i) / 200.0;
        const double u = before.knots().domain_start() +
                         t * (before.knots().domain_end() - before.knots().domain_start());
        INFO("u = " << u);
        REQUIRE((before.evaluate(u) - after.evaluate(u)).norm() < tolerance);
    }
}

void check_same_surface(const NurbsSurface& before, const NurbsSurface& after, double tolerance) {
    for (int i = 0; i <= 30; ++i) {
        for (int j = 0; j <= 30; ++j) {
            const double u = static_cast<double>(i) / 30.0;
            const double v = static_cast<double>(j) / 30.0;
            INFO("u = " << u << ", v = " << v);
            REQUIRE((before.evaluate(u, v) - after.evaluate(u, v)).norm() < tolerance);
        }
    }
}

} // namespace

TEST_CASE("inserting a knot into a curve does not move it", "[nurbs][insertion]") {
    const KnotVector knots{3, {0, 0, 0, 0, 0.3, 0.7, 1, 1, 1, 1}};
    const std::vector<Eigen::Vector3d> points{
        {0, 0, 0}, {1, 2, 0}, {2, -1, 1}, {3, 0.5, 2}, {4, 1, -1}, {5, 0, 0}};
    const std::vector<double> weights{1.0, 2.0, 0.5, 1.5, 1.0, 3.0};
    const NurbsCurve curve{knots, points, weights};

    SECTION("a new knot value") {
        const NurbsCurve refined = n2s::insert_knot(curve, 0.5);

        CHECK(refined.num_control_points() == curve.num_control_points() + 1);
        CHECK(refined.knots().multiplicity(0.5) == 1);
        CHECK(refined.degree() == 3);
        check_same_curve(curve, refined, 1e-12);
    }

    SECTION("raising the multiplicity of an existing knot") {
        const NurbsCurve refined = n2s::insert_knot(curve, 0.3, 2);

        CHECK(refined.num_control_points() == curve.num_control_points() + 2);
        CHECK(refined.knots().multiplicity(0.3) == 3);
        check_same_curve(curve, refined, 1e-12);
    }

    SECTION("several insertions in a row") {
        NurbsCurve refined = curve;
        for (const double u : {0.1, 0.45, 0.62, 0.9}) {
            refined = n2s::insert_knot(refined, u);
        }

        CHECK(refined.num_control_points() == curve.num_control_points() + 4);
        check_same_curve(curve, refined, 1e-12);
    }
}

TEST_CASE("insertion preserves an exact rational circle", "[nurbs][insertion]") {
    // Insertion is a sequence of affine combinations of the *weighted* control
    // points. Doing it on the Cartesian points instead would leave the circle
    // visibly dented, so this is the test that catches that mistake.
    const double radius = 1.0;
    const NurbsCurve circle = n2s::testing::full_circle(radius);

    NurbsCurve refined = circle;
    for (const double u : {0.1, 0.35, 0.6, 0.85}) {
        refined = n2s::insert_knot(refined, u);
    }

    double worst = 0.0;
    for (int i = 0; i <= 500; ++i) {
        const double u = static_cast<double>(i) / 500.0;
        worst = std::max(worst, std::abs(refined.evaluate(u).norm() - radius));
    }

    INFO("worst radial error after insertion " << worst);
    CHECK(worst < 1e-12);
    CHECK(refined.is_rational());
}

TEST_CASE("raising a knot to full multiplicity puts a control point on the curve",
          "[nurbs][insertion]") {
    // Bezier extraction: once an interior knot reaches multiplicity p, the
    // curve interpolates the corresponding control point.
    const KnotVector knots{3, {0, 0, 0, 0, 0.4, 1, 1, 1, 1}};
    const std::vector<Eigen::Vector3d> points{
        {0, 0, 0}, {1, 2, 0}, {2, -1, 1}, {3, 0.5, 2}, {4, 1, -1}};
    const NurbsCurve curve = NurbsCurve::bspline(knots, points);

    const NurbsCurve refined = n2s::insert_knot(curve, 0.4, 2);
    REQUIRE(refined.knots().multiplicity(0.4) == 3);

    const Eigen::Vector3d on_curve = curve.evaluate(0.4);
    double closest = std::numeric_limits<double>::infinity();
    for (const Eigen::Vector3d& p : refined.control_points()) {
        closest = std::min(closest, (p - on_curve).norm());
    }

    INFO("closest control point is " << closest << " away from C(0.4)");
    CHECK(closest < 1e-12);
}

TEST_CASE("inserting a knot into a surface does not move it", "[nurbs][insertion]") {
    const NurbsSurface patch = n2s::testing::bicubic_bezier(n2s::testing::wavy_bicubic_net());

    SECTION("in the u direction") {
        const NurbsSurface refined = n2s::insert_knot(patch, Direction::U, 0.5);

        CHECK(refined.num_u() == patch.num_u() + 1);
        CHECK(refined.num_v() == patch.num_v());
        CHECK(refined.knots_u().multiplicity(0.5) == 1);
        check_same_surface(patch, refined, 1e-12);
    }

    SECTION("in the v direction") {
        const NurbsSurface refined = n2s::insert_knot(patch, Direction::V, 0.25, 2);

        CHECK(refined.num_u() == patch.num_u());
        CHECK(refined.num_v() == patch.num_v() + 2);
        CHECK(refined.knots_v().multiplicity(0.25) == 2);
        check_same_surface(patch, refined, 1e-12);
    }

    SECTION("in both directions") {
        NurbsSurface refined = n2s::insert_knot(patch, Direction::U, 0.35);
        refined = n2s::insert_knot(refined, Direction::V, 0.6);
        refined = n2s::insert_knot(refined, Direction::U, 0.8);

        CHECK(refined.num_u() == patch.num_u() + 2);
        CHECK(refined.num_v() == patch.num_v() + 1);
        check_same_surface(patch, refined, 1e-12);
    }
}

TEST_CASE("insertion preserves an exact rational cylinder", "[nurbs][insertion]") {
    const double radius = 2.5;
    const NurbsSurface cyl = n2s::testing::cylinder(radius, 1.0);

    NurbsSurface refined = n2s::insert_knot(cyl, Direction::U, 0.125);
    refined = n2s::insert_knot(refined, Direction::V, 0.5);

    double worst = 0.0;
    for (int i = 0; i <= 120; ++i) {
        for (int j = 0; j <= 8; ++j) {
            const double u = static_cast<double>(i) / 120.0;
            const double v = static_cast<double>(j) / 8.0;
            const Eigen::Vector3d p = refined.evaluate(u, v);
            worst = std::max(worst, std::abs(std::hypot(p.x(), p.y()) - radius));
        }
    }

    INFO("worst radial error " << worst);
    CHECK(worst < 1e-12);
}

TEST_CASE("invalid insertions are rejected", "[nurbs][insertion]") {
    const KnotVector knots{3, {0, 0, 0, 0, 0.4, 1, 1, 1, 1}};
    const std::vector<Eigen::Vector3d> points{
        {0, 0, 0}, {1, 2, 0}, {2, -1, 1}, {3, 0.5, 2}, {4, 1, -1}};
    const NurbsCurve curve = NurbsCurve::bspline(knots, points);

    SECTION("outside the domain") {
        CHECK_THROWS_AS(n2s::insert_knot(curve, 1.5), std::invalid_argument);
        CHECK_THROWS_AS(n2s::insert_knot(curve, -0.2), std::invalid_argument);
    }

    SECTION("beyond the degree") {
        // 0.4 already has multiplicity 1, so three more would make it 4 > p.
        CHECK_THROWS_AS(n2s::insert_knot(curve, 0.4, 3), std::invalid_argument);
    }

    SECTION("a non-positive count") {
        CHECK_THROWS_AS(n2s::insert_knot(curve, 0.6, 0), std::invalid_argument);
    }

    SECTION("at a clamped end, where the multiplicity is already p + 1") {
        CHECK_THROWS_AS(n2s::insert_knot(curve, 0.0), std::invalid_argument);
        CHECK_THROWS_AS(n2s::insert_knot(curve, 1.0), std::invalid_argument);
    }
}
