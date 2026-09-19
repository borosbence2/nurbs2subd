#include "n2s/nurbs/projection.hpp"
#include "n2s/nurbs/surface.hpp"

#include "support/analytic_nurbs.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

using Catch::Approx;
using n2s::KnotVector;
using n2s::NurbsSurface;
using n2s::ProjectionStatus;
using n2s::SurfaceProjection;

TEST_CASE("projecting a point of the surface recovers its parameters", "[nurbs][projection]") {
    // The round-trip oracle from the plan: project S(u,v) and get (u,v) back.
    const NurbsSurface patch = n2s::testing::bicubic_bezier(n2s::testing::wavy_bicubic_net());

    for (int i = 1; i < 8; ++i) {
        for (int j = 1; j < 8; ++j) {
            const double u = static_cast<double>(i) / 8.0;
            const double v = static_cast<double>(j) / 8.0;

            const SurfaceProjection p = n2s::project_to_surface(patch, patch.evaluate(u, v));

            INFO("u = " << u << ", v = " << v << ", status " << static_cast<int>(p.status)
                        << ", iterations " << p.iterations);
            REQUIRE(p.converged());
            CHECK(p.u == Approx(u).margin(1e-9));
            CHECK(p.v == Approx(v).margin(1e-9));
            CHECK(p.distance < 1e-10);
            CHECK_FALSE(p.on_boundary);
        }
    }
}

TEST_CASE("projecting onto a plane drops the perpendicular", "[nurbs][projection]") {
    const KnotVector knots{1, {0, 0, 1, 1}};
    const std::vector<Eigen::Vector3d> net{{0, 0, 0}, {0, 4, 0}, {6, 0, 0}, {6, 4, 0}};
    const NurbsSurface plane = NurbsSurface::bspline(knots, knots, net);

    const Eigen::Vector3d target{3.0, 1.0, 7.5};
    const SurfaceProjection p = n2s::project_to_surface(plane, target);

    REQUIRE(p.converged());
    CHECK(p.status == ProjectionStatus::Perpendicular);
    CHECK(p.point.x() == Approx(3.0).margin(1e-10));
    CHECK(p.point.y() == Approx(1.0).margin(1e-10));
    CHECK(p.point.z() == Approx(0.0).margin(1e-10));
    CHECK(p.distance == Approx(7.5).margin(1e-10));
}

TEST_CASE("projecting onto a sphere lands on the radial line", "[nurbs][projection]") {
    const double radius = 2.0;
    const NurbsSurface ball = n2s::testing::sphere(radius);

    // A point outside the sphere: the closest surface point is along the ray
    // from the centre, at distance |target| - radius.
    const Eigen::Vector3d target{3.0, 1.0, 2.0};
    const SurfaceProjection p = n2s::project_to_surface(ball, target);

    REQUIRE(p.converged());
    CHECK(p.distance == Approx(target.norm() - radius).epsilon(1e-9));
    CHECK((p.point - radius * target.normalized()).norm() < 1e-9);
}

TEST_CASE("a target beyond the patch edge clamps to the boundary", "[nurbs][projection]") {
    // The closest point of this patch to a target well past the u = 1 edge is
    // on that edge, where the residual cannot be perpendicular to the surface.
    const KnotVector knots{1, {0, 0, 1, 1}};
    const std::vector<Eigen::Vector3d> net{{0, 0, 0}, {0, 1, 0}, {1, 0, 0}, {1, 1, 0}};
    const NurbsSurface patch = NurbsSurface::bspline(knots, knots, net);

    const Eigen::Vector3d target{5.0, 0.5, 0.0};
    const SurfaceProjection p = n2s::project_to_surface(patch, target);

    CHECK(p.on_boundary);
    CHECK(p.u == Approx(1.0).margin(1e-12));
    CHECK(p.v == Approx(0.5).margin(1e-9));
    CHECK(p.distance == Approx(4.0).margin(1e-9));
    CHECK(p.status == ProjectionStatus::ClampedToBoundary);
    CHECK(p.converged());
}

TEST_CASE("projection reports the status rather than guessing", "[nurbs][projection]") {
    const NurbsSurface patch = n2s::testing::bicubic_bezier(n2s::testing::wavy_bicubic_net());

    SECTION("an exhausted iteration budget is reported, not hidden") {
        n2s::ProjectionOptions options;
        options.max_iterations = 1;
        options.samples_per_span_u = 1;
        options.samples_per_span_v = 1;
        options.distance_tolerance = 1e-16;
        options.cosine_tolerance = 1e-16;

        const SurfaceProjection p =
            n2s::project_to_surface(patch, Eigen::Vector3d{0.4, 0.6, 5.0}, options);

        CHECK(p.status == ProjectionStatus::IterationLimit);
        CHECK_FALSE(p.converged());
        CHECK(p.iterations == 1);
    }

    SECTION("a bad iteration limit is rejected outright") {
        n2s::ProjectionOptions options;
        options.max_iterations = 0;
        CHECK_THROWS_AS(n2s::project_to_surface(patch, Eigen::Vector3d::Zero(), options),
                        std::invalid_argument);
    }
}

TEST_CASE("the reported point always matches the reported parameters", "[nurbs][projection]") {
    // Whatever the status, (u, v) and point must agree, so that a caller can
    // use the best iterate without re-evaluating or branching on status.
    const NurbsSurface ball = n2s::testing::sphere(1.0);

    const std::vector<Eigen::Vector3d> targets{
        {0.0, 0.0, 0.0}, {5.0, 0.0, 0.0}, {0.1, -0.2, 0.95}, {0.0, 0.0, 3.0}};

    for (const Eigen::Vector3d& target : targets) {
        const SurfaceProjection p = n2s::project_to_surface(ball, target);

        INFO("target " << target.transpose() << ", status " << static_cast<int>(p.status));
        CHECK((p.point - ball.evaluate(p.u, p.v)).norm() < 1e-14);
        CHECK(p.distance == Approx((target - p.point).norm()).margin(1e-14));
    }
}
