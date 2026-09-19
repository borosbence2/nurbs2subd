#include "n2s/nurbs/projection.hpp"

#include <Eigen/LU> // Matrix2d::determinant, Matrix2d::inverse
#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace n2s {

namespace {

double clamp_to(double value, double low, double high) {
    return std::min(std::max(value, low), high);
}

struct GridSeed {
    double u;
    double v;
    double squared_distance;
};

/// Coarse search for the basin Newton should start in. Newton on the
/// perpendicularity conditions converges to whatever stationary point it is
/// nearest, including maxima and saddles, so the starting iterate decides the
/// answer on anything other than a convex patch.
GridSeed grid_search(const NurbsSurface& surface,
                     const Eigen::Vector3d& target,
                     const ProjectionOptions& options) {
    const double u0 = surface.knots_u().domain_start();
    const double u1 = surface.knots_u().domain_end();
    const double v0 = surface.knots_v().domain_start();
    const double v1 = surface.knots_v().domain_end();

    const auto steps_u =
        static_cast<int>(surface.knots_u().num_spans()) * std::max(1, options.samples_per_span_u);
    const auto steps_v =
        static_cast<int>(surface.knots_v().num_spans()) * std::max(1, options.samples_per_span_v);

    GridSeed best{u0, v0, std::numeric_limits<double>::infinity()};

    for (int i = 0; i <= steps_u; ++i) {
        const double u = u0 + (u1 - u0) * static_cast<double>(i) / static_cast<double>(steps_u);
        for (int j = 0; j <= steps_v; ++j) {
            const double v = v0 + (v1 - v0) * static_cast<double>(j) / static_cast<double>(steps_v);
            const double d2 = (surface.evaluate(u, v) - target).squaredNorm();
            if (d2 < best.squared_distance) {
                best = GridSeed{u, v, d2};
            }
        }
    }

    return best;
}

} // namespace

SurfaceProjection project_to_surface(const NurbsSurface& surface,
                                     const Eigen::Vector3d& target,
                                     const ProjectionOptions& options) {
    if (options.max_iterations < 1) {
        throw std::invalid_argument(
            fmt::format("max_iterations must be at least 1, got {}", options.max_iterations));
    }

    const double u0 = surface.knots_u().domain_start();
    const double u1 = surface.knots_u().domain_end();
    const double v0 = surface.knots_v().domain_start();
    const double v1 = surface.knots_v().domain_end();

    const GridSeed seed = grid_search(surface, target, options);
    double u = seed.u;
    double v = seed.v;

    auto make_result = [&](ProjectionStatus status, int iterations) {
        const Eigen::Vector3d point = surface.evaluate(u, v);
        const bool on_boundary = u <= u0 || u >= u1 || v <= v0 || v >= v1;
        return SurfaceProjection{
            .u = u,
            .v = v,
            .point = point,
            .distance = (target - point).norm(),
            .status = status,
            .iterations = iterations,
            .on_boundary = on_boundary,
        };
    };

    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        const SurfaceDerivatives d = surface.derivatives(u, v, 2);
        const Eigen::Vector3d residual = d.position() - target;
        const double residual_norm = residual.norm();

        // Criterion 1: the target is on the surface. Checked first because the
        // cosine criterion below divides by the residual norm.
        if (residual_norm <= options.distance_tolerance) {
            return make_result(ProjectionStatus::PointCoincident, iteration);
        }

        const double du_norm = d.du().norm();
        const double dv_norm = d.dv().norm();
        if (du_norm < tol::kDegenerateDerivative || dv_norm < tol::kDegenerateDerivative) {
            return make_result(ProjectionStatus::DegenerateJacobian, iteration);
        }

        // Criterion 2: the residual is perpendicular to the tangent plane, so
        // this is a genuine foot of the perpendicular.
        const double cos_u = std::abs(d.du().dot(residual)) / (du_norm * residual_norm);
        const double cos_v = std::abs(d.dv().dot(residual)) / (dv_norm * residual_norm);
        if (cos_u <= options.cosine_tolerance && cos_v <= options.cosine_tolerance) {
            return make_result(ProjectionStatus::Perpendicular, iteration);
        }

        // Newton step on f = r.Su, g = r.Sv. The Jacobian is the Hessian of
        // half the squared distance, so it is symmetric.
        const double f = residual.dot(d.du());
        const double g = residual.dot(d.dv());

        Eigen::Matrix2d jacobian;
        jacobian(0, 0) = d.du().squaredNorm() + residual.dot(d.duu());
        jacobian(0, 1) = d.du().dot(d.dv()) + residual.dot(d.duv());
        jacobian(1, 0) = jacobian(0, 1);
        jacobian(1, 1) = d.dv().squaredNorm() + residual.dot(d.dvv());

        const double determinant = jacobian.determinant();
        if (std::abs(determinant) < tol::kDegenerateDerivative) {
            return make_result(ProjectionStatus::DegenerateJacobian, iteration);
        }

        const Eigen::Vector2d rhs{-f, -g};
        const Eigen::Vector2d step = jacobian.inverse() * rhs;

        const double next_u = clamp_to(u + step.x(), u0, u1);
        const double next_v = clamp_to(v + step.y(), v0, v1);

        // Criterion 3: the iterate has stopped moving in model space. Measured
        // in model space rather than parameter space so that the test means the
        // same thing whatever the parameterisation does.
        const Eigen::Vector3d movement = (next_u - u) * d.du() + (next_v - v) * d.dv();
        const bool stalled = movement.norm() <= options.distance_tolerance;

        const bool clamped = (next_u != u + step.x()) || (next_v != v + step.y());

        u = next_u;
        v = next_v;

        if (stalled) {
            // A stall against the clamp is the constrained minimum on the edge
            // of the domain; a stall in the interior is not an answer, it is a
            // tolerance that this surface cannot meet.
            return make_result(clamped ? ProjectionStatus::ClampedToBoundary
                                       : ProjectionStatus::Stalled,
                               iteration);
        }
    }

    return make_result(ProjectionStatus::IterationLimit, options.max_iterations);
}

} // namespace n2s
