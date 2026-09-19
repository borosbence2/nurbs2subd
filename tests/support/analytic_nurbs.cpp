#include "analytic_nurbs.hpp"

#include "n2s/nurbs/knot_vector.hpp"

#include <cmath>
#include <cstddef>

namespace n2s::testing {

namespace {

const double kHalfSqrt2 = std::sqrt(2.0) / 2.0;

/// Control points and weights of the unit circle in the xy-plane, as four
/// rational quadratic Bezier segments. The points are the corners and edge
/// midpoints of the circumscribed square; the corner weights of sqrt(2)/2 are
/// what makes the representation exact.
struct UnitCircle {
    std::vector<Eigen::Vector2d> points{
        {1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}};
    std::vector<double> weights{
        1.0, kHalfSqrt2, 1.0, kHalfSqrt2, 1.0, kHalfSqrt2, 1.0, kHalfSqrt2, 1.0};
    KnotVector knots{2, {0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1}};
};

double bernstein(int j, int p, double t) {
    double binom = 1.0;
    for (int i = 0; i < j; ++i) {
        binom = binom * static_cast<double>(p - i) / static_cast<double>(i + 1);
    }
    return binom * std::pow(t, j) * std::pow(1.0 - t, p - j);
}

} // namespace

NurbsCurve quarter_circle(double radius) {
    const KnotVector knots{2, {0, 0, 0, 1, 1, 1}};
    const std::vector<Eigen::Vector3d> points{
        {radius, 0.0, 0.0}, {radius, radius, 0.0}, {0.0, radius, 0.0}};
    const std::vector<double> weights{1.0, kHalfSqrt2, 1.0};
    return NurbsCurve{knots, points, weights};
}

NurbsCurve full_circle(double radius) {
    const UnitCircle circle;
    std::vector<Eigen::Vector3d> points;
    points.reserve(circle.points.size());
    for (const Eigen::Vector2d& p : circle.points) {
        points.emplace_back(radius * p.x(), radius * p.y(), 0.0);
    }
    return NurbsCurve{circle.knots, points, circle.weights};
}

NurbsSurface cylinder(double radius, double height) {
    const UnitCircle circle;
    const KnotVector knots_v{1, {0, 0, 1, 1}};

    const std::size_t num_u = circle.points.size();
    std::vector<Eigen::Vector3d> points;
    std::vector<double> weights;
    points.reserve(num_u * 2);
    weights.reserve(num_u * 2);

    // u is the slow index, so each circle control point contributes a bottom
    // and a top entry in sequence.
    for (std::size_t i = 0; i < num_u; ++i) {
        const Eigen::Vector2d& p = circle.points[i];
        points.emplace_back(radius * p.x(), radius * p.y(), 0.0);
        points.emplace_back(radius * p.x(), radius * p.y(), height);
        weights.push_back(circle.weights[i]);
        weights.push_back(circle.weights[i]);
    }

    return NurbsSurface{circle.knots, knots_v, points, weights};
}

NurbsSurface sphere(double radius) {
    const UnitCircle circle;

    // Semicircular profile in the xz-plane, from the north pole down to the
    // south pole: two rational quadratic segments. Stored as (x_radius, z).
    const KnotVector knots_v{2, {0, 0, 0, 0.5, 0.5, 1, 1, 1}};
    const std::vector<Eigen::Vector2d> profile{
        {0.0, radius}, {radius, radius}, {radius, 0.0}, {radius, -radius}, {0.0, -radius}};
    const std::vector<double> profile_weights{1.0, kHalfSqrt2, 1.0, kHalfSqrt2, 1.0};

    // Surface of revolution about the z-axis (Piegl & Tiller A8.1): the
    // profile radius scales the circle control points and the weights
    // multiply.
    std::vector<Eigen::Vector3d> points;
    std::vector<double> weights;
    points.reserve(circle.points.size() * profile.size());
    weights.reserve(circle.points.size() * profile.size());

    for (std::size_t i = 0; i < circle.points.size(); ++i) {
        for (std::size_t j = 0; j < profile.size(); ++j) {
            const double r = profile[j].x();
            points.emplace_back(r * circle.points[i].x(), r * circle.points[i].y(), profile[j].y());
            weights.push_back(circle.weights[i] * profile_weights[j]);
        }
    }

    return NurbsSurface{circle.knots, knots_v, points, weights};
}

NurbsSurface bicubic_bezier(const std::vector<Eigen::Vector3d>& control_net) {
    const KnotVector knots{3, {0, 0, 0, 0, 1, 1, 1, 1}};
    return NurbsSurface::bspline(knots, knots, control_net);
}

std::vector<Eigen::Vector3d> wavy_bicubic_net() {
    std::vector<Eigen::Vector3d> net;
    net.reserve(16);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            const double x = static_cast<double>(i) / 3.0;
            const double y = static_cast<double>(j) / 3.0;
            // Deliberately asymmetric, so that a transposed index would show up.
            const double z = 0.4 * std::sin(2.3 * x) + 0.25 * std::cos(1.7 * y) + 0.15 * x * y;
            net.emplace_back(x, y, z);
        }
    }
    return net;
}

Eigen::Vector3d
bernstein_patch(const std::vector<Eigen::Vector3d>& control_net, double u, double v) {
    Eigen::Vector3d result = Eigen::Vector3d::Zero();
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result += bernstein(i, 3, u) * bernstein(j, 3, v) *
                      control_net[static_cast<std::size_t>(i * 4 + j)];
        }
    }
    return result;
}

} // namespace n2s::testing
