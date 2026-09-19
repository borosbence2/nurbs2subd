#include "n2s/trim/cases.hpp"

#include "n2s/nurbs/knot_vector.hpp"

#include <fmt/format.h>

#include <array>
#include <cmath>
#include <numbers>
#include <string>
#include <utility>

namespace n2s::cases {

namespace {

const double kHalfSqrt2 = std::sqrt(2.0) / 2.0;

/// Bezier coefficients of the identity over [low, high] at degree 2: a linear
/// function is exact at any degree, and the coefficients of the identity are
/// evenly spaced.
std::array<double, 3> linear_coefficients(double low, double high) {
    return {low, 0.5 * (low + high), high};
}

/// Degree-2 Bezier coefficients of `f(t)^2` where `f` maps [0,1] linearly onto
/// [low, high].
///
/// For a quadratic `q(t) = a + b t + c t^2` the degree-2 Bezier coefficients
/// are `q(0)`, `q(0) + b/2`, `q(1)`. Writing the square out gives
/// `a = low^2`, `b = 2 low (high - low)`, `c = (high - low)^2`.
std::array<double, 3> square_coefficients(double low, double high) {
    const double span = high - low;
    const double a = low * low;
    const double b = 2.0 * low * span;
    return {a, a + 0.5 * b, high * high};
}

/// Straight segment as a degree-1 NURBS curve.
NurbsCurve2 segment(const Eigen::Vector2d& from, const Eigen::Vector2d& to) {
    return NurbsCurve2::bspline(KnotVector{1, {0, 0, 1, 1}},
                                std::vector<Eigen::Vector2d>{from, to});
}

/// Circular arc of exactly 90 degrees as a rational quadratic, going
/// counter-clockwise from `start_angle`.
NurbsCurve2 quarter_arc(const Eigen::Vector2d& centre, double radius, double start_angle) {
    const double c = std::cos(start_angle);
    const double s = std::sin(start_angle);

    // Rotate the canonical arc (1,0) -> (1,1) -> (0,1) into place. The middle
    // control point is the intersection of the end tangents, at distance
    // radius*sqrt(2) from the centre, and the weight of sqrt(2)/2 there is what
    // makes the arc exactly circular.
    const auto rotate = [&](double x, double y) {
        return Eigen::Vector2d{centre.x() + radius * (c * x - s * y),
                               centre.y() + radius * (s * x + c * y)};
    };

    return NurbsCurve2{KnotVector{2, {0, 0, 0, 1, 1, 1}},
                       {rotate(1.0, 0.0), rotate(1.0, 1.0), rotate(0.0, 1.0)},
                       {1.0, kHalfSqrt2, 1.0}};
}

/// Biquadratic Bezier patch from a height function that is the sum of a
/// u-dependent and a v-dependent quadratic. Both the paraboloid and the saddle
/// have exactly that shape, which is why both are exact at degree 2.
NurbsSurface separable_quadratic_patch(double half_width, double curvature, double v_sign) {
    const KnotVector knots{2, {0, 0, 0, 1, 1, 1}};

    const std::array<double, 3> position = linear_coefficients(-half_width, half_width);
    const std::array<double, 3> squares = square_coefficients(-half_width, half_width);

    std::vector<Eigen::Vector3d> net;
    net.reserve(9);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const auto ui = static_cast<std::size_t>(i);
            const auto vj = static_cast<std::size_t>(j);
            const double z = curvature * (squares[ui] + v_sign * squares[vj]);
            net.emplace_back(position[ui], position[vj], z);
        }
    }

    return NurbsSurface::bspline(knots, knots, net);
}

/// Control points and weights of the unit circle as four rational quadratic
/// segments, counter-clockwise from (1, 0).
struct UnitCircle {
    std::vector<Eigen::Vector2d> points{
        {1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}};
    std::vector<double> weights{
        1.0, kHalfSqrt2, 1.0, kHalfSqrt2, 1.0, kHalfSqrt2, 1.0, kHalfSqrt2, 1.0};
    KnotVector knots{2, {0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1}};
};

} // namespace

NurbsSurface plane(double width, double height) {
    const KnotVector knots{1, {0, 0, 1, 1}};
    const std::vector<Eigen::Vector3d> net{
        {0.0, 0.0, 0.0}, {0.0, height, 0.0}, {width, 0.0, 0.0}, {width, height, 0.0}};
    return NurbsSurface::bspline(knots, knots, net);
}

NurbsSurface paraboloid(double half_width, double curvature) {
    return separable_quadratic_patch(half_width, curvature, +1.0);
}

NurbsSurface saddle(double half_width, double curvature) {
    return separable_quadratic_patch(half_width, curvature, -1.0);
}

NurbsSurface cylinder(double radius, double height) {
    const UnitCircle circle;
    const KnotVector knots_v{1, {0, 0, 1, 1}};

    std::vector<Eigen::Vector3d> points;
    std::vector<double> weights;
    points.reserve(circle.points.size() * 2);
    weights.reserve(circle.points.size() * 2);

    // u is the slow index, so each circle control point contributes its bottom
    // and top entry in sequence.
    for (std::size_t i = 0; i < circle.points.size(); ++i) {
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

    // Semicircular profile in the xz-plane from the north pole to the south,
    // stored as (radius from the axis, height). Two rational quadratics.
    const KnotVector knots_v{2, {0, 0, 0, 0.5, 0.5, 1, 1, 1}};
    const std::vector<Eigen::Vector2d> profile{
        {0.0, radius}, {radius, radius}, {radius, 0.0}, {radius, -radius}, {0.0, -radius}};
    const std::vector<double> profile_weights{1.0, kHalfSqrt2, 1.0, kHalfSqrt2, 1.0};

    // Revolution about z (Piegl & Tiller A8.1): the profile radius scales the
    // circle control points and the weights multiply.
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

TrimLoop rectangle_loop(const Eigen::Vector2d& low, const Eigen::Vector2d& high) {
    const Eigen::Vector2d a{low.x(), low.y()};
    const Eigen::Vector2d b{high.x(), low.y()};
    const Eigen::Vector2d c{high.x(), high.y()};
    const Eigen::Vector2d d{low.x(), high.y()};

    return TrimLoop{{segment(a, b), segment(b, c), segment(c, d), segment(d, a)}};
}

TrimLoop circle_loop(const Eigen::Vector2d& centre, double radius) {
    const UnitCircle circle;
    const KnotVector arc_knots{2, {0, 0, 0, 1, 1, 1}};

    // Split the nine-point control polygon into its four arcs rather than
    // rotating a canonical arc into place four times. Consecutive arcs then
    // share an endpoint that is the *same expression* and so the same bits,
    // instead of one that agrees only to within cos(pi/2) being 6e-17 rather
    // than 0. The loop closes exactly, and the validator has nothing to report.
    std::vector<NurbsCurve2> arcs;
    arcs.reserve(4);
    for (std::size_t k = 0; k < 4; ++k) {
        std::vector<Eigen::Vector2d> points;
        std::vector<double> weights;
        points.reserve(3);
        weights.reserve(3);
        for (std::size_t i = 2 * k; i <= 2 * k + 2; ++i) {
            points.push_back(centre + radius * circle.points[i]);
            weights.push_back(circle.weights[i]);
        }
        arcs.emplace_back(arc_knots, std::move(points), std::move(weights));
    }
    return TrimLoop{std::move(arcs)};
}

TrimLoop l_shape_loop(const Eigen::Vector2d& low,
                      const Eigen::Vector2d& high,
                      const Eigen::Vector2d& notch) {
    // Counter-clockwise, cutting the block out of the high-u, high-v corner.
    const Eigen::Vector2d a{low.x(), low.y()};
    const Eigen::Vector2d b{high.x(), low.y()};
    const Eigen::Vector2d c{high.x(), notch.y()};
    const Eigen::Vector2d d{notch.x(), notch.y()};
    const Eigen::Vector2d e{notch.x(), high.y()};
    const Eigen::Vector2d f{low.x(), high.y()};

    return TrimLoop{
        {segment(a, b), segment(b, c), segment(c, d), segment(d, e), segment(e, f), segment(f, a)}};
}

TrimLoop
quarter_arc_corner_loop(const Eigen::Vector2d& low, const Eigen::Vector2d& high, double radius) {
    // The arc is centred on the low corner and runs from (low.x + r, low.y) to
    // (low.x, low.y + r). Traversed backwards so that the whole loop stays
    // counter-clockwise while the disc is removed rather than added.
    const Eigen::Vector2d arc_start{low.x() + radius, low.y()};
    const Eigen::Vector2d arc_end{low.x(), low.y() + radius};

    const Eigen::Vector2d b{high.x(), low.y()};
    const Eigen::Vector2d c{high.x(), high.y()};
    const Eigen::Vector2d d{low.x(), high.y()};

    return TrimLoop{{segment(arc_start, b),
                     segment(b, c),
                     segment(c, d),
                     segment(d, arc_end),
                     reversed(quarter_arc(low, radius, 0.0))}};
}

TrimmedCase square_with_circular_hole(const NurbsSurface& surface, double radius) {
    TrimLoop hole = circle_loop({0.5, 0.5}, radius);
    hole.reverse(); // holes run clockwise

    return TrimmedCase{
        .name = fmt::format("square_with_circular_hole_r{:g}", radius),
        .surface = surface,
        .region = TrimRegion{rectangle_loop({0.0, 0.0}, {1.0, 1.0}), {std::move(hole)}},
        .analytic_domain_area = 1.0 - std::numbers::pi * radius * radius,
    };
}

TrimmedCase quarter_arc_corner(const NurbsSurface& surface, double radius) {
    return TrimmedCase{
        .name = fmt::format("quarter_arc_corner_r{:g}", radius),
        .surface = surface,
        .region = TrimRegion{quarter_arc_corner_loop({0.0, 0.0}, {1.0, 1.0}, radius)},
        .analytic_domain_area = 1.0 - 0.25 * std::numbers::pi * radius * radius,
    };
}

TrimmedCase l_shape(const NurbsSurface& surface, const Eigen::Vector2d& notch) {
    const double removed = (1.0 - notch.x()) * (1.0 - notch.y());
    return TrimmedCase{
        .name = fmt::format("l_shape_{:g}x{:g}", notch.x(), notch.y()),
        .surface = surface,
        .region = TrimRegion{l_shape_loop({0.0, 0.0}, {1.0, 1.0}, notch)},
        .analytic_domain_area = 1.0 - removed,
    };
}

std::vector<TrimmedCase> all_synthetic_cases() {
    const std::vector<std::pair<std::string, NurbsSurface>> surfaces{
        {"plane", plane()},
        {"paraboloid", paraboloid()},
        {"saddle", saddle()},
        {"cylinder", cylinder()},
    };

    std::vector<TrimmedCase> cases;
    cases.reserve(surfaces.size() * 3);

    for (const auto& [surface_name, surface] : surfaces) {
        std::vector<TrimmedCase> generated;
        generated.push_back(square_with_circular_hole(surface));
        generated.push_back(quarter_arc_corner(surface));
        generated.push_back(l_shape(surface));

        for (TrimmedCase& one : generated) {
            one.name = fmt::format("{}_{}", surface_name, one.name);
            cases.push_back(std::move(one));
        }
    }

    return cases;
}

} // namespace n2s::cases
