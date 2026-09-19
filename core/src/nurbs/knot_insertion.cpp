#include "n2s/nurbs/knot_insertion.hpp"

#include "n2s/tolerances.hpp"

#include <fmt/format.h>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace n2s {

namespace {

std::size_t as_index(int i) {
    return static_cast<std::size_t>(i);
}

/// Homogeneous control point (w*x, w*y, w*z, w). The insertion recurrence is a
/// sequence of affine combinations, which is only geometry-preserving for a
/// rational curve when it is carried out on the weighted points.
Eigen::Vector4d to_homogeneous(const Eigen::Vector3d& point, double weight) {
    return Eigen::Vector4d{weight * point.x(), weight * point.y(), weight * point.z(), weight};
}

struct InsertionPlan {
    std::size_t span;          ///< k in Piegl & Tiller A5.1.
    int multiplicity;          ///< s, the existing multiplicity of the knot.
    int times;                 ///< r, how many times to insert.
    std::vector<double> knots; ///< The knot vector of the result.
};

InsertionPlan plan_insertion(const KnotVector& kv, double value, int times) {
    if (times < 1) {
        throw std::invalid_argument(
            fmt::format("knot insertion count must be at least 1, got {}", times));
    }
    if (value < kv.domain_start() - tol::kKnot || value > kv.domain_end() + tol::kKnot) {
        throw std::invalid_argument(fmt::format("cannot insert knot {} outside the domain [{}, {}]",
                                                value,
                                                kv.domain_start(),
                                                kv.domain_end()));
    }

    const int multiplicity = kv.multiplicity(value);
    const int degree = kv.degree();

    // Multiplicity above the degree makes the curve discontinuous there, which
    // KnotVector rejects on construction anyway; failing here gives the caller
    // a message that names the real problem.
    if (multiplicity + times > degree) {
        throw std::invalid_argument(fmt::format(
            "inserting knot {} {} times would raise its multiplicity to {}, above the degree {}",
            value,
            times,
            multiplicity + times,
            degree));
    }

    const std::size_t span = kv.find_span(value);

    std::vector<double> knots;
    knots.reserve(kv.size() + as_index(times));
    for (std::size_t i = 0; i <= span; ++i) {
        knots.push_back(kv[i]);
    }
    for (int i = 0; i < times; ++i) {
        knots.push_back(value);
    }
    for (std::size_t i = span + 1; i < kv.size(); ++i) {
        knots.push_back(kv[i]);
    }

    return InsertionPlan{span, multiplicity, times, std::move(knots)};
}

/// Piegl & Tiller A5.1, on homogeneous control points.
std::vector<Eigen::Vector4d> insert_into(const KnotVector& kv,
                                         const std::vector<Eigen::Vector4d>& points,
                                         double value,
                                         const InsertionPlan& plan) {
    const int p = kv.degree();
    const int k = static_cast<int>(plan.span);
    const int s = plan.multiplicity;
    const int r = plan.times;
    const auto np = static_cast<int>(points.size()) - 1;

    std::vector<Eigen::Vector4d> result(points.size() + as_index(r), Eigen::Vector4d::Zero());

    // Control points before and after the affected window are unchanged.
    for (int i = 0; i <= k - p; ++i) {
        result[as_index(i)] = points[as_index(i)];
    }
    for (int i = k - s; i <= np; ++i) {
        result[as_index(i + r)] = points[as_index(i)];
    }

    // The window itself, refined r times by de Boor's recurrence.
    std::vector<Eigen::Vector4d> window(as_index(p - s) + 1, Eigen::Vector4d::Zero());
    for (int i = 0; i <= p - s; ++i) {
        window[as_index(i)] = points[as_index(k - p + i)];
    }

    int last = 0;
    for (int j = 1; j <= r; ++j) {
        const int l = k - p + j;
        for (int i = 0; i <= p - j - s; ++i) {
            const double lower = kv[as_index(l + i)];
            const double upper = kv[as_index(i + k + 1)];
            const double alpha = (value - lower) / (upper - lower);
            window[as_index(i)] =
                alpha * window[as_index(i + 1)] + (1.0 - alpha) * window[as_index(i)];
        }
        result[as_index(l)] = window[0];
        result[as_index(k + r - j - s)] = window[as_index(p - j - s)];
        last = l;
    }

    for (int i = last + 1; i < k - s; ++i) {
        result[as_index(i)] = window[as_index(i - last)];
    }

    return result;
}

std::vector<Eigen::Vector4d> homogenise(const std::vector<Eigen::Vector3d>& points,
                                        const std::vector<double>& weights) {
    std::vector<Eigen::Vector4d> result;
    result.reserve(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        result.push_back(to_homogeneous(points[i], weights[i]));
    }
    return result;
}

void dehomogenise(const std::vector<Eigen::Vector4d>& homogeneous,
                  std::vector<Eigen::Vector3d>& points,
                  std::vector<double>& weights) {
    points.clear();
    weights.clear();
    points.reserve(homogeneous.size());
    weights.reserve(homogeneous.size());
    for (const Eigen::Vector4d& h : homogeneous) {
        points.emplace_back(h.x() / h.w(), h.y() / h.w(), h.z() / h.w());
        weights.push_back(h.w());
    }
}

} // namespace

NurbsCurve insert_knot(const NurbsCurve& curve, double u, int times) {
    const InsertionPlan plan = plan_insertion(curve.knots(), u, times);

    const std::vector<Eigen::Vector4d> inserted =
        insert_into(curve.knots(), homogenise(curve.control_points(), curve.weights()), u, plan);

    std::vector<Eigen::Vector3d> points;
    std::vector<double> weights;
    dehomogenise(inserted, points, weights);

    return NurbsCurve{
        KnotVector{curve.degree(), plan.knots}, std::move(points), std::move(weights)};
}

NurbsSurface
insert_knot(const NurbsSurface& surface, Direction direction, double value, int times) {
    const KnotVector& target = direction == Direction::U ? surface.knots_u() : surface.knots_v();
    const InsertionPlan plan = plan_insertion(target, value, times);

    const std::size_t num_u = surface.num_u();
    const std::size_t num_v = surface.num_v();
    const std::size_t new_count =
        (direction == Direction::U ? num_u : num_v) + as_index(plan.times);

    const std::size_t rows = direction == Direction::U ? num_v : num_u;
    const std::size_t new_num_u = direction == Direction::U ? new_count : num_u;
    const std::size_t new_num_v = direction == Direction::U ? num_v : new_count;

    std::vector<Eigen::Vector3d> points(new_num_u * new_num_v, Eigen::Vector3d::Zero());
    std::vector<double> weights(new_num_u * new_num_v, 1.0);

    // Insert into every line of the control net that runs along `direction`:
    // a tensor product surface is a curve in each direction, so the curve
    // algorithm applied line by line is exactly the surface algorithm.
    for (std::size_t line = 0; line < rows; ++line) {
        std::vector<Eigen::Vector4d> homogeneous;
        homogeneous.reserve(direction == Direction::U ? num_u : num_v);

        if (direction == Direction::U) {
            for (std::size_t i = 0; i < num_u; ++i) {
                homogeneous.push_back(
                    to_homogeneous(surface.control_point(i, line), surface.weight(i, line)));
            }
        } else {
            for (std::size_t j = 0; j < num_v; ++j) {
                homogeneous.push_back(
                    to_homogeneous(surface.control_point(line, j), surface.weight(line, j)));
            }
        }

        const std::vector<Eigen::Vector4d> inserted = insert_into(target, homogeneous, value, plan);

        for (std::size_t n = 0; n < inserted.size(); ++n) {
            const Eigen::Vector4d& h = inserted[n];
            const std::size_t flat =
                direction == Direction::U ? n * new_num_v + line : line * new_num_v + n;
            points[flat] = Eigen::Vector3d{h.x() / h.w(), h.y() / h.w(), h.z() / h.w()};
            weights[flat] = h.w();
        }
    }

    KnotVector knots_u =
        direction == Direction::U ? KnotVector{surface.degree_u(), plan.knots} : surface.knots_u();
    KnotVector knots_v =
        direction == Direction::V ? KnotVector{surface.degree_v(), plan.knots} : surface.knots_v();

    return NurbsSurface{
        std::move(knots_u), std::move(knots_v), std::move(points), std::move(weights)};
}

} // namespace n2s
