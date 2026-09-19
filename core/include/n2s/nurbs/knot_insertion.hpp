#pragma once

#include "n2s/nurbs/curve.hpp"
#include "n2s/nurbs/surface.hpp"

#include <vector>

namespace n2s {

/// Which parametric direction of a surface an operation applies to.
enum class Direction { U, V };

namespace detail {

/// Everything `insert_knot` needs to know before touching any control point:
/// where the knot goes, how many times, and what the resulting knot vector
/// looks like. Computed once per insertion and reused for every line of a
/// surface control net.
struct InsertionPlan {
    std::size_t span;          ///< k in Piegl & Tiller A5.1.
    int multiplicity;          ///< s, the existing multiplicity of the knot.
    int times;                 ///< r, how many times to insert.
    std::vector<double> knots; ///< The knot vector of the result.
};

InsertionPlan plan_insertion(const KnotVector& kv, double value, int times);

/// Piegl & Tiller A5.1 on homogeneous control points of any dimension. The
/// recurrence is a sequence of affine combinations, which is only
/// geometry-preserving for a rational curve when carried out on the *weighted*
/// points, hence the homogeneous representation.
template<int Dim>
std::vector<Eigen::Matrix<double, Dim + 1, 1>>
insert_into(const KnotVector& kv,
            const std::vector<Eigen::Matrix<double, Dim + 1, 1>>& points,
            double value,
            const InsertionPlan& plan);

} // namespace detail

/// Inserts the knot `u` into the curve `times` times, returning a curve with
/// the same geometry and `times` extra control points.
///
/// Knot insertion is exact: the returned curve evaluates identically to the
/// input at every parameter, to rounding. It is the mechanism behind Bezier
/// extraction, behind subdividing a curve at a parameter, and behind matching
/// the knot structure of two curves that have to share a seam.
///
/// Throws `std::invalid_argument` if `u` is outside the parametric domain or if
/// the insertion would raise the multiplicity of an interior knot above the
/// degree, which would disconnect the curve.
template<int Dim>
NurbsCurveT<Dim> insert_knot(const NurbsCurveT<Dim>& curve, double u, int times = 1);

/// Surface counterpart, inserting into one direction. The insertion is applied
/// to every row (or column) of the control net, which is what makes the result
/// identical to the input surface.
NurbsSurface
insert_knot(const NurbsSurface& surface, Direction direction, double value, int times = 1);

// ---------------------------------------------------------------------------
// Implementation.
// ---------------------------------------------------------------------------

namespace detail {

template<int Dim>
Eigen::Matrix<double, Dim + 1, 1> to_homogeneous(const Eigen::Matrix<double, Dim, 1>& point,
                                                 double weight) {
    Eigen::Matrix<double, Dim + 1, 1> result;
    result.template head<Dim>() = weight * point;
    result[Dim] = weight;
    return result;
}

template<int Dim>
std::vector<Eigen::Matrix<double, Dim + 1, 1>>
insert_into(const KnotVector& kv,
            const std::vector<Eigen::Matrix<double, Dim + 1, 1>>& points,
            double value,
            const InsertionPlan& plan) {
    using Homogeneous = Eigen::Matrix<double, Dim + 1, 1>;

    const int p = kv.degree();
    const int k = static_cast<int>(plan.span);
    const int s = plan.multiplicity;
    const int r = plan.times;
    const auto np = static_cast<int>(points.size()) - 1;

    const auto index = [](int i) { return static_cast<std::size_t>(i); };

    std::vector<Homogeneous> result(points.size() + index(r), Homogeneous::Zero());

    // Control points before and after the affected window are unchanged.
    for (int i = 0; i <= k - p; ++i) {
        result[index(i)] = points[index(i)];
    }
    for (int i = k - s; i <= np; ++i) {
        result[index(i + r)] = points[index(i)];
    }

    // The window itself, refined r times by de Boor's recurrence.
    std::vector<Homogeneous> window(index(p - s) + 1, Homogeneous::Zero());
    for (int i = 0; i <= p - s; ++i) {
        window[index(i)] = points[index(k - p + i)];
    }

    int last = 0;
    for (int j = 1; j <= r; ++j) {
        const int l = k - p + j;
        for (int i = 0; i <= p - j - s; ++i) {
            const double lower = kv[index(l + i)];
            const double upper = kv[index(i + k + 1)];
            const double alpha = (value - lower) / (upper - lower);
            window[index(i)] = alpha * window[index(i + 1)] + (1.0 - alpha) * window[index(i)];
        }
        result[index(l)] = window[0];
        result[index(k + r - j - s)] = window[index(p - j - s)];
        last = l;
    }

    for (int i = last + 1; i < k - s; ++i) {
        result[index(i)] = window[index(i - last)];
    }

    return result;
}

} // namespace detail

template<int Dim>
NurbsCurveT<Dim> insert_knot(const NurbsCurveT<Dim>& curve, double u, int times) {
    using Point = typename NurbsCurveT<Dim>::Point;
    using Homogeneous = Eigen::Matrix<double, Dim + 1, 1>;

    const detail::InsertionPlan plan = detail::plan_insertion(curve.knots(), u, times);

    std::vector<Homogeneous> homogeneous;
    homogeneous.reserve(curve.num_control_points());
    for (std::size_t i = 0; i < curve.num_control_points(); ++i) {
        homogeneous.push_back(
            detail::to_homogeneous<Dim>(curve.control_points()[i], curve.weights()[i]));
    }

    const std::vector<Homogeneous> inserted =
        detail::insert_into<Dim>(curve.knots(), homogeneous, u, plan);

    std::vector<Point> points;
    std::vector<double> weights;
    points.reserve(inserted.size());
    weights.reserve(inserted.size());
    for (const Homogeneous& h : inserted) {
        points.push_back(h.template head<Dim>() / h[Dim]);
        weights.push_back(h[Dim]);
    }

    return NurbsCurveT<Dim>{
        KnotVector{curve.degree(), plan.knots}, std::move(points), std::move(weights)};
}

} // namespace n2s
