#pragma once

#include "n2s/nurbs/basis.hpp"
#include "n2s/nurbs/detail/binomial.hpp"
#include "n2s/nurbs/knot_vector.hpp"
#include "n2s/tolerances.hpp"

#include <Eigen/Core>
#include <fmt/format.h>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace n2s {

/// A NURBS curve of dimension `Dim`: control points, positive weights and a
/// knot vector.
///
/// Two dimensions are in use. `NurbsCurve` (3D) is model-space geometry.
/// `NurbsCurve2` (2D) lives in the parametric `(u, v)` domain of a surface and
/// is what trim curves are made of. The evaluation mathematics is identical, so
/// it lives here once rather than being copied per dimension.
///
/// Control points are stored in Cartesian form with the weights alongside them,
/// which is what files carry and what a user edits. The evaluation routines
/// convert to homogeneous coordinates internally.
///
/// Validity is a class invariant; the constructor throws
/// `std::invalid_argument` on a control point count that disagrees with the
/// knot vector, a mismatched weight count, or a non-positive weight.
template<int Dim>
class NurbsCurveT {
public:
    static_assert(Dim == 2 || Dim == 3, "only 2D domain curves and 3D model curves are used");

    using Point = Eigen::Matrix<double, Dim, 1>;

    NurbsCurveT(KnotVector knots, std::vector<Point> control_points, std::vector<double> weights);

    /// Non-rational curve: every weight is 1.
    static NurbsCurveT bspline(KnotVector knots, std::vector<Point> control_points);

    const KnotVector& knots() const noexcept { return knots_; }

    const std::vector<Point>& control_points() const noexcept { return points_; }

    const std::vector<double>& weights() const noexcept { return weights_; }

    int degree() const noexcept { return knots_.degree(); }

    std::size_t num_control_points() const noexcept { return points_.size(); }

    /// True if any weight differs from the first one. A curve with uniform
    /// weights is polynomial, and several algorithms can take a shorter path.
    bool is_rational() const noexcept;

    /// Point on the curve. `u` is clamped to the parametric domain.
    Point evaluate(double u) const;

    /// Derivatives of orders 0 through `max_order`, so `result[k]` is
    /// `d^k C / du^k`. Piegl & Tiller A3.2 for the homogeneous derivatives and
    /// A4.2 for the rational quotient rule.
    std::vector<Point> derivatives(double u, int max_order) const;

private:
    KnotVector knots_;
    std::vector<Point> points_;
    std::vector<double> weights_;
};

/// Model-space curve.
using NurbsCurve = NurbsCurveT<3>;

/// Curve in the parametric domain of a surface, i.e. a trim curve.
using NurbsCurve2 = NurbsCurveT<2>;

/// The same curve traversed backwards: control points and weights reversed, and
/// the knot vector mirrored within its own domain. The image is identical, the
/// parameterisation runs the other way. Needed to flip a trim loop whose
/// orientation is wrong without disturbing its geometry.
template<int Dim>
NurbsCurveT<Dim> reversed(const NurbsCurveT<Dim>& curve);

// ---------------------------------------------------------------------------
// Implementation.
// ---------------------------------------------------------------------------

template<int Dim>
NurbsCurveT<Dim>::NurbsCurveT(KnotVector knots,
                              std::vector<Point> control_points,
                              std::vector<double> weights)
    : knots_(std::move(knots)),
      points_(std::move(control_points)),
      weights_(std::move(weights)) {
    if (points_.size() != knots_.num_control_points()) {
        throw std::invalid_argument(fmt::format("knot vector implies {} control points, got {}",
                                                knots_.num_control_points(),
                                                points_.size()));
    }
    if (weights_.size() != points_.size()) {
        throw std::invalid_argument(
            fmt::format("got {} control points but {} weights", points_.size(), weights_.size()));
    }
    for (std::size_t i = 0; i < weights_.size(); ++i) {
        if (!(weights_[i] > tol::kMinWeight)) {
            throw std::invalid_argument(fmt::format(
                "weight {} is {}, but weights must be positive; a non-positive weight breaks "
                "both the convex hull property and the rational derivative formulas",
                i,
                weights_[i]));
        }
    }
}

template<int Dim>
NurbsCurveT<Dim> NurbsCurveT<Dim>::bspline(KnotVector knots, std::vector<Point> control_points) {
    std::vector<double> weights(control_points.size(), 1.0);
    return NurbsCurveT{std::move(knots), std::move(control_points), std::move(weights)};
}

template<int Dim>
bool NurbsCurveT<Dim>::is_rational() const noexcept {
    for (const double w : weights_) {
        if (std::abs(w - weights_.front()) > tol::kMinWeight) {
            return true;
        }
    }
    return false;
}

template<int Dim>
typename NurbsCurveT<Dim>::Point NurbsCurveT<Dim>::evaluate(double u) const {
    const std::size_t span = knots_.find_span(u);
    const std::vector<double> n = bspline::basis_functions(knots_, span, u);

    // Accumulate in homogeneous coordinates and divide once, which is both
    // cheaper and better conditioned than a weighted average of the points.
    Point numerator = Point::Zero();
    double denominator = 0.0;

    const std::size_t first = span - static_cast<std::size_t>(degree());
    for (std::size_t j = 0; j < n.size(); ++j) {
        const double weighted = n[j] * weights_[first + j];
        numerator += weighted * points_[first + j];
        denominator += weighted;
    }

    return numerator / denominator;
}

template<int Dim>
std::vector<typename NurbsCurveT<Dim>::Point> NurbsCurveT<Dim>::derivatives(double u,
                                                                            int max_order) const {
    if (max_order < 0) {
        throw std::invalid_argument(
            fmt::format("derivative order must not be negative, got {}", max_order));
    }

    const std::size_t span = knots_.find_span(u);
    const bspline::BasisDerivatives ders = bspline::basis_derivatives(knots_, span, u, max_order);
    const std::size_t first = span - static_cast<std::size_t>(degree());

    // Homogeneous derivatives: A^(k) is the weighted point numerator and w^(k)
    // the weight denominator. Piegl & Tiller A3.2 applied to (w*P, w).
    const auto num_orders = static_cast<std::size_t>(max_order) + 1;
    std::vector<Point> a(num_orders, Point::Zero());
    std::vector<double> w(num_orders, 0.0);

    for (std::size_t k = 0; k < num_orders; ++k) {
        for (std::size_t j = 0; j < ders.num_functions(); ++j) {
            const double weighted = ders(k, j) * weights_[first + j];
            a[k] += weighted * points_[first + j];
            w[k] += weighted;
        }
    }

    // Rational quotient rule, Piegl & Tiller A4.2:
    //   C^(k) = ( A^(k) - sum_{i=1..k} C(k,i) w^(i) C^(k-i) ) / w
    std::vector<Point> result(num_orders, Point::Zero());
    for (std::size_t k = 0; k < num_orders; ++k) {
        Point value = a[k];
        for (std::size_t i = 1; i <= k; ++i) {
            const double c = detail::binomial(static_cast<int>(k), static_cast<int>(i));
            value -= c * w[i] * result[k - i];
        }
        result[k] = value / w[0];
    }

    return result;
}

template<int Dim>
NurbsCurveT<Dim> reversed(const NurbsCurveT<Dim>& curve) {
    const KnotVector& kv = curve.knots();
    const double a = kv[0];
    const double b = kv[kv.size() - 1];

    // U'[i] = a + b - U[m - i]. Mirroring about the midpoint of the knot
    // vector keeps the domain where it was, so the reversed curve is still
    // parameterised over the same interval.
    std::vector<double> knots;
    knots.reserve(kv.size());
    for (std::size_t i = 0; i < kv.size(); ++i) {
        knots.push_back(a + b - kv[kv.size() - 1 - i]);
    }

    std::vector<typename NurbsCurveT<Dim>::Point> points(curve.control_points().rbegin(),
                                                         curve.control_points().rend());
    std::vector<double> weights(curve.weights().rbegin(), curve.weights().rend());

    return NurbsCurveT<Dim>{
        KnotVector{curve.degree(), std::move(knots)}, std::move(points), std::move(weights)};
}

} // namespace n2s
