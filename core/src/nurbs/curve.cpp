#include "n2s/nurbs/curve.hpp"

#include "n2s/nurbs/basis.hpp"
#include "n2s/tolerances.hpp"

#include "nurbs/detail/binomial.hpp"

#include <fmt/format.h>

#include <cmath>
#include <stdexcept>
#include <utility>

namespace n2s {

namespace {

std::size_t as_index(int i) {
    return static_cast<std::size_t>(i);
}

} // namespace

NurbsCurve::NurbsCurve(KnotVector knots,
                       std::vector<Eigen::Vector3d> control_points,
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

NurbsCurve NurbsCurve::bspline(KnotVector knots, std::vector<Eigen::Vector3d> control_points) {
    std::vector<double> weights(control_points.size(), 1.0);
    return NurbsCurve{std::move(knots), std::move(control_points), std::move(weights)};
}

bool NurbsCurve::is_rational() const noexcept {
    for (const double w : weights_) {
        if (std::abs(w - weights_.front()) > tol::kMinWeight) {
            return true;
        }
    }
    return false;
}

Eigen::Vector3d NurbsCurve::evaluate(double u) const {
    const std::size_t span = knots_.find_span(u);
    const std::vector<double> n = bspline::basis_functions(knots_, span, u);

    // Accumulate in homogeneous coordinates and divide once, which is both
    // cheaper and better conditioned than a weighted average of the points.
    Eigen::Vector3d numerator = Eigen::Vector3d::Zero();
    double denominator = 0.0;

    const std::size_t first = span - as_index(degree());
    for (std::size_t j = 0; j < n.size(); ++j) {
        const double weighted = n[j] * weights_[first + j];
        numerator += weighted * points_[first + j];
        denominator += weighted;
    }

    return numerator / denominator;
}

std::vector<Eigen::Vector3d> NurbsCurve::derivatives(double u, int max_order) const {
    if (max_order < 0) {
        throw std::invalid_argument(
            fmt::format("derivative order must not be negative, got {}", max_order));
    }

    const std::size_t span = knots_.find_span(u);
    const bspline::BasisDerivatives ders = bspline::basis_derivatives(knots_, span, u, max_order);
    const std::size_t first = span - as_index(degree());

    // Homogeneous derivatives: A^(k) is the weighted point numerator and w^(k)
    // the weight denominator. Piegl & Tiller A3.2 applied to (w*P, w).
    const std::size_t num_orders = as_index(max_order) + 1;
    std::vector<Eigen::Vector3d> a(num_orders, Eigen::Vector3d::Zero());
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
    std::vector<Eigen::Vector3d> result(num_orders, Eigen::Vector3d::Zero());
    for (std::size_t k = 0; k < num_orders; ++k) {
        Eigen::Vector3d value = a[k];
        for (std::size_t i = 1; i <= k; ++i) {
            const double c = detail::binomial(static_cast<int>(k), static_cast<int>(i));
            value -= c * w[i] * result[k - i];
        }
        result[k] = value / w[0];
    }

    return result;
}

} // namespace n2s
