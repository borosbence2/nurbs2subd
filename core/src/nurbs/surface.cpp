#include "n2s/nurbs/surface.hpp"

#include "n2s/nurbs/basis.hpp"
#include "n2s/nurbs/detail/binomial.hpp"
#include "n2s/tolerances.hpp"

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

SurfaceDerivatives::SurfaceDerivatives(int max_order)
    : max_order_(max_order),
      values_(as_index((max_order + 1) * (max_order + 1)), Eigen::Vector3d::Zero()) {
    if (max_order < 0) {
        throw std::invalid_argument(
            fmt::format("derivative order must not be negative, got {}", max_order));
    }
}

const Eigen::Vector3d& SurfaceDerivatives::operator()(int du_order, int dv_order) const {
    return values_[as_index(du_order * (max_order_ + 1) + dv_order)];
}

Eigen::Vector3d& SurfaceDerivatives::operator()(int du_order, int dv_order) {
    return values_[as_index(du_order * (max_order_ + 1) + dv_order)];
}

NurbsSurface::NurbsSurface(KnotVector knots_u,
                           KnotVector knots_v,
                           std::vector<Eigen::Vector3d> control_points,
                           std::vector<double> weights)
    : knots_u_(std::move(knots_u)),
      knots_v_(std::move(knots_v)),
      points_(std::move(control_points)),
      weights_(std::move(weights)) {
    const std::size_t expected = knots_u_.num_control_points() * knots_v_.num_control_points();
    if (points_.size() != expected) {
        throw std::invalid_argument(
            fmt::format("knot vectors imply a {} x {} control net ({} points), got {}",
                        knots_u_.num_control_points(),
                        knots_v_.num_control_points(),
                        expected,
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

NurbsSurface NurbsSurface::bspline(KnotVector knots_u,
                                   KnotVector knots_v,
                                   std::vector<Eigen::Vector3d> control_points) {
    std::vector<double> weights(control_points.size(), 1.0);
    return NurbsSurface{
        std::move(knots_u), std::move(knots_v), std::move(control_points), std::move(weights)};
}

bool NurbsSurface::is_rational() const noexcept {
    for (const double w : weights_) {
        if (std::abs(w - weights_.front()) > tol::kMinWeight) {
            return true;
        }
    }
    return false;
}

Eigen::Vector3d NurbsSurface::evaluate(double u, double v) const {
    const std::size_t span_u = knots_u_.find_span(u);
    const std::size_t span_v = knots_v_.find_span(v);
    const std::vector<double> nu = bspline::basis_functions(knots_u_, span_u, u);
    const std::vector<double> nv = bspline::basis_functions(knots_v_, span_v, v);

    const std::size_t first_u = span_u - as_index(degree_u());
    const std::size_t first_v = span_v - as_index(degree_v());

    Eigen::Vector3d numerator = Eigen::Vector3d::Zero();
    double denominator = 0.0;

    for (std::size_t i = 0; i < nu.size(); ++i) {
        for (std::size_t j = 0; j < nv.size(); ++j) {
            const double weighted = nu[i] * nv[j] * weight(first_u + i, first_v + j);
            numerator += weighted * control_point(first_u + i, first_v + j);
            denominator += weighted;
        }
    }

    return numerator / denominator;
}

SurfaceDerivatives NurbsSurface::derivatives(double u, double v, int max_order) const {
    if (max_order < 0) {
        throw std::invalid_argument(
            fmt::format("derivative order must not be negative, got {}", max_order));
    }

    const std::size_t span_u = knots_u_.find_span(u);
    const std::size_t span_v = knots_v_.find_span(v);
    const bspline::BasisDerivatives du = bspline::basis_derivatives(knots_u_, span_u, u, max_order);
    const bspline::BasisDerivatives dv = bspline::basis_derivatives(knots_v_, span_v, v, max_order);

    const std::size_t first_u = span_u - as_index(degree_u());
    const std::size_t first_v = span_v - as_index(degree_v());

    // Homogeneous derivatives A^(k,l) (weighted points) and w^(k,l) (weights).
    SurfaceDerivatives a(max_order);
    std::vector<double> w(as_index((max_order + 1) * (max_order + 1)), 0.0);
    auto weight_deriv = [&](int k, int l) -> double& {
        return w[as_index(k * (max_order + 1) + l)];
    };

    for (int k = 0; k <= max_order; ++k) {
        for (int l = 0; k + l <= max_order; ++l) {
            Eigen::Vector3d point_sum = Eigen::Vector3d::Zero();
            double weight_sum = 0.0;

            for (std::size_t i = 0; i < du.num_functions(); ++i) {
                for (std::size_t j = 0; j < dv.num_functions(); ++j) {
                    const double basis = du(as_index(k), i) * dv(as_index(l), j);
                    const double weighted = basis * weight(first_u + i, first_v + j);
                    point_sum += weighted * control_point(first_u + i, first_v + j);
                    weight_sum += weighted;
                }
            }

            a(k, l) = point_sum;
            weight_deriv(k, l) = weight_sum;
        }
    }

    // Rational quotient rule, Piegl & Tiller A4.4. The three correction sums
    // subtract, in order, the u-derivatives of the weight, the v-derivatives,
    // and the mixed ones.
    SurfaceDerivatives result(max_order);
    const double w00 = weight_deriv(0, 0);

    for (int k = 0; k <= max_order; ++k) {
        for (int l = 0; k + l <= max_order; ++l) {
            Eigen::Vector3d value = a(k, l);

            for (int i = 1; i <= k; ++i) {
                value -= detail::binomial(k, i) * weight_deriv(i, 0) * result(k - i, l);
            }
            for (int j = 1; j <= l; ++j) {
                value -= detail::binomial(l, j) * weight_deriv(0, j) * result(k, l - j);
            }
            for (int i = 1; i <= k; ++i) {
                Eigen::Vector3d inner = Eigen::Vector3d::Zero();
                for (int j = 1; j <= l; ++j) {
                    inner += detail::binomial(l, j) * weight_deriv(i, j) * result(k - i, l - j);
                }
                value -= detail::binomial(k, i) * inner;
            }

            result(k, l) = value / w00;
        }
    }

    return result;
}

} // namespace n2s
