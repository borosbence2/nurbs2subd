#pragma once

#include "n2s/nurbs/knot_vector.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace n2s::bspline {

/// The `max_order + 1` by `degree + 1` table produced by
/// `basis_derivatives`: `operator()(k, j)` is the k-th derivative with respect
/// to u of the basis function `N_{span - degree + j, degree}`.
class BasisDerivatives {
public:
    BasisDerivatives(std::size_t num_orders, std::size_t num_functions)
        : num_functions_(num_functions),
          values_(num_orders * num_functions, 0.0) {}

    double operator()(std::size_t order, std::size_t function) const {
        return values_[order * num_functions_ + function];
    }

    double& operator()(std::size_t order, std::size_t function) {
        return values_[order * num_functions_ + function];
    }

    std::size_t num_orders() const noexcept { return values_.size() / num_functions_; }

    std::size_t num_functions() const noexcept { return num_functions_; }

private:
    std::size_t num_functions_;
    std::vector<double> values_;
};

/// The `degree + 1` basis functions that do not vanish at `u`, in order
/// `N_{span - degree}` .. `N_{span}`. Piegl & Tiller A2.2.
///
/// `span` must come from `kv.find_span(u)`. `out` must have exactly
/// `degree + 1` elements.
void basis_functions(const KnotVector& kv, std::size_t span, double u, std::span<double> out);

/// Allocating convenience overload of the above.
std::vector<double> basis_functions(const KnotVector& kv, std::size_t span, double u);

/// Derivatives of the non-vanishing basis functions, orders 0 through
/// `max_order`. Piegl & Tiller A2.3.
///
/// Orders above the degree are identically zero and are returned as such.
BasisDerivatives basis_derivatives(const KnotVector& kv, std::size_t span, double u, int max_order);

} // namespace n2s::bspline
