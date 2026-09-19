#pragma once

#include <cstddef>

namespace n2s::detail {

/// Binomial coefficient C(n, k) for the small n used by the rational
/// derivative formulas (Piegl & Tiller A4.2 and A4.4). Computed rather than
/// tabulated so that the derivative order is not silently capped.
constexpr double binomial(int n, int k) {
    if (k < 0 || k > n) {
        return 0.0;
    }
    double result = 1.0;
    for (int i = 0; i < k; ++i) {
        result = result * static_cast<double>(n - i) / static_cast<double>(i + 1);
    }
    return result;
}

} // namespace n2s::detail
