#include "n2s/nurbs/basis.hpp"

#include <fmt/format.h>

#include <stdexcept>
#include <vector>

namespace n2s::bspline {

namespace {

std::size_t as_index(int i) {
    return static_cast<std::size_t>(i);
}

/// Row-major scratch matrix with signed indices. Piegl & Tiller A2.3 indexes
/// with quantities that go transiently negative (`rk = r - k`), so the
/// arithmetic is kept in `int` and converted once, here.
class Scratch {
public:
    Scratch(int rows, int cols)
        : cols_(cols),
          values_(as_index(rows * cols), 0.0) {}

    double& operator()(int i, int j) { return values_[as_index(i * cols_ + j)]; }

    double operator()(int i, int j) const { return values_[as_index(i * cols_ + j)]; }

private:
    int cols_;
    std::vector<double> values_;
};

} // namespace

void basis_functions(const KnotVector& kv, std::size_t span, double u, std::span<double> out) {
    const int p = kv.degree();
    if (out.size() != as_index(p) + 1) {
        throw std::invalid_argument(fmt::format(
            "basis_functions needs {} output slots for degree {}, got {}", p + 1, p, out.size()));
    }

    // Piegl & Tiller A2.2. `left` and `right` are the textbook temporaries; the
    // recurrence never divides by zero because span is a non-empty span.
    std::vector<double> left(as_index(p) + 1, 0.0);
    std::vector<double> right(as_index(p) + 1, 0.0);

    out[0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[as_index(j)] = u - kv[span + 1 - as_index(j)];
        right[as_index(j)] = kv[span + as_index(j)] - u;

        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            const double temp = out[as_index(r)] / (right[as_index(r + 1)] + left[as_index(j - r)]);
            out[as_index(r)] = saved + right[as_index(r + 1)] * temp;
            saved = left[as_index(j - r)] * temp;
        }
        out[as_index(j)] = saved;
    }
}

std::vector<double> basis_functions(const KnotVector& kv, std::size_t span, double u) {
    std::vector<double> result(as_index(kv.degree()) + 1, 0.0);
    basis_functions(kv, span, u, result);
    return result;
}

BasisDerivatives
basis_derivatives(const KnotVector& kv, std::size_t span, double u, int max_order) {
    if (max_order < 0) {
        throw std::invalid_argument(
            fmt::format("derivative order must not be negative, got {}", max_order));
    }

    const int p = kv.degree();
    BasisDerivatives ders(as_index(max_order) + 1, as_index(p) + 1);

    // Piegl & Tiller A2.3. `ndu` holds the basis functions in its upper
    // triangle and the knot differences in its lower triangle; `a` alternates
    // between two rows of coefficients for the derivative recurrence.
    Scratch ndu(p + 1, p + 1);
    Scratch a(2, p + 1);
    std::vector<double> left(as_index(p) + 1, 0.0);
    std::vector<double> right(as_index(p) + 1, 0.0);

    ndu(0, 0) = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[as_index(j)] = u - kv[span + 1 - as_index(j)];
        right[as_index(j)] = kv[span + as_index(j)] - u;

        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            ndu(j, r) = right[as_index(r + 1)] + left[as_index(j - r)];
            const double temp = ndu(r, j - 1) / ndu(j, r);

            ndu(r, j) = saved + right[as_index(r + 1)] * temp;
            saved = left[as_index(j - r)] * temp;
        }
        ndu(j, j) = saved;
    }

    for (int j = 0; j <= p; ++j) {
        ders(0, as_index(j)) = ndu(j, p);
    }

    // Orders above the degree are identically zero; the recurrence below would
    // index outside `ndu` if it were allowed to run that far.
    const int order_limit = max_order < p ? max_order : p;

    for (int r = 0; r <= p; ++r) {
        int s1 = 0;
        int s2 = 1;
        a(0, 0) = 1.0;

        for (int k = 1; k <= order_limit; ++k) {
            double d = 0.0;
            const int rk = r - k;
            const int pk = p - k;

            if (r >= k) {
                a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
                d = a(s2, 0) * ndu(rk, pk);
            }

            const int j1 = (rk >= -1) ? 1 : -rk;
            const int j2 = (r - 1 <= pk) ? k - 1 : p - r;

            for (int j = j1; j <= j2; ++j) {
                a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
                d += a(s2, j) * ndu(rk + j, pk);
            }

            if (r <= pk) {
                a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
                d += a(s2, k) * ndu(r, pk);
            }

            ders(as_index(k), as_index(r)) = d;

            const int swap = s1;
            s1 = s2;
            s2 = swap;
        }
    }

    // The recurrence above computes the derivatives up to the factor
    // p! / (p - k)!, applied here.
    int factor = p;
    for (int k = 1; k <= order_limit; ++k) {
        for (int j = 0; j <= p; ++j) {
            ders(as_index(k), as_index(j)) *= static_cast<double>(factor);
        }
        factor *= (p - k);
    }

    return ders;
}

} // namespace n2s::bspline
