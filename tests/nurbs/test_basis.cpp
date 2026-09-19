#include "n2s/nurbs/basis.hpp"
#include "n2s/nurbs/knot_vector.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

using Catch::Approx;
using n2s::KnotVector;
namespace bspline = n2s::bspline;

namespace {

/// Binomial coefficient, small arguments only. Test oracle.
double binomial(int n, int k) {
    double result = 1.0;
    for (int i = 0; i < k; ++i) {
        result = result * static_cast<double>(n - i) / static_cast<double>(i + 1);
    }
    return result;
}

/// Bernstein polynomial B_{j,p}(u). Test oracle: on a Bezier knot vector the
/// B-spline basis functions reduce to exactly these.
double bernstein(int j, int p, double u) {
    return binomial(p, j) * std::pow(u, j) * std::pow(1.0 - u, p - j);
}

/// Central-difference derivative of the basis function with local index `j`.
/// Used as the oracle for the analytic derivatives.
double finite_difference(const KnotVector& kv, double u, std::size_t j, int order, double h) {
    auto value_at = [&](double t) {
        const std::size_t span = kv.find_span(t);
        const std::vector<double> n = bspline::basis_functions(kv, span, t);
        // The non-vanishing functions are N_{span-p} .. N_{span}; convert the
        // global index j into that local window, or return 0 outside it.
        const auto first = static_cast<std::ptrdiff_t>(span) - kv.degree();
        const auto local = static_cast<std::ptrdiff_t>(j) - first;
        if (local < 0 || local >= static_cast<std::ptrdiff_t>(n.size())) {
            return 0.0;
        }
        return n[static_cast<std::size_t>(local)];
    };

    if (order == 1) {
        return (value_at(u + h) - value_at(u - h)) / (2.0 * h);
    }
    if (order == 2) {
        return (value_at(u + h) - 2.0 * value_at(u) + value_at(u - h)) / (h * h);
    }
    FAIL("unsupported finite difference order");
    return 0.0;
}

} // namespace

TEST_CASE("basis functions form a partition of unity", "[nurbs][basis]") {
    const std::vector<KnotVector> vectors{
        KnotVector{1, {0, 0, 0.3, 0.6, 1, 1}},
        KnotVector{2, {0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}},
        KnotVector{3, {0, 0, 0, 0, 0.25, 0.5, 0.5, 0.75, 1, 1, 1, 1}},
        KnotVector::uniform_clamped(5, 9),
    };

    for (const KnotVector& kv : vectors) {
        const double a = kv.domain_start();
        const double b = kv.domain_end();

        for (int i = 0; i <= 100; ++i) {
            const double u = a + (b - a) * static_cast<double>(i) / 100.0;
            const std::size_t span = kv.find_span(u);
            const std::vector<double> n = bspline::basis_functions(kv, span, u);

            REQUIRE(n.size() == static_cast<std::size_t>(kv.degree()) + 1);

            double sum = 0.0;
            for (const double value : n) {
                CHECK(value >= -1e-15); // basis functions are non-negative
                sum += value;
            }
            CHECK(sum == Approx(1.0).margin(1e-14));
        }
    }
}

TEST_CASE("on a Bezier knot vector the basis reduces to Bernstein polynomials", "[nurbs][basis]") {
    for (int p = 1; p <= 5; ++p) {
        const KnotVector kv = KnotVector::uniform_clamped(p, static_cast<std::size_t>(p) + 1);
        REQUIRE(kv.num_spans() == 1);

        for (int i = 0; i <= 20; ++i) {
            const double u = static_cast<double>(i) / 20.0;
            const std::size_t span = kv.find_span(u);
            const std::vector<double> n = bspline::basis_functions(kv, span, u);

            for (int j = 0; j <= p; ++j) {
                INFO("p = " << p << ", j = " << j << ", u = " << u);
                CHECK(n[static_cast<std::size_t>(j)] == Approx(bernstein(j, p, u)).margin(1e-13));
            }
        }
    }
}

TEST_CASE("the zeroth derivative row equals the basis functions", "[nurbs][basis]") {
    const KnotVector kv{3, {0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1}};

    for (int i = 0; i <= 40; ++i) {
        const double u = static_cast<double>(i) / 40.0;
        const std::size_t span = kv.find_span(u);

        const std::vector<double> n = bspline::basis_functions(kv, span, u);
        const bspline::BasisDerivatives ders = bspline::basis_derivatives(kv, span, u, 2);

        for (std::size_t j = 0; j < n.size(); ++j) {
            CHECK(ders(0, j) == Approx(n[j]).margin(1e-15));
        }
    }
}

TEST_CASE("basis derivatives match central finite differences", "[nurbs][basis]") {
    const std::vector<KnotVector> vectors{
        KnotVector{2, {0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}},
        KnotVector{3, {0, 0, 0, 0, 0.25, 0.5, 0.5, 0.75, 1, 1, 1, 1}},
        KnotVector::uniform_clamped(4, 8),
    };

    for (const KnotVector& kv : vectors) {
        const double a = kv.domain_start();
        const double b = kv.domain_end();
        const double h = 1e-5 * (b - a);

        // Sample strictly inside spans: at a knot of multiplicity p the second
        // derivative is genuinely discontinuous and a central difference there
        // compares two different polynomials.
        for (int i = 1; i < 60; ++i) {
            const double u = a + (b - a) * (static_cast<double>(i) + 0.37) / 60.0;
            if (kv.multiplicity(u) > 0) {
                continue;
            }

            const std::size_t span = kv.find_span(u);
            const bspline::BasisDerivatives ders = bspline::basis_derivatives(kv, span, u, 2);
            const auto first =
                static_cast<std::size_t>(static_cast<std::ptrdiff_t>(span) - kv.degree());

            for (std::size_t j = 0; j <= static_cast<std::size_t>(kv.degree()); ++j) {
                INFO("degree " << kv.degree() << ", u = " << u << ", j = " << j);
                CHECK(ders(1, j) == Approx(finite_difference(kv, u, first + j, 1, h)).margin(1e-6));
                CHECK(ders(2, j) == Approx(finite_difference(kv, u, first + j, 2, h)).margin(1e-4));
            }
        }
    }
}

TEST_CASE("derivatives of the basis sum to zero", "[nurbs][basis]") {
    // Differentiating the partition of unity: sum_j N_j(u) == 1 for all u, so
    // every derivative order above zero must sum to exactly zero.
    const KnotVector kv{3, {0, 0, 0, 0, 0.3, 0.55, 0.8, 1, 1, 1, 1}};

    for (int i = 0; i <= 50; ++i) {
        const double u = static_cast<double>(i) / 50.0;
        const std::size_t span = kv.find_span(u);
        const bspline::BasisDerivatives ders = bspline::basis_derivatives(kv, span, u, 3);

        for (std::size_t k = 1; k <= 3; ++k) {
            double sum = 0.0;
            for (std::size_t j = 0; j <= static_cast<std::size_t>(kv.degree()); ++j) {
                sum += ders(k, j);
            }
            INFO("order " << k << " at u = " << u);
            CHECK(sum == Approx(0.0).margin(1e-9));
        }
    }
}

TEST_CASE("derivative orders above the degree vanish", "[nurbs][basis]") {
    const KnotVector kv = KnotVector::uniform_clamped(2, 5);
    const double u = 0.42;
    const std::size_t span = kv.find_span(u);

    const bspline::BasisDerivatives ders = bspline::basis_derivatives(kv, span, u, 4);

    for (std::size_t k = 3; k <= 4; ++k) {
        for (std::size_t j = 0; j <= 2; ++j) {
            CHECK(ders(k, j) == Approx(0.0).margin(1e-15));
        }
    }
}
