#include "n2s/nurbs/knot_vector.hpp"

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

KnotVector::KnotVector(int degree, std::vector<double> knots)
    : degree_(degree),
      knots_(std::move(knots)) {
    if (degree_ < 1) {
        throw std::invalid_argument(
            fmt::format("knot vector degree must be at least 1, got {}", degree_));
    }

    const std::size_t minimum = 2 * as_index(degree_ + 1);
    if (knots_.size() < minimum) {
        throw std::invalid_argument(
            fmt::format("a degree {} knot vector needs at least {} knots, got {}",
                        degree_,
                        minimum,
                        knots_.size()));
    }

    for (std::size_t i = 1; i < knots_.size(); ++i) {
        if (knots_[i] < knots_[i - 1]) {
            throw std::invalid_argument(
                fmt::format("knot vector must be non-decreasing: U[{}] = {} < U[{}] = {}",
                            i,
                            knots_[i],
                            i - 1,
                            knots_[i - 1]));
        }
    }

    if (domain_end() - domain_start() < tol::kMinParametricSpan) {
        throw std::invalid_argument(fmt::format(
            "knot vector has a degenerate domain [{}, {}]", domain_start(), domain_end()));
    }

    // Interior knots only. Multiplicity p + 1 there would split the curve into
    // disconnected pieces, which nothing downstream of here expects.
    const std::size_t first_interior = as_index(degree_) + 1;
    const std::size_t last_interior = knots_.size() - as_index(degree_) - 1;
    for (std::size_t i = first_interior; i < last_interior;) {
        std::size_t j = i;
        while (j < last_interior && std::abs(knots_[j] - knots_[i]) <= tol::kKnot) {
            ++j;
        }
        const auto mult = static_cast<int>(j - i);
        if (mult > degree_) {
            throw std::invalid_argument(
                fmt::format("interior knot {} has multiplicity {}, which exceeds the degree {}",
                            knots_[i],
                            mult,
                            degree_));
        }
        i = j;
    }
}

KnotVector KnotVector::uniform_clamped(int degree, std::size_t num_control_points) {
    if (degree < 1) {
        throw std::invalid_argument(
            fmt::format("knot vector degree must be at least 1, got {}", degree));
    }
    if (num_control_points < as_index(degree + 1)) {
        throw std::invalid_argument(
            fmt::format("a degree {} curve needs at least {} control points, got {}",
                        degree,
                        degree + 1,
                        num_control_points));
    }

    // m + 1 = (n + 1) + (p + 1) knots, of which the first and last p + 1 are
    // repeated and the rest are evenly spaced strictly inside (0, 1).
    const std::size_t num_interior = num_control_points - as_index(degree + 1);
    std::vector<double> knots;
    knots.reserve(num_control_points + as_index(degree + 1));

    knots.insert(knots.end(), as_index(degree + 1), 0.0);
    for (std::size_t i = 1; i <= num_interior; ++i) {
        knots.push_back(static_cast<double>(i) / static_cast<double>(num_interior + 1));
    }
    knots.insert(knots.end(), as_index(degree + 1), 1.0);

    return KnotVector{degree, std::move(knots)};
}

std::size_t KnotVector::num_control_points() const noexcept {
    return knots_.size() - as_index(degree_) - 1;
}

std::size_t KnotVector::num_spans() const noexcept {
    std::size_t spans = 0;
    const std::size_t last = knots_.size() - as_index(degree_) - 1;
    for (std::size_t i = as_index(degree_); i < last; ++i) {
        if (knots_[i + 1] - knots_[i] > tol::kKnot) {
            ++spans;
        }
    }
    return spans;
}

double KnotVector::domain_start() const noexcept {
    return knots_[as_index(degree_)];
}

double KnotVector::domain_end() const noexcept {
    return knots_[knots_.size() - as_index(degree_) - 1];
}

std::size_t KnotVector::find_span(double u) const {
    const std::size_t n = num_control_points() - 1;

    // Clamping rather than asserting: callers routinely arrive with parameters
    // that sit outside the domain by a rounding error, and both sampling and
    // closest-point projection walk right up to the boundary.
    if (u <= domain_start()) {
        // Step over any empty spans at the start of an unclamped vector.
        std::size_t span = as_index(degree_);
        while (span < n && knots_[span + 1] - knots_[span] <= tol::kKnot) {
            ++span;
        }
        return span;
    }
    if (u >= domain_end()) {
        // The textbook special case: u == U[m-p] belongs to the last non-empty
        // span, not to the empty one that follows it.
        std::size_t span = n;
        while (span > as_index(degree_) && knots_[span + 1] - knots_[span] <= tol::kKnot) {
            --span;
        }
        return span;
    }

    // Binary search, Piegl & Tiller A2.1.
    std::size_t low = as_index(degree_);
    std::size_t high = n + 1;
    std::size_t mid = (low + high) / 2;
    while (u < knots_[mid] || u >= knots_[mid + 1]) {
        if (u < knots_[mid]) {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    return mid;
}

int KnotVector::multiplicity(double u) const noexcept {
    int count = 0;
    for (const double knot : knots_) {
        if (std::abs(knot - u) <= tol::kKnot) {
            ++count;
        }
    }
    return count;
}

bool KnotVector::is_clamped() const noexcept {
    const std::size_t p = as_index(degree_);
    for (std::size_t i = 1; i <= p; ++i) {
        if (std::abs(knots_[i] - knots_[0]) > tol::kKnot) {
            return false;
        }
        if (std::abs(knots_[knots_.size() - 1 - i] - knots_.back()) > tol::kKnot) {
            return false;
        }
    }
    return true;
}

} // namespace n2s
