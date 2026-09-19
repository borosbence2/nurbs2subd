#pragma once

#include <cstddef>
#include <vector>

namespace n2s {

/// A B-spline knot vector together with the degree it belongs to.
///
/// Piegl & Tiller notation is used throughout this file and its implementation:
/// degree `p`, knots `U[0..m]`, and `n + 1 = m - p` control points.
///
/// Validity is a class invariant: an instance that exists is usable. The
/// constructor throws `std::invalid_argument` rather than producing an object
/// that every downstream algorithm would have to re-check.
class KnotVector {
public:
    /// Throws `std::invalid_argument` if `knots` is not a valid knot vector for
    /// `degree`: degree below 1, fewer than `2 * (degree + 1)` knots, a
    /// decreasing step, a degenerate domain, or an interior knot repeated more
    /// than `degree` times.
    KnotVector(int degree, std::vector<double> knots);

    /// Clamped knot vector on [0, 1] with evenly spaced interior knots.
    static KnotVector uniform_clamped(int degree, std::size_t num_control_points);

    int degree() const noexcept { return degree_; }

    /// Number of knots, `m + 1`.
    std::size_t size() const noexcept { return knots_.size(); }

    /// Number of control points the vector implies, `n + 1 = m - p`.
    std::size_t num_control_points() const noexcept;

    /// Number of non-empty knot spans inside the domain.
    std::size_t num_spans() const noexcept;

    const std::vector<double>& knots() const noexcept { return knots_; }

    double operator[](std::size_t i) const { return knots_[i]; }

    /// Start of the parametric domain, `U[p]`.
    double domain_start() const noexcept;

    /// End of the parametric domain, `U[m - p]`.
    double domain_end() const noexcept;

    /// Index `i` of the knot span with `U[i] <= u < U[i+1]`, following
    /// Piegl & Tiller A2.1. `u` is clamped to the domain first, and
    /// `u == domain_end()` resolves to the last non-empty span so that the
    /// result always indexes a valid set of basis functions.
    std::size_t find_span(double u) const;

    /// How many times `u` appears in the knot vector, compared with
    /// `tol::kKnot`. Zero if `u` is not a knot.
    int multiplicity(double u) const noexcept;

    /// True if the first and last `degree + 1` knots are each equal, i.e. the
    /// curve interpolates its first and last control points.
    bool is_clamped() const noexcept;

private:
    int degree_;
    std::vector<double> knots_;
};

} // namespace n2s
