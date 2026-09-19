#pragma once

#include "n2s/nurbs/knot_vector.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace n2s {

/// A NURBS curve in 3D: control points, positive weights and a knot vector.
///
/// Control points are stored in Cartesian form with the weights alongside them,
/// which is what files carry and what a user edits. The evaluation routines
/// convert to homogeneous coordinates internally.
///
/// Validity is a class invariant; the constructor throws `std::invalid_argument`
/// on a control point count that disagrees with the knot vector, a mismatched
/// weight count, or a non-positive weight.
class NurbsCurve {
public:
    NurbsCurve(KnotVector knots,
               std::vector<Eigen::Vector3d> control_points,
               std::vector<double> weights);

    /// Non-rational curve: every weight is 1.
    static NurbsCurve bspline(KnotVector knots, std::vector<Eigen::Vector3d> control_points);

    const KnotVector& knots() const noexcept { return knots_; }

    const std::vector<Eigen::Vector3d>& control_points() const noexcept { return points_; }

    const std::vector<double>& weights() const noexcept { return weights_; }

    int degree() const noexcept { return knots_.degree(); }

    std::size_t num_control_points() const noexcept { return points_.size(); }

    /// True if any weight differs from the first one. A curve with uniform
    /// weights is polynomial, and several algorithms can take a shorter path.
    bool is_rational() const noexcept;

    /// Point on the curve. `u` is clamped to the parametric domain.
    Eigen::Vector3d evaluate(double u) const;

    /// Derivatives of orders 0 through `max_order`, so `result[k]` is
    /// `d^k C / du^k`. Piegl & Tiller A3.2 for the homogeneous derivatives and
    /// A4.2 for the rational quotient rule.
    std::vector<Eigen::Vector3d> derivatives(double u, int max_order) const;

private:
    KnotVector knots_;
    std::vector<Eigen::Vector3d> points_;
    std::vector<double> weights_;
};

} // namespace n2s
