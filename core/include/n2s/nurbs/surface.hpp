#pragma once

#include "n2s/nurbs/knot_vector.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace n2s {

/// Table of mixed partial derivatives of a surface, `operator()(k, l)` being
/// `d^(k+l) S / du^k dv^l`. Entries with `k + l > max_order` are not computed
/// and read back as zero.
class SurfaceDerivatives {
public:
    explicit SurfaceDerivatives(int max_order);

    const Eigen::Vector3d& operator()(int du_order, int dv_order) const;
    Eigen::Vector3d& operator()(int du_order, int dv_order);

    int max_order() const noexcept { return max_order_; }

    const Eigen::Vector3d& position() const { return (*this)(0, 0); }

    const Eigen::Vector3d& du() const { return (*this)(1, 0); }

    const Eigen::Vector3d& dv() const { return (*this)(0, 1); }

    const Eigen::Vector3d& duu() const { return (*this)(2, 0); }

    const Eigen::Vector3d& duv() const { return (*this)(1, 1); }

    const Eigen::Vector3d& dvv() const { return (*this)(0, 2); }

private:
    int max_order_;
    std::vector<Eigen::Vector3d> values_;
};

/// A NURBS surface in 3D.
///
/// The control net is stored row-major with u as the slow index: the control
/// point at `(i, j)` sits at `i * num_v() + j`, where `i` runs over the u
/// direction and `j` over v.
///
/// Validity is a class invariant; the constructor throws
/// `std::invalid_argument` on a net size that disagrees with the knot vectors,
/// a mismatched weight count, or a non-positive weight.
class NurbsSurface {
public:
    NurbsSurface(KnotVector knots_u,
                 KnotVector knots_v,
                 std::vector<Eigen::Vector3d> control_points,
                 std::vector<double> weights);

    /// Non-rational surface: every weight is 1.
    static NurbsSurface
    bspline(KnotVector knots_u, KnotVector knots_v, std::vector<Eigen::Vector3d> control_points);

    const KnotVector& knots_u() const noexcept { return knots_u_; }

    const KnotVector& knots_v() const noexcept { return knots_v_; }

    int degree_u() const noexcept { return knots_u_.degree(); }

    int degree_v() const noexcept { return knots_v_.degree(); }

    std::size_t num_u() const noexcept { return knots_u_.num_control_points(); }

    std::size_t num_v() const noexcept { return knots_v_.num_control_points(); }

    const Eigen::Vector3d& control_point(std::size_t i, std::size_t j) const {
        return points_[i * num_v() + j];
    }

    double weight(std::size_t i, std::size_t j) const { return weights_[i * num_v() + j]; }

    const std::vector<Eigen::Vector3d>& control_points() const noexcept { return points_; }

    const std::vector<double>& weights() const noexcept { return weights_; }

    /// True if any weight differs from the first one.
    bool is_rational() const noexcept;

    /// Point on the surface. `(u, v)` is clamped to the parametric domain.
    Eigen::Vector3d evaluate(double u, double v) const;

    /// All mixed partials with `k + l <= max_order`. Piegl & Tiller A3.6 for
    /// the homogeneous derivatives and A4.4 for the rational quotient rule.
    SurfaceDerivatives derivatives(double u, double v, int max_order) const;

private:
    KnotVector knots_u_;
    KnotVector knots_v_;
    std::vector<Eigen::Vector3d> points_;
    std::vector<double> weights_;
};

} // namespace n2s
