#pragma once

#include "n2s/nurbs/curve.hpp"
#include "n2s/nurbs/surface.hpp"

#include <Eigen/Core>

#include <vector>

/// Exact NURBS representations of shapes whose geometry is known in closed
/// form. These are test oracles: the conics below are *exact* under the
/// rational representation, so a correct evaluator reproduces them to machine
/// precision rather than approximately.
namespace n2s::testing {

/// Quarter circle of the given radius in the xy-plane, from (r, 0, 0) to
/// (0, r, 0). Rational quadratic, one Bezier segment.
NurbsCurve quarter_circle(double radius);

/// Full circle of the given radius in the xy-plane, as four rational quadratic
/// segments (the standard nine control point representation).
NurbsCurve full_circle(double radius);

/// Circular cylinder of the given radius and height, axis along z, base at
/// z = 0. Rational quadratic in u (around), linear in v (along the axis).
NurbsSurface cylinder(double radius, double height);

/// Full sphere of the given radius, centred at the origin, built as a surface
/// of revolution of a semicircular profile about the z-axis.
///
/// Note that v = 0 and v = 1 are the poles, where the control net degenerates
/// and the surface normal is undefined. Tests must stay away from them.
NurbsSurface sphere(double radius);

/// Bicubic Bezier patch from a 4x4 control net, row-major with u as the slow
/// index. Non-rational.
NurbsSurface bicubic_bezier(const std::vector<Eigen::Vector3d>& control_net);

/// A 4x4 control net with no symmetry, used wherever an arbitrary but
/// reproducible bicubic patch is needed.
std::vector<Eigen::Vector3d> wavy_bicubic_net();

/// Closed-form evaluation of a bicubic Bezier patch from its control net.
/// Test oracle for `NurbsSurface::evaluate` on a Bezier knot vector.
Eigen::Vector3d
bernstein_patch(const std::vector<Eigen::Vector3d>& control_net, double u, double v);

} // namespace n2s::testing
