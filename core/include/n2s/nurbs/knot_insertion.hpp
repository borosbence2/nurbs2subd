#pragma once

#include "n2s/nurbs/curve.hpp"
#include "n2s/nurbs/surface.hpp"

namespace n2s {

/// Which parametric direction of a surface an operation applies to.
enum class Direction { U, V };

/// Inserts the knot `u` into the curve `times` times, returning a curve with
/// the same geometry and `times` extra control points.
///
/// Knot insertion is exact: the returned curve evaluates identically to the
/// input at every parameter, to rounding. It is the mechanism behind Bezier
/// extraction, behind subdividing a curve at a parameter, and behind matching
/// the knot structure of two curves that have to share a seam.
///
/// Throws `std::invalid_argument` if `u` is outside the parametric domain or if
/// the insertion would raise the multiplicity of an interior knot above the
/// degree, which would disconnect the curve.
NurbsCurve insert_knot(const NurbsCurve& curve, double u, int times = 1);

/// Surface counterpart, inserting into one direction. The insertion is applied
/// to every row (or column) of the control net, which is what makes the result
/// identical to the input surface.
NurbsSurface
insert_knot(const NurbsSurface& surface, Direction direction, double value, int times = 1);

} // namespace n2s
