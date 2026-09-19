#pragma once

namespace n2s {

/// Named tolerances. CLAUDE.md forbids magic numbers in geometry code: every
/// comparison against a small quantity must name one of these (or a documented
/// per-case tolerance read from the case file).
///
/// All values are absolute unless the name says otherwise. Model-space
/// tolerances assume inputs normalised to roughly unit size; anything that
/// scales with the model takes a relative tolerance and multiplies by the
/// bounding box diagonal at the call site.
namespace tol {

/// Two knot values closer than this are the same knot. Knot vectors come from
/// files and from knot insertion, so exact equality is not reliable, but knots
/// are never "almost" distinct in practice either.
inline constexpr double kKnot = 1e-12;

/// Below this a weight is treated as degenerate and the NURBS is rejected.
/// Negative or zero weights break the convex hull property and the rational
/// derivative formulas.
inline constexpr double kMinWeight = 1e-12;

/// A parametric interval shorter than this is empty. Used to reject degenerate
/// knot vectors and zero-area domains.
inline constexpr double kMinParametricSpan = 1e-12;

/// Length below which a surface partial derivative counts as degenerate, so the
/// surface normal and the fundamental forms are undefined at that point
/// (a pole of a sphere, a collapsed patch edge).
inline constexpr double kDegenerateDerivative = 1e-10;

/// Closest-point projection: the target counts as reached once the residual
/// `|S(u,v) - P|` falls below this. Model-space, so callers working at a scale
/// far from unity should scale it by the bounding box diagonal.
inline constexpr double kProjectionDistance = 1e-12;

/// Closest-point projection: the residual counts as perpendicular to the
/// surface once `|Su.r| / (|Su| |r|)` and its v counterpart fall below this.
/// This is the criterion that actually terminates the iteration when the target
/// is off the surface, where the residual never reaches zero.
inline constexpr double kProjectionCosine = 1e-12;

} // namespace tol

} // namespace n2s
