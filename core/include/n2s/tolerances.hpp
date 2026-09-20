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

/// Two consecutive trim curves in a loop count as joined when the gap between
/// one's end and the next one's start falls below this. Parametric-domain
/// units, so it is relative to the `(u, v)` box, not to model size.
inline constexpr double kTrimClosure = 1e-9;

/// A closure gap at or below this is what exact arithmetic would have produced
/// had doubles been exact -- evaluating two rational arcs at a shared endpoint
/// lands a few ulps apart. Such a gap is still snapped shut, but silently: a
/// repair log full of 1e-17 entries trains the reader to ignore it, which is
/// how the 1e-3 entry gets missed.
inline constexpr double kNegligibleClosureGap = 1e-15;

/// Below this a polygon's signed area is too small to give a reliable
/// orientation, which means the loop is degenerate rather than merely small.
inline constexpr double kMinLoopArea = 1e-14;

/// Closest-point projection: the target counts as reached once the residual
/// `|S(u,v) - P|` falls below this. Model-space, so callers working at a scale
/// far from unity should scale it by the bounding box diagonal.
inline constexpr double kProjectionDistance = 1e-12;

/// Closest-point projection: the residual counts as perpendicular to the
/// surface once `|Su.r| / (|Su| |r|)` and its v counterpart fall below this.
/// This is the criterion that actually terminates the iteration when the target
/// is off the surface, where the residual never reaches zero.
inline constexpr double kProjectionCosine = 1e-12;

/// How far outside a layout quad's local `[0,1]^2` an inverted domain point may
/// land and still count as inside that quad. Parametric-domain units, like
/// `kTrimClosure`. A sample sitting exactly on an edge shared by two quads
/// inverts to `u` or `v` of 0 or 1 in both, and rounding decides which of the
/// two it misses by an ulp; without slack it would belong to neither.
inline constexpr double kLayoutContainment = 1e-9;

/// The bilinear inversion treats a quad as a parallelogram, and solves a linear
/// system instead of a quadratic, when the quadratic coefficient falls below
/// this fraction of the quad's own size. Relative rather than absolute, because
/// an absolute area threshold means something different at every model scale --
/// the defect the trim tolerances already exposed once on real data.
inline constexpr double kRelativeBilinearDegeneracy = 1e-12;

} // namespace tol

} // namespace n2s
