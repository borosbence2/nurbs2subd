#pragma once

#include "n2s/nurbs/curve.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace n2s {

/// Orientation of a closed loop in the parametric domain.
enum class Orientation { CounterClockwise, Clockwise };

/// Default polygonisation density for the cheap 2D tests (orientation,
/// self-intersection, point containment). These never need the surface, so they
/// deliberately do not use the 3D-aware adaptive sampler: they only have to be
/// fine enough to get the topology right.
inline constexpr int kLoopTestSamplesPerCurve = 64;

/// A closed loop of trim curves in the `(u, v)` domain of a surface.
///
/// The curves are ordered head to tail: the end of curve `i` is the start of
/// curve `i + 1`, and the end of the last is the start of the first. That the
/// joins actually close to within tolerance is *not* a class invariant --
/// closure gaps are exactly what real CAD input gets wrong, so they are
/// measured and repaired by `validate_and_repair` rather than rejected here.
class TrimLoop {
public:
    /// Throws `std::invalid_argument` if `curves` is empty.
    explicit TrimLoop(std::vector<NurbsCurve2> curves);

    const std::vector<NurbsCurve2>& curves() const noexcept { return curves_; }

    std::size_t size() const noexcept { return curves_.size(); }

    /// Mutable access, for the repair pass. Ordinary code should not need it.
    std::vector<NurbsCurve2>& mutable_curves() noexcept { return curves_; }

    /// Largest distance between the end of a curve and the start of the next,
    /// taken cyclically. Zero for a properly closed loop.
    double worst_closure_gap() const;

    /// Polyline through the loop, `samples_per_curve` points from each curve
    /// plus the closing point. Used by every cheap 2D test.
    std::vector<Eigen::Vector2d> polygonise(int samples_per_curve = kLoopTestSamplesPerCurve) const;

    /// Signed area of the polygonised loop by the shoelace formula. Positive
    /// means counter-clockwise. The value converges to the true signed area as
    /// the sampling is refined, but the *sign* is stable long before the
    /// magnitude is, which is all the orientation check needs.
    double signed_area(int samples_per_curve = kLoopTestSamplesPerCurve) const;

    Orientation orientation(int samples_per_curve = kLoopTestSamplesPerCurve) const;

    /// Reverses the traversal direction: the curve order and each curve's own
    /// parameterisation. The point set of the loop is unchanged.
    void reverse();

private:
    std::vector<NurbsCurve2> curves_;
};

/// A trimmed parametric domain: one outer boundary and any number of holes.
///
/// By convention, enforced by `validate_and_repair`, the outer loop runs
/// counter-clockwise and holes run clockwise. That convention is what makes
/// "inside" well defined for the triangulator and for point containment.
class TrimRegion {
public:
    /// Throws `std::invalid_argument` only for structurally impossible input.
    /// Everything that can be repaired is left to `validate_and_repair`.
    explicit TrimRegion(TrimLoop outer, std::vector<TrimLoop> holes = {});

    const TrimLoop& outer() const noexcept { return outer_; }

    const std::vector<TrimLoop>& holes() const noexcept { return holes_; }

    TrimLoop& mutable_outer() noexcept { return outer_; }

    std::vector<TrimLoop>& mutable_holes() noexcept { return holes_; }

    /// True if `point` lies inside the outer loop and outside every hole,
    /// by the even-odd rule on the polygonised loops.
    bool contains(const Eigen::Vector2d& point,
                  int samples_per_curve = kLoopTestSamplesPerCurve) const;

    /// Area of the trimmed region: outer area less the hole areas, using the
    /// polygonised loops. Approximate by construction, converging as the
    /// sampling is refined.
    double area(int samples_per_curve = kLoopTestSamplesPerCurve) const;

private:
    TrimLoop outer_;
    std::vector<TrimLoop> holes_;
};

/// Even-odd point-in-polygon test on an explicit polyline. Exposed because the
/// triangulator needs it on already-sampled boundaries and should not have to
/// re-polygonise.
bool point_in_polygon(const std::vector<Eigen::Vector2d>& polygon, const Eigen::Vector2d& point);

/// Shoelace signed area of an explicit polyline, positive for counter-clockwise.
double polygon_signed_area(const std::vector<Eigen::Vector2d>& polygon);

} // namespace n2s
