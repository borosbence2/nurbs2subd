#pragma once

#include "n2s/nurbs/surface.hpp"
#include "n2s/tolerances.hpp"

#include <Eigen/Core>

namespace n2s {

/// Why a projection stopped. Reported rather than swallowed: an unconverged
/// projection silently treated as converged is exactly how an error metric ends
/// up quietly optimistic.
enum class ProjectionStatus {
    /// The residual fell below the distance tolerance: the target lies on the
    /// surface to within tolerance.
    PointCoincident,
    /// The residual became perpendicular to both partials: a genuine foot of
    /// the perpendicular.
    Perpendicular,
    /// The closest point lies on the edge of the parametric domain and the
    /// clamp holds it there. This is a real answer -- the constrained minimum --
    /// even though the residual is not perpendicular to the surface.
    ClampedToBoundary,
    /// The iterate stopped moving in the interior of the domain without either
    /// tolerance being met. Usually means the tolerances are tighter than the
    /// conditioning of this surface supports.
    Stalled,
    /// The iteration limit was hit while still moving.
    IterationLimit,
    /// The Newton system became singular, which means the surface is degenerate
    /// at the iterate (a pole, a collapsed edge).
    DegenerateJacobian,
};

/// Result of projecting a point onto a surface. `u`, `v` and `point` always
/// hold the best iterate found, whatever the status, so a caller that only
/// wants an approximate footpoint can use them without branching.
struct SurfaceProjection {
    double u;
    double v;
    Eigen::Vector3d point; ///< S(u, v) at the reported parameters.
    double distance;       ///< |target - point|.
    ProjectionStatus status;
    int iterations;
    bool on_boundary; ///< The parameters sit on the edge of the domain.

    /// True when the reported point is a genuine local minimum of the distance,
    /// including one constrained to the edge of the domain. Code that needs
    /// strict perpendicularity -- a normal-deviation metric, say -- should test
    /// `status` instead, since a boundary minimum has no perpendicularity to
    /// measure.
    bool converged() const {
        return status == ProjectionStatus::PointCoincident ||
               status == ProjectionStatus::Perpendicular ||
               status == ProjectionStatus::ClampedToBoundary;
    }
};

struct ProjectionOptions {
    /// Initial search grid density, per non-empty knot span in each direction.
    /// Scaling with the spans rather than using a fixed grid keeps the
    /// initialisation honest on surfaces with many knots, where a fixed grid
    /// would miss features and seed Newton in the wrong basin.
    int samples_per_span_u = 8;
    int samples_per_span_v = 8;

    int max_iterations = 64;
    double distance_tolerance = tol::kProjectionDistance;
    double cosine_tolerance = tol::kProjectionCosine;
};

/// Closest point on a surface to `target`: a grid search for the starting
/// iterate followed by Newton iteration on the two perpendicularity conditions
/// `r.Su = 0`, `r.Sv = 0` with `r = S(u,v) - target`, clamping the parameters
/// to the domain at every step (Piegl & Tiller section 6.1).
///
/// The grid search makes this robust rather than fast, and it finds *a* local
/// minimum, not provably the global one: a fine enough grid is what makes the
/// difference on wrinkled surfaces.
SurfaceProjection project_to_surface(const NurbsSurface& surface,
                                     const Eigen::Vector3d& target,
                                     const ProjectionOptions& options = {});

} // namespace n2s
