#pragma once

#include "n2s/metrics/statistics.hpp"
#include "n2s/nurbs/projection.hpp"
#include "n2s/nurbs/surface.hpp"
#include "n2s/subd/subdivision.hpp"
#include "n2s/trim/trim_loop.hpp"

#include <Eigen/Core>

#include <functional>
#include <vector>

namespace n2s::metrics {

/// Maps a point of the subdivision control mesh, `(face, u, v)`, to the point
/// of the NURBS parametric domain it is meant to approximate.
///
/// This correspondence is what makes *parametric* error meaningful: without it
/// one can only ask how far the limit surface is from the NURBS as a set, not
/// whether the right part of one sits over the right part of the other. Two
/// surfaces can be geometrically close while badly reparameterised, and only
/// the parametric error sees that.
using DomainMap = std::function<Eigen::Vector2d(const LimitLocation&)>;

/// The correspondence for a control mesh that is an `rows x columns` grid laid
/// over `[0,1]^2` the way `ControlMesh::grid` builds one: face `(i, j)` with
/// local `(u, v)` covers the domain rectangle `[i, i+1] x [j, j+1]` scaled down
/// to the unit square.
///
/// Throws `std::invalid_argument` if the grid is too small to have faces.
DomainMap grid_domain_map(int rows, int columns);

struct SurfaceErrorOptions {
    /// Limit samples per control mesh face edge. The sample count grows as the
    /// square of this.
    int samples_per_face = 8;

    /// Samples per side of the grid used for the NURBS-to-limit direction of
    /// the two-sided distance.
    int nurbs_samples_per_side = 40;

    ProjectionOptions projection;

    /// Skip samples whose domain point falls outside the trimmed region. The
    /// limit surface covers the whole layout, and a layout may overhang the
    /// trim; measuring error there compares against surface that was never
    /// meant to be approximated.
    bool restrict_to_trimmed_region = true;
};

/// Everything measured between a limit surface and the NURBS it approximates.
struct SurfaceError {
    /// `|L(face,u,v) - S(phi(face,u,v))|`, using the domain correspondence.
    /// Empty (count 0) when no correspondence was supplied.
    Stats parametric;

    /// `|L - closest point on the NURBS|`. Independent of any correspondence.
    Stats geometric;

    /// The other direction: NURBS samples to their closest point on the limit
    /// surface. Needed because the one-sided distance above can be small while
    /// the limit surface fails to cover part of the NURBS entirely.
    Stats reverse_geometric;

    /// `max(geometric.max, reverse_geometric.max)`, the two-sided Hausdorff
    /// distance over the samples taken.
    double hausdorff = 0.0;

    /// Angle between the limit normal and the NURBS normal at corresponding
    /// points, in degrees. Needs the domain correspondence.
    Stats normal_degrees;

    /// `|H_limit - H_nurbs|` and `|K_limit - K_nurbs|` at corresponding points.
    Stats mean_curvature;
    Stats gaussian_curvature;
};

/// Measures a limit surface against the NURBS it approximates.
///
/// `domain_map` may be empty, in which case the parametric, normal and
/// curvature statistics come back with `count == 0` rather than being invented
/// from a guessed correspondence.
SurfaceError measure_surface_error(const SubdivisionSurface& limit,
                                   const NurbsSurface& nurbs,
                                   const TrimRegion& region,
                                   const DomainMap& domain_map = {},
                                   const SurfaceErrorOptions& options = {});

/// Per-sample record, for writing out a samples.csv and for colouring a figure.
struct ErrorSample {
    LimitLocation location;
    Eigen::Vector2d domain_point{0.0, 0.0};
    Eigen::Vector3d limit_point{0.0, 0.0, 0.0};
    double parametric = 0.0;
    double geometric = 0.0;
    double normal_degrees = 0.0;
    bool measured = true;

    /// `||H_limit| - |H_nurbs||` and `|K_limit - K_nurbs|`, the curvature
    /// deviations. Magnitudes for the mean because the two surfaces may carry
    /// opposite normal orientation, which flips its sign without any
    /// geometric difference; Gaussian curvature is orientation-free.
    double mean_curvature = 0.0;
    double gaussian_curvature = 0.0;

    /// Whether the two above hold a measurement. Separate from `measured`
    /// because curvature needs second derivatives on both surfaces and can
    /// fail on its own, at a degenerate point where position and normal are
    /// perfectly fine. A sample with `measured == false` never has curvature.
    bool curvature_measured = false;
};

/// The same measurement, keeping every sample. Use when the samples themselves
/// are wanted; `measure_surface_error` is the summary-only path.
///
/// Curvature is filled in only when `domain_map` is supplied, since it is a
/// deviation from the NURBS at the corresponding point and there is no
/// corresponding point without a correspondence.
std::vector<ErrorSample> sample_surface_error(const SubdivisionSurface& limit,
                                              const NurbsSurface& nurbs,
                                              const TrimRegion& region,
                                              const DomainMap& domain_map = {},
                                              const SurfaceErrorOptions& options = {});

/// Deviation of the limit surface's boundary from the trim curve's image on the
/// NURBS.
///
/// This is the metric that says whether the converted patch stops where the
/// trim said it should. It is measured against the *image on the surface* of
/// the trim curve, not against the domain curve, because the domain distance
/// means nothing physical.
struct BoundaryErrorOptions {
    /// Samples along each boundary edge of the layout.
    int samples_per_edge = 32;

    /// Density of the polyline the trim curve is approximated by when finding
    /// the closest point on it. Finer costs time and bounds the error of the
    /// measurement itself.
    int trim_samples_per_curve = 256;
};

Stats measure_boundary_error(const SubdivisionSurface& limit,
                             const NurbsSurface& nurbs,
                             const TrimRegion& region,
                             const BoundaryErrorOptions& options = {});

} // namespace n2s::metrics
