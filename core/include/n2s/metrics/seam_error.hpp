#pragma once

#include "n2s/metrics/statistics.hpp"
#include "n2s/subd/subdivision.hpp"

#include <Eigen/Core>

#include <vector>

namespace n2s::metrics {

/// One edge of a quad face, named by which side of the face's local `(u, v)`
/// square it is.
enum class FaceEdge {
    VMin, ///< v = 0, walked with increasing u.
    UMax, ///< u = 1, walked with increasing v.
    VMax, ///< v = 1, walked with increasing u.
    UMin, ///< u = 0, walked with increasing v.
};

/// One side of a seam: a face and the edge of it that lies on the seam.
struct SeamEdge {
    int face = 0;
    FaceEdge edge = FaceEdge::VMin;

    /// Set when this side's natural walk runs opposite to the other side's, so
    /// that the two are compared point for corresponding point rather than end
    /// to end.
    bool reversed = false;
};

/// How two patches meet along a shared edge.
struct SeamError {
    /// `|L_a(t) - L_b(t)|` along the seam. This is *the* number R3 exists to
    /// drive to zero -- not small, zero, to the last bit the representation
    /// allows.
    Stats gap;

    /// Angle between the two limit normals across the seam, in degrees. Shared
    /// boundary control points make the gap vanish but say nothing about
    /// tangency, so this is generally *not* zero, and the plan treats closing
    /// it as a separate, optional improvement.
    Stats normal_degrees;
};

/// Converts an edge parameter into the face-local `(u, v)` it names.
Eigen::Vector2d seam_parameter(const SeamEdge& edge, double t);

/// Measures how two patches meet along a shared edge, sampling `samples + 1`
/// points across it.
///
/// The two patches are sampled at corresponding parameters, which is what makes
/// the gap meaningful: comparing the two edges as point *sets* would report
/// zero for two curves that trace the same path at different speeds, and that
/// is not what watertightness means when the two patches are later refined.
SeamError measure_seam(const SubdivisionSurface& a,
                       const SeamEdge& edge_a,
                       const SubdivisionSurface& b,
                       const SeamEdge& edge_b,
                       int samples = 64);

/// Per-sample seam record, for writing out and for plotting the gap along the
/// seam rather than only its maximum.
struct SeamSample {
    double t = 0.0;
    Eigen::Vector3d point_a{0.0, 0.0, 0.0};
    Eigen::Vector3d point_b{0.0, 0.0, 0.0};
    double gap = 0.0;
    double normal_degrees = 0.0;
};

std::vector<SeamSample> sample_seam(const SubdivisionSurface& a,
                                    const SeamEdge& edge_a,
                                    const SubdivisionSurface& b,
                                    const SeamEdge& edge_b,
                                    int samples = 64);

} // namespace n2s::metrics
