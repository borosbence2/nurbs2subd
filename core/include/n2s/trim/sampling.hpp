#pragma once

#include "n2s/nurbs/surface.hpp"
#include "n2s/trim/trim_loop.hpp"

#include <Eigen/Core>

#include <vector>

namespace n2s {

enum class SamplingMode {
    /// Subdivide until both the 3D chord length and the 3D sagitta are within
    /// tolerance. This is the mode everything should use.
    Adaptive,
    /// Evenly spaced parameters per knot span, ignoring how the curve behaves
    /// on the surface. Present only so that experiments can quantify what the
    /// adaptive mode buys, which is a defect of the thesis implementation that
    /// this project sets out to correct.
    Uniform,
};

struct SamplingOptions {
    SamplingMode mode = SamplingMode::Adaptive;

    /// Upper bound on the model-space length of one polyline segment. This is
    /// the arc-length term: it bounds how far apart consecutive samples are on
    /// the *surface*, not in the parametric domain, so a trim curve crossing a
    /// stretched region of the parameterisation gets more samples there.
    double max_segment_length = 0.05;

    /// Upper bound on the model-space distance between the true curve on the
    /// surface and the chord standing in for it. This is the curvature term:
    /// it forces samples into tight turns and allows long chords on flat runs.
    double max_sagitta = 0.005;

    /// Safety net on the recursion depth. A degenerate parameterisation can
    /// otherwise subdivide forever without either bound converging.
    int max_subdivisions = 14;

    /// Uniform mode only: how many samples each knot span receives.
    int samples_per_span = 8;
};

/// Parameters at which to sample one trim curve, ascending, from the start of
/// its domain to the end inclusive.
///
/// Every distinct knot of the curve is always present, in both modes. A knot is
/// where the curve's continuity drops, so a sampling that steps over one can
/// smooth a corner away and produce a trimmed region with the wrong shape --
/// no tolerance on chord length or sagitta will catch that, because both are
/// measured between samples.
std::vector<double> sample_parameters(const NurbsCurve2& curve,
                                      const NurbsSurface& surface,
                                      const SamplingOptions& options = {});

/// The `(u, v)` polyline for one trim curve, including both endpoints.
std::vector<Eigen::Vector2d> sample_curve(const NurbsCurve2& curve,
                                          const NurbsSurface& surface,
                                          const SamplingOptions& options = {});

/// The `(u, v)` polyline for a closed loop. Each curve contributes its start
/// but not its end, so the result is a closed polygon with no duplicated
/// vertex: the last point joins back to the first.
std::vector<Eigen::Vector2d>
sample_loop(const TrimLoop& loop, const NurbsSurface& surface, const SamplingOptions& options = {});

/// Model-space length of the polyline that `sample_loop` produces, mapped onto
/// the surface. Used by the tests and by the sampling-density sweeps.
double sampled_loop_length(const TrimLoop& loop,
                           const NurbsSurface& surface,
                           const SamplingOptions& options = {});

} // namespace n2s
