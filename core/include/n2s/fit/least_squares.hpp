#pragma once

#include "n2s/fit/domain_layout.hpp"
#include "n2s/fit/interpolate.hpp"
#include "n2s/nurbs/surface.hpp"
#include "n2s/subd/control_mesh.hpp"
#include "n2s/subd/subdivision.hpp"
#include "n2s/trim/trim_loop.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace n2s::fit {

/// Which solver produced the answer.
enum class SolverUsed {
    None,
    /// `SimplicialLDLT` on the normal equations, the fast path.
    NormalEquations,
    /// `SparseQR` on the stacked system, used when the normal equations look
    /// too ill-conditioned to trust.
    Qr,
};

std::string to_string(SolverUsed solver);

struct LeastSquaresOptions {
    /// Samples per side of the regular grid laid over the layout's domain
    /// bounding box. The usable sample count is lower, since samples outside
    /// the trim region or outside the layout are discarded. Grows as the
    /// square, and is one of the axes R2 sweeps.
    int samples_per_side = 64;

    /// Weight of the fairness term. Zero is plain least squares.
    ///
    /// Not scale free: the data term is a squared length in model units and
    /// the fairness term is a squared umbrella offset, also in model units, so
    /// lambda is dimensionless but its useful magnitude still depends on how
    /// many samples there are relative to control points. `normalise_lambda`
    /// takes that dependence out.
    double lambda = 0.0;

    /// Divide the data term by the sample count and the fairness term by the
    /// row count before weighting, so that lambda means the same thing across
    /// a sample-density sweep. Without it, raising the density alone weakens
    /// the fairness term, and a sweep over one axis silently moves the other.
    bool normalise_lambda = true;

    /// Above this estimated condition number the normal equations are
    /// abandoned for `SparseQR` on the stacked system, which works with the
    /// square root of that conditioning.
    double max_condition = 1e10;

    /// Iterations used by the condition estimate. It is a lower bound that
    /// tightens with more of them.
    int condition_iterations = 40;
};

struct LeastSquaresReport {
    /// Residual, convergence and notes, shared with the other fit methods.
    /// `max_interpolation_error` here is `max |A V - P|` over the samples.
    FitReport fit;

    std::size_t samples_requested = 0;
    std::size_t samples_used = 0;
    std::size_t samples_outside_trim = 0;
    std::size_t samples_outside_layout = 0;

    std::size_t free_vertices = 0;
    std::size_t fixed_vertices = 0;

    /// Lower-bound estimate of the 2-norm condition number of the matrix that
    /// was factorised. Lower bound, not the value: it comes from power and
    /// inverse-power iteration, which approach the extreme eigenvalues from
    /// the inside. A large estimate is therefore trustworthy; a small one only
    /// means nothing bad was found in the iterations allowed.
    double condition_estimate = 0.0;

    SolverUsed solver = SolverUsed::None;

    /// `||L V||^2` of the fitted mesh.
    double fairness_energy = 0.0;
};

/// Least-squares fit of a layout to a surface, with a fairness term.
///
/// Minimises `||A V - P||^2 + lambda ||L V||^2` over the *interior* control
/// points, where `A` is the limit stencil matrix at the sample locations, `P`
/// the NURBS points at the corresponding domain points, and `L` the umbrella
/// operator of `fairness.hpp`.
///
/// **The boundary, and why it is constrained rather than solved.** The plan
/// specifies constrained boundary control points, and they are constrained to
/// the values R1's interpolating solve gives them: the ones making the limit
/// boundary pass through the layout's boundary vertices, which lie on the trim.
/// Three things follow, and they are the reason this is the right constraint
/// rather than merely a legal one:
///
///  - R2 and R1 then have *identical* boundaries, so comparing them on the
///    same layout isolates what the interior fit does. A boundary that moved
///    too would mix the two effects with no way to separate them afterwards;
///  - R3 needs the seam control points fixed before either patch is fitted,
///    and this is exactly that operation performed on every boundary at once;
///  - it costs nothing to obtain. With EDGE_AND_CORNER the limit position of a
///    boundary vertex depends only on boundary control points, so the R1 square
///    system already decouples, and its boundary block *is* the boundary
///    sub-solve. No approximation is involved in reusing it.
///
/// The constraint is applied by elimination rather than by KKT multipliers: the
/// fixed columns move to the right-hand side, leaving a smaller system that is
/// still symmetric positive definite, which is what lets `SimplicialLDLT` be
/// the fast path at all.
///
/// **What this does not do.** Samples are weighted equally. A uniform grid in
/// the domain is not uniform in `(u, v)` once the layout's quads are not
/// rectangles -- see `LayoutLocator` -- so regions covered by stretched quads
/// carry more samples and pull harder. R2's sample-density sweep is what
/// measures whether that matters; until it has, the equal weighting is a stated
/// choice rather than a justified one.
ControlMesh solve_least_squares(const DomainLayout& layout,
                                const NurbsSurface& surface,
                                const TrimRegion& region,
                                LeastSquaresReport& report,
                                const LeastSquaresOptions& options = {},
                                const SubdivisionOptions& subdivision = {});

} // namespace n2s::fit
