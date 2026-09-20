#pragma once

#include "n2s/fit/domain_layout.hpp"
#include "n2s/nurbs/surface.hpp"
#include "n2s/subd/control_mesh.hpp"
#include "n2s/subd/subdivision.hpp"

#include <Eigen/Core>

#include <string>
#include <vector>

namespace n2s::fit {

/// Where a fit is asked to put the limit surface: one target point per layout
/// vertex, together with the `(face, u, v)` at which the limit surface is
/// evaluated to reach it.
///
/// The correspondence used is the layout's own, so the target for a vertex is
/// the NURBS point it maps to. For a vertex on the layout boundary that point
/// lies on the trim curve, which is what makes the boundary condition here
/// "the limit boundary interpolates the trim" rather than anything extra: the
/// boundary is not special-cased, it simply has targets that sit on the trim.
struct InterpolationTargets {
    std::vector<LimitLocation> locations;
    std::vector<Eigen::Vector3d> points;
    std::vector<Eigen::Vector2d> domain_points;
};

/// Builds one target per layout vertex.
InterpolationTargets interpolation_targets(const DomainLayout& layout, const NurbsSurface& surface);

/// How a fit turned out. Reported rather than assumed: a solve that did not
/// converge still returns a control mesh, and using it without checking is how
/// a bad fit becomes a plausible-looking error number.
struct FitReport {
    /// Largest `|limit(vertex) - target|` after fitting. For the direct solve
    /// this is the residual of the linear system and should be at solver
    /// precision; for an iterative method it is how far it got.
    double max_interpolation_error = 0.0;
    double rms_interpolation_error = 0.0;

    bool converged = false;
    int iterations = 0;

    /// Free text for anything the caller should know: a singular system, an
    /// iteration limit, a suspiciously large residual.
    std::vector<std::string> notes;
};

/// Solves for the control points that make the limit surface pass through
/// every target.
///
/// This is the square system of Halstead, Kass and DeRose (1993): one equation
/// per layout vertex, one unknown per control point, the matrix being the limit
/// stencil matrix `A` restricted to the vertex locations. Boundary and interior
/// are solved together. They decouple anyway -- with EDGE_AND_CORNER the limit
/// position of a boundary vertex depends only on boundary control points -- so
/// solving them as one system costs nothing and avoids a second code path that
/// could disagree with the first.
///
/// `A` is not symmetric, so this uses `SparseLU` rather than the `SimplicialLDLT`
/// the plan names for the least-squares normal equations of R2.
///
/// What this does *not* claim: the limit boundary curve interpolates the trim
/// **at the layout's boundary vertices**, not everywhere between them. Between
/// two of them it is the cubic B-spline through the solved control points,
/// which approximates the trim. `metrics::measure_boundary_error` measures
/// exactly that residual, and it is the number R3 will care about.
ControlMesh solve_interpolation(const DomainLayout& layout,
                                const NurbsSurface& surface,
                                FitReport& report,
                                const SubdivisionOptions& options = {});

struct PiaOptions {
    int max_iterations = 500;

    /// Stop when the largest control point update falls below this. Relative
    /// to the model: an absolute threshold would mean something different on
    /// every case, which is the defect real data already exposed once in the
    /// trim tolerances.
    double relative_update_tolerance = 1e-12;
};

/// Progressive iterative approximation: start from the targets as control
/// points and repeatedly add the current residual.
///
/// `V <- V + (targets - A V)`, which converges because the limit stencil matrix
/// of a Catmull-Clark surface has its spectrum inside the unit disc about 1.
///
/// The point of having this *and* `solve_interpolation` is defect 5 of the
/// thesis: the two were compared as though they were different methods. They
/// solve the same square system, so PIA converges to the direct solution and
/// any difference between them is PIA not having converged yet.
/// `pia_error_history` exists to demonstrate that rather than assert it.
ControlMesh solve_pia(const DomainLayout& layout,
                      const NurbsSurface& surface,
                      FitReport& report,
                      const PiaOptions& pia = {},
                      const SubdivisionOptions& options = {});

/// One PIA run, keeping the error at every iteration.
struct PiaHistory {
    ControlMesh mesh;
    FitReport report;

    /// Max interpolation error after each iteration, starting with iteration 1.
    std::vector<double> max_error;

    /// Largest control point movement at each iteration.
    std::vector<double> max_update;

    /// Max `|V_pia - V_direct|` after each iteration: the distance to the
    /// direct solution, which is the quantity defect 5 is about.
    std::vector<double> distance_to_direct;
};

/// Runs PIA while recording its convergence, and compares it against the
/// direct solve at every iteration.
PiaHistory pia_error_history(const DomainLayout& layout,
                             const NurbsSurface& surface,
                             const PiaOptions& pia = {},
                             const SubdivisionOptions& options = {});

} // namespace n2s::fit
