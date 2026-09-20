#include "n2s/fit/interpolate.hpp"

#include "n2s/fit/layout.hpp"

#include <Eigen/SparseLU>
#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace n2s::fit {

namespace {

/// Bounding box diagonal of the targets, the scale tolerances are measured
/// against.
double model_scale(const std::vector<Eigen::Vector3d>& points) {
    if (points.empty()) {
        return 1.0;
    }
    Eigen::Vector3d low = points.front();
    Eigen::Vector3d high = points.front();
    for (const Eigen::Vector3d& p : points) {
        low = low.cwiseMin(p);
        high = high.cwiseMax(p);
    }
    return std::max((high - low).norm(), 1e-12);
}

Eigen::MatrixXd as_rows(const std::vector<Eigen::Vector3d>& points) {
    Eigen::MatrixXd rows(static_cast<Eigen::Index>(points.size()), 3);
    for (std::size_t i = 0; i < points.size(); ++i) {
        rows.row(static_cast<Eigen::Index>(i)) = points[i];
    }
    return rows;
}

/// max |A V - targets| over the rows.
double max_row_norm(const Eigen::MatrixXd& residual) {
    double worst = 0.0;
    for (Eigen::Index i = 0; i < residual.rows(); ++i) {
        worst = std::max(worst, residual.row(i).norm());
    }
    return worst;
}

void fill_report(FitReport& report, const Eigen::MatrixXd& residual) {
    report.max_interpolation_error = max_row_norm(residual);
    report.rms_interpolation_error = residual.rows() == 0
                                         ? 0.0
                                         : std::sqrt(residual.rowwise().squaredNorm().sum() /
                                                     static_cast<double>(residual.rows()));
}

} // namespace

InterpolationTargets interpolation_targets(const DomainLayout& layout,
                                           const NurbsSurface& surface) {
    if (layout.quads.empty()) {
        throw std::invalid_argument("cannot build targets for an empty layout");
    }

    InterpolationTargets targets;
    targets.locations.assign(layout.num_vertices(), LimitLocation{-1, 0.0, 0.0});

    // Each vertex is addressed through one of the faces that uses it, at that
    // face's corner. Which face is immaterial: the limit point of a vertex is a
    // property of the vertex.
    static constexpr double kCorners[4][2] = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    for (std::size_t f = 0; f < layout.quads.size(); ++f) {
        for (std::size_t corner = 0; corner < 4; ++corner) {
            const auto vertex = static_cast<std::size_t>(layout.quads[f][corner]);
            if (targets.locations[vertex].face < 0) {
                targets.locations[vertex] =
                    LimitLocation{static_cast<int>(f), kCorners[corner][0], kCorners[corner][1]};
            }
        }
    }

    targets.domain_points = layout.vertices;
    targets.points.reserve(layout.num_vertices());
    for (std::size_t i = 0; i < layout.num_vertices(); ++i) {
        if (targets.locations[i].face < 0) {
            throw std::invalid_argument(fmt::format("layout vertex {} belongs to no quad", i));
        }
        const Eigen::Vector2d& p = layout.vertices[i];
        targets.points.push_back(surface.evaluate(p.x(), p.y()));
    }

    return targets;
}

ControlMesh solve_interpolation(const DomainLayout& layout,
                                const NurbsSurface& surface,
                                FitReport& report,
                                const SubdivisionOptions& options) {
    const InterpolationTargets targets = interpolation_targets(layout, surface);

    // The topology is what matters here; the initial positions only have to be
    // valid, because the solve replaces them.
    const ControlMesh initial = lift(layout, surface);
    const SubdivisionSurface subdivision{initial, options};

    const LimitMatrices matrices =
        subdivision.build_limit_matrices(targets.locations, DerivativeOrder::None);

    Eigen::SparseMatrix<double> a = matrices.position;
    a.makeCompressed();

    if (a.rows() != a.cols()) {
        throw std::runtime_error(fmt::format(
            "the interpolation system is {} x {}, not square; there must be exactly one "
            "target per control point",
            a.rows(),
            a.cols()));
    }

    // SparseLU, not SimplicialLDLT: the limit stencil matrix is not symmetric.
    // The plan's LDLT belongs to R2's normal equations, which are.
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(a);
    solver.factorize(a);

    ControlMesh solved = initial;
    if (solver.info() != Eigen::Success) {
        report.converged = false;
        report.notes.push_back(fmt::format(
            "the interpolation system could not be factorised ({}). The control mesh returned "
            "is the unfitted starting point, not a fit.",
            solver.lastErrorMessage()));
        return solved;
    }

    const Eigen::MatrixXd solution = solver.solve(as_rows(targets.points));
    if (solver.info() != Eigen::Success) {
        report.converged = false;
        report.notes.emplace_back("the interpolation system factorised but the solve failed");
        return solved;
    }

    solved.set_vertices(to_points(solution));

    const Eigen::MatrixXd residual = a * solution - as_rows(targets.points);
    fill_report(report, residual);
    report.iterations = 1;

    // A direct solve either works to near machine precision or has hit a
    // conditioning problem worth knowing about; there is no useful middle.
    const double scale = model_scale(targets.points);
    report.converged = report.max_interpolation_error <= 1e-9 * scale;
    if (!report.converged) {
        report.notes.push_back(fmt::format(
            "the direct solve left a residual of {:.3e} on a model of size {:.3e}, which is "
            "far larger than a well-conditioned factorisation should. Treat the fit as "
            "suspect.",
            report.max_interpolation_error,
            scale));
    }

    return solved;
}

namespace {

/// Shared by `solve_pia` and `pia_error_history`.
struct PiaState {
    Eigen::SparseMatrix<double> a;
    Eigen::MatrixXd targets;
    Eigen::MatrixXd vertices;
    double scale = 1.0;
};

PiaState make_pia_state(const SubdivisionOptions& options,
                        const ControlMesh& initial,
                        const InterpolationTargets& targets) {
    const SubdivisionSurface subdivision{initial, options};
    const LimitMatrices matrices =
        subdivision.build_limit_matrices(targets.locations, DerivativeOrder::None);

    PiaState state;
    state.a = matrices.position;
    state.a.makeCompressed();
    state.targets = as_rows(targets.points);
    // PIA's usual start: take the targets themselves as the control points.
    state.vertices = state.targets;
    state.scale = model_scale(targets.points);
    return state;
}

} // namespace

ControlMesh solve_pia(const DomainLayout& layout,
                      const NurbsSurface& surface,
                      FitReport& report,
                      const PiaOptions& pia,
                      const SubdivisionOptions& options) {
    if (pia.max_iterations < 1) {
        throw std::invalid_argument(
            fmt::format("PIA needs at least one iteration, got {}", pia.max_iterations));
    }

    const InterpolationTargets targets = interpolation_targets(layout, surface);
    ControlMesh mesh = lift(layout, surface);
    PiaState state = make_pia_state(options, mesh, targets);

    const double tolerance = pia.relative_update_tolerance * state.scale;

    for (int iteration = 1; iteration <= pia.max_iterations; ++iteration) {
        const Eigen::MatrixXd residual = state.targets - state.a * state.vertices;
        state.vertices += residual;
        report.iterations = iteration;

        if (max_row_norm(residual) <= tolerance) {
            report.converged = true;
            break;
        }
    }

    const Eigen::MatrixXd residual = state.a * state.vertices - state.targets;
    fill_report(report, residual);
    if (!report.converged) {
        report.notes.push_back(fmt::format(
            "PIA stopped at the iteration limit of {} with a residual of {:.3e}. Comparing an "
            "unconverged PIA against the direct solve measures the shortfall, not a difference "
            "between two methods.",
            pia.max_iterations,
            report.max_interpolation_error));
    }

    mesh.set_vertices(to_points(state.vertices));
    return mesh;
}

PiaHistory pia_error_history(const DomainLayout& layout,
                             const NurbsSurface& surface,
                             const PiaOptions& pia,
                             const SubdivisionOptions& options) {
    const InterpolationTargets targets = interpolation_targets(layout, surface);

    FitReport direct_report;
    const ControlMesh direct = solve_interpolation(layout, surface, direct_report, options);
    const Eigen::MatrixXd direct_vertices = control_point_matrix(direct);

    // ControlMesh has no default constructor -- a mesh without topology is not
    // a thing this codebase allows -- so the history is built around one.
    PiaHistory history{lift(layout, surface), FitReport{}, {}, {}, {}};
    PiaState state = make_pia_state(options, history.mesh, targets);

    const double tolerance = pia.relative_update_tolerance * state.scale;

    for (int iteration = 1; iteration <= pia.max_iterations; ++iteration) {
        const Eigen::MatrixXd residual = state.targets - state.a * state.vertices;
        state.vertices += residual;

        history.max_update.push_back(max_row_norm(residual));
        history.max_error.push_back(max_row_norm(state.a * state.vertices - state.targets));
        history.distance_to_direct.push_back(max_row_norm(state.vertices - direct_vertices));
        history.report.iterations = iteration;

        if (history.max_update.back() <= tolerance) {
            history.report.converged = true;
            break;
        }
    }

    fill_report(history.report, state.a * state.vertices - state.targets);
    history.mesh.set_vertices(to_points(state.vertices));

    if (!direct_report.converged) {
        history.report.notes.emplace_back(
            "the direct solve this run compares against did not converge, so the "
            "distance-to-direct column is measured against an unreliable reference");
    }

    return history;
}

} // namespace n2s::fit
