#include "n2s/fit/least_squares.hpp"

#include "n2s/fit/fairness.hpp"
#include "n2s/fit/layout.hpp"
#include "n2s/fit/locate.hpp"
#include "n2s/tolerances.hpp"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace n2s::fit {

namespace {

/// One usable sample: where to evaluate the limit surface, and what it should
/// land on.
struct Sample {
    LimitLocation location;
    Eigen::Vector3d target;
};

/// Samples on a regular grid over the layout's domain bounding box, keeping
/// those that fall inside both the trim region and the layout.
///
/// Both rejections are counted rather than merely skipped. A configuration that
/// discards most of its samples still produces a fit, and the only sign that it
/// was built from a handful of points is the count.
std::vector<Sample> gather_samples(const DomainLayout& layout,
                                   const NurbsSurface& surface,
                                   const TrimRegion& region,
                                   const LayoutLocator& locator,
                                   int samples_per_side,
                                   LeastSquaresReport& report) {
    Eigen::Vector2d low = layout.vertices.front();
    Eigen::Vector2d high = layout.vertices.front();
    for (const Eigen::Vector2d& v : layout.vertices) {
        low = low.cwiseMin(v);
        high = high.cwiseMax(v);
    }

    // Polygonised once. TrimRegion::contains re-evaluates every curve per call,
    // which at these sample counts dominates everything else in the fit.
    const std::vector<Eigen::Vector2d> outer = region.outer().polygonise();
    std::vector<std::vector<Eigen::Vector2d>> holes;
    holes.reserve(region.holes().size());
    for (const TrimLoop& hole : region.holes()) {
        holes.push_back(hole.polygonise());
    }
    const auto inside_trim = [&](const Eigen::Vector2d& p) {
        if (!point_in_polygon(outer, p)) {
            return false;
        }
        return std::none_of(
            holes.begin(), holes.end(), [&](const std::vector<Eigen::Vector2d>& hole) {
                return point_in_polygon(hole, p);
            });
    };

    std::vector<Sample> samples;
    const auto side = static_cast<double>(samples_per_side - 1);

    for (int i = 0; i < samples_per_side; ++i) {
        for (int j = 0; j < samples_per_side; ++j) {
            const Eigen::Vector2d p{low.x() + (high.x() - low.x()) * static_cast<double>(i) / side,
                                    low.y() + (high.y() - low.y()) * static_cast<double>(j) / side};

            ++report.samples_requested;

            if (!inside_trim(p)) {
                ++report.samples_outside_trim;
                continue;
            }

            const std::optional<LimitLocation> location = locator.locate(p);
            if (!location.has_value()) {
                ++report.samples_outside_layout;
                continue;
            }

            samples.push_back(Sample{*location, surface.evaluate(p.x(), p.y())});
        }
    }

    return samples;
}

/// Lower-bound estimate of the 2-norm condition number of a symmetric positive
/// definite matrix, from power iteration for the largest eigenvalue and inverse
/// power iteration -- reusing the factorisation -- for the smallest.
///
/// Both iterations approach their eigenvalue from the inside, so the ratio
/// understates the true condition number. That asymmetry is the useful one: a
/// large estimate is evidence, a small one only says nothing was found yet.
template<typename Solver>
double estimate_condition(const Eigen::SparseMatrix<double>& matrix,
                          const Solver& solver,
                          int iterations) {
    const Eigen::Index n = matrix.rows();
    if (n == 0) {
        return 0.0;
    }

    // A fixed, non-random start, so the estimate is reproducible run to run.
    // Not the constant vector: it is orthogonal to too many eigenvectors of the
    // operators that show up here, and the iteration would start on a subspace.
    Eigen::VectorXd v = Eigen::VectorXd::LinSpaced(n, 1.0, 2.0).normalized();

    double largest = 0.0;
    for (int k = 0; k < iterations; ++k) {
        Eigen::VectorXd next = matrix * v;
        const double norm = next.norm();
        if (!(norm > 0.0)) {
            break;
        }
        v = next / norm;
        largest = norm;
    }

    Eigen::VectorXd w = Eigen::VectorXd::LinSpaced(n, 1.0, 2.0).normalized();
    double smallest_inverse = 0.0;
    for (int k = 0; k < iterations; ++k) {
        Eigen::VectorXd next = solver.solve(w);
        const double norm = next.norm();
        if (!(norm > 0.0) || !next.allFinite()) {
            break;
        }
        w = next / norm;
        smallest_inverse = norm;
    }

    if (!(largest > 0.0) || !(smallest_inverse > 0.0)) {
        return 0.0;
    }
    return largest * smallest_inverse;
}

double max_row_norm(const Eigen::MatrixXd& rows) {
    double worst = 0.0;
    for (Eigen::Index i = 0; i < rows.rows(); ++i) {
        worst = std::max(worst, rows.row(i).norm());
    }
    return worst;
}

} // namespace

std::string to_string(SolverUsed solver) {
    switch (solver) {
    case SolverUsed::NormalEquations:
        return "normal_equations";
    case SolverUsed::Qr:
        return "qr";
    case SolverUsed::None:
        break;
    }
    return "none";
}

ControlMesh solve_least_squares(const DomainLayout& layout,
                                const NurbsSurface& surface,
                                const TrimRegion& region,
                                LeastSquaresReport& report,
                                const LeastSquaresOptions& options,
                                const SubdivisionOptions& subdivision) {
    if (layout.quads.empty()) {
        throw std::invalid_argument("cannot fit an empty layout");
    }
    if (options.samples_per_side < 2) {
        throw std::invalid_argument(
            fmt::format("samples_per_side must be at least 2, got {}", options.samples_per_side));
    }
    if (options.lambda < 0.0) {
        throw std::invalid_argument(
            fmt::format("lambda must not be negative, got {}", options.lambda));
    }

    report = LeastSquaresReport{};

    // The boundary, taken from the R1 solve. The square system decouples on the
    // boundary, so its boundary block is the boundary sub-solve exactly; see
    // the header for why that is the constraint.
    FitReport boundary_report;
    const ControlMesh interpolated =
        solve_interpolation(layout, surface, boundary_report, subdivision);
    for (const std::string& note : boundary_report.notes) {
        report.fit.notes.push_back("boundary: " + note);
    }
    if (!boundary_report.converged) {
        report.fit.notes.emplace_back(
            "the boundary solve did not converge, so the constraint this fit is built on is "
            "not the one it claims");
    }

    ControlMesh mesh = interpolated;
    const SubdivisionSurface limit{mesh, subdivision};
    const LayoutLocator locator{layout};

    const std::vector<Sample> samples =
        gather_samples(layout, surface, region, locator, options.samples_per_side, report);
    report.samples_used = samples.size();

    // Partition the vertices. `free` are the interior ones the system solves
    // for; the rest keep the values the boundary solve gave them.
    const auto n = static_cast<Eigen::Index>(mesh.num_vertices());
    std::vector<Eigen::Index> column_of(static_cast<std::size_t>(n), -1);
    std::vector<int> free_vertices;
    for (Eigen::Index v = 0; v < n; ++v) {
        if (!mesh.is_boundary_vertex(static_cast<int>(v))) {
            column_of[static_cast<std::size_t>(v)] =
                static_cast<Eigen::Index>(free_vertices.size());
            free_vertices.push_back(static_cast<int>(v));
        }
    }
    const auto k = static_cast<Eigen::Index>(free_vertices.size());
    report.free_vertices = free_vertices.size();
    report.fixed_vertices = static_cast<std::size_t>(n) - free_vertices.size();

    if (k == 0 || samples.empty()) {
        report.fit.notes.push_back(fmt::format(
            "nothing to solve: {} free control points and {} usable samples. The mesh is the "
            "boundary solve alone.",
            k,
            samples.size()));
        report.fit.converged = samples.empty() && k == 0;
        report.solver = SolverUsed::None;
        const FairnessOperator fairness = umbrella_operator(mesh);
        report.fairness_energy = fairness_energy(fairness, mesh);
        return mesh;
    }

    // A, at the sample locations, split into its free and fixed columns.
    std::vector<LimitLocation> locations;
    locations.reserve(samples.size());
    Eigen::MatrixXd targets(static_cast<Eigen::Index>(samples.size()), 3);
    for (std::size_t s = 0; s < samples.size(); ++s) {
        locations.push_back(samples[s].location);
        targets.row(static_cast<Eigen::Index>(s)) = samples[s].target.transpose();
    }

    const Eigen::SparseMatrix<double, Eigen::RowMajor> full =
        limit.build_limit_matrices(locations, DerivativeOrder::None).position;

    const auto m = static_cast<Eigen::Index>(samples.size());
    Eigen::SparseMatrix<double> free_block(m, k);
    Eigen::MatrixXd fixed_contribution = Eigen::MatrixXd::Zero(m, 3);
    {
        const Eigen::MatrixXd positions = vertex_matrix(mesh);
        std::vector<Eigen::Triplet<double>> triplets;
        for (Eigen::Index row = 0; row < full.outerSize(); ++row) {
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(full, row); it;
                 ++it) {
                const Eigen::Index column = column_of[static_cast<std::size_t>(it.col())];
                if (column >= 0) {
                    triplets.emplace_back(
                        static_cast<int>(row), static_cast<int>(column), it.value());
                } else {
                    fixed_contribution.row(row) += it.value() * positions.row(it.col());
                }
            }
        }
        free_block.setFromTriplets(triplets.begin(), triplets.end());
    }

    // The data the free columns have to account for: what the fixed columns do
    // not already supply.
    const Eigen::MatrixXd data_rhs = targets - fixed_contribution;

    // The fairness term, split the same way. `L_free x = -L_fixed V_fixed` is
    // the fairness residual expressed in the free unknowns alone.
    const FairnessOperator fairness = umbrella_operator(mesh);
    Eigen::SparseMatrix<double> fair_free(fairness.rows(), k);
    Eigen::MatrixXd fair_rhs = Eigen::MatrixXd::Zero(fairness.rows(), 3);
    if (fairness.rows() > 0) {
        const Eigen::MatrixXd positions = vertex_matrix(mesh);
        std::vector<Eigen::Triplet<double>> triplets;
        for (Eigen::Index column = 0; column < fairness.matrix.outerSize(); ++column) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(fairness.matrix, column); it; ++it) {
                const Eigen::Index mapped = column_of[static_cast<std::size_t>(it.col())];
                if (mapped >= 0) {
                    triplets.emplace_back(
                        static_cast<int>(it.row()), static_cast<int>(mapped), it.value());
                } else {
                    fair_rhs.row(it.row()) -= it.value() * positions.row(it.col());
                }
            }
        }
        fair_free.setFromTriplets(triplets.begin(), triplets.end());
    }

    // Weights. Normalising by the row counts makes lambda mean the same thing
    // whatever the sample density, which is what lets the two sweeps vary one
    // axis at a time.
    double data_weight = 1.0;
    double fair_weight = options.lambda;
    if (options.normalise_lambda) {
        data_weight = 1.0 / static_cast<double>(m);
        fair_weight =
            fairness.rows() > 0 ? options.lambda / static_cast<double>(fairness.rows()) : 0.0;
    }

    const Eigen::SparseMatrix<double> data_t = free_block.transpose();
    Eigen::SparseMatrix<double> normal = data_weight * (data_t * free_block);
    Eigen::MatrixXd rhs = data_weight * (data_t * data_rhs);
    if (fairness.rows() > 0 && fair_weight > 0.0) {
        const Eigen::SparseMatrix<double> fair_t = fair_free.transpose();
        normal += fair_weight * (fair_t * fair_free);
        rhs += fair_weight * (fair_t * fair_rhs);
    }
    normal.makeCompressed();

    Eigen::MatrixXd solution;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
    ldlt.compute(normal);

    bool use_qr = ldlt.info() != Eigen::Success;
    if (!use_qr) {
        report.condition_estimate = estimate_condition(normal, ldlt, options.condition_iterations);
        use_qr = report.condition_estimate > options.max_condition;
        if (use_qr) {
            report.fit.notes.push_back(fmt::format(
                "the normal equations estimate a condition number of at least {:.3g}, above the "
                "{:.3g} allowed, so the stacked system was solved with SparseQR instead. QR works "
                "with the square root of that conditioning.",
                report.condition_estimate,
                options.max_condition));
        }
    } else {
        report.fit.notes.emplace_back(
            "the normal equations were not positive definite, so SparseQR solved the stacked "
            "system instead");
    }

    if (!use_qr) {
        solution = ldlt.solve(rhs);
        report.solver = SolverUsed::NormalEquations;
        if (ldlt.info() != Eigen::Success) {
            use_qr = true;
            report.fit.notes.emplace_back(
                "the LDLT back-substitution failed, so SparseQR solved the stacked system");
        }
    }

    if (use_qr) {
        // Stacked, not normalised: QR on [sqrt(w) A ; sqrt(lambda') L] avoids
        // forming A^T A at all, and it is the squaring in that product that
        // costs the digits worth rescuing here.
        const double data_root = std::sqrt(data_weight);
        const double fair_root = std::sqrt(fair_weight);
        const bool with_fairness = fairness.rows() > 0 && fair_weight > 0.0;
        const Eigen::Index stacked_rows = m + (with_fairness ? fairness.rows() : 0);

        std::vector<Eigen::Triplet<double>> triplets;
        for (Eigen::Index column = 0; column < free_block.outerSize(); ++column) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(free_block, column); it; ++it) {
                triplets.emplace_back(
                    static_cast<int>(it.row()), static_cast<int>(it.col()), data_root * it.value());
            }
        }
        if (with_fairness) {
            for (Eigen::Index column = 0; column < fair_free.outerSize(); ++column) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(fair_free, column); it; ++it) {
                    triplets.emplace_back(static_cast<int>(m + it.row()),
                                          static_cast<int>(it.col()),
                                          fair_root * it.value());
                }
            }
        }

        Eigen::SparseMatrix<double> stacked(stacked_rows, k);
        stacked.setFromTriplets(triplets.begin(), triplets.end());
        stacked.makeCompressed();

        Eigen::MatrixXd stacked_rhs(stacked_rows, 3);
        stacked_rhs.topRows(m) = data_root * data_rhs;
        if (with_fairness) {
            stacked_rhs.bottomRows(fairness.rows()) = fair_root * fair_rhs;
        }

        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr;
        qr.compute(stacked);
        if (qr.info() != Eigen::Success) {
            report.fit.converged = false;
            report.solver = SolverUsed::Qr;
            report.fit.notes.emplace_back(
                "SparseQR could not factorise the stacked system either; the mesh returned is "
                "the boundary solve with the layout lifted inside it");
            return mesh;
        }
        solution = qr.solve(stacked_rhs);
        report.solver = SolverUsed::Qr;
    }

    if (!solution.allFinite()) {
        report.fit.converged = false;
        report.fit.notes.emplace_back(
            "the solve produced non-finite control points, so the result was discarded");
        return mesh;
    }

    std::vector<Eigen::Vector3d> fitted = mesh.vertices();
    for (std::size_t i = 0; i < free_vertices.size(); ++i) {
        fitted[static_cast<std::size_t>(free_vertices[i])] =
            solution.row(static_cast<Eigen::Index>(i)).transpose();
    }
    mesh.set_vertices(std::move(fitted));

    // The residual, measured on the mesh actually returned rather than inferred
    // from the solve, so a mistake in writing the solution back shows up here.
    const Eigen::MatrixXd residual = full * vertex_matrix(mesh) - targets;
    report.fit.max_interpolation_error = max_row_norm(residual);
    report.fit.rms_interpolation_error =
        std::sqrt(residual.rowwise().squaredNorm().sum() / static_cast<double>(m));
    report.fit.iterations = 1;
    report.fit.converged = true;
    report.fairness_energy = fairness_energy(umbrella_operator(mesh), mesh);

    return mesh;
}

} // namespace n2s::fit
