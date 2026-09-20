#include "n2s/fit/interpolate.hpp"
#include "n2s/fit/layout.hpp"
#include "n2s/fit/refine.hpp"
#include "n2s/io/case_json.hpp"
#include "n2s/metrics/surface_error.hpp"
#include "n2s/trim/cases.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <vector>

using Catch::Approx;
using n2s::ControlMesh;
using n2s::NurbsSurface;
using n2s::SubdivisionSurface;
using n2s::fit::DomainLayout;
using n2s::fit::FitReport;

namespace {

/// Largest distance between a fitted limit surface and the points it was asked
/// to pass through. This is what "interpolation" has to mean; the residual the
/// solver reports is the same quantity computed from the linear algebra, so
/// checking it independently here guards against the matrix and the evaluator
/// disagreeing.
double interpolation_error(const ControlMesh& mesh,
                           const DomainLayout& layout,
                           const NurbsSurface& surface) {
    const n2s::fit::InterpolationTargets targets = n2s::fit::interpolation_targets(layout, surface);
    const SubdivisionSurface subdivision{mesh};
    const std::vector<n2s::LimitSample> samples =
        subdivision.evaluate_limit(targets.locations, n2s::DerivativeOrder::None);

    double worst = 0.0;
    for (std::size_t i = 0; i < samples.size(); ++i) {
        worst = std::max(worst, (samples[i].position - targets.points[i]).norm());
    }
    return worst;
}

DomainLayout refined_grid(int rows, int columns, int factor = 1) {
    return n2s::fit::refine_quads(n2s::fit::grid_domain_layout(rows, columns), factor);
}

} // namespace

// ---------------------------------------------------------------------------
// Refinement.
// ---------------------------------------------------------------------------

TEST_CASE("thirding multiplies the faces by nine", "[fit][refine]") {
    const DomainLayout base = n2s::fit::grid_domain_layout(3, 3); // 4 quads
    const DomainLayout thirded = n2s::fit::third_edges(base);

    CHECK(base.num_quads() == 4);
    CHECK(thirded.num_quads() == 36);

    // A 3x3 grid thirded is a 7x7 grid of vertices.
    CHECK(thirded.num_vertices() == 49);
}

TEST_CASE("refinement shares the vertices along a common edge", "[fit][refine]") {
    // A duplicated vertex on a shared edge is a crack in the control mesh, and
    // Catmull-Clark would treat the two sides as separate boundaries. The
    // vertex count is the cheapest way to detect it: an unshared refinement of
    // a 2x1 layout would produce 32 vertices rather than 28.
    DomainLayout two_quads;
    two_quads.vertices = {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {2, 0}, {2, 1}};
    two_quads.quads = {{0, 1, 2, 3}, {1, 4, 5, 2}};

    const DomainLayout thirded = n2s::fit::third_edges(two_quads);
    CHECK(thirded.num_quads() == 18);
    CHECK(thirded.num_vertices() == 4 * 7); // a 4-by-7 grid of vertices

    // Every position appears exactly once.
    for (std::size_t i = 0; i < thirded.num_vertices(); ++i) {
        for (std::size_t j = i + 1; j < thirded.num_vertices(); ++j) {
            INFO("vertices " << i << " and " << j);
            REQUIRE((thirded.vertices[i] - thirded.vertices[j]).norm() > 1e-12);
        }
    }
}

TEST_CASE("refinement preserves the region the layout covers", "[fit][refine]") {
    const DomainLayout base = n2s::fit::grid_domain_layout(3, 4);
    const DomainLayout thirded = n2s::fit::third_edges(base);

    const auto area = [](const DomainLayout& layout) {
        double total = 0.0;
        for (const std::array<int, 4>& q : layout.quads) {
            const Eigen::Vector2d& a = layout.vertices[static_cast<std::size_t>(q[0])];
            const Eigen::Vector2d& b = layout.vertices[static_cast<std::size_t>(q[1])];
            const Eigen::Vector2d& c = layout.vertices[static_cast<std::size_t>(q[2])];
            const Eigen::Vector2d& d = layout.vertices[static_cast<std::size_t>(q[3])];
            const auto cross = [](const Eigen::Vector2d& p, const Eigen::Vector2d& q2) {
                return p.x() * q2.y() - p.y() * q2.x();
            };
            total += 0.5 * std::abs(cross(b - a, c - a)) + 0.5 * std::abs(cross(c - a, d - a));
        }
        return total;
    };

    CHECK(area(thirded) == Approx(area(base)).epsilon(1e-12));
}

TEST_CASE("a factor of one is the identity, and a bad factor is rejected", "[fit][refine]") {
    const DomainLayout base = n2s::fit::grid_domain_layout(3, 3);
    CHECK(n2s::fit::refine_quads(base, 1).num_quads() == base.num_quads());
    CHECK_THROWS_AS(n2s::fit::refine_quads(base, 0), std::invalid_argument);
    CHECK_THROWS_AS(n2s::fit::refine_quads(DomainLayout{}, 3), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// The direct solve.
// ---------------------------------------------------------------------------

TEST_CASE("the fitted limit surface passes through its targets", "[fit][interpolate]") {
    // The defining property. Before fitting, the limit surface of a lifted
    // layout misses its own vertices by the usual Catmull-Clark shrinkage;
    // afterwards it passes through them to solver precision.
    const NurbsSurface surface = n2s::cases::saddle(1.0, 0.6);
    const DomainLayout layout = refined_grid(4, 4, 3);

    const double before = interpolation_error(n2s::fit::lift(layout, surface), layout, surface);

    FitReport report;
    const ControlMesh fitted = n2s::fit::solve_interpolation(layout, surface, report);
    const double after = interpolation_error(fitted, layout, surface);

    INFO("before " << before << ", after " << after << ", reported residual "
                   << report.max_interpolation_error);
    CHECK(before > 1e-3);
    CHECK(after < 1e-10);
    CHECK(report.converged);

    // The solver's own residual and an independent evaluation agree, so the
    // matrix and the limit evaluator are consistent.
    CHECK(report.max_interpolation_error == Approx(after).margin(1e-11));
}

TEST_CASE("interpolation is exact on a plane", "[fit][interpolate]") {
    // A planar layout on a plane is already its own solution, so the fit must
    // leave it where it is rather than drifting.
    const NurbsSurface plane = n2s::cases::plane();
    const DomainLayout layout = refined_grid(4, 4, 2);

    FitReport report;
    const ControlMesh fitted = n2s::fit::solve_interpolation(layout, plane, report);

    CHECK(report.converged);
    for (const Eigen::Vector3d& v : fitted.vertices()) {
        CHECK(std::abs(v.z()) < 1e-12);
    }
}

TEST_CASE("fitting lowers the measured error against the NURBS", "[fit][interpolate]") {
    // Interpolating the layout vertices is not the same as approximating the
    // surface well everywhere, so this is worth measuring rather than assuming.
    const NurbsSurface surface = n2s::cases::saddle(1.0, 0.6);
    const n2s::TrimRegion region{n2s::cases::rectangle_loop({0.05, 0.05}, {0.95, 0.95})};
    const DomainLayout layout = refined_grid(4, 4, 2);

    n2s::metrics::SurfaceErrorOptions options;
    options.samples_per_face = 2;
    options.nurbs_samples_per_side = 12;
    options.restrict_to_trimmed_region = false;

    const n2s::metrics::DomainMap map = n2s::fit::bilinear_domain_map(layout);

    const auto measure = [&](const ControlMesh& mesh) {
        return n2s::metrics::measure_surface_error(
                   SubdivisionSurface{mesh}, surface, region, map, options)
            .geometric.rms;
    };

    FitReport report;
    const double unfitted = measure(n2s::fit::lift(layout, surface));
    const double fitted = measure(n2s::fit::solve_interpolation(layout, surface, report));

    INFO("unfitted rms " << unfitted << ", fitted rms " << fitted);
    CHECK(fitted < unfitted);
}

// ---------------------------------------------------------------------------
// PIA, and defect 5 of the thesis.
// ---------------------------------------------------------------------------

TEST_CASE("PIA converges to the direct solution", "[fit][pia][defect5]") {
    // Defect 5: the thesis compared progressive iteration against the direct
    // solve as though they were different methods. They are the same square
    // system, so PIA converges *to* the direct solution and any difference
    // between them is PIA not having converged yet.
    //
    // This is the test that settles it. If PIA ever converged to something
    // else, the distance below would plateau at a non-zero value.
    const NurbsSurface surface = n2s::cases::saddle(1.0, 0.6);
    const DomainLayout layout = refined_grid(4, 4, 2);

    n2s::fit::PiaOptions pia;
    pia.max_iterations = 400;

    const n2s::fit::PiaHistory history = n2s::fit::pia_error_history(layout, surface, pia);

    REQUIRE(history.distance_to_direct.size() > 5);
    INFO("after " << history.report.iterations << " iterations: error "
                  << history.report.max_interpolation_error << ", distance to direct "
                  << history.distance_to_direct.back());

    // It gets there.
    CHECK(history.report.converged);
    CHECK(history.distance_to_direct.back() < 1e-9);

    // And it approaches monotonically, rather than wandering.
    for (std::size_t i = 1; i < history.distance_to_direct.size(); ++i) {
        INFO("iteration " << i);
        REQUIRE(history.distance_to_direct[i] <= history.distance_to_direct[i - 1] * 1.0001);
    }

    // The early iterations are genuinely far off, so the comparison is not
    // vacuous: a difference between the two methods really is visible before
    // convergence, which is precisely what the thesis measured and reported as
    // a difference in method.
    CHECK(history.distance_to_direct.front() > 1e-4);
}

TEST_CASE("an unconverged PIA says so", "[fit][pia]") {
    const NurbsSurface surface = n2s::cases::saddle(1.0, 0.6);
    const DomainLayout layout = refined_grid(4, 4, 2);

    n2s::fit::PiaOptions pia;
    pia.max_iterations = 2;

    FitReport report;
    n2s::fit::solve_pia(layout, surface, report, pia);

    CHECK_FALSE(report.converged);
    CHECK(report.iterations == 2);
    REQUIRE_FALSE(report.notes.empty());
    CHECK(report.notes.front().find("iteration limit") != std::string::npos);
}

TEST_CASE("PIA and the direct solve agree once converged", "[fit][pia][defect5]") {
    const NurbsSurface surface = n2s::cases::paraboloid(1.0, 0.5);
    const DomainLayout layout = refined_grid(3, 4, 3);

    FitReport direct_report;
    const ControlMesh direct = n2s::fit::solve_interpolation(layout, surface, direct_report);

    n2s::fit::PiaOptions pia;
    pia.max_iterations = 600;
    FitReport pia_report;
    const ControlMesh iterated = n2s::fit::solve_pia(layout, surface, pia_report, pia);

    REQUIRE(direct_report.converged);
    REQUIRE(pia_report.converged);

    double worst = 0.0;
    for (std::size_t i = 0; i < direct.num_vertices(); ++i) {
        worst = std::max(worst, (direct.vertices()[i] - iterated.vertices()[i]).norm());
    }

    INFO("worst control point difference " << worst << " after " << pia_report.iterations
                                           << " PIA iterations");
    CHECK(worst < 1e-9);
}

// ---------------------------------------------------------------------------
// The thesis case.
// ---------------------------------------------------------------------------

TEST_CASE("the thesis layout fits once it has interior freedom", "[fit][interpolate]") {
    // The hand-authored layout has 26 vertices of which 22 are on the boundary,
    // leaving four free interior control points. Thirding is what gives the
    // interior enough freedom for a fit to mean anything, and it is the
    // refinement the thesis itself used.
    const std::filesystem::path path = std::filesystem::path{N2S_DATA_DIR} / "thesis_doublevb.json";
    if (!std::filesystem::exists(path)) {
        SKIP("the thesis case file is not present");
    }

    const n2s::io::Case loaded = n2s::io::read_case(path);
    REQUIRE(loaded.layout.has_value());
    CHECK(loaded.layout->num_quads() == 14);
    CHECK(loaded.layout->num_vertices() == 26);

    const DomainLayout thirded = n2s::fit::third_edges(*loaded.layout);
    CHECK(thirded.num_quads() == 126);

    FitReport report;
    const ControlMesh fitted = n2s::fit::solve_interpolation(thirded, loaded.surface, report);

    INFO("residual " << report.max_interpolation_error);
    for (const std::string& note : report.notes) {
        INFO(note);
    }
    CHECK(report.converged);
    CHECK(interpolation_error(fitted, thirded, loaded.surface) < 1e-9);
}
