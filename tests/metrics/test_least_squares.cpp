#include "n2s/fit/fairness.hpp"
#include "n2s/fit/layout.hpp"
#include "n2s/fit/least_squares.hpp"
#include "n2s/fit/refine.hpp"
#include "n2s/io/case_json.hpp"
#include "n2s/metrics/surface_error.hpp"
#include "n2s/trim/cases.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <vector>

using Catch::Approx;
using n2s::ControlMesh;
using n2s::NurbsSurface;
using n2s::TrimRegion;
using n2s::fit::DomainLayout;
using n2s::fit::LeastSquaresOptions;
using n2s::fit::LeastSquaresReport;

namespace {

TrimRegion unit_square() {
    return TrimRegion{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0})};
}

LeastSquaresOptions cheap() {
    LeastSquaresOptions options;
    options.samples_per_side = 24;
    return options;
}

/// Max distance from the limit surface to the NURBS over the layout, measured
/// through the shared correspondence rather than by projection, so the number
/// is comparable between two fits of the same layout.
double
parametric_max(const ControlMesh& mesh, const DomainLayout& layout, const NurbsSurface& surface) {
    const n2s::SubdivisionSurface limit{mesh};
    const n2s::metrics::DomainMap map = n2s::fit::bilinear_domain_map(layout);

    std::vector<n2s::LimitLocation> locations;
    for (std::size_t f = 0; f < layout.num_quads(); ++f) {
        for (const double u : {0.125, 0.375, 0.625, 0.875}) {
            for (const double v : {0.125, 0.375, 0.625, 0.875}) {
                locations.push_back(n2s::LimitLocation{static_cast<int>(f), u, v});
            }
        }
    }

    const std::vector<n2s::LimitSample> samples =
        limit.evaluate_limit(locations, n2s::DerivativeOrder::None);

    double worst = 0.0;
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const Eigen::Vector2d domain = map(locations[i]);
        const Eigen::Vector3d expected = surface.evaluate(domain.x(), domain.y());
        worst = std::max(worst, (samples[i].position - expected).norm());
    }
    return worst;
}

DomainLayout thesis_layout() {
    const std::filesystem::path path = std::filesystem::path{N2S_DATA_DIR} / "thesis_doublevb.json";
    const n2s::io::Case loaded = n2s::io::read_case(path);
    REQUIRE(loaded.layout.has_value());
    return *loaded.layout;
}

} // namespace

TEST_CASE("a plane is fitted exactly", "[fit][leastsquares][oracle]") {
    // The sharpest oracle available: a plane is in the span of the limit
    // surfaces of a planar control net, so the least-squares residual has to be
    // zero to solver precision. Any non-zero answer is the assembly, not the
    // approximation.
    const NurbsSurface plane = n2s::cases::plane(1.0, 1.0);
    const DomainLayout layout = n2s::fit::grid_domain_layout(6, 6);

    LeastSquaresReport report;
    const ControlMesh fitted =
        n2s::fit::solve_least_squares(layout, plane, unit_square(), report, cheap());

    INFO("solver " << n2s::fit::to_string(report.solver) << ", samples " << report.samples_used);
    REQUIRE(report.fit.converged);
    CHECK(report.samples_used > 0);
    CHECK(report.fit.max_interpolation_error < 1e-10);
    CHECK(parametric_max(fitted, layout, plane) < 1e-10);
}

TEST_CASE("the boundary is left exactly where the interpolating solve put it",
          "[fit][leastsquares]") {
    // The constraint this fit is built on. R2 and R1 must agree on the boundary
    // vertex for vertex, or comparing them on one layout compares two things at
    // once.
    const NurbsSurface saddle = n2s::cases::saddle(1.0, 0.5);
    const DomainLayout layout = n2s::fit::grid_domain_layout(6, 6);

    n2s::fit::FitReport interpolation;
    const ControlMesh interpolated = n2s::fit::solve_interpolation(layout, saddle, interpolation);

    LeastSquaresOptions options = cheap();
    options.lambda = 1e-3;
    LeastSquaresReport report;
    const ControlMesh fitted =
        n2s::fit::solve_least_squares(layout, saddle, unit_square(), report, options);

    REQUIRE(fitted.num_vertices() == interpolated.num_vertices());
    std::size_t boundary = 0;
    for (std::size_t v = 0; v < fitted.num_vertices(); ++v) {
        if (!fitted.is_boundary_vertex(static_cast<int>(v))) {
            continue;
        }
        ++boundary;
        INFO("boundary vertex " << v);
        CHECK((fitted.vertices()[v] - interpolated.vertices()[v]).norm() ==
              Approx(0.0).margin(0.0));
    }
    REQUIRE(boundary > 0);
    CHECK(report.fixed_vertices == boundary);
    CHECK(report.free_vertices == fitted.num_vertices() - boundary);
}

TEST_CASE("least squares beats interpolation away from the vertices",
          "[fit][leastsquares][oracle]") {
    // Interpolation pins the surface at the layout vertices and lets it do what
    // it likes between them; least squares spends that freedom on the samples
    // in between. On a curved surface sampled densely, the second has to win on
    // the samples the first never saw.
    const NurbsSurface saddle = n2s::cases::saddle(1.0, 0.5);
    const DomainLayout layout = n2s::fit::grid_domain_layout(6, 6);

    n2s::fit::FitReport interpolation;
    const ControlMesh interpolated = n2s::fit::solve_interpolation(layout, saddle, interpolation);

    LeastSquaresReport report;
    const ControlMesh fitted =
        n2s::fit::solve_least_squares(layout, saddle, unit_square(), report, cheap());

    const double interpolating = parametric_max(interpolated, layout, saddle);
    const double least_squares = parametric_max(fitted, layout, saddle);

    INFO("interpolating " << interpolating << ", least squares " << least_squares);
    REQUIRE(report.fit.converged);
    CHECK(least_squares < interpolating);
}

TEST_CASE("raising lambda lowers the fairness energy and raises the residual",
          "[fit][leastsquares][oracle]") {
    // The trade-off the parameter exists to make. Both directions are checked:
    // a lambda that only ever lowered the energy without costing anything would
    // mean the data term was not connected to the solve.
    const NurbsSurface saddle = n2s::cases::saddle(1.0, 0.5);
    const DomainLayout layout = n2s::fit::grid_domain_layout(7, 7);

    double previous_energy = std::numeric_limits<double>::infinity();
    double previous_residual = 0.0;

    for (const double lambda : {0.0, 1e-6, 1e-3, 1.0}) {
        LeastSquaresOptions options = cheap();
        options.lambda = lambda;

        LeastSquaresReport report;
        const ControlMesh fitted =
            n2s::fit::solve_least_squares(layout, saddle, unit_square(), report, options);

        INFO("lambda " << lambda << " energy " << report.fairness_energy << " residual "
                       << report.fit.rms_interpolation_error);
        REQUIRE(report.fit.converged);

        CHECK(report.fairness_energy <= previous_energy * (1.0 + 1e-9));
        CHECK(report.fit.rms_interpolation_error >= previous_residual * (1.0 - 1e-9));

        previous_energy = report.fairness_energy;
        previous_residual = report.fit.rms_interpolation_error;
    }
}

TEST_CASE("a fairness-only solve reaches a near-zero energy", "[fit][leastsquares][oracle]") {
    // Driving lambda up until the data term is negligible leaves the umbrella
    // as the only thing being minimised, and its minimum subject to a fixed
    // boundary is the discrete thin-plate solution, whose energy is far below
    // the fit's. Checked as a large ratio rather than an absolute value,
    // because the boundary keeps it from reaching zero.
    const NurbsSurface saddle = n2s::cases::saddle(1.0, 0.5);
    const DomainLayout layout = n2s::fit::grid_domain_layout(7, 7);

    LeastSquaresOptions plain = cheap();
    LeastSquaresReport unfaired;
    n2s::fit::solve_least_squares(layout, saddle, unit_square(), unfaired, plain);

    LeastSquaresOptions heavy = cheap();
    heavy.lambda = 1e8;
    LeastSquaresReport faired;
    n2s::fit::solve_least_squares(layout, saddle, unit_square(), faired, heavy);

    INFO("unfaired " << unfaired.fairness_energy << ", faired " << faired.fairness_energy);
    REQUIRE(unfaired.fairness_energy > 0.0);
    CHECK(faired.fairness_energy < unfaired.fairness_energy * 1e-3);
}

TEST_CASE("samples outside the trim and outside the layout are counted", "[fit][leastsquares]") {
    // Both rejections have to be visible. A configuration that throws most of
    // its samples away still returns a fit, and the count is the only thing
    // that says it was built from a handful of points.
    const NurbsSurface plane = n2s::cases::plane(1.0, 1.0);
    const DomainLayout layout = n2s::fit::grid_domain_layout(5, 5);

    // A trim covering a quarter of the layout, so most samples fall outside it.
    const TrimRegion corner{n2s::cases::rectangle_loop({0.0, 0.0}, {0.5, 0.5})};

    LeastSquaresReport report;
    n2s::fit::solve_least_squares(layout, plane, corner, report, cheap());

    CHECK(report.samples_requested == 24 * 24);
    CHECK(report.samples_outside_trim > 0);
    CHECK(report.samples_used > 0);
    CHECK(report.samples_used + report.samples_outside_trim + report.samples_outside_layout ==
          report.samples_requested);
}

TEST_CASE("the condition estimate is reported and the solver named", "[fit][leastsquares]") {
    const NurbsSurface saddle = n2s::cases::saddle(1.0, 0.5);
    const DomainLayout layout = n2s::fit::grid_domain_layout(6, 6);

    LeastSquaresReport report;
    n2s::fit::solve_least_squares(layout, saddle, unit_square(), report, cheap());

    CHECK(report.solver == n2s::fit::SolverUsed::NormalEquations);
    CHECK(report.condition_estimate >= 1.0);
    CHECK(std::isfinite(report.condition_estimate));
}

TEST_CASE("an impossible condition limit forces the QR fallback", "[fit][leastsquares]") {
    // The fallback path has to be exercised by something, and a threshold below
    // 1 is the cheapest way to demand it. The answer must be the same fit, not
    // merely a different code path that terminates.
    const NurbsSurface saddle = n2s::cases::saddle(1.0, 0.5);
    const DomainLayout layout = n2s::fit::grid_domain_layout(6, 6);

    LeastSquaresOptions normal = cheap();
    LeastSquaresReport via_normal;
    const ControlMesh a =
        n2s::fit::solve_least_squares(layout, saddle, unit_square(), via_normal, normal);

    LeastSquaresOptions forced = cheap();
    forced.max_condition = 0.5;
    LeastSquaresReport via_qr;
    const ControlMesh b =
        n2s::fit::solve_least_squares(layout, saddle, unit_square(), via_qr, forced);

    REQUIRE(via_normal.solver == n2s::fit::SolverUsed::NormalEquations);
    REQUIRE(via_qr.solver == n2s::fit::SolverUsed::Qr);
    REQUIRE(via_qr.fit.converged);
    CHECK_FALSE(via_qr.fit.notes.empty());

    REQUIRE(a.num_vertices() == b.num_vertices());
    double worst = 0.0;
    for (std::size_t v = 0; v < a.num_vertices(); ++v) {
        worst = std::max(worst, (a.vertices()[v] - b.vertices()[v]).norm());
    }
    INFO("largest disagreement between the two solvers " << worst);
    CHECK(worst < 1e-9);
}

TEST_CASE("the thesis layout fits, and the boundary is R1's", "[fit][leastsquares]") {
    // The case the whole milestone is about. Thirded, as R1 measured it.
    const DomainLayout layout = n2s::fit::refine_quads(thesis_layout(), 3);
    const std::filesystem::path path = std::filesystem::path{N2S_DATA_DIR} / "thesis_doublevb.json";
    const n2s::io::Case loaded = n2s::io::read_case(path);

    LeastSquaresOptions options;
    options.samples_per_side = 48;
    options.lambda = 1e-4;

    LeastSquaresReport report;
    const ControlMesh fitted =
        n2s::fit::solve_least_squares(layout, loaded.surface, loaded.region, report, options);

    INFO("samples " << report.samples_used << " of " << report.samples_requested << ", free "
                    << report.free_vertices << ", condition " << report.condition_estimate);
    REQUIRE(report.fit.converged);
    CHECK(report.samples_used > report.free_vertices); // over-determined, as intended
    CHECK(report.free_vertices > 0);
    CHECK(std::isfinite(report.fit.rms_interpolation_error));

    n2s::fit::FitReport interpolation;
    const ControlMesh interpolated =
        n2s::fit::solve_interpolation(layout, loaded.surface, interpolation);
    for (std::size_t v = 0; v < fitted.num_vertices(); ++v) {
        if (fitted.is_boundary_vertex(static_cast<int>(v))) {
            INFO("boundary vertex " << v);
            CHECK((fitted.vertices()[v] - interpolated.vertices()[v]).norm() ==
                  Approx(0.0).margin(0.0));
        }
    }
}

TEST_CASE("bad options are rejected by name", "[fit][leastsquares]") {
    const NurbsSurface plane = n2s::cases::plane(1.0, 1.0);
    const DomainLayout layout = n2s::fit::grid_domain_layout(4, 4);
    LeastSquaresReport report;

    LeastSquaresOptions too_few = cheap();
    too_few.samples_per_side = 1;
    CHECK_THROWS_AS(n2s::fit::solve_least_squares(layout, plane, unit_square(), report, too_few),
                    std::invalid_argument);

    LeastSquaresOptions negative = cheap();
    negative.lambda = -1.0;
    CHECK_THROWS_AS(n2s::fit::solve_least_squares(layout, plane, unit_square(), report, negative),
                    std::invalid_argument);
}
