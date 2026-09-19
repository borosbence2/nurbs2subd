#include "n2s/trim/cases.hpp"
#include "n2s/trim/sampling.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

using Catch::Approx;
using n2s::NurbsCurve2;
using n2s::NurbsSurface;
using n2s::SamplingMode;
using n2s::SamplingOptions;
using n2s::TrimLoop;

namespace {

/// Longest model-space gap between consecutive samples of a loop.
double worst_segment_length(const TrimLoop& loop,
                            const NurbsSurface& surface,
                            const SamplingOptions& options) {
    const std::vector<Eigen::Vector2d> points = n2s::sample_loop(loop, surface, options);
    double worst = 0.0;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const Eigen::Vector2d& a = points[i];
        const Eigen::Vector2d& b = points[(i + 1) % points.size()];
        worst = std::max(worst,
                         (surface.evaluate(b.x(), b.y()) - surface.evaluate(a.x(), a.y())).norm());
    }
    return worst;
}

/// Worst model-space deviation of the true curve from the sampled polyline,
/// measured at the midpoint of every segment. This is the quantity
/// `max_sagitta` is supposed to bound.
double worst_midpoint_deviation(const NurbsCurve2& curve,
                                const NurbsSurface& surface,
                                const SamplingOptions& options) {
    const std::vector<double> parameters = n2s::sample_parameters(curve, surface, options);

    const auto lift = [&](double t) {
        const Eigen::Vector2d p = curve.evaluate(t);
        return surface.evaluate(p.x(), p.y());
    };

    double worst = 0.0;
    for (std::size_t i = 0; i + 1 < parameters.size(); ++i) {
        const double a = parameters[i];
        const double b = parameters[i + 1];
        const Eigen::Vector3d pa = lift(a);
        const Eigen::Vector3d pb = lift(b);
        const Eigen::Vector3d pm = lift(0.5 * (a + b));
        worst = std::max(worst, (pm - 0.5 * (pa + pb)).norm());
    }
    return worst;
}

} // namespace

TEST_CASE("every knot survives sampling, in both modes", "[trim][sampling]") {
    // A knot is where continuity drops. Stepping over one can round off a
    // corner, and no chord or sagitta tolerance would notice, because both are
    // measured between samples that have already skipped it.
    const NurbsSurface surface = n2s::cases::plane();
    const NurbsCurve2 curve{
        n2s::KnotVector{2, {0, 0, 0, 0.25, 0.25, 0.6, 1, 1, 1}},
        {{0.0, 0.0}, {0.2, 0.5}, {0.5, 0.1}, {0.7, 0.8}, {0.9, 0.2}, {1.0, 0.6}},
        {1.0, 1.0, 2.0, 1.0, 1.0, 1.0}};

    for (const SamplingMode mode : {SamplingMode::Adaptive, SamplingMode::Uniform}) {
        SamplingOptions options;
        options.mode = mode;
        options.max_segment_length = 0.5; // deliberately loose
        options.max_sagitta = 0.5;

        const std::vector<double> parameters = n2s::sample_parameters(curve, surface, options);

        for (const double knot : {0.0, 0.25, 0.6, 1.0}) {
            const bool present = std::any_of(parameters.begin(), parameters.end(), [&](double t) {
                return std::abs(t - knot) < 1e-12;
            });
            INFO("mode " << static_cast<int>(mode) << ", knot " << knot);
            CHECK(present);
        }
    }
}

TEST_CASE("adaptive sampling respects the model-space segment length bound", "[trim][sampling]") {
    // The cylinder is the interesting case: equal steps in the parametric
    // domain are equal steps in arc length around it, but the domain is the
    // unit square whatever the radius, so only a 3D-aware sampler gets the
    // density right.
    const NurbsSurface surface = n2s::cases::cylinder(3.0, 2.0);
    const TrimLoop loop = n2s::cases::rectangle_loop({0.1, 0.1}, {0.9, 0.9});

    for (const double bound : {0.5, 0.2, 0.05}) {
        SamplingOptions options;
        options.max_segment_length = bound;
        options.max_sagitta = bound; // keep the sagitta from being the binding constraint

        const double worst = worst_segment_length(loop, surface, options);
        INFO("bound " << bound << ", worst segment " << worst);
        CHECK(worst <= bound * 1.0001);
    }
}

TEST_CASE("adaptive sampling respects the sagitta bound", "[trim][sampling]") {
    const NurbsSurface surface = n2s::cases::cylinder(2.0, 1.0);
    // A trim curve running all the way around the cylinder: maximally curved.
    const NurbsCurve2 around = n2s::cases::rectangle_loop({0.0, 0.2}, {1.0, 0.8}).curves().front();

    double previous = 1.0;
    for (const double bound : {0.05, 0.01, 0.002}) {
        SamplingOptions options;
        options.max_segment_length = 10.0; // leave the sagitta in charge
        options.max_sagitta = bound;

        const double worst = worst_midpoint_deviation(around, surface, options);
        INFO("sagitta bound " << bound << ", worst deviation " << worst);
        CHECK(worst <= bound);
        CHECK(worst < previous);
        previous = worst;
    }
}

TEST_CASE("adaptive sampling puts points where the curvature is", "[trim][sampling]") {
    // On a cylinder a trim line along v is a straight ruling in 3D, while one
    // along u wraps around the circumference. With a curvature-driven sampler
    // the curved one must collect far more samples; with the uniform sampler
    // the thesis used, both get the same number regardless.
    const NurbsSurface surface = n2s::cases::cylinder(2.0, 2.0);
    const TrimLoop rect = n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0});

    SamplingOptions adaptive;
    adaptive.max_segment_length = 1.0;
    adaptive.max_sagitta = 0.002;

    const NurbsCurve2& along_u = rect.curves()[0]; // (0,0) -> (1,0), around
    const NurbsCurve2& along_v = rect.curves()[1]; // (1,0) -> (1,1), straight ruling

    const std::size_t curved = n2s::sample_parameters(along_u, surface, adaptive).size();
    const std::size_t straight = n2s::sample_parameters(along_v, surface, adaptive).size();

    INFO("curved edge " << curved << " samples, straight edge " << straight);
    CHECK(curved > 4 * straight);

    SamplingOptions uniform;
    uniform.mode = SamplingMode::Uniform;
    uniform.samples_per_span = 8;
    CHECK(n2s::sample_parameters(along_u, surface, uniform).size() ==
          n2s::sample_parameters(along_v, surface, uniform).size());
}

TEST_CASE("the sampled loop length converges to the true length", "[trim][sampling]") {
    // A trim circle of radius r in the domain of a plane of unit extent is a
    // circle of radius r in space, so its length is known exactly.
    const NurbsSurface surface = n2s::cases::plane(1.0, 1.0);
    const double radius = 0.3;
    const TrimLoop circle = n2s::cases::circle_loop({0.5, 0.5}, radius);
    const double exact = 2.0 * std::numbers::pi * radius;

    double previous_error = 1.0;
    for (const double bound : {0.02, 0.005, 0.001}) {
        SamplingOptions options;
        options.max_segment_length = 10.0;
        options.max_sagitta = bound;

        const double length = n2s::sampled_loop_length(circle, surface, options);
        const double error = std::abs(length - exact);
        INFO("sagitta bound " << bound << ", length " << length << ", error " << error);
        CHECK(length < exact); // an inscribed polyline is always short
        CHECK(error < previous_error);

        // The length deficit of a circle sampled to a sagitta bound s is
        // linear in s, not quadratic: a chord of half-angle theta has sagitta
        // r(1 - cos theta) ~ r theta^2 / 2 and loses r theta^3 / 3 of length,
        // and there are pi / theta of them, giving a total deficit of about
        // 2 pi s / 3 ~ 2.1 s. Bisection stops as soon as the bound is met, so
        // the realised sagitta is usually under it and the true error lands
        // below that estimate.
        CHECK(error < 3.0 * bound);
        previous_error = error;
    }
}

TEST_CASE("a closed loop samples without a duplicated vertex", "[trim][sampling]") {
    const NurbsSurface surface = n2s::cases::plane();
    const TrimLoop loop = n2s::cases::circle_loop({0.5, 0.5}, 0.25);

    const std::vector<Eigen::Vector2d> points = n2s::sample_loop(loop, surface);
    REQUIRE(points.size() > 8);

    // Consecutive vertices, wrapping, are all distinct: a repeated point would
    // become a zero-length constraint edge in the triangulation.
    for (std::size_t i = 0; i < points.size(); ++i) {
        const Eigen::Vector2d& a = points[i];
        const Eigen::Vector2d& b = points[(i + 1) % points.size()];
        INFO("vertex " << i);
        CHECK((a - b).norm() > 1e-12);
    }
}

TEST_CASE("nonsensical sampling tolerances are rejected", "[trim][sampling]") {
    const NurbsSurface surface = n2s::cases::plane();
    const NurbsCurve2 curve = n2s::cases::rectangle_loop({0, 0}, {1, 1}).curves().front();

    SamplingOptions options;
    options.max_segment_length = 0.0;
    CHECK_THROWS_AS(n2s::sample_parameters(curve, surface, options), std::invalid_argument);

    options.max_segment_length = 0.1;
    options.max_sagitta = -1.0;
    CHECK_THROWS_AS(n2s::sample_parameters(curve, surface, options), std::invalid_argument);
}
