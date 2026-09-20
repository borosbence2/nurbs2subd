#include "n2s/metrics/surface_error.hpp"
#include "n2s/trim/cases.hpp"
#include "n2s/trim/validate.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using Catch::Approx;
using n2s::ControlMesh;
using n2s::NurbsSurface;
using n2s::SubdivisionSurface;
using n2s::TrimRegion;
using n2s::metrics::Stats;
using n2s::metrics::SurfaceError;
using n2s::metrics::SurfaceErrorOptions;

namespace {

/// A grid control mesh lying exactly on the unit plane, so that its limit
/// surface *is* that plane. Optionally lifted by `height`, which moves the
/// limit surface to a parallel plane and makes every distance analytic.
ControlMesh planar_grid(int size, double height = 0.0) {
    ControlMesh mesh = ControlMesh::grid(size, size);
    if (height == 0.0) {
        return mesh;
    }
    std::vector<Eigen::Vector3d> vertices = mesh.vertices();
    for (Eigen::Vector3d& v : vertices) {
        v.z() += height;
    }
    mesh.set_vertices(std::move(vertices));
    return mesh;
}

SurfaceErrorOptions cheap_options() {
    SurfaceErrorOptions options;
    options.samples_per_face = 2;
    options.nurbs_samples_per_side = 12;
    options.restrict_to_trimmed_region = false;
    return options;
}

} // namespace

// ---------------------------------------------------------------------------
// Statistics.
// ---------------------------------------------------------------------------

TEST_CASE("the accumulator reports max, mean and RMS", "[metrics][statistics]") {
    n2s::metrics::Accumulator accumulator;
    for (const double value : {3.0, 0.0, 4.0}) {
        accumulator.add(value);
    }

    const Stats stats = accumulator.result();
    CHECK(stats.count == 3);
    CHECK(stats.max == Approx(4.0));
    CHECK(stats.mean == Approx(7.0 / 3.0));
    CHECK(stats.rms == Approx(std::sqrt(25.0 / 3.0)));
    CHECK(stats.complete());
}

TEST_CASE("unmeasured samples are excluded, not counted as zero", "[metrics][statistics]") {
    // Folding a failed measurement in as zero flatters the result; folding it
    // in as something large invents data. Neither is acceptable in a number
    // that goes into a paper, so they are excluded and counted.
    n2s::metrics::Accumulator accumulator;
    accumulator.add(2.0);
    accumulator.add_unmeasured();
    accumulator.add(4.0);

    const Stats stats = accumulator.result();
    CHECK(stats.count == 2);
    CHECK(stats.unmeasured == 1);
    CHECK(stats.mean == Approx(3.0));
    CHECK_FALSE(stats.complete());
}

TEST_CASE("a non-finite sample is treated as unmeasured", "[metrics][statistics]") {
    // One NaN admitted into the sums turns every statistic into NaN and hides
    // which sample was at fault.
    n2s::metrics::Accumulator accumulator;
    accumulator.add(1.0);
    accumulator.add(std::numeric_limits<double>::quiet_NaN());
    accumulator.add(std::numeric_limits<double>::infinity());

    const Stats stats = accumulator.result();
    CHECK(stats.count == 1);
    CHECK(stats.unmeasured == 2);
    CHECK(stats.mean == Approx(1.0));
}

TEST_CASE("an empty accumulator reports zeros, not NaN", "[metrics][statistics]") {
    const Stats stats = n2s::metrics::Accumulator{}.result();
    CHECK(stats.count == 0);
    CHECK(stats.max == 0.0);
    CHECK(stats.mean == 0.0);
    CHECK(stats.rms == 0.0);
}

// ---------------------------------------------------------------------------
// The domain correspondence.
// ---------------------------------------------------------------------------

TEST_CASE("the grid domain map covers the unit domain", "[metrics][domain]") {
    const int size = 5;
    const n2s::metrics::DomainMap map = n2s::metrics::grid_domain_map(size, size);

    // Face 0 corner (0,0) is the domain origin; the last face's far corner is
    // (1,1). Between them the map is affine per face.
    CHECK((map(n2s::LimitLocation{0, 0.0, 0.0}) - Eigen::Vector2d{0.0, 0.0}).norm() < 1e-15);

    const int faces = (size - 1) * (size - 1);
    CHECK((map(n2s::LimitLocation{faces - 1, 1.0, 1.0}) - Eigen::Vector2d{1.0, 1.0}).norm() <
          1e-15);

    // Face (i, j) starts at (i/(size-1), j/(size-1)).
    CHECK((map(n2s::LimitLocation{1 * (size - 1) + 2, 0.0, 0.0}) - Eigen::Vector2d{0.25, 0.5})
              .norm() < 1e-15);

    CHECK_THROWS_AS(n2s::metrics::grid_domain_map(1, 5), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Exact oracles: a planar control net whose limit surface is the plane.
// ---------------------------------------------------------------------------

TEST_CASE("a limit surface that equals the NURBS reports zero error", "[metrics][oracle]") {
    // The sharpest oracle available: a planar control net has a planar limit
    // surface, and the plane is exactly representable as a NURBS. Every error
    // must therefore come back at machine precision, and any that does not is
    // measuring something other than what it claims.
    const int size = 6;
    const NurbsSurface plane = n2s::cases::plane(1.0, 1.0);
    const SubdivisionSurface limit{planar_grid(size)};
    const TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0})};

    const SurfaceError error = n2s::metrics::measure_surface_error(
        limit, plane, region, n2s::metrics::grid_domain_map(size, size), cheap_options());

    INFO("geometric max " << error.geometric.max << ", parametric max " << error.parametric.max
                          << ", normal max " << error.normal_degrees.max);

    REQUIRE(error.geometric.count > 0);
    CHECK(error.geometric.max < 1e-12);
    CHECK(error.hausdorff < 1e-9);
    CHECK(error.normal_degrees.max < 1e-9);
    CHECK(error.mean_curvature.max < 1e-9);
    CHECK(error.gaussian_curvature.max < 1e-9);
    CHECK(error.geometric.complete());
}

TEST_CASE("the parametric correspondence is exact on interior faces", "[metrics][oracle]") {
    // A cubic B-spline reproduces linear functions exactly, so on an interior
    // face the limit surface's parameterisation of an evenly spaced planar grid
    // agrees with the affine grid domain map to machine precision.
    //
    // Not so on the boundary faces: the corner rule reparameterises the last
    // span, so the same domain map is slightly off there. That is a real
    // property of the surface rather than a defect of the metric, and it is
    // why the parametric error is reported separately from the geometric one.
    const int size = 6;
    const NurbsSurface plane = n2s::cases::plane(1.0, 1.0);
    const SubdivisionSurface limit{planar_grid(size)};
    const n2s::metrics::DomainMap map = n2s::metrics::grid_domain_map(size, size);

    const int faces_per_side = size - 1;
    double worst_interior = 0.0;
    for (int i = 1; i + 1 < faces_per_side; ++i) {
        for (int j = 1; j + 1 < faces_per_side; ++j) {
            const int face = i * faces_per_side + j;
            for (const double u : {0.0, 0.5, 1.0}) {
                for (const double v : {0.0, 0.5, 1.0}) {
                    const n2s::LimitLocation location{face, u, v};
                    const Eigen::Vector2d domain_point = map(location);
                    const Eigen::Vector3d expected =
                        plane.evaluate(domain_point.x(), domain_point.y());
                    const Eigen::Vector3d actual = limit.evaluate_limit(location).position;
                    worst_interior = std::max(worst_interior, (expected - actual).norm());
                }
            }
        }
    }

    INFO("worst interior parametric deviation " << worst_interior);
    CHECK(worst_interior < 1e-12);
}

TEST_CASE("a translated control mesh reports exactly the translation", "[metrics][oracle]") {
    // Lifting the planar control net by d puts the limit surface on the plane
    // z = d, so every distance to the original plane is exactly d, and the
    // normals stay parallel. Anything other than d means the projection or the
    // accumulation is wrong.
    const int size = 6;
    const double height = 0.125;
    const NurbsSurface plane = n2s::cases::plane(1.0, 1.0);
    const SubdivisionSurface limit{planar_grid(size, height)};
    const TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0})};

    const SurfaceError error = n2s::metrics::measure_surface_error(
        limit, plane, region, n2s::metrics::grid_domain_map(size, size), cheap_options());

    INFO("geometric mean " << error.geometric.mean << " vs " << height);
    CHECK(error.geometric.max == Approx(height).epsilon(1e-9));
    CHECK(error.geometric.mean == Approx(height).epsilon(1e-9));
    CHECK(error.geometric.rms == Approx(height).epsilon(1e-9));
    CHECK(error.parametric.max == Approx(height).epsilon(1e-9));

    // Parallel planes: the tangent planes agree exactly.
    CHECK(error.normal_degrees.max < 1e-9);
}

TEST_CASE("the reverse direction catches a limit surface that does not cover",
          "[metrics][oracle]") {
    // A layout covering only part of the domain is close to the NURBS
    // everywhere it exists, so the one-sided distance stays small. Only the
    // reverse direction notices that half the NURBS has nothing above it --
    // which is exactly why the plan asks for a two-sided measure.
    const NurbsSurface plane = n2s::cases::plane(1.0, 1.0);
    const TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0})};

    ControlMesh half = ControlMesh::grid(5, 5);
    std::vector<Eigen::Vector3d> vertices = half.vertices();
    for (Eigen::Vector3d& v : vertices) {
        v.x() *= 0.5; // squeeze the layout into the left half of the domain
    }
    half.set_vertices(std::move(vertices));

    const SubdivisionSurface limit{half};
    const SurfaceError error =
        n2s::metrics::measure_surface_error(limit, plane, region, {}, cheap_options());

    INFO("forward max " << error.geometric.max << ", reverse max " << error.reverse_geometric.max);
    CHECK(error.geometric.max < 1e-9);        // everything it does cover is on the plane
    CHECK(error.reverse_geometric.max > 0.4); // but it covers only half of it
    CHECK(error.hausdorff == Approx(error.reverse_geometric.max));
}

TEST_CASE("samples outside the trimmed region are excluded", "[metrics][trim]") {
    // The layout covers the whole domain but the trim cuts a hole in it.
    // Measuring error over the hole would compare against surface that was
    // never meant to be approximated.
    const int size = 5;
    const NurbsSurface plane = n2s::cases::plane(1.0, 1.0);
    const SubdivisionSurface limit{planar_grid(size)};

    n2s::TrimLoop hole = n2s::cases::circle_loop({0.5, 0.5}, 0.3);
    hole.reverse();
    TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0}), {std::move(hole)}};
    REQUIRE(n2s::validate_and_repair(region).ok());

    SurfaceErrorOptions options = cheap_options();
    options.restrict_to_trimmed_region = true;

    const SurfaceError error = n2s::metrics::measure_surface_error(
        limit, plane, region, n2s::metrics::grid_domain_map(size, size), options);

    INFO(error.geometric.count << " measured, " << error.geometric.unmeasured << " excluded");
    CHECK(error.geometric.unmeasured > 0);
    CHECK(error.geometric.count > 0);
    CHECK_FALSE(error.geometric.complete());
}

TEST_CASE("without a domain map only the correspondence-free metrics are reported",
          "[metrics][oracle]") {
    // Reported as absent rather than invented from a guessed correspondence.
    const NurbsSurface plane = n2s::cases::plane(1.0, 1.0);
    const SubdivisionSurface limit{planar_grid(5)};
    const TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0})};

    const SurfaceError error =
        n2s::metrics::measure_surface_error(limit, plane, region, {}, cheap_options());

    CHECK(error.geometric.count > 0);
    CHECK(error.parametric.count == 0);
    CHECK(error.normal_degrees.count == 0);
    CHECK(error.mean_curvature.count == 0);
}

TEST_CASE("every sample carries the curvature the summary reports", "[metrics][curvature]") {
    // The summary used to recompute curvature in a second pass of its own.
    // Now it reads what the samples carry, so the two must still agree -- if
    // they drift, a figure drawn from samples.csv stops describing the numbers
    // in metrics.json, and nothing in either would say so.
    const int size = 6;
    const NurbsSurface saddle = n2s::cases::saddle(1.0, 0.5);
    const TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0})};
    const n2s::metrics::DomainMap map = n2s::metrics::grid_domain_map(size, size);

    // A control net that is *not* the surface, so the deviation is non-zero
    // and an accidental agreement at zero cannot pass for a real one.
    const SubdivisionSurface limit{planar_grid(size)};

    const SurfaceError summary =
        n2s::metrics::measure_surface_error(limit, saddle, region, map, cheap_options());
    const std::vector<n2s::metrics::ErrorSample> samples =
        n2s::metrics::sample_surface_error(limit, saddle, region, map, cheap_options());

    double mean_max = 0.0;
    double gaussian_max = 0.0;
    std::size_t measured = 0;
    for (const n2s::metrics::ErrorSample& sample : samples) {
        // Curvature must never outlive the sample carrying it. Curvature is
        // computed before the projection that can fail, so this is the
        // invariant that keeps the curvature statistics over the same set of
        // samples as every statistic printed next to them.
        REQUIRE((!sample.curvature_measured || sample.measured));
        if (!sample.curvature_measured) {
            continue;
        }
        ++measured;
        mean_max = std::max(mean_max, sample.mean_curvature);
        gaussian_max = std::max(gaussian_max, sample.gaussian_curvature);
    }

    CHECK(summary.mean_curvature.count <= summary.geometric.count);

    REQUIRE(measured > 0);
    CHECK(summary.mean_curvature.count == measured);
    CHECK(summary.mean_curvature.max == Approx(mean_max));
    CHECK(summary.gaussian_curvature.max == Approx(gaussian_max));

    // The saddle really is curved, so this is a measurement and not a zero.
    CHECK(mean_max > 1e-6);
}

TEST_CASE("a sample with no correspondence carries no curvature", "[metrics][curvature]") {
    // Curvature is a deviation from the NURBS at the corresponding point.
    // Without a correspondence there is no such point, and reporting zero
    // would read as a perfect curvature match.
    const NurbsSurface saddle = n2s::cases::saddle(1.0, 0.5);
    const SubdivisionSurface limit{planar_grid(5)};
    const TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0})};

    const std::vector<n2s::metrics::ErrorSample> samples =
        n2s::metrics::sample_surface_error(limit, saddle, region, {}, cheap_options());

    REQUIRE_FALSE(samples.empty());
    for (const n2s::metrics::ErrorSample& sample : samples) {
        REQUIRE_FALSE(sample.curvature_measured);
    }
}

TEST_CASE("an unmeasured sample never claims a curvature", "[metrics][curvature]") {
    // A sample dropped for falling outside the trim must not carry curvature
    // either, or a map drawn over the trimmed region would show values from
    // outside it.
    const int size = 5;
    const NurbsSurface saddle = n2s::cases::saddle(1.0, 0.5);
    const SubdivisionSurface limit{planar_grid(size)};

    // A trim covering only a corner, so most samples fall outside it.
    const TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {0.4, 0.4})};

    SurfaceErrorOptions options = cheap_options();
    options.restrict_to_trimmed_region = true;

    const std::vector<n2s::metrics::ErrorSample> samples = n2s::metrics::sample_surface_error(
        limit, saddle, region, n2s::metrics::grid_domain_map(size, size), options);

    std::size_t outside = 0;
    for (const n2s::metrics::ErrorSample& sample : samples) {
        if (!sample.measured) {
            ++outside;
            CHECK_FALSE(sample.curvature_measured);
        }
    }
    REQUIRE(outside > 0);
}

// ---------------------------------------------------------------------------
// Boundary deviation.
// ---------------------------------------------------------------------------

TEST_CASE("a layout ending on the trim reports zero boundary deviation", "[metrics][boundary]") {
    // The planar grid's limit boundary is the unit square's edge, which is
    // exactly the image of the rectangular trim loop.
    const NurbsSurface plane = n2s::cases::plane(1.0, 1.0);
    const SubdivisionSurface limit{planar_grid(6)};
    const TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0})};

    const Stats boundary = n2s::metrics::measure_boundary_error(limit, plane, region);

    INFO("boundary max " << boundary.max << " over " << boundary.count << " samples");
    REQUIRE(boundary.count > 0);
    CHECK(boundary.max < 1e-9);
}

TEST_CASE("a layout that stops short reports the shortfall", "[metrics][boundary]") {
    const NurbsSurface plane = n2s::cases::plane(1.0, 1.0);
    const TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0})};

    // Shrink the layout toward the origin: its boundary now sits well inside
    // the trim curve.
    ControlMesh shrunk = ControlMesh::grid(6, 6);
    std::vector<Eigen::Vector3d> vertices = shrunk.vertices();
    for (Eigen::Vector3d& v : vertices) {
        v.x() = 0.1 + 0.8 * v.x();
        v.y() = 0.1 + 0.8 * v.y();
    }
    shrunk.set_vertices(std::move(vertices));

    const Stats boundary =
        n2s::metrics::measure_boundary_error(SubdivisionSurface{shrunk}, plane, region);

    INFO("boundary max " << boundary.max);
    CHECK(boundary.max == Approx(0.1).epsilon(1e-3));
}

TEST_CASE("a closed control mesh reports no boundary samples", "[metrics][boundary]") {
    // Zero samples, not a zero deviation: a cube has no boundary, and reporting
    // 0.0 would read as a perfect result.
    const Stats boundary = n2s::metrics::measure_boundary_error(
        SubdivisionSurface{ControlMesh::cube()},
        n2s::cases::plane(),
        TrimRegion{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0})});

    CHECK(boundary.count == 0);
}
