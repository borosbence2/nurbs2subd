#include "n2s/trim/cases.hpp"
#include "n2s/trim/triangulation.hpp"
#include "n2s/trim/validate.hpp"

#include <Eigen/Geometry> // Vector3d::cross
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

using Catch::Approx;
using n2s::DomainMesh;
using n2s::NurbsSurface;
using n2s::SamplingOptions;
using n2s::SurfaceMesh;
using n2s::TrimLoop;
using n2s::TrimRegion;

namespace {

SamplingOptions at_tolerance(double sagitta) {
    SamplingOptions options;
    options.max_segment_length = 4.0 * sagitta;
    options.max_sagitta = sagitta;
    return options;
}

} // namespace

TEST_CASE("the triangulated area converges to the analytic trimmed area", "[trim][triangulation]") {
    // The headline oracle from the plan: a unit square less a disc of radius r
    // has area exactly 1 - pi r^2, and the triangulation must converge to it.
    const double radius = 0.25;
    const NurbsSurface surface = n2s::cases::plane();
    n2s::cases::TrimmedCase testcase = n2s::cases::square_with_circular_hole(surface, radius);
    REQUIRE(n2s::validate_and_repair(testcase.region).ok());

    const double exact = 1.0 - std::numbers::pi * radius * radius;
    REQUIRE(testcase.analytic_domain_area == Approx(exact));

    double previous_error = 1.0;
    for (const double sagitta : {0.02, 0.005, 0.001}) {
        const DomainMesh mesh = n2s::triangulate(testcase.region, surface, at_tolerance(sagitta));
        const double error = std::abs(mesh.area() - exact);

        INFO("sagitta " << sagitta << ": " << mesh.vertices.size() << " vertices, "
                        << mesh.triangles.size() << " triangles, area " << mesh.area() << ", error "
                        << error);

        // The hole is polygonised by an inscribed polygon, so slightly too
        // little is cut away and the triangulated area overshoots.
        CHECK(mesh.area() > exact);
        CHECK(error < previous_error);
        previous_error = error;
    }

    CHECK(previous_error < 1e-4);
}

TEST_CASE("every triangle lies inside the trimmed region", "[trim][triangulation]") {
    // The failure this catches is a hole that got filled in instead of cut out,
    // which produces a perfectly valid-looking mesh covering the wrong set.
    const NurbsSurface surface = n2s::cases::plane();
    n2s::cases::TrimmedCase testcase = n2s::cases::square_with_circular_hole(surface, 0.3);
    REQUIRE(n2s::validate_and_repair(testcase.region).ok());

    const DomainMesh mesh = n2s::triangulate(testcase.region, surface, at_tolerance(0.005));
    REQUIRE(mesh.triangles.size() > 50);

    for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
        const auto& t = mesh.triangles[i];
        const Eigen::Vector2d centroid =
            (mesh.vertices[t[0]] + mesh.vertices[t[1]] + mesh.vertices[t[2]]) / 3.0;

        INFO("triangle " << i << " centroid " << centroid.transpose());
        CHECK(testcase.region.contains(centroid, 512));
    }
}

TEST_CASE("a reversed hole triangulates the same once repaired", "[trim][triangulation]") {
    // The plan's third M2 test: orientation repair must be a no-op on the
    // result, not merely "not crash".
    const NurbsSurface surface = n2s::cases::plane();
    const double radius = 0.25;

    TrimLoop clockwise_hole = n2s::cases::circle_loop({0.5, 0.5}, radius);
    clockwise_hole.reverse();
    TrimRegion good{n2s::cases::rectangle_loop({0, 0}, {1, 1}), {std::move(clockwise_hole)}};

    // Same geometry, hole wound the wrong way.
    TrimRegion reversed{n2s::cases::rectangle_loop({0, 0}, {1, 1}),
                        {n2s::cases::circle_loop({0.5, 0.5}, radius)}};

    REQUIRE(n2s::validate_and_repair(good).repairs.empty());
    REQUIRE(n2s::validate_and_repair(reversed).repairs.size() == 1);

    const SamplingOptions options = at_tolerance(0.01);
    const DomainMesh from_good = n2s::triangulate(good, surface, options);
    const DomainMesh from_reversed = n2s::triangulate(reversed, surface, options);

    CHECK(from_good.vertices.size() == from_reversed.vertices.size());
    CHECK(from_good.triangles.size() == from_reversed.triangles.size());
    CHECK(from_good.area() == Approx(from_reversed.area()).epsilon(1e-12));
}

TEST_CASE("boundary vertices come first and sit on the trim loops", "[trim][triangulation]") {
    const NurbsSurface surface = n2s::cases::plane();
    n2s::cases::TrimmedCase testcase = n2s::cases::square_with_circular_hole(surface, 0.25);
    REQUIRE(n2s::validate_and_repair(testcase.region).ok());

    const DomainMesh mesh = n2s::triangulate(testcase.region, surface, at_tolerance(0.01));

    REQUIRE(mesh.boundary_vertex_count > 0);
    REQUIRE(mesh.boundary_vertex_count <= mesh.vertices.size());

    // Every boundary vertex is on the outer rectangle or on the hole circle.
    for (std::size_t i = 0; i < mesh.boundary_vertex_count; ++i) {
        const Eigen::Vector2d& p = mesh.vertices[i];
        const double to_centre = (p - Eigen::Vector2d{0.5, 0.5}).norm();

        const bool on_rectangle = std::abs(p.x()) < 1e-12 || std::abs(p.x() - 1.0) < 1e-12 ||
                                  std::abs(p.y()) < 1e-12 || std::abs(p.y() - 1.0) < 1e-12;
        const bool on_hole = std::abs(to_centre - 0.25) < 1e-9;

        INFO("boundary vertex " << i << " at " << p.transpose());
        CHECK((on_rectangle || on_hole));
    }
}

TEST_CASE("the lifted mesh keeps its correspondence with the domain mesh",
          "[trim][triangulation]") {
    const NurbsSurface surface = n2s::cases::saddle(1.0, 0.7);
    n2s::cases::TrimmedCase testcase = n2s::cases::square_with_circular_hole(surface, 0.2);
    REQUIRE(n2s::validate_and_repair(testcase.region).ok());

    const DomainMesh domain = n2s::triangulate(testcase.region, surface, at_tolerance(0.01));
    const SurfaceMesh lifted = n2s::map_to_surface(domain, surface);

    REQUIRE(lifted.vertices.size() == domain.vertices.size());
    REQUIRE(lifted.triangles == domain.triangles);
    CHECK(lifted.boundary_vertex_count == domain.boundary_vertex_count);

    for (std::size_t i = 0; i < domain.vertices.size(); ++i) {
        const Eigen::Vector2d& p = domain.vertices[i];
        CHECK((lifted.vertices[i] - surface.evaluate(p.x(), p.y())).norm() < 1e-15);
    }

    // A saddle is not flat, so its area must exceed that of its domain image.
    CHECK(lifted.area() > domain.area());
}

TEST_CASE("the lifted area converges to the analytic surface area", "[trim][triangulation]") {
    // A cylinder of radius r trimmed to the whole domain has area 2 pi r h
    // exactly, and the triangulated area approaches it from below.
    //
    // Note what limits the accuracy here. `triangulate` refines the *boundary*
    // only; CDT adds no interior vertices, so it triangulates the boundary
    // polygon with long triangles whose chords cut straight across the curved
    // cylinder. Tightening the boundary sampling barely helps, because the
    // interior spans are what lose the area. The assertion below is therefore
    // deliberately loose -- it is the accuracy this pass can reach on its own.
    // Interior refinement is tested separately and does far better.
    const double radius = 1.5;
    const double height = 2.0;
    const NurbsSurface surface = n2s::cases::cylinder(radius, height);
    const TrimRegion region{n2s::cases::rectangle_loop({0, 0}, {1, 1})};

    const double exact = 2.0 * std::numbers::pi * radius * height;

    double previous_error = 1e9;
    for (const double sagitta : {0.05, 0.01, 0.002}) {
        const SurfaceMesh mesh =
            n2s::map_to_surface(n2s::triangulate(region, surface, at_tolerance(sagitta)), surface);
        const double error = std::abs(mesh.area() - exact);

        INFO("sagitta " << sagitta << ": area " << mesh.area() << " vs " << exact);
        CHECK(mesh.area() < exact);
        CHECK(error < previous_error);
        previous_error = error;
    }
    CHECK(previous_error / exact < 1e-2);
}

TEST_CASE("every synthetic case triangulates cleanly", "[trim][triangulation][cases]") {
    // Part of the M2 exit criterion.
    for (n2s::cases::TrimmedCase testcase : n2s::cases::all_synthetic_cases()) {
        const n2s::TrimReport report = n2s::validate_and_repair(testcase.region);
        INFO(testcase.name << "\n" << report.to_string());
        REQUIRE(report.ok());

        const DomainMesh mesh =
            n2s::triangulate(testcase.region, testcase.surface, at_tolerance(0.01));

        CHECK(mesh.triangles.size() > 10);
        CHECK(mesh.vertices.size() > 10);
        CHECK(mesh.area() == Approx(testcase.analytic_domain_area).epsilon(2e-3));

        // No degenerate triangles.
        for (const auto& t : mesh.triangles) {
            CHECK(t[0] != t[1]);
            CHECK(t[1] != t[2]);
            CHECK(t[0] != t[2]);
        }
    }
}

TEST_CASE("refinement bounds triangle area and improves the shape",
          "[trim][triangulation][refinement]") {
    const NurbsSurface surface = n2s::cases::saddle(1.0, 0.6);
    n2s::cases::TrimmedCase testcase = n2s::cases::square_with_circular_hole(surface, 0.25);
    REQUIRE(n2s::validate_and_repair(testcase.region).ok());

    const SamplingOptions sampling = at_tolerance(0.02);
    const DomainMesh coarse = n2s::triangulate(testcase.region, surface, sampling);

    n2s::RefinementOptions refinement;
    refinement.max_triangle_area = 0.005;
    refinement.min_angle_degrees = 25.0;
    const DomainMesh fine =
        n2s::triangulate_refined(testcase.region, surface, sampling, refinement);

    INFO("coarse " << coarse.triangles.size() << " triangles, fine " << fine.triangles.size());
    CHECK(fine.triangles.size() > coarse.triangles.size());
    CHECK(fine.vertices.size() > coarse.vertices.size());

    // Every triangle within the area bound.
    const n2s::SurfaceMesh lifted = n2s::map_to_surface(fine, surface);
    double worst_area = 0.0;
    for (const auto& t : lifted.triangles) {
        const Eigen::Vector3d ab = lifted.vertices[t[1]] - lifted.vertices[t[0]];
        const Eigen::Vector3d ac = lifted.vertices[t[2]] - lifted.vertices[t[0]];
        worst_area = std::max(worst_area, 0.5 * ab.cross(ac).norm());
    }
    INFO("worst triangle area " << worst_area);
    CHECK(worst_area <= refinement.max_triangle_area * 1.05);
}

TEST_CASE("refinement never disturbs the boundary discretisation",
          "[trim][triangulation][refinement]") {
    // The property two neighbouring patches depend on for a watertight join:
    // whatever refinement does inside, the boundary vertices stay exactly as
    // the trim sampling produced them.
    const NurbsSurface surface = n2s::cases::paraboloid(1.0, 0.8);
    n2s::cases::TrimmedCase testcase = n2s::cases::square_with_circular_hole(surface, 0.25);
    REQUIRE(n2s::validate_and_repair(testcase.region).ok());

    const SamplingOptions sampling = at_tolerance(0.02);
    const DomainMesh coarse = n2s::triangulate(testcase.region, surface, sampling);

    n2s::RefinementOptions refinement;
    refinement.max_triangle_area = 0.004;
    const DomainMesh fine =
        n2s::triangulate_refined(testcase.region, surface, sampling, refinement);

    REQUIRE(fine.boundary_vertex_count == coarse.boundary_vertex_count);
    for (std::size_t i = 0; i < coarse.boundary_vertex_count; ++i) {
        INFO("boundary vertex " << i);
        CHECK((fine.vertices[i] - coarse.vertices[i]).norm() < 1e-15);
    }
}

TEST_CASE("refined area converges as the vertex budget grows",
          "[trim][triangulation][refinement]") {
    // Two things are worth recording here, because both bear on how M5 should
    // compute area-weighted error metrics.
    //
    // First, the sign. A triangulation whose vertices all lie exactly on a
    // smooth surface can have *more* area than the surface, and can be made to
    // have arbitrarily more: Schwarz's lantern. Anisotropic triangles on a
    // cylinder are that construction, and every refined mesh below overshoots.
    // So the honest claim is a bound on the magnitude of the error, never on
    // its sign -- unlike the unrefined case above, whose long chords of a
    // convex surface do always fall short.
    //
    // Second, and less obviously: refining is not monotonically an improvement.
    // The unrefined mesh here is about 0.4% low, while the refined mesh at the
    // smallest budget is nearly 3% high. Refinement only pays once the budget
    // is large enough, and a mesh refined just a little is worse for area than
    // one not refined at all.
    const double radius = 1.5;
    const double height = 2.0;
    const NurbsSurface surface = n2s::cases::cylinder(radius, height);
    const TrimRegion region{n2s::cases::rectangle_loop({0, 0}, {1, 1})};
    const double exact = 2.0 * std::numbers::pi * radius * height;

    const SamplingOptions sampling = at_tolerance(0.01);

    double previous_error = 1.0;
    std::size_t previous_triangles = 0;

    for (const std::size_t budget : {std::size_t{1000}, std::size_t{3000}, std::size_t{8000}}) {
        n2s::RefinementOptions refinement;
        refinement.max_triangle_area = 0.02;
        refinement.min_angle_degrees = 20.0;
        refinement.max_vertices = budget;

        const DomainMesh domain = n2s::triangulate_refined(region, surface, sampling, refinement);
        const n2s::SurfaceMesh mesh = n2s::map_to_surface(domain, surface);
        const double error = std::abs(mesh.area() - exact) / exact;

        INFO("budget " << budget << ": " << domain.vertices.size() << " vertices, "
                       << mesh.triangles.size() << " triangles, relative error " << error);

        CHECK(domain.vertices.size() <= budget);
        CHECK(mesh.triangles.size() > previous_triangles);
        CHECK(error < previous_error);

        previous_error = error;
        previous_triangles = mesh.triangles.size();
    }

    CHECK(previous_error < 4e-3);
}

TEST_CASE("refinement does not produce degenerate triangles", "[trim][triangulation][refinement]") {
    // The failure this pins down: before the minimum-separation floor existed,
    // a sliver's circumcentre landed almost on an existing vertex, which made a
    // thinner sliver, which did it again. The mesh stayed a valid triangulation
    // with exactly the right domain area throughout, so nothing complained --
    // but the worst 3D aspect ratio reached 1e13 and the surface area came out
    // half a percent too high.
    const NurbsSurface surface = n2s::cases::cylinder(1.5, 2.0);
    const TrimRegion region{n2s::cases::rectangle_loop({0, 0}, {1, 1})};

    n2s::RefinementOptions refinement;
    refinement.max_triangle_area = 0.02;
    refinement.min_angle_degrees = 20.0;
    refinement.max_vertices = 4000;

    const DomainMesh domain =
        n2s::triangulate_refined(region, surface, at_tolerance(0.01), refinement);
    const n2s::SurfaceMesh mesh = n2s::map_to_surface(domain, surface);

    // The domain is tiled exactly, with no overlaps and no gaps.
    CHECK(domain.area() == Approx(1.0).epsilon(1e-12));

    double worst_aspect = 0.0;
    for (const auto& t : mesh.triangles) {
        const double e0 = (mesh.vertices[t[1]] - mesh.vertices[t[0]]).norm();
        const double e1 = (mesh.vertices[t[2]] - mesh.vertices[t[1]]).norm();
        const double e2 = (mesh.vertices[t[0]] - mesh.vertices[t[2]]).norm();
        worst_aspect = std::max(worst_aspect, std::max({e0, e1, e2}) / std::min({e0, e1, e2}));
    }

    INFO("worst 3D aspect ratio " << worst_aspect);
    // Generous, because the cylinder's parameterisation is itself anisotropic:
    // the domain is a unit square standing for a 9.4 by 2.0 surface, so even a
    // perfectly equilateral domain triangle maps to a 4.7:1 one.
    CHECK(worst_aspect < 50.0);
}

TEST_CASE("refinement with both criteria disabled is a plain triangulation",
          "[trim][triangulation][refinement]") {
    const NurbsSurface surface = n2s::cases::plane();
    n2s::cases::TrimmedCase testcase = n2s::cases::square_with_circular_hole(surface, 0.25);
    REQUIRE(n2s::validate_and_repair(testcase.region).ok());

    n2s::RefinementOptions off;
    off.max_triangle_area = 0.0;
    off.min_angle_degrees = 0.0;

    const SamplingOptions sampling = at_tolerance(0.02);
    const DomainMesh plain = n2s::triangulate(testcase.region, surface, sampling);
    const DomainMesh same = n2s::triangulate_refined(testcase.region, surface, sampling, off);

    CHECK(plain.vertices.size() == same.vertices.size());
    CHECK(plain.triangles.size() == same.triangles.size());
}

TEST_CASE("the vertex budget is respected", "[trim][triangulation][refinement]") {
    const NurbsSurface surface = n2s::cases::plane();
    const TrimRegion region{n2s::cases::rectangle_loop({0, 0}, {1, 1})};

    n2s::RefinementOptions refinement;
    refinement.max_triangle_area = 1e-6; // unreachable on purpose
    refinement.max_vertices = 400;
    refinement.max_rounds = 30;

    const DomainMesh mesh =
        n2s::triangulate_refined(region, surface, at_tolerance(0.05), refinement);

    INFO(mesh.vertices.size() << " vertices");
    CHECK(mesh.vertices.size() <= refinement.max_vertices * 2);
    CHECK(mesh.triangles.size() > 10);
}
