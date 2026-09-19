#include "n2s/metrics/seam_error.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

using Catch::Approx;
using n2s::ControlMesh;
using n2s::SubdivisionSurface;
using n2s::metrics::FaceEdge;
using n2s::metrics::SeamEdge;
using n2s::metrics::SeamError;

namespace {

constexpr int kSize = 6; // vertices per side
constexpr int kFaces = kSize - 1;

double height(int i, int j) {
    return 0.25 * std::sin(0.8 * i) + 0.18 * std::cos(1.1 * j);
}

/// A grid patch spanning x in [x0, x0 + 1], y in [0, 1].
///
/// The heights depend on the *global* column index, so two patches built with
/// adjacent offsets agree exactly along the column they share -- which is what
/// makes them candidates for a watertight join.
ControlMesh patch(int column_offset, double dz = 0.0) {
    ControlMesh mesh = ControlMesh::grid(kSize, kSize);
    std::vector<Eigen::Vector3d> vertices = mesh.vertices();

    for (int i = 0; i < kSize; ++i) {
        for (int j = 0; j < kSize; ++j) {
            const int global_i = column_offset + i;
            Eigen::Vector3d& v = vertices[static_cast<std::size_t>(i * kSize + j)];
            v.x() = static_cast<double>(global_i) / static_cast<double>(kSize - 1);
            v.y() = static_cast<double>(j) / static_cast<double>(kSize - 1);
            v.z() = height(global_i, j) + dz;
        }
    }

    mesh.set_vertices(std::move(vertices));
    return mesh;
}

/// The seam runs along the shared column: patch A's u = 1 edge against patch
/// B's u = 0 edge, face row by face row.
std::vector<std::pair<SeamEdge, SeamEdge>> seam_pairs() {
    std::vector<std::pair<SeamEdge, SeamEdge>> pairs;
    for (int j = 0; j < kFaces; ++j) {
        const SeamEdge left{(kFaces - 1) * kFaces + j, FaceEdge::UMax, false};
        const SeamEdge right{0 * kFaces + j, FaceEdge::UMin, false};
        pairs.emplace_back(left, right);
    }
    return pairs;
}

} // namespace

TEST_CASE("seam_parameter walks the named edge", "[metrics][seam]") {
    CHECK((n2s::metrics::seam_parameter({0, FaceEdge::VMin, false}, 0.25) -
           Eigen::Vector2d{0.25, 0.0})
              .norm() < 1e-15);
    CHECK((n2s::metrics::seam_parameter({0, FaceEdge::UMax, false}, 0.25) -
           Eigen::Vector2d{1.0, 0.25})
              .norm() < 1e-15);
    CHECK((n2s::metrics::seam_parameter({0, FaceEdge::VMax, false}, 0.25) -
           Eigen::Vector2d{0.25, 1.0})
              .norm() < 1e-15);
    CHECK((n2s::metrics::seam_parameter({0, FaceEdge::UMin, false}, 0.25) -
           Eigen::Vector2d{0.0, 0.25})
              .norm() < 1e-15);

    // Reversed walks the same edge the other way.
    CHECK(
        (n2s::metrics::seam_parameter({0, FaceEdge::VMin, true}, 0.25) - Eigen::Vector2d{0.75, 0.0})
            .norm() < 1e-15);
}

TEST_CASE("patches sharing boundary control points join with no gap at all",
          "[metrics][seam][watertight]") {
    // The property the whole watertightness argument rests on, measured rather
    // than assumed. M3 established it on the stencil weights: the limit
    // boundary curve is a function of the boundary control points alone. Two
    // patches whose shared column of control points is identical therefore
    // produce the *same* boundary curve, and the gap is not small -- it is
    // zero to the last bit the arithmetic carries.
    //
    // This is the metric R3 will drive, so it is worth knowing now that it can
    // actually read zero rather than bottoming out at some tolerance.
    const SubdivisionSurface left{patch(0)};
    const SubdivisionSurface right{patch(kSize - 1)};

    double worst_gap = 0.0;
    for (const auto& [edge_left, edge_right] : seam_pairs()) {
        const SeamError error = n2s::metrics::measure_seam(left, edge_left, right, edge_right);

        INFO("faces " << edge_left.face << " and " << edge_right.face << ", gap max "
                      << error.gap.max);
        REQUIRE(error.gap.count > 0);
        CHECK(error.gap.max < 1e-14);
        worst_gap = std::max(worst_gap, error.gap.max);
    }

    INFO("worst gap across the whole seam " << worst_gap);
    CHECK(worst_gap < 1e-14);
}

TEST_CASE("the seam normals do not agree just because the gap is zero",
          "[metrics][seam][watertight]") {
    // Shared boundary control points give position continuity and nothing
    // more. Recording this now matters: a zero gap is easy to mistake for a
    // smooth join, and closing the normal angle is a separate piece of work
    // that the plan treats as optional.
    const SubdivisionSurface left{patch(0)};
    const SubdivisionSurface right{patch(kSize - 1)};

    const auto pairs = seam_pairs();
    const SeamError error =
        n2s::metrics::measure_seam(left, pairs[1].first, right, pairs[1].second);

    CHECK(error.gap.max < 1e-14);
    INFO("normal deviation across the seam: max " << error.normal_degrees.max << " degrees");
    CHECK(error.normal_degrees.count > 0);
    CHECK(error.normal_degrees.max > 1e-6); // G0, not G1
    CHECK(error.normal_degrees.max < 90.0);
}

TEST_CASE("a displaced patch reports the displacement as the gap", "[metrics][seam]") {
    const double displacement = 0.05;
    const SubdivisionSurface left{patch(0)};
    const SubdivisionSurface right{patch(kSize - 1, displacement)};

    const auto pairs = seam_pairs();
    const SeamError error =
        n2s::metrics::measure_seam(left, pairs[2].first, right, pairs[2].second);

    INFO("gap max " << error.gap.max << ", mean " << error.gap.mean);
    CHECK(error.gap.max == Approx(displacement).epsilon(1e-9));
    CHECK(error.gap.mean == Approx(displacement).epsilon(1e-9));

    // A rigid shift in z leaves the two surfaces parallel, so the normal
    // deviation is unchanged by it.
    const SeamError undisplaced = n2s::metrics::measure_seam(
        left, pairs[2].first, SubdivisionSurface{patch(kSize - 1)}, pairs[2].second);
    CHECK(error.normal_degrees.max == Approx(undisplaced.normal_degrees.max).epsilon(1e-9));
}

TEST_CASE("a patch measured against itself has a zero seam", "[metrics][seam]") {
    const SubdivisionSurface surface{patch(0)};
    const SeamEdge edge{3, FaceEdge::UMax, false};

    const SeamError error = n2s::metrics::measure_seam(surface, edge, surface, edge);
    CHECK(error.gap.max == 0.0);
    CHECK(error.normal_degrees.max < 1e-12);
}

TEST_CASE("comparing a seam end-to-end instead of point-for-point is caught", "[metrics][seam]") {
    // Two edges tracing the same curve in opposite directions have the same
    // point set but do not correspond. The reversed flag exists for that, and
    // without it the gap is large even though the curves coincide -- which is
    // the correct answer, because the two patches would not refine compatibly.
    const SubdivisionSurface left{patch(0)};
    const SubdivisionSurface right{patch(kSize - 1)};

    const auto pairs = seam_pairs();
    SeamEdge flipped = pairs[0].second;
    flipped.reversed = true;

    const SeamError misaligned = n2s::metrics::measure_seam(left, pairs[0].first, right, flipped);
    CHECK(misaligned.gap.max > 1e-3);
}

TEST_CASE("a nonsensical sample count is rejected", "[metrics][seam]") {
    const SubdivisionSurface surface{patch(0)};
    CHECK_THROWS_AS(
        n2s::metrics::measure_seam(
            surface, {0, FaceEdge::VMin, false}, surface, {0, FaceEdge::VMin, false}, 0),
        std::invalid_argument);
}
