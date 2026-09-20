#include "n2s/io/case_json.hpp"
#include "n2s/trim/cases.hpp"
#include "n2s/trim/triangulation.hpp"
#include "n2s/trim/validate.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <filesystem>
#include <stdexcept>
#include <string>

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using n2s::ControlMesh;
using n2s::TrimLoop;
using n2s::TrimRegion;
using n2s::io::Case;

namespace {

/// Scratch path that removes itself, so a failing assertion cannot leave a
/// stale file for the next run to read.
class ScratchFile {
public:
    explicit ScratchFile(const std::string& name)
        : path_(std::filesystem::temp_directory_path() / ("n2s_test_" + name)) {
        std::filesystem::remove(path_);
    }

    ~ScratchFile() { std::filesystem::remove(path_); }

    ScratchFile(const ScratchFile&) = delete;
    ScratchFile& operator=(const ScratchFile&) = delete;

    const std::filesystem::path& path() const { return path_; }

private:
    std::filesystem::path path_;
};

Case sample_case() {
    n2s::cases::TrimmedCase generated =
        n2s::cases::square_with_circular_hole(n2s::cases::saddle(1.0, 0.6), 0.25);
    REQUIRE(n2s::validate_and_repair(generated.region).ok());

    return Case{
        .name = generated.name,
        .surface = generated.surface,
        .region = std::move(generated.region),
        .layout = std::nullopt,
        .control_mesh = ControlMesh::grid(4, 4),
        .camera = n2s::io::CameraPose{{4.0, -3.0, 2.5}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}},
        .color_ranges = {{"error", {0.0, 0.01}}, {"mean_curvature", {-2.0, 2.0}}},
    };
}

void check_same_surface(const n2s::NurbsSurface& a, const n2s::NurbsSurface& b) {
    REQUIRE(a.num_u() == b.num_u());
    REQUIRE(a.num_v() == b.num_v());
    for (int i = 0; i <= 8; ++i) {
        for (int j = 0; j <= 8; ++j) {
            const double u = static_cast<double>(i) / 8.0;
            const double v = static_cast<double>(j) / 8.0;
            REQUIRE((a.evaluate(u, v) - b.evaluate(u, v)).norm() < 1e-15);
        }
    }
}

void check_same_loop(const TrimLoop& a, const TrimLoop& b) {
    REQUIRE(a.size() == b.size());
    for (std::size_t c = 0; c < a.size(); ++c) {
        const n2s::NurbsCurve2& first = a.curves()[c];
        const n2s::NurbsCurve2& second = b.curves()[c];
        for (int i = 0; i <= 16; ++i) {
            const double t = first.knots().domain_start() +
                             (first.knots().domain_end() - first.knots().domain_start()) *
                                 static_cast<double>(i) / 16.0;
            REQUIRE((first.evaluate(t) - second.evaluate(t)).norm() < 1e-15);
        }
    }
}

} // namespace

TEST_CASE("a whole case survives a JSON round trip", "[io][case]") {
    const Case original = sample_case();
    const Case restored = n2s::io::case_from_json(n2s::io::to_json(original));

    CHECK(restored.name == original.name);
    check_same_surface(original.surface, restored.surface);

    check_same_loop(original.region.outer(), restored.region.outer());
    REQUIRE(restored.region.holes().size() == original.region.holes().size());
    for (std::size_t h = 0; h < original.region.holes().size(); ++h) {
        check_same_loop(original.region.holes()[h], restored.region.holes()[h]);
    }

    REQUIRE(restored.control_mesh.has_value());
    CHECK(restored.control_mesh->num_vertices() == original.control_mesh->num_vertices());
    CHECK(restored.control_mesh->num_quads() == original.control_mesh->num_quads());

    REQUIRE(restored.camera.has_value());
    CHECK((restored.camera->eye - original.camera->eye).norm() < 1e-15);
    CHECK((restored.camera->target - original.camera->target).norm() < 1e-15);

    REQUIRE(restored.color_ranges.size() == 2);
    CHECK(restored.color_ranges.at("error").high == Approx(0.01));
    CHECK(restored.color_ranges.at("mean_curvature").low == Approx(-2.0));
}

TEST_CASE("a case survives a file round trip", "[io][case]") {
    const ScratchFile file("case_roundtrip.json");
    const Case original = sample_case();

    n2s::io::write_case(file.path(), original);
    REQUIRE(std::filesystem::exists(file.path()));

    const Case restored = n2s::io::read_case(file.path());
    check_same_surface(original.surface, restored.surface);
    CHECK(restored.name == original.name);
}

TEST_CASE("the optional parts really are optional", "[io][case]") {
    n2s::cases::TrimmedCase generated = n2s::cases::quarter_arc_corner(n2s::cases::plane(), 0.5);
    REQUIRE(n2s::validate_and_repair(generated.region).ok());

    const Case minimal{
        .name = "minimal",
        .surface = generated.surface,
        .region = std::move(generated.region),
        .layout = std::nullopt,
        .control_mesh = std::nullopt,
        .camera = std::nullopt,
        .color_ranges = {},
    };

    const nlohmann::json document = n2s::io::to_json(minimal);
    CHECK(document.find("control_mesh") == document.end());
    CHECK(document.find("camera") == document.end());
    CHECK(document.find("color_ranges") == document.end());

    const Case restored = n2s::io::case_from_json(document);
    CHECK_FALSE(restored.control_mesh.has_value());
    CHECK_FALSE(restored.camera.has_value());
    CHECK(restored.color_ranges.empty());
    CHECK(restored.region.holes().empty());
}

TEST_CASE("crease and corner tags survive the round trip", "[io][case]") {
    ControlMesh mesh = ControlMesh::grid(3, 3);
    mesh.set_crease(0, 3, 2.5);
    mesh.set_corner(0, 10.0);

    n2s::cases::TrimmedCase generated = n2s::cases::l_shape(n2s::cases::plane(), {0.5, 0.5});
    REQUIRE(n2s::validate_and_repair(generated.region).ok());

    const Case original{
        .name = "tagged",
        .surface = generated.surface,
        .region = std::move(generated.region),
        .layout = std::nullopt,
        .control_mesh = std::move(mesh),
        .camera = std::nullopt,
        .color_ranges = {},
    };

    const Case restored = n2s::io::case_from_json(n2s::io::to_json(original));
    REQUIRE(restored.control_mesh.has_value());

    REQUIRE(restored.control_mesh->creases().size() == 1);
    CHECK(restored.control_mesh->creases()[0].sharpness == Approx(2.5));
    REQUIRE(restored.control_mesh->corners().size() == 1);
    CHECK(restored.control_mesh->corners()[0].vertex == 0);
    CHECK(restored.control_mesh->corners()[0].sharpness == Approx(10.0));
}

TEST_CASE("an empty or inverted colour range is rejected", "[io][case]") {
    // Defect 6 of the thesis was curvature plots with no fixed colour scale. A
    // range that is empty or backwards is the same failure wearing a disguise:
    // it produces a picture that looks authoritative and compares nothing.
    nlohmann::json document = n2s::io::to_json(sample_case());

    SECTION("inverted") {
        document["color_ranges"]["error"] = nlohmann::json::array({1.0, 0.0});
        CHECK_THROWS_WITH(n2s::io::case_from_json(document),
                          ContainsSubstring("empty or inverted"));
    }

    SECTION("degenerate") {
        document["color_ranges"]["error"] = nlohmann::json::array({0.5, 0.5});
        CHECK_THROWS_WITH(n2s::io::case_from_json(document),
                          ContainsSubstring("empty or inverted"));
    }

    SECTION("not a pair") {
        document["color_ranges"]["error"] = nlohmann::json::array({0.0});
        CHECK_THROWS_WITH(n2s::io::case_from_json(document), ContainsSubstring("low, high"));
    }
}

TEST_CASE("malformed case documents name the offending field", "[io][case]") {
    SECTION("wrong format tag") {
        nlohmann::json document = n2s::io::to_json(sample_case());
        document["format"] = "n2s-surface";
        CHECK_THROWS_WITH(n2s::io::case_from_json(document), ContainsSubstring("n2s-case"));
    }

    SECTION("unknown version") {
        nlohmann::json document = n2s::io::to_json(sample_case());
        document["version"] = 42;
        CHECK_THROWS_WITH(n2s::io::case_from_json(document), ContainsSubstring("version"));
    }

    SECTION("missing trim") {
        nlohmann::json document = n2s::io::to_json(sample_case());
        document.erase("trim");
        CHECK_THROWS_WITH(n2s::io::case_from_json(document), ContainsSubstring("trim"));
    }

    SECTION("a model-space curve where a domain curve belongs") {
        // The separate "n2s-curve2" tag exists precisely so this fails here
        // rather than producing a trim loop in the wrong space.
        nlohmann::json document = n2s::io::to_json(sample_case());
        document["trim"]["outer"][0]["format"] = "n2s-curve";
        CHECK_THROWS_WITH(n2s::io::case_from_json(document), ContainsSubstring("n2s-curve2"));
    }

    SECTION("the failing curve is identified by index") {
        nlohmann::json document = n2s::io::to_json(sample_case());
        document["trim"]["outer"][2].erase("knots");
        CHECK_THROWS_WITH(n2s::io::case_from_json(document), ContainsSubstring("curve 2"));
    }
}

TEST_CASE("file errors name the path", "[io][case]") {
    const std::filesystem::path missing =
        std::filesystem::temp_directory_path() / "n2s_test_no_such_case.json";
    std::filesystem::remove(missing);

    CHECK_THROWS_WITH(n2s::io::read_case(missing), ContainsSubstring("no_such_case"));
}

TEST_CASE("a round-tripped case still triangulates identically", "[io][case]") {
    // The round trip has to preserve the geometry well enough that the whole
    // downstream pipeline is unchanged, not merely well enough that the numbers
    // look similar.
    const Case original = sample_case();
    const Case restored = n2s::io::case_from_json(n2s::io::to_json(original));

    const n2s::DomainMesh from_original = n2s::triangulate(original.region, original.surface);
    const n2s::DomainMesh from_restored = n2s::triangulate(restored.region, restored.surface);

    REQUIRE(from_original.vertices.size() == from_restored.vertices.size());
    CHECK(from_original.triangles.size() == from_restored.triangles.size());
    CHECK(from_original.area() == Approx(from_restored.area()).epsilon(1e-14));
}
