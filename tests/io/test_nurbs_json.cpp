#include "n2s/io/nurbs_json.hpp"
#include "n2s/trim/cases.hpp"

#include "support/analytic_nurbs.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using n2s::KnotVector;
using n2s::NurbsCurve;
using n2s::NurbsSurface;

namespace {

/// Unique scratch path that removes itself, so a failing assertion cannot leave
/// a stale file behind for the next run to read.
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

void check_same_surface(const NurbsSurface& a, const NurbsSurface& b) {
    REQUIRE(a.num_u() == b.num_u());
    REQUIRE(a.num_v() == b.num_v());
    REQUIRE(a.degree_u() == b.degree_u());
    REQUIRE(a.degree_v() == b.degree_v());

    for (int i = 0; i <= 12; ++i) {
        for (int j = 0; j <= 12; ++j) {
            const double u = static_cast<double>(i) / 12.0;
            const double v = static_cast<double>(j) / 12.0;
            INFO("u = " << u << ", v = " << v);
            REQUIRE((a.evaluate(u, v) - b.evaluate(u, v)).norm() < 1e-15);
        }
    }
}

n2s::NurbsCurve circle_for_dimension_check() {
    return n2s::testing::full_circle(1.0);
}

} // namespace

TEST_CASE("a rational curve survives a JSON round trip", "[io][json]") {
    const NurbsCurve circle = n2s::testing::full_circle(2.0);
    const NurbsCurve restored = n2s::io::curve_from_json(n2s::io::to_json(circle));

    REQUIRE(restored.num_control_points() == circle.num_control_points());
    REQUIRE(restored.degree() == circle.degree());
    CHECK(restored.is_rational());

    for (int i = 0; i <= 200; ++i) {
        const double u = static_cast<double>(i) / 200.0;
        INFO("u = " << u);
        REQUIRE((circle.evaluate(u) - restored.evaluate(u)).norm() < 1e-15);
    }
}

TEST_CASE("a rational surface survives a JSON round trip", "[io][json]") {
    const NurbsSurface ball = n2s::testing::sphere(1.25);
    check_same_surface(ball, n2s::io::surface_from_json(n2s::io::to_json(ball)));
}

TEST_CASE("geometry survives a file round trip", "[io][json]") {
    const ScratchFile file("surface_roundtrip.json");
    const NurbsSurface patch = n2s::testing::bicubic_bezier(n2s::testing::wavy_bicubic_net());

    n2s::io::write_surface(file.path(), patch);
    REQUIRE(std::filesystem::exists(file.path()));

    check_same_surface(patch, n2s::io::read_surface(file.path()));
}

TEST_CASE("weights are optional and default to a polynomial surface", "[io][json]") {
    nlohmann::json document =
        n2s::io::to_json(n2s::testing::bicubic_bezier(n2s::testing::wavy_bicubic_net()));
    document.erase("weights");

    const NurbsSurface restored = n2s::io::surface_from_json(document);
    CHECK_FALSE(restored.is_rational());
    for (const double w : restored.weights()) {
        CHECK(w == Approx(1.0));
    }
}

TEST_CASE("malformed documents are rejected with the field named", "[io][json]") {
    const NurbsCurve circle = n2s::testing::full_circle(1.0);

    SECTION("missing field") {
        nlohmann::json document = n2s::io::to_json(circle);
        document.erase("knots");
        CHECK_THROWS_WITH(n2s::io::curve_from_json(document), ContainsSubstring("knots"));
    }

    SECTION("wrong format tag") {
        nlohmann::json document = n2s::io::to_json(circle);
        document["format"] = "n2s-surface";
        CHECK_THROWS_WITH(n2s::io::curve_from_json(document), ContainsSubstring("n2s-curve"));
    }

    SECTION("unknown schema version") {
        nlohmann::json document = n2s::io::to_json(circle);
        document["version"] = 99;
        CHECK_THROWS_WITH(n2s::io::curve_from_json(document), ContainsSubstring("version"));
    }

    SECTION("weight count disagrees with the control points") {
        nlohmann::json document = n2s::io::to_json(circle);
        document["weights"] = std::vector<double>{1.0, 1.0};
        CHECK_THROWS_WITH(n2s::io::curve_from_json(document), ContainsSubstring("weights"));
    }

    SECTION("a control point of the wrong width") {
        // The message names the width it wanted, because the same reader now
        // serves 2D domain curves and 3D model curves and "wrong number of
        // components" is ambiguous between them.
        nlohmann::json document = n2s::io::to_json(circle);
        document["control_points"][2] = {1.0, 2.0};
        CHECK_THROWS_WITH(n2s::io::curve_from_json(document), ContainsSubstring("3 numbers"));
    }
}

TEST_CASE("domain curves round-trip and are kept distinct from model curves", "[io][json]") {
    const n2s::NurbsCurve2 arc = n2s::cases::circle_loop({0.5, 0.5}, 0.25).curves().front();
    const nlohmann::json document = n2s::io::to_json(arc);
    CHECK(document.at("format") == "n2s-curve2");

    const n2s::NurbsCurve2 restored = n2s::io::curve2_from_json(document);
    for (int i = 0; i <= 32; ++i) {
        const double t = static_cast<double>(i) / 32.0;
        INFO("t = " << t);
        REQUIRE((arc.evaluate(t) - restored.evaluate(t)).norm() < 1e-15);
    }

    // A 2D curve must not read as a 3D one, or a trim loop would silently
    // become model geometry.
    CHECK_THROWS_WITH(n2s::io::curve_from_json(document), ContainsSubstring("n2s-curve"));
    CHECK_THROWS_WITH(n2s::io::curve2_from_json(n2s::io::to_json(circle_for_dimension_check())),
                      ContainsSubstring("n2s-curve2"));
}

TEST_CASE("a transposed control net is caught by the redundant shape fields", "[io][json]") {
    // A 4x6 net and a 6x4 net hold the same number of points, so the only thing
    // standing between a transposed file and a plausible but wrong surface is
    // the explicit num_u / num_v pair.
    const KnotVector knots_u = KnotVector::uniform_clamped(3, 4);
    const KnotVector knots_v = KnotVector::uniform_clamped(3, 6);
    const std::vector<Eigen::Vector3d> net(24, Eigen::Vector3d::Zero());
    const NurbsSurface patch = NurbsSurface::bspline(knots_u, knots_v, net);

    nlohmann::json document = n2s::io::to_json(patch);
    std::swap(document["num_u"], document["num_v"]);

    CHECK_THROWS_WITH(n2s::io::surface_from_json(document), ContainsSubstring("6 x 4"));
}

TEST_CASE("file errors name the path", "[io][json]") {
    const std::filesystem::path missing =
        std::filesystem::temp_directory_path() / "n2s_test_does_not_exist.json";
    std::filesystem::remove(missing);

    CHECK_THROWS_WITH(n2s::io::read_surface(missing), ContainsSubstring("does_not_exist"));
}

TEST_CASE("invalid but well-formed documents fail in the geometry constructor", "[io][json]") {
    nlohmann::json document = n2s::io::to_json(n2s::testing::full_circle(1.0));
    document["weights"][3] = -1.0;

    CHECK_THROWS_AS(n2s::io::curve_from_json(document), std::invalid_argument);
}
