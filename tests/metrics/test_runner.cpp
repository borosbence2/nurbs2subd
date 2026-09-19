#include "n2s/experiment/runner.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using n2s::experiment::RunConfig;
using n2s::experiment::RunResult;

namespace {

/// Scratch results root that cleans up after itself.
class ScratchResults {
public:
    explicit ScratchResults(const std::string& name)
        : path_(std::filesystem::temp_directory_path() / ("n2s_results_" + name)) {
        std::filesystem::remove_all(path_);
    }

    ~ScratchResults() {
        std::error_code ignored;
        std::filesystem::remove_all(path_, ignored);
    }

    ScratchResults(const ScratchResults&) = delete;
    ScratchResults& operator=(const ScratchResults&) = delete;

    const std::filesystem::path& path() const { return path_; }

private:
    std::filesystem::path path_;
};

std::filesystem::path case_file(const char* name) {
    return std::filesystem::path{N2S_DATA_DIR} / name;
}

/// The reference configuration. Deliberately small, and deliberately free of
/// the interior refinement pass, whose vertex insertion order is the part of
/// the pipeline most sensitive to floating-point differences between
/// platforms. What it exercises is the metric arithmetic, which is what the
/// regression is guarding.
RunConfig reference_config() {
    RunConfig config;
    config.name = "regression";
    config.case_path = case_file("saddle_square_with_circular_hole_r0.25.json");
    config.layout_rows = 5;
    config.layout_columns = 5;

    config.sampling.mode = n2s::SamplingMode::Uniform;
    config.sampling.samples_per_span = 8;

    config.refinement.max_triangle_area = 0.0; // refinement off
    config.refinement.min_angle_degrees = 0.0;

    config.error.samples_per_face = 2;
    config.error.nurbs_samples_per_side = 8;
    config.error.restrict_to_trimmed_region = true;

    config.boundary.samples_per_edge = 8;
    config.boundary.trim_samples_per_curve = 64;

    config.export_level = 1;
    config.write_samples = true;
    return config;
}

} // namespace

TEST_CASE("a run writes every artefact the plan requires", "[experiment][runner]") {
    const ScratchResults results("artefacts");
    const RunResult result = n2s::experiment::run(reference_config(), results.path());

    const std::filesystem::path directory = results.path() / result.run_id;
    REQUIRE(std::filesystem::is_directory(directory));

    for (const char* name : {"config.json",
                             "metrics.json",
                             "meta.json",
                             "samples.csv",
                             "control_mesh.obj",
                             "limit_surface.obj",
                             "trimmed_nurbs.obj"}) {
        INFO(name);
        CHECK(std::filesystem::exists(directory / name));
        CHECK(std::filesystem::file_size(directory / name) > 0);
    }
}

TEST_CASE("meta.json records the provenance the hard rule demands", "[experiment][runner]") {
    // CLAUDE.md: every experiment records the config, the git hash and a
    // timestamp. Checked here rather than trusted, because a result that
    // cannot be traced to a commit is not usable in a write-up.
    const ScratchResults results("meta");
    const RunResult result = n2s::experiment::run(reference_config(), results.path());

    std::ifstream stream(results.path() / result.run_id / "meta.json");
    REQUIRE(stream);
    nlohmann::json meta;
    stream >> meta;

    CHECK(meta.at("run_id").get<std::string>() == result.run_id);
    CHECK_FALSE(meta.at("git_sha").get<std::string>().empty());
    CHECK_FALSE(meta.at("timestamp_utc").get<std::string>().empty());
    CHECK(meta.at("timestamp_utc").get<std::string>().back() == 'Z');
    CHECK(meta.contains("timings_ms"));
    CHECK(meta.at("timings_ms").contains("surface_error"));

    // The config travels with the result, not just its name.
    std::ifstream config_stream(results.path() / result.run_id / "config.json");
    nlohmann::json written_config;
    config_stream >> written_config;
    CHECK(written_config.at("layout_rows").get<int>() == 5);
}

TEST_CASE("a run config round-trips through JSON", "[experiment][runner]") {
    const RunConfig original = reference_config();
    const RunConfig restored =
        n2s::experiment::run_config_from_json(n2s::experiment::to_json(original), {});

    CHECK(restored.name == original.name);
    CHECK(restored.case_path == original.case_path);
    CHECK(restored.layout_rows == original.layout_rows);
    CHECK(restored.error.samples_per_face == original.error.samples_per_face);
    CHECK(restored.sampling.mode == n2s::SamplingMode::Uniform);
    CHECK(restored.refinement.max_triangle_area == Approx(0.0));
}

TEST_CASE("a config path is resolved once, not twice", "[experiment][runner]") {
    // The bug this pins down: a sweep serialises its base config, which already
    // holds a resolved case path, then re-reads it. Resolving against the base
    // directory a second time produced a doubled-up path and the sweep could
    // not open its own case.
    const nlohmann::json document{{"case", "sub/case.json"}};

    const RunConfig relative = n2s::experiment::run_config_from_json(document, "base");
    CHECK(relative.case_path == std::filesystem::path{"base"} / "sub/case.json");

    const RunConfig already_resolved = n2s::experiment::run_config_from_json(document, {});
    CHECK(already_resolved.case_path == std::filesystem::path{"sub/case.json"});
}

TEST_CASE("a run over an unusable case fails loudly", "[experiment][runner]") {
    const ScratchResults results("missing");
    RunConfig config = reference_config();
    config.case_path = case_file("there_is_no_such_case.json");

    CHECK_THROWS_WITH(n2s::experiment::run(config, results.path()),
                      ContainsSubstring("there_is_no_such_case"));
}

TEST_CASE("a sweep varies every axis it is given", "[experiment][sweep]") {
    const ScratchResults results("sweep");

    n2s::experiment::SweepConfig config;
    config.name = "axes";
    config.base = reference_config();
    config.base.write_samples = false;
    config.axes["layout_rows"] = {5, 7};
    config.axes["layout_columns"] = {5, 7};

    const n2s::experiment::SweepResult result = n2s::experiment::sweep(config, results.path());

    REQUIRE(result.runs.size() == 4);
    CHECK(std::filesystem::exists(results.path() / result.sweep_id / "sweep.json"));

    // Every combination produces its own control vertex count, which is the
    // cheapest proof that both axes reached the run rather than one silently
    // keeping its default.
    std::vector<std::size_t> counts;
    for (const RunResult& run : result.runs) {
        counts.push_back(run.control_vertices);
    }
    std::sort(counts.begin(), counts.end());
    CHECK(counts == std::vector<std::size_t>{25, 35, 35, 49});
}

TEST_CASE("a sweep needs axes to sweep", "[experiment][sweep]") {
    const nlohmann::json document{{"base", {{"case", "x.json"}}},
                                  {"axes", nlohmann::json::object()}};
    CHECK_THROWS_WITH(n2s::experiment::sweep_config_from_json(document), ContainsSubstring("axes"));
}

// ---------------------------------------------------------------------------
// The regression run.
// ---------------------------------------------------------------------------

TEST_CASE("the reference run's metrics stay within tolerance", "[experiment][regression]") {
    // A fixed, small run whose numbers are pinned. Its job is to notice when a
    // change somewhere in the pipeline moves the measured error, whether or not
    // that change was meant to.
    //
    // The tolerance is relative and loose enough to survive a different
    // compiler or standard library, and tight enough that anything which
    // actually alters the geometry or the metric arithmetic will trip it. If
    // one of these fails after a deliberate improvement, update the number
    // *and say so in the Progress notes* -- an unexplained edit here is how a
    // regression test stops meaning anything.
    const ScratchResults results("regression");
    const RunResult result = n2s::experiment::run(reference_config(), results.path());

    INFO("geometric max " << result.surface.geometric.max << ", rms "
                          << result.surface.geometric.rms << ", n "
                          << result.surface.geometric.count);
    INFO("boundary max " << result.boundary.max);
    INFO("normal max " << result.surface.normal_degrees.max);

    CHECK(result.control_vertices == 25);
    CHECK_FALSE(result.layout_from_case);

    // Measured 2026-09-19 on GCC 14.2 / Windows, Debug, at the config above.
    CHECK(result.surface.geometric.max == Approx(0.0296162).epsilon(1e-3));
    CHECK(result.surface.geometric.rms == Approx(0.0124712).epsilon(1e-2));
    CHECK(result.surface.hausdorff == Approx(0.0416667).epsilon(1e-3));
    CHECK(result.surface.normal_degrees.max == Approx(8.13010).epsilon(1e-2));
    CHECK(result.boundary.max == Approx(0.0416616).epsilon(1e-2));

    // Structural facts that must not drift silently.
    CHECK(result.surface.geometric.count > 0);
    CHECK(result.surface.parametric.count == result.surface.geometric.count);
    CHECK(result.surface.reverse_geometric.count > 0);
}
