// Headless experiment runner.
//
// Every subcommand that produces numbers writes them to a result directory
// alongside the config that produced them and a meta.json recording the git
// hash and timestamp. That is a hard rule in CLAUDE.md, and doing it here
// rather than by convention means an untraceable result cannot be produced
// through this tool at all.

#include "n2s/build_info.hpp"
#include "n2s/experiment/runner.hpp"
#include "n2s/io/case_json.hpp"
#include "n2s/trim/cases.hpp"
#include "n2s/trim/validate.hpp"

#include <CLI/CLI.hpp>
#include <fmt/core.h>
#include <nlohmann/json.hpp>

#include <exception>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <system_error>

namespace {

nlohmann::json read_json(const std::filesystem::path& path) {
    std::ifstream stream(path);
    if (!stream) {
        throw std::runtime_error(fmt::format("cannot open {} for reading", path.string()));
    }
    return nlohmann::json::parse(stream);
}

void print_stats(const char* label, const n2s::metrics::Stats& stats) {
    if (stats.count == 0) {
        fmt::print("  {:<22} (not measured)\n", label);
        return;
    }
    fmt::print("  {:<22} max {:<12.6g} rms {:<12.6g} mean {:<12.6g} n={}",
               label,
               stats.max,
               stats.rms,
               stats.mean,
               stats.count);
    if (stats.unmeasured > 0) {
        fmt::print("  [{} unmeasured]", stats.unmeasured);
    }
    fmt::print("\n");
}

void report(const n2s::experiment::RunResult& result, const std::filesystem::path& results_root) {
    fmt::print("\n{} ({})\n", result.run_id, result.case_name);
    fmt::print("  layout {}, {} control vertices, {} domain triangles\n",
               result.layout_from_case ? "from case" : "generated",
               result.control_vertices,
               result.domain_triangles);

    print_stats("parametric error", result.surface.parametric);
    print_stats("geometric error", result.surface.geometric);
    print_stats("reverse geometric", result.surface.reverse_geometric);
    fmt::print("  {:<22} {:.6g}\n", "hausdorff", result.surface.hausdorff);
    print_stats("normal (degrees)", result.surface.normal_degrees);
    print_stats("mean curvature", result.surface.mean_curvature);
    print_stats("gaussian curvature", result.surface.gaussian_curvature);
    print_stats("boundary deviation", result.boundary);

    for (const std::string& note : result.notes) {
        fmt::print("  note: {}\n", note);
    }
    fmt::print("  written to {}\n", (results_root / result.run_id).string());
}

int run_command(const std::filesystem::path& config_path,
                const std::filesystem::path& results_root) {
    const nlohmann::json document = read_json(config_path);
    const n2s::experiment::RunConfig config =
        n2s::experiment::run_config_from_json(document, config_path.parent_path());

    const n2s::experiment::RunResult result = n2s::experiment::run(config, results_root);
    report(result, results_root);
    return 0;
}

int sweep_command(const std::filesystem::path& config_path,
                  const std::filesystem::path& results_root) {
    const nlohmann::json document = read_json(config_path);
    const n2s::experiment::SweepConfig config =
        n2s::experiment::sweep_config_from_json(document, config_path.parent_path());

    const n2s::experiment::SweepResult result = n2s::experiment::sweep(config, results_root);

    fmt::print("\nsweep {}: {} runs\n", result.sweep_id, result.runs.size());
    for (std::size_t i = 0; i < result.runs.size(); ++i) {
        std::string point;
        for (const auto& [key, value] : result.points[i]) {
            point += fmt::format("{}={} ", key, value.dump());
        }
        fmt::print("  {:<44} geometric max {:.6g}\n", point, result.runs[i].surface.geometric.max);
    }
    fmt::print("  index written to {}\n", (results_root / result.sweep_id / "sweep.json").string());
    return 0;
}

/// Writes every synthetic case to JSON, so that `data/` is generated from the
/// same code the tests use rather than hand-maintained. Regenerating it after a
/// change to the case definitions is one command, and the diff shows exactly
/// what moved.
int export_cases_command(const std::filesystem::path& directory) {
    std::error_code error;
    std::filesystem::create_directories(directory, error);
    if (error) {
        fmt::print(stderr, "cannot create {}: {}\n", directory.string(), error.message());
        return 1;
    }

    int written = 0;
    for (n2s::cases::TrimmedCase& generated : n2s::cases::all_synthetic_cases()) {
        const n2s::TrimReport report_for_case = n2s::validate_and_repair(generated.region);
        if (!report_for_case.ok()) {
            fmt::print(
                stderr, "{} did not validate:\n{}", generated.name, report_for_case.to_string());
            return 1;
        }

        const n2s::io::Case test_case{
            .name = generated.name,
            .surface = generated.surface,
            .region = std::move(generated.region),
            .control_mesh = std::nullopt,
            .camera = std::nullopt,
            // Ranges have to come from somewhere, and a file with none would
            // make every figure drawn from it auto-ranged. These are starting
            // points to be tightened per experiment, not measurements.
            .color_ranges = {{"mean_curvature", {-2.0, 2.0}},
                             {"gaussian_curvature", {-4.0, 4.0}},
                             {"error", {0.0, 0.05}}},
        };

        const std::filesystem::path path = directory / (generated.name + ".json");
        n2s::io::write_case(path, test_case);
        fmt::print("{}\n", path.string());
        ++written;
    }

    fmt::print("wrote {} cases to {}\n", written, directory.string());
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    CLI::App app{"nurbs2subd - trimmed NURBS to Catmull-Clark conversion", "nurbs2subd"};
    app.require_subcommand(0, 1);

    bool show_version = false;
    app.add_flag("-v,--version", show_version, "Print version and build provenance");

    std::filesystem::path results_root{"results"};
    app.add_option("--results", results_root, "Directory to write result directories into")
        ->capture_default_str();

    std::filesystem::path run_config;
    CLI::App* run = app.add_subcommand("run", "Run one experiment from a config file");
    run->add_option("config", run_config, "Path to the experiment config JSON")
        ->required()
        ->check(CLI::ExistingFile);

    std::filesystem::path sweep_config;
    CLI::App* sweep = app.add_subcommand("sweep", "Run a parameter sweep from a config file");
    sweep->add_option("config", sweep_config, "Path to the sweep config JSON")
        ->required()
        ->check(CLI::ExistingFile);

    std::filesystem::path export_directory{"data"};
    CLI::App* export_cases =
        app.add_subcommand("export-cases", "Write the synthetic test cases to JSON");
    export_cases->add_option("directory", export_directory, "Where to write the case files")
        ->capture_default_str();

    CLI11_PARSE(app, argc, argv);

    if (show_version) {
        fmt::print("{}\n", n2s::build_info_string());
        return 0;
    }

    try {
        if (*run) {
            return run_command(run_config, results_root);
        }
        if (*sweep) {
            return sweep_command(sweep_config, results_root);
        }
        if (*export_cases) {
            return export_cases_command(export_directory);
        }
    } catch (const std::exception& error) {
        // Reported rather than allowed to terminate: a failed run must say why,
        // because the alternative is a missing result directory and no clue.
        fmt::print(stderr, "error: {}\n", error.what());
        return 1;
    }

    fmt::print("{}", app.help());
    return 0;
}
