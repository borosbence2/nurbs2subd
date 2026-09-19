// Headless experiment runner. The pipeline stages it drives arrive with M5; for
// now it exists so that the build, the argument surface and the CI smoke test
// are in place.

#include "n2s/build_info.hpp"
#include "n2s/io/case_json.hpp"
#include "n2s/trim/cases.hpp"
#include "n2s/trim/validate.hpp"

#include <CLI/CLI.hpp>
#include <fmt/core.h>

#include <filesystem>
#include <string>

namespace {

int run_command(const std::filesystem::path& config) {
    fmt::print("run: {}\n", config.string());
    fmt::print("The conversion pipeline is not implemented yet (plan milestone M5).\n");
    return 1;
}

int sweep_command(const std::filesystem::path& config) {
    fmt::print("sweep: {}\n", config.string());
    fmt::print("Parameter sweeps are not implemented yet (plan milestone M5).\n");
    return 1;
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
        const n2s::TrimReport report = n2s::validate_and_repair(generated.region);
        if (!report.ok()) {
            fmt::print(stderr, "{} did not validate:\n{}", generated.name, report.to_string());
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

    if (*run) {
        return run_command(run_config);
    }
    if (*sweep) {
        return sweep_command(sweep_config);
    }
    if (*export_cases) {
        return export_cases_command(export_directory);
    }

    fmt::print("{}", app.help());
    return 0;
}
