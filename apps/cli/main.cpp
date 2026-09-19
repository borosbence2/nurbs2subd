// Headless experiment runner. The pipeline stages it drives arrive with M5; for
// now it exists so that the build, the argument surface and the CI smoke test
// are in place.

#include "n2s/build_info.hpp"

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

    fmt::print("{}", app.help());
    return 0;
}
