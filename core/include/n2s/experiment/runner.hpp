#pragma once

#include "n2s/io/case_json.hpp"
#include "n2s/metrics/surface_error.hpp"
#include "n2s/trim/sampling.hpp"
#include "n2s/trim/triangulation.hpp"

#include <nlohmann/json.hpp>

#include <filesystem>
#include <map>
#include <string>
#include <vector>

/// The experiment harness: one config in, one reproducible result directory
/// out.
///
/// CLAUDE.md requires every experiment to be reproducible from a config file,
/// with the results recording the config, the git hash and a timestamp. That is
/// enforced here rather than left to discipline: `run` copies the config into
/// the result directory and writes a `meta.json` next to it, so a result that
/// cannot be traced back to its inputs cannot be produced by this code path at
/// all.
namespace n2s::experiment {

struct RunConfig {
    /// Name of the run. The result directory is `<name>-<timestamp>`.
    std::string name = "run";

    /// Case file to load, relative to the config file's own directory.
    std::filesystem::path case_path;

    /// Layout size, used only when the case ships no control mesh of its own.
    int layout_rows = 6;
    int layout_columns = 6;

    SamplingOptions sampling;
    RefinementOptions refinement;
    metrics::SurfaceErrorOptions error;
    metrics::BoundaryErrorOptions boundary;

    /// Uniform refinement level for the exported OBJ of the limit surface.
    int export_level = 2;

    /// Write the per-sample CSV. Large, so optional.
    bool write_samples = true;
};

/// Parses a run config. Relative paths inside it are resolved against
/// `base_directory`, so a config is portable between machines.
RunConfig run_config_from_json(const nlohmann::json& document,
                               const std::filesystem::path& base_directory = ".");

nlohmann::json to_json(const RunConfig& config);

struct RunResult {
    std::string run_id;
    std::string case_name;

    metrics::SurfaceError surface;
    metrics::Stats boundary;

    std::size_t control_vertices = 0;
    std::size_t domain_triangles = 0;
    bool layout_from_case = false;

    /// Wall-clock milliseconds per stage, in execution order.
    std::vector<std::pair<std::string, double>> timings_ms;

    /// Anything the run wants the reader to know: trim repairs, unmeasured
    /// samples, a layout that had to be generated.
    std::vector<std::string> notes;
};

nlohmann::json to_json(const RunResult& result);

/// Runs one experiment and writes
/// `results/<run-id>/{config.json, metrics.json, samples.csv, *.obj, meta.json}`.
///
/// `meta.json` carries the git hash, the timestamp and the timings. Throws
/// `std::runtime_error` if the case will not load or the trim region does not
/// validate, because a run over a broken region would produce numbers that
/// describe nothing.
RunResult run(const RunConfig& config, const std::filesystem::path& results_root = "results");

/// A parameter sweep: a base config plus a list of overrides to apply to it.
///
/// Each point of the grid is a full run with its own result directory, and the
/// sweep writes an index listing them so the analysis scripts have one file to
/// start from.
struct SweepConfig {
    std::string name = "sweep";
    RunConfig base;

    /// Values to sweep, by dotted config key, e.g. "layout_rows" or
    /// "error.samples_per_face". The grid is the Cartesian product.
    std::map<std::string, std::vector<nlohmann::json>> axes;
};

SweepConfig sweep_config_from_json(const nlohmann::json& document,
                                   const std::filesystem::path& base_directory = ".");

struct SweepResult {
    std::string sweep_id;
    std::vector<RunResult> runs;
    std::vector<std::map<std::string, nlohmann::json>> points;
};

SweepResult sweep(const SweepConfig& config, const std::filesystem::path& results_root = "results");

} // namespace n2s::experiment
