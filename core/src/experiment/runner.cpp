#include "n2s/experiment/runner.hpp"

#include "n2s/build_info.hpp"
#include "n2s/fit/layout.hpp"
#include "n2s/fit/refine.hpp"
#include "n2s/subd/subdivision.hpp"
#include "n2s/trim/validate.hpp"

#include <fmt/format.h>

#include <chrono>
#include <ctime>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace n2s::experiment {

namespace {

using Clock = std::chrono::steady_clock;

/// Times a stage and records it, so that `meta.json` reports where a run spent
/// its time without anyone having to instrument by hand.
class StageTimer {
public:
    explicit StageTimer(std::vector<std::pair<std::string, double>>& into)
        : into_(into) {}

    template<typename Fn>
    auto run(const std::string& name, Fn&& body) -> decltype(body()) {
        const Clock::time_point start = Clock::now();
        if constexpr (std::is_void_v<decltype(body())>) {
            body();
            record(name, start);
        } else {
            auto result = body();
            record(name, start);
            return result;
        }
    }

private:
    void record(const std::string& name, Clock::time_point start) {
        const std::chrono::duration<double, std::milli> elapsed = Clock::now() - start;
        into_.emplace_back(name, elapsed.count());
    }

    std::vector<std::pair<std::string, double>>& into_;
};

/// UTC, ISO 8601, second resolution. Used in the run id, so it must be
/// filesystem-safe: colons are not.
std::string timestamp(bool filesystem_safe) {
    const std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &now);
#else
    gmtime_r(&now, &utc);
#endif

    // Two separate calls rather than a ternary between the format strings:
    // fmt checks the format at compile time, and a runtime-selected literal is
    // not a constant expression.
    if (filesystem_safe) {
        return fmt::format("{:04d}{:02d}{:02d}-{:02d}{:02d}{:02d}",
                           utc.tm_year + 1900,
                           utc.tm_mon + 1,
                           utc.tm_mday,
                           utc.tm_hour,
                           utc.tm_min,
                           utc.tm_sec);
    }
    return fmt::format("{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}Z",
                       utc.tm_year + 1900,
                       utc.tm_mon + 1,
                       utc.tm_mday,
                       utc.tm_hour,
                       utc.tm_min,
                       utc.tm_sec);
}

template<typename T>
T value_or(const nlohmann::json& document, const char* key, T fallback) {
    const auto found = document.find(key);
    return found == document.end() ? fallback : found->get<T>();
}

nlohmann::json stats_to_json(const metrics::Stats& stats) {
    return nlohmann::json{
        {"max", stats.max},
        {"mean", stats.mean},
        {"rms", stats.rms},
        {"count", stats.count},
        {"unmeasured", stats.unmeasured},
    };
}

void write_json(const std::filesystem::path& path, const nlohmann::json& document) {
    std::ofstream stream(path);
    if (!stream) {
        throw std::runtime_error(fmt::format("cannot open {} for writing", path.string()));
    }
    stream << document.dump(2) << '\n';
}

/// Applies a dotted key like "error.samples_per_face" to a config document.
void apply_override(nlohmann::json& document,
                    const std::string& dotted_key,
                    const nlohmann::json& value) {
    nlohmann::json* current = &document;
    std::size_t start = 0;
    while (true) {
        const std::size_t dot = dotted_key.find('.', start);
        const std::string part = dotted_key.substr(start, dot - start);
        if (dot == std::string::npos) {
            (*current)[part] = value;
            return;
        }
        current = &(*current)[part];
        start = dot + 1;
    }
}

} // namespace

std::string to_string(FitMethod method) {
    switch (method) {
    case FitMethod::Interpolate:
        return "interpolate";
    case FitMethod::Pia:
        return "pia";
    case FitMethod::None:
        break;
    }
    return "none";
}

FitMethod fit_method_from_string(const std::string& name) {
    if (name == "interpolate") {
        return FitMethod::Interpolate;
    }
    if (name == "pia") {
        return FitMethod::Pia;
    }
    if (name == "none") {
        return FitMethod::None;
    }
    throw std::runtime_error(
        fmt::format("unknown fit method \"{}\"; expected one of none, interpolate, pia", name));
}

RunConfig run_config_from_json(const nlohmann::json& document,
                               const std::filesystem::path& base_directory) {
    RunConfig config;
    config.name = value_or<std::string>(document, "name", "run");

    const auto case_path = document.find("case");
    if (case_path == document.end()) {
        throw std::runtime_error("run config is missing the required field \"case\"");
    }
    const std::filesystem::path given = case_path->get<std::string>();
    // Relative paths resolve against the config file, so a config and the case
    // it names travel together. An empty base means the path is already
    // resolved -- which is the case when a sweep re-reads a config it just
    // serialised, and resolving twice there produced a doubled-up path.
    config.case_path =
        (given.is_absolute() || base_directory.empty()) ? given : base_directory / given;

    config.layout_rows = value_or(document, "layout_rows", config.layout_rows);
    config.layout_columns = value_or(document, "layout_columns", config.layout_columns);
    config.export_level = value_or(document, "export_level", config.export_level);
    config.write_samples = value_or(document, "write_samples", config.write_samples);

    if (const auto found = document.find("fit"); found != document.end()) {
        config.fit.method = fit_method_from_string(value_or<std::string>(*found, "method", "none"));
        config.fit.layout_refinement =
            value_or(*found, "layout_refinement", config.fit.layout_refinement);
        config.fit.write_convergence =
            value_or(*found, "write_convergence", config.fit.write_convergence);

        if (const auto pia = found->find("pia"); pia != found->end()) {
            config.fit.pia.max_iterations =
                value_or(*pia, "max_iterations", config.fit.pia.max_iterations);
            config.fit.pia.relative_update_tolerance = value_or(
                *pia, "relative_update_tolerance", config.fit.pia.relative_update_tolerance);
        }
    }

    if (const auto found = document.find("sampling"); found != document.end()) {
        config.sampling.max_segment_length =
            value_or(*found, "max_segment_length", config.sampling.max_segment_length);
        config.sampling.max_sagitta = value_or(*found, "max_sagitta", config.sampling.max_sagitta);
        config.sampling.samples_per_span =
            value_or(*found, "samples_per_span", config.sampling.samples_per_span);
        if (value_or<std::string>(*found, "mode", "adaptive") == "uniform") {
            config.sampling.mode = SamplingMode::Uniform;
        }
    }

    if (const auto found = document.find("refinement"); found != document.end()) {
        config.refinement.max_triangle_area =
            value_or(*found, "max_triangle_area", config.refinement.max_triangle_area);
        config.refinement.min_angle_degrees =
            value_or(*found, "min_angle_degrees", config.refinement.min_angle_degrees);
        config.refinement.max_vertices =
            value_or(*found, "max_vertices", config.refinement.max_vertices);
    }

    if (const auto found = document.find("error"); found != document.end()) {
        config.error.samples_per_face =
            value_or(*found, "samples_per_face", config.error.samples_per_face);
        config.error.nurbs_samples_per_side =
            value_or(*found, "nurbs_samples_per_side", config.error.nurbs_samples_per_side);
        config.error.restrict_to_trimmed_region =
            value_or(*found, "restrict_to_trimmed_region", config.error.restrict_to_trimmed_region);
    }

    if (const auto found = document.find("boundary"); found != document.end()) {
        config.boundary.samples_per_edge =
            value_or(*found, "samples_per_edge", config.boundary.samples_per_edge);
        config.boundary.trim_samples_per_curve =
            value_or(*found, "trim_samples_per_curve", config.boundary.trim_samples_per_curve);
    }

    return config;
}

nlohmann::json to_json(const RunConfig& config) {
    return nlohmann::json{
        {"name", config.name},
        {"case", config.case_path.generic_string()},
        {"layout_rows", config.layout_rows},
        {"layout_columns", config.layout_columns},
        {"export_level", config.export_level},
        {"write_samples", config.write_samples},
        {"fit",
         {{"method", to_string(config.fit.method)},
          {"layout_refinement", config.fit.layout_refinement},
          {"write_convergence", config.fit.write_convergence},
          {"pia",
           {{"max_iterations", config.fit.pia.max_iterations},
            {"relative_update_tolerance", config.fit.pia.relative_update_tolerance}}}}},
        {"sampling",
         {{"mode", config.sampling.mode == SamplingMode::Uniform ? "uniform" : "adaptive"},
          {"max_segment_length", config.sampling.max_segment_length},
          {"max_sagitta", config.sampling.max_sagitta},
          {"samples_per_span", config.sampling.samples_per_span}}},
        {"refinement",
         {{"max_triangle_area", config.refinement.max_triangle_area},
          {"min_angle_degrees", config.refinement.min_angle_degrees},
          {"max_vertices", config.refinement.max_vertices}}},
        {"error",
         {{"samples_per_face", config.error.samples_per_face},
          {"nurbs_samples_per_side", config.error.nurbs_samples_per_side},
          {"restrict_to_trimmed_region", config.error.restrict_to_trimmed_region}}},
        {"boundary",
         {{"samples_per_edge", config.boundary.samples_per_edge},
          {"trim_samples_per_curve", config.boundary.trim_samples_per_curve}}},
    };
}

nlohmann::json to_json(const RunResult& result) {
    nlohmann::json timings = nlohmann::json::object();
    for (const auto& [stage, milliseconds] : result.timings_ms) {
        timings[stage] = milliseconds;
    }

    return nlohmann::json{
        {"run_id", result.run_id},
        {"case", result.case_name},
        {"layout_from_case", result.layout_from_case},
        {"control_vertices", result.control_vertices},
        {"fit",
         {{"method", to_string(result.fit_method)},
          {"layout_refinement", result.layout_refinement},
          {"converged", result.fit_report.converged},
          {"iterations", result.fit_report.iterations},
          {"max_interpolation_error", result.fit_report.max_interpolation_error},
          {"rms_interpolation_error", result.fit_report.rms_interpolation_error},
          {"notes", result.fit_report.notes}}},
        {"domain_triangles", result.domain_triangles},
        {"surface_error",
         {{"parametric", stats_to_json(result.surface.parametric)},
          {"geometric", stats_to_json(result.surface.geometric)},
          {"reverse_geometric", stats_to_json(result.surface.reverse_geometric)},
          {"hausdorff", result.surface.hausdorff},
          {"normal_degrees", stats_to_json(result.surface.normal_degrees)},
          {"mean_curvature", stats_to_json(result.surface.mean_curvature)},
          {"gaussian_curvature", stats_to_json(result.surface.gaussian_curvature)}}},
        {"boundary_error", stats_to_json(result.boundary)},
        {"timings_ms", std::move(timings)},
        {"notes", result.notes},
    };
}

RunResult run(const RunConfig& config, const std::filesystem::path& results_root) {
    RunResult result;
    result.run_id = fmt::format("{}-{}", config.name, timestamp(true));

    StageTimer timer(result.timings_ms);

    io::Case loaded = timer.run("load_case", [&] { return io::read_case(config.case_path); });
    result.case_name = loaded.name;

    const TrimReport report =
        timer.run("validate_trim", [&] { return validate_and_repair(loaded.region); });
    for (const std::string& repair : report.repairs) {
        result.notes.push_back("trim repaired: " + repair);
    }
    if (!report.ok()) {
        throw std::runtime_error(
            fmt::format("the trim region of {} did not validate, so there is nothing meaningful to "
                        "measure:\n{}",
                        config.case_path.string(),
                        report.to_string()));
    }

    const DomainMesh domain = timer.run("triangulate", [&] {
        return triangulate_refined(
            loaded.region, loaded.surface, config.sampling, config.refinement);
    });
    result.domain_triangles = domain.triangles.size();

    // The layout. A domain layout is preferred wherever one exists, because it
    // is the only form that can be refined and fitted: refinement subdivides
    // quads in the domain, and a fit needs the correspondence to know what
    // each control point is supposed to approximate.
    result.layout_from_case = loaded.layout.has_value() || loaded.control_mesh.has_value();
    result.layout_refinement = config.fit.layout_refinement;
    result.fit_method = config.fit.method;

    std::optional<fit::DomainLayout> domain_layout;
    if (loaded.layout.has_value()) {
        domain_layout = *loaded.layout;
    } else if (!loaded.control_mesh.has_value()) {
        domain_layout = fit::grid_domain_layout(config.layout_rows, config.layout_columns);
        result.notes.push_back(
            fmt::format("the case ships no layout, so a {} x {} grid was generated.",
                        config.layout_rows,
                        config.layout_columns));
    }

    if (domain_layout.has_value() && config.fit.layout_refinement > 1) {
        const std::size_t before = domain_layout->num_quads();
        domain_layout = timer.run("refine_layout", [&] {
            return fit::refine_quads(*domain_layout, config.fit.layout_refinement);
        });
        result.notes.push_back(fmt::format("layout refined {}x per side: {} quads -> {}",
                                           config.fit.layout_refinement,
                                           before,
                                           domain_layout->num_quads()));
    }

    metrics::DomainMap domain_map;
    ControlMesh layout = [&] {
        if (!domain_layout.has_value()) {
            // A bare control mesh: usable, but neither refinable nor fittable,
            // and it carries no correspondence.
            result.notes.emplace_back(
                "the case supplied bare control points rather than a domain layout, so it "
                "cannot be refined or fitted and no correspondence is known. The parametric, "
                "normal and curvature statistics are absent rather than guessed.");
            return *loaded.control_mesh;
        }

        domain_map = fit::bilinear_domain_map(*domain_layout);

        switch (config.fit.method) {
        case FitMethod::Interpolate:
            return timer.run("fit", [&] {
                return fit::solve_interpolation(*domain_layout, loaded.surface, result.fit_report);
            });
        case FitMethod::Pia:
            return timer.run("fit", [&] {
                return fit::solve_pia(
                    *domain_layout, loaded.surface, result.fit_report, config.fit.pia);
            });
        case FitMethod::None:
            break;
        }

        result.notes.emplace_back(
            "no fit was applied, so the control points are simply the layout lifted onto the "
            "surface. Catmull-Clark pulls its limit surface inside the control net, and that "
            "shrinkage dominates the error below.");
        return fit::lift(*domain_layout, loaded.surface);
    }();

    result.control_vertices = layout.num_vertices();

    if (config.fit.method != FitMethod::None && !result.fit_report.converged) {
        result.notes.emplace_back(
            "the fit did not converge, so every error below describes a surface that is not "
            "the one the method was supposed to produce");
    }
    for (const std::string& note : result.fit_report.notes) {
        result.notes.push_back("fit: " + note);
    }

    const SubdivisionSurface limit{std::move(layout)};

    result.surface = timer.run("surface_error", [&] {
        return metrics::measure_surface_error(
            limit, loaded.surface, loaded.region, domain_map, config.error);
    });
    result.boundary = timer.run("boundary_error", [&] {
        return metrics::measure_boundary_error(
            limit, loaded.surface, loaded.region, config.boundary);
    });

    if (!result.surface.geometric.complete()) {
        result.notes.push_back(
            fmt::format("{} of {} surface samples could not be measured and are excluded",
                        result.surface.geometric.unmeasured,
                        result.surface.geometric.unmeasured + result.surface.geometric.count));
    }

    // ---- Write the result directory ----------------------------------------
    const std::filesystem::path directory = results_root / result.run_id;
    timer.run("write_results", [&] {
        std::filesystem::create_directories(directory);

        write_json(directory / "config.json", to_json(config));
        write_json(directory / "metrics.json", to_json(result));

        limit.control_mesh().write_obj(directory / "control_mesh.obj");
        limit.refine_uniform(config.export_level).write_obj(directory / "limit_surface.obj");
        map_to_surface(domain, loaded.surface).write_obj(directory / "trimmed_nurbs.obj");

        // The per-iteration history behind the defect 5 comparison. Written
        // whenever PIA ran, so a result that reports the finding can also
        // show it.
        if (config.fit.method == FitMethod::Pia && config.fit.write_convergence &&
            domain_layout.has_value()) {
            const fit::PiaHistory history =
                fit::pia_error_history(*domain_layout, loaded.surface, config.fit.pia);

            std::ofstream convergence(directory / "pia_convergence.csv");
            if (!convergence) {
                throw std::runtime_error("cannot open pia_convergence.csv for writing");
            }
            convergence << "iteration,max_error,max_update,distance_to_direct\n";
            for (std::size_t i = 0; i < history.max_error.size(); ++i) {
                convergence << fmt::format("{},{:.17g},{:.17g},{:.17g}\n",
                                           i + 1,
                                           history.max_error[i],
                                           history.max_update[i],
                                           history.distance_to_direct[i]);
            }
        }

        if (config.write_samples) {
            std::ofstream csv(directory / "samples.csv");
            if (!csv) {
                throw std::runtime_error("cannot open samples.csv for writing");
            }
            csv << "face,u,v,domain_u,domain_v,x,y,z,parametric,geometric,normal_degrees,"
                   "measured,mean_curvature,gaussian_curvature,curvature_measured\n";
            for (const metrics::ErrorSample& sample : metrics::sample_surface_error(
                     limit, loaded.surface, loaded.region, domain_map, config.error)) {
                csv << fmt::format("{},{:.17g},{:.17g},{:.17g},{:.17g},{:.17g},{:.17g},{:.17g},"
                                   "{:.17g},{:.17g},{:.17g},{},{:.17g},{:.17g},{}\n",
                                   sample.location.face,
                                   sample.location.u,
                                   sample.location.v,
                                   sample.domain_point.x(),
                                   sample.domain_point.y(),
                                   sample.limit_point.x(),
                                   sample.limit_point.y(),
                                   sample.limit_point.z(),
                                   sample.parametric,
                                   sample.geometric,
                                   sample.normal_degrees,
                                   sample.measured ? 1 : 0,
                                   sample.mean_curvature,
                                   sample.gaussian_curvature,
                                   sample.curvature_measured ? 1 : 0);
            }
        }
    });

    // meta.json last, so that the timings it records include everything above.
    const BuildInfo& build = build_info();
    write_json(directory / "meta.json",
               nlohmann::json{
                   {"run_id", result.run_id},
                   {"timestamp_utc", timestamp(false)},
                   {"git_sha", std::string{build.git_sha}},
                   {"git_describe", std::string{build.git_describe}},
                   {"git_dirty", std::string{build.git_dirty}},
                   {"version", std::string{build.version}},
                   {"build_type", std::string{build.build_type}},
                   {"compiler", std::string{build.compiler}},
                   {"timings_ms", to_json(result)["timings_ms"]},
               });

    return result;
}

SweepConfig sweep_config_from_json(const nlohmann::json& document,
                                   const std::filesystem::path& base_directory) {
    SweepConfig config;
    config.name = value_or<std::string>(document, "name", "sweep");

    const auto base = document.find("base");
    if (base == document.end()) {
        throw std::runtime_error("sweep config is missing the required field \"base\"");
    }
    config.base = run_config_from_json(*base, base_directory);

    const auto axes = document.find("axes");
    if (axes == document.end() || axes->empty()) {
        throw std::runtime_error(
            "sweep config needs a non-empty \"axes\" object mapping config keys to value lists");
    }
    for (const auto& [key, values] : axes->items()) {
        if (!values.is_array() || values.empty()) {
            throw std::runtime_error(
                fmt::format("axis \"{}\" must be a non-empty array of values", key));
        }
        config.axes[key] = values.get<std::vector<nlohmann::json>>();
    }

    return config;
}

SweepResult sweep(const SweepConfig& config, const std::filesystem::path& results_root) {
    SweepResult result;
    result.sweep_id = fmt::format("{}-{}", config.name, timestamp(true));

    // Cartesian product of the axes, built iteratively so the axis count is not
    // fixed at compile time.
    std::vector<std::map<std::string, nlohmann::json>> points{{}};
    for (const auto& [key, values] : config.axes) {
        std::vector<std::map<std::string, nlohmann::json>> expanded;
        expanded.reserve(points.size() * values.size());
        for (const auto& point : points) {
            for (const nlohmann::json& value : values) {
                auto extended = point;
                extended[key] = value;
                expanded.push_back(std::move(extended));
            }
        }
        points = std::move(expanded);
    }

    const std::filesystem::path directory = results_root / result.sweep_id;
    std::filesystem::create_directories(directory);

    nlohmann::json index = nlohmann::json::array();

    for (std::size_t i = 0; i < points.size(); ++i) {
        nlohmann::json document = to_json(config.base);
        for (const auto& [key, value] : points[i]) {
            apply_override(document, key, value);
        }
        document["name"] = fmt::format("{}-{:03d}", config.name, i);

        // The serialised base already holds a resolved case path, so this
        // re-read must not resolve it again.
        const RunConfig point_config = run_config_from_json(document, {});
        const RunResult run_result = run(point_config, results_root);

        nlohmann::json entry{{"run_id", run_result.run_id}, {"point", points[i]}};
        entry["metrics"] = to_json(run_result);
        index.push_back(std::move(entry));

        result.runs.push_back(run_result);
        result.points.push_back(points[i]);
    }

    write_json(directory / "sweep.json",
               nlohmann::json{{"sweep_id", result.sweep_id},
                              {"base", to_json(config.base)},
                              {"runs", std::move(index)}});

    return result;
}

} // namespace n2s::experiment
