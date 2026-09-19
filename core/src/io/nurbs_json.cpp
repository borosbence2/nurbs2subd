#include "n2s/io/nurbs_json.hpp"

#include <fmt/format.h>

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace n2s::io {

namespace {

// Takes the name as const char* rather than std::string or std::string_view:
// both would be a class-type temporary at every call site, and a function
// returning a reference while taking one trips -Wdangling-reference.
const nlohmann::json& require_field(const nlohmann::json& document, const char* name) {
    const auto it = document.find(name);
    if (it == document.end()) {
        throw std::runtime_error(fmt::format("missing required field \"{}\"", name));
    }
    return *it;
}

void require_format(const nlohmann::json& document, const char* expected) {
    const nlohmann::json& format = require_field(document, "format");
    if (!format.is_string() || format.get<std::string>() != expected) {
        throw std::runtime_error(
            fmt::format("expected \"format\": \"{}\", got {}", expected, format.dump()));
    }

    const nlohmann::json& version = require_field(document, "version");
    if (!version.is_number_integer() || version.get<int>() != kNurbsJsonVersion) {
        throw std::runtime_error(fmt::format(
            "unsupported schema version {}; this build reads version {}. Refusing to guess "
            "at the meaning of a format it does not know.",
            version.dump(),
            kNurbsJsonVersion));
    }
}

std::vector<double> read_doubles(const nlohmann::json& document, const char* name) {
    const nlohmann::json& value = require_field(document, name);
    if (!value.is_array()) {
        throw std::runtime_error(fmt::format("field \"{}\" must be an array", name));
    }
    return value.get<std::vector<double>>();
}

std::vector<Eigen::Vector3d> read_points(const nlohmann::json& document, const char* name) {
    const nlohmann::json& value = require_field(document, name);
    if (!value.is_array()) {
        throw std::runtime_error(fmt::format("field \"{}\" must be an array", name));
    }

    std::vector<Eigen::Vector3d> points;
    points.reserve(value.size());
    for (std::size_t i = 0; i < value.size(); ++i) {
        const nlohmann::json& entry = value[i];
        if (!entry.is_array() || entry.size() != 3) {
            throw std::runtime_error(fmt::format(
                "\"{}\"[{}] must be an array of three numbers, got {}", name, i, entry.dump()));
        }
        points.emplace_back(entry[0].get<double>(), entry[1].get<double>(), entry[2].get<double>());
    }
    return points;
}

/// Weights are optional: a document without them describes a polynomial curve
/// or surface. A document *with* them must size them correctly.
std::vector<double> read_weights(const nlohmann::json& document, std::size_t expected) {
    const auto it = document.find("weights");
    if (it == document.end()) {
        return std::vector<double>(expected, 1.0);
    }
    if (!it->is_array()) {
        throw std::runtime_error("field \"weights\" must be an array");
    }
    auto weights = it->get<std::vector<double>>();
    if (weights.size() != expected) {
        throw std::runtime_error(
            fmt::format("\"weights\" has {} entries but there are {} control points",
                        weights.size(),
                        expected));
    }
    return weights;
}

int read_degree(const nlohmann::json& document, const char* name) {
    const nlohmann::json& value = require_field(document, name);
    if (!value.is_number_integer()) {
        throw std::runtime_error(fmt::format("field \"{}\" must be an integer", name));
    }
    return value.get<int>();
}

nlohmann::json points_to_json(const std::vector<Eigen::Vector3d>& points) {
    nlohmann::json array = nlohmann::json::array();
    for (const Eigen::Vector3d& p : points) {
        array.push_back({p.x(), p.y(), p.z()});
    }
    return array;
}

nlohmann::json read_document(const std::filesystem::path& path) {
    std::ifstream stream(path);
    if (!stream) {
        throw std::runtime_error(fmt::format("cannot open {} for reading", path.string()));
    }

    try {
        return nlohmann::json::parse(stream);
    } catch (const nlohmann::json::parse_error& error) {
        throw std::runtime_error(fmt::format("{}: {}", path.string(), error.what()));
    }
}

void write_document(const std::filesystem::path& path, const nlohmann::json& document) {
    std::ofstream stream(path);
    if (!stream) {
        throw std::runtime_error(fmt::format("cannot open {} for writing", path.string()));
    }
    stream << document.dump(2) << '\n';
}

} // namespace

nlohmann::json to_json(const NurbsCurve& curve) {
    return nlohmann::json{
        {"format", "n2s-curve"},
        {"version", kNurbsJsonVersion},
        {"degree", curve.degree()},
        {"knots", curve.knots().knots()},
        {"control_points", points_to_json(curve.control_points())},
        {"weights", curve.weights()},
    };
}

nlohmann::json to_json(const NurbsSurface& surface) {
    return nlohmann::json{
        {"format", "n2s-surface"},
        {"version", kNurbsJsonVersion},
        {"degree_u", surface.degree_u()},
        {"degree_v", surface.degree_v()},
        {"knots_u", surface.knots_u().knots()},
        {"knots_v", surface.knots_v().knots()},
        {"num_u", surface.num_u()},
        {"num_v", surface.num_v()},
        {"control_points", points_to_json(surface.control_points())},
        {"weights", surface.weights()},
    };
}

NurbsCurve curve_from_json(const nlohmann::json& document) {
    require_format(document, "n2s-curve");

    const int degree = read_degree(document, "degree");
    std::vector<double> knots = read_doubles(document, "knots");
    std::vector<Eigen::Vector3d> points = read_points(document, "control_points");
    std::vector<double> weights = read_weights(document, points.size());

    return NurbsCurve{KnotVector{degree, std::move(knots)}, std::move(points), std::move(weights)};
}

NurbsSurface surface_from_json(const nlohmann::json& document) {
    require_format(document, "n2s-surface");

    const int degree_u = read_degree(document, "degree_u");
    const int degree_v = read_degree(document, "degree_v");
    std::vector<double> knots_u = read_doubles(document, "knots_u");
    std::vector<double> knots_v = read_doubles(document, "knots_v");
    std::vector<Eigen::Vector3d> points = read_points(document, "control_points");
    std::vector<double> weights = read_weights(document, points.size());

    KnotVector kv_u{degree_u, std::move(knots_u)};
    KnotVector kv_v{degree_v, std::move(knots_v)};

    // num_u and num_v are redundant with the knot vectors, and that is the
    // point: a transposed control net still has the right total size, so the
    // only way to catch it is to state the shape twice and compare.
    const auto num_u = require_field(document, "num_u").get<std::size_t>();
    const auto num_v = require_field(document, "num_v").get<std::size_t>();
    if (num_u != kv_u.num_control_points() || num_v != kv_v.num_control_points()) {
        throw std::runtime_error(fmt::format(
            "the document says the control net is {} x {}, but the knot vectors imply {} x {}",
            num_u,
            num_v,
            kv_u.num_control_points(),
            kv_v.num_control_points()));
    }

    return NurbsSurface{std::move(kv_u), std::move(kv_v), std::move(points), std::move(weights)};
}

void write_curve(const std::filesystem::path& path, const NurbsCurve& curve) {
    write_document(path, to_json(curve));
}

void write_surface(const std::filesystem::path& path, const NurbsSurface& surface) {
    write_document(path, to_json(surface));
}

NurbsCurve read_curve(const std::filesystem::path& path) {
    return curve_from_json(read_document(path));
}

NurbsSurface read_surface(const std::filesystem::path& path) {
    return surface_from_json(read_document(path));
}

} // namespace n2s::io
