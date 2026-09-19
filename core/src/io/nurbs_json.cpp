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

template<int Dim>
std::vector<Eigen::Matrix<double, Dim, 1>> read_points(const nlohmann::json& document,
                                                       const char* name) {
    const nlohmann::json& value = require_field(document, name);
    if (!value.is_array()) {
        throw std::runtime_error(fmt::format("field \"{}\" must be an array", name));
    }

    std::vector<Eigen::Matrix<double, Dim, 1>> points;
    points.reserve(value.size());
    for (std::size_t i = 0; i < value.size(); ++i) {
        const nlohmann::json& entry = value[i];
        if (!entry.is_array() || entry.size() != Dim) {
            throw std::runtime_error(fmt::format(
                "\"{}\"[{}] must be an array of {} numbers, got {}", name, i, Dim, entry.dump()));
        }
        Eigen::Matrix<double, Dim, 1> point;
        for (int k = 0; k < Dim; ++k) {
            point[k] = entry[static_cast<std::size_t>(k)].get<double>();
        }
        points.push_back(point);
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

template<int Dim>
nlohmann::json points_to_json(const std::vector<Eigen::Matrix<double, Dim, 1>>& points) {
    nlohmann::json array = nlohmann::json::array();
    for (const Eigen::Matrix<double, Dim, 1>& p : points) {
        nlohmann::json entry = nlohmann::json::array();
        for (int k = 0; k < Dim; ++k) {
            entry.push_back(p[k]);
        }
        array.push_back(std::move(entry));
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

namespace {

template<int Dim>
const char* curve_tag() {
    return Dim == 2 ? "n2s-curve2" : "n2s-curve";
}

template<int Dim>
nlohmann::json curve_to_json(const NurbsCurveT<Dim>& curve) {
    return nlohmann::json{
        {"format", curve_tag<Dim>()},
        {"version", kNurbsJsonVersion},
        {"degree", curve.degree()},
        {"knots", curve.knots().knots()},
        {"control_points", points_to_json<Dim>(curve.control_points())},
        {"weights", curve.weights()},
    };
}

template<int Dim>
NurbsCurveT<Dim> curve_from_document(const nlohmann::json& document) {
    require_format(document, curve_tag<Dim>());

    const int degree = read_degree(document, "degree");
    std::vector<double> knots = read_doubles(document, "knots");
    std::vector<Eigen::Matrix<double, Dim, 1>> points =
        read_points<Dim>(document, "control_points");
    std::vector<double> weights = read_weights(document, points.size());

    return NurbsCurveT<Dim>{
        KnotVector{degree, std::move(knots)}, std::move(points), std::move(weights)};
}

} // namespace

nlohmann::json to_json(const NurbsCurve& curve) {
    return curve_to_json<3>(curve);
}

nlohmann::json to_json(const NurbsCurve2& curve) {
    return curve_to_json<2>(curve);
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
        {"control_points", points_to_json<3>(surface.control_points())},
        {"weights", surface.weights()},
    };
}

NurbsCurve curve_from_json(const nlohmann::json& document) {
    return curve_from_document<3>(document);
}

NurbsCurve2 curve2_from_json(const nlohmann::json& document) {
    return curve_from_document<2>(document);
}

NurbsSurface surface_from_json(const nlohmann::json& document) {
    require_format(document, "n2s-surface");

    const int degree_u = read_degree(document, "degree_u");
    const int degree_v = read_degree(document, "degree_v");
    std::vector<double> knots_u = read_doubles(document, "knots_u");
    std::vector<double> knots_v = read_doubles(document, "knots_v");
    std::vector<Eigen::Vector3d> points = read_points<3>(document, "control_points");
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
