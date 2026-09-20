#include "n2s/io/case_json.hpp"

#include "n2s/io/nurbs_json.hpp"

#include <fmt/format.h>

#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace n2s::io {

namespace {

const nlohmann::json& require(const nlohmann::json& document, const char* name) {
    const auto it = document.find(name);
    if (it == document.end()) {
        throw std::runtime_error(fmt::format("missing required field \"{}\"", name));
    }
    return *it;
}

nlohmann::json vector_to_json(const Eigen::Vector3d& v) {
    return nlohmann::json::array({v.x(), v.y(), v.z()});
}

Eigen::Vector3d vector_from_json(const nlohmann::json& value, const char* name) {
    if (!value.is_array() || value.size() != 3) {
        throw std::runtime_error(
            fmt::format("\"{}\" must be an array of three numbers, got {}", name, value.dump()));
    }
    return Eigen::Vector3d{value[0].get<double>(), value[1].get<double>(), value[2].get<double>()};
}

nlohmann::json loop_to_json(const TrimLoop& loop) {
    nlohmann::json curves = nlohmann::json::array();
    for (const NurbsCurve2& curve : loop.curves()) {
        curves.push_back(to_json(curve));
    }
    return curves;
}

TrimLoop loop_from_json(const nlohmann::json& value, const std::string& where) {
    if (!value.is_array() || value.empty()) {
        throw std::runtime_error(fmt::format("{} must be a non-empty array of trim curves", where));
    }

    std::vector<NurbsCurve2> curves;
    curves.reserve(value.size());
    for (std::size_t i = 0; i < value.size(); ++i) {
        try {
            curves.push_back(curve2_from_json(value[i]));
        } catch (const std::exception& error) {
            throw std::runtime_error(fmt::format("{}, curve {}: {}", where, i, error.what()));
        }
    }
    return TrimLoop{std::move(curves)};
}

nlohmann::json mesh_to_json(const ControlMesh& mesh) {
    nlohmann::json vertices = nlohmann::json::array();
    for (const Eigen::Vector3d& v : mesh.vertices()) {
        vertices.push_back(vector_to_json(v));
    }

    nlohmann::json quads = nlohmann::json::array();
    for (const std::array<int, 4>& q : mesh.quads()) {
        quads.push_back(nlohmann::json::array({q[0], q[1], q[2], q[3]}));
    }

    nlohmann::json creases = nlohmann::json::array();
    for (const ControlMesh::Crease& crease : mesh.creases()) {
        creases.push_back(
            nlohmann::json{{"v0", crease.v0}, {"v1", crease.v1}, {"sharpness", crease.sharpness}});
    }

    nlohmann::json corners = nlohmann::json::array();
    for (const ControlMesh::Corner& corner : mesh.corners()) {
        corners.push_back(
            nlohmann::json{{"vertex", corner.vertex}, {"sharpness", corner.sharpness}});
    }

    return nlohmann::json{{"vertices", std::move(vertices)},
                          {"quads", std::move(quads)},
                          {"creases", std::move(creases)},
                          {"corners", std::move(corners)}};
}

ControlMesh mesh_from_json(const nlohmann::json& value) {
    const nlohmann::json& vertex_array = require(value, "vertices");
    std::vector<Eigen::Vector3d> vertices;
    vertices.reserve(vertex_array.size());
    for (std::size_t i = 0; i < vertex_array.size(); ++i) {
        vertices.push_back(
            vector_from_json(vertex_array[i], fmt::format("vertices[{}]", i).c_str()));
    }

    const nlohmann::json& quad_array = require(value, "quads");
    std::vector<std::array<int, 4>> quads;
    quads.reserve(quad_array.size());
    for (std::size_t i = 0; i < quad_array.size(); ++i) {
        const nlohmann::json& entry = quad_array[i];
        if (!entry.is_array() || entry.size() != 4) {
            throw std::runtime_error(fmt::format(
                "quads[{}] must be an array of four vertex indices, got {}", i, entry.dump()));
        }
        quads.push_back(
            {entry[0].get<int>(), entry[1].get<int>(), entry[2].get<int>(), entry[3].get<int>()});
    }

    ControlMesh mesh{std::move(vertices), std::move(quads)};

    if (const auto creases = value.find("creases"); creases != value.end()) {
        for (const nlohmann::json& crease : *creases) {
            mesh.set_crease(require(crease, "v0").get<int>(),
                            require(crease, "v1").get<int>(),
                            require(crease, "sharpness").get<double>());
        }
    }
    if (const auto corners = value.find("corners"); corners != value.end()) {
        for (const nlohmann::json& corner : *corners) {
            mesh.set_corner(require(corner, "vertex").get<int>(),
                            require(corner, "sharpness").get<double>());
        }
    }

    return mesh;
}

nlohmann::json layout_to_json(const fit::DomainLayout& layout) {
    nlohmann::json vertices = nlohmann::json::array();
    for (const Eigen::Vector2d& v : layout.vertices) {
        vertices.push_back(nlohmann::json::array({v.x(), v.y()}));
    }

    nlohmann::json quads = nlohmann::json::array();
    for (const std::array<int, 4>& q : layout.quads) {
        quads.push_back(nlohmann::json::array({q[0], q[1], q[2], q[3]}));
    }

    return nlohmann::json{{"vertices", std::move(vertices)}, {"quads", std::move(quads)}};
}

fit::DomainLayout layout_from_json(const nlohmann::json& value) {
    fit::DomainLayout layout;

    const nlohmann::json& vertices = require(value, "vertices");
    layout.vertices.reserve(vertices.size());
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const nlohmann::json& entry = vertices[i];
        if (!entry.is_array() || entry.size() != 2) {
            throw std::runtime_error(
                fmt::format("layout.vertices[{}] must be a (u, v) pair, got {}", i, entry.dump()));
        }
        layout.vertices.emplace_back(entry[0].get<double>(), entry[1].get<double>());
    }

    const nlohmann::json& quads = require(value, "quads");
    layout.quads.reserve(quads.size());
    for (std::size_t i = 0; i < quads.size(); ++i) {
        const nlohmann::json& entry = quads[i];
        if (!entry.is_array() || entry.size() != 4) {
            throw std::runtime_error(fmt::format(
                "layout.quads[{}] must be four vertex indices, got {}", i, entry.dump()));
        }
        layout.quads.push_back(
            {entry[0].get<int>(), entry[1].get<int>(), entry[2].get<int>(), entry[3].get<int>()});
    }

    return layout;
}

} // namespace

nlohmann::json to_json(const Case& test_case) {
    nlohmann::json holes = nlohmann::json::array();
    for (const TrimLoop& hole : test_case.region.holes()) {
        holes.push_back(loop_to_json(hole));
    }

    nlohmann::json document{
        {"format", "n2s-case"},
        {"version", kCaseJsonVersion},
        {"name", test_case.name},
        {"surface", to_json(test_case.surface)},
        {"trim",
         nlohmann::json{{"outer", loop_to_json(test_case.region.outer())},
                        {"holes", std::move(holes)}}},
    };

    if (test_case.layout.has_value()) {
        document["layout"] = layout_to_json(*test_case.layout);
    }

    if (test_case.control_mesh.has_value()) {
        document["control_mesh"] = mesh_to_json(*test_case.control_mesh);
    }

    if (test_case.camera.has_value()) {
        document["camera"] = nlohmann::json{{"eye", vector_to_json(test_case.camera->eye)},
                                            {"target", vector_to_json(test_case.camera->target)},
                                            {"up", vector_to_json(test_case.camera->up)}};
    }

    if (!test_case.color_ranges.empty()) {
        nlohmann::json ranges = nlohmann::json::object();
        for (const auto& [name, range] : test_case.color_ranges) {
            ranges[name] = nlohmann::json::array({range.low, range.high});
        }
        document["color_ranges"] = std::move(ranges);
    }

    return document;
}

Case case_from_json(const nlohmann::json& document) {
    const nlohmann::json& format = require(document, "format");
    if (!format.is_string() || format.get<std::string>() != "n2s-case") {
        throw std::runtime_error(
            fmt::format("expected \"format\": \"n2s-case\", got {}", format.dump()));
    }
    const nlohmann::json& version = require(document, "version");
    if (!version.is_number_integer() || version.get<int>() != kCaseJsonVersion) {
        throw std::runtime_error(
            fmt::format("unsupported case schema version {}; this build reads version {}",
                        version.dump(),
                        kCaseJsonVersion));
    }

    const nlohmann::json& trim = require(document, "trim");

    std::vector<TrimLoop> holes;
    if (const auto found = trim.find("holes"); found != trim.end()) {
        for (std::size_t i = 0; i < found->size(); ++i) {
            holes.push_back(loop_from_json((*found)[i], fmt::format("trim.holes[{}]", i)));
        }
    }

    Case test_case{
        .name = require(document, "name").get<std::string>(),
        .surface = surface_from_json(require(document, "surface")),
        .region =
            TrimRegion{loop_from_json(require(trim, "outer"), "trim.outer"), std::move(holes)},
        .layout = std::nullopt,
        .control_mesh = std::nullopt,
        .camera = std::nullopt,
        .color_ranges = {},
    };

    if (const auto found = document.find("layout"); found != document.end()) {
        test_case.layout = layout_from_json(*found);
    }

    if (const auto found = document.find("control_mesh"); found != document.end()) {
        test_case.control_mesh = mesh_from_json(*found);
    }

    if (const auto found = document.find("camera"); found != document.end()) {
        CameraPose camera;
        camera.eye = vector_from_json(require(*found, "eye"), "camera.eye");
        camera.target = vector_from_json(require(*found, "target"), "camera.target");
        camera.up = vector_from_json(require(*found, "up"), "camera.up");
        test_case.camera = camera;
    }

    if (const auto found = document.find("color_ranges"); found != document.end()) {
        for (const auto& [name, value] : found->items()) {
            if (!value.is_array() || value.size() != 2) {
                throw std::runtime_error(
                    fmt::format("color_ranges.{} must be [low, high], got {}", name, value.dump()));
            }
            const ScalarRange range{value[0].get<double>(), value[1].get<double>()};
            if (!(range.low < range.high)) {
                throw std::runtime_error(
                    fmt::format("color_ranges.{} is [{}, {}], which is empty or inverted; a "
                                "colour range has to be a real interval for a figure to mean "
                                "anything",
                                name,
                                range.low,
                                range.high));
            }
            test_case.color_ranges.emplace(name, range);
        }
    }

    return test_case;
}

void write_case(const std::filesystem::path& path, const Case& test_case) {
    std::ofstream stream(path);
    if (!stream) {
        throw std::runtime_error(fmt::format("cannot open {} for writing", path.string()));
    }
    stream << to_json(test_case).dump(2) << '\n';
}

Case read_case(const std::filesystem::path& path) {
    std::ifstream stream(path);
    if (!stream) {
        throw std::runtime_error(fmt::format("cannot open {} for reading", path.string()));
    }

    try {
        return case_from_json(nlohmann::json::parse(stream));
    } catch (const nlohmann::json::parse_error& error) {
        throw std::runtime_error(fmt::format("{}: {}", path.string(), error.what()));
    } catch (const std::runtime_error& error) {
        throw std::runtime_error(fmt::format("{}: {}", path.string(), error.what()));
    }
}

} // namespace n2s::io
