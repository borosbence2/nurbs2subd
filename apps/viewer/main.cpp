// Polyscope viewer for nurbs2subd.
//
// The pipeline stages and every layer the plan asks for are here; the geometry
// behind them lives in scene.cpp, which knows nothing about Polyscope so that
// it can be exercised headlessly. `--selftest` does exactly that and exits,
// which is what CI runs, because nothing else about a viewer can be checked
// without a human looking at it.

#include "n2s/build_info.hpp"
#include "n2s/fit/layout.hpp"
#include "n2s/io/case_json.hpp"
#include "n2s/trim/cases.hpp"
#include "n2s/trim/validate.hpp"

#include "scene.hpp"

#include <fmt/core.h>
#include <imgui.h>
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/view.h>

#include <array>
#include <exception>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace {

using n2s::viewer::Scene;
using n2s::viewer::SceneOptions;

/// Everything the ImGui panel owns.
struct Viewer {
    std::string case_path;
    std::optional<n2s::io::Case> loaded;
    Scene scene;
    SceneOptions options;

    std::string status = "No case loaded.";
    std::string screenshot_path = "figure.png";
    bool have_scene = false;
};

Viewer g_viewer;

std::vector<glm::vec3> to_glm(const std::vector<Eigen::Vector3d>& points) {
    std::vector<glm::vec3> result;
    result.reserve(points.size());
    for (const Eigen::Vector3d& p : points) {
        result.emplace_back(
            static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()));
    }
    return result;
}

std::vector<std::array<size_t, 4>> to_faces(const std::vector<std::array<int, 4>>& quads) {
    std::vector<std::array<size_t, 4>> result;
    result.reserve(quads.size());
    for (const std::array<int, 4>& q : quads) {
        result.push_back({static_cast<size_t>(q[0]),
                          static_cast<size_t>(q[1]),
                          static_cast<size_t>(q[2]),
                          static_cast<size_t>(q[3])});
    }
    return result;
}

std::vector<std::array<size_t, 3>> to_faces(const std::vector<std::array<std::size_t, 3>>& tris) {
    std::vector<std::array<size_t, 3>> result;
    result.reserve(tris.size());
    for (const std::array<std::size_t, 3>& t : tris) {
        result.push_back({t[0], t[1], t[2]});
    }
    return result;
}

/// Attaches a scalar with its range pinned. Never auto-ranged: `setMapRange` is
/// always called, with the range the scene decided, and the numbers are shown
/// in the panel so a reader of the figure can see them.
template<typename Structure>
void attach_scalars(Structure* structure, const std::vector<n2s::viewer::ScalarField>& fields) {
    for (const n2s::viewer::ScalarField& field : fields) {
        std::vector<double> values = field.values;
        auto* quantity = structure->addVertexScalarQuantity(field.name, values);
        quantity->setMapRange({field.range.low, field.range.high});
        quantity->setColorMap(field.name == "isophote" ? "phase" : "coolwarm");
    }
}

void register_scene(const Scene& scene) {
    polyscope::removeAllStructures();

    if (!scene.nurbs.vertices.empty()) {
        auto* nurbs = polyscope::registerSurfaceMesh(
            "NURBS surface", to_glm(scene.nurbs.vertices), to_faces(scene.nurbs.quads));
        attach_scalars(nurbs, scene.nurbs_scalars);
        nurbs->setEnabled(true);
    }

    if (!scene.control_net_nodes.empty()) {
        std::vector<std::array<size_t, 2>> edges;
        edges.reserve(scene.control_net_edges.size());
        for (const std::array<int, 2>& e : scene.control_net_edges) {
            edges.push_back({static_cast<size_t>(e[0]), static_cast<size_t>(e[1])});
        }
        auto* net = polyscope::registerCurveNetwork(
            "NURBS control net", to_glm(scene.control_net_nodes), edges);
        net->setEnabled(false);
    }

    if (!scene.trimmed.vertices.empty()) {
        auto* trimmed = polyscope::registerSurfaceMesh(
            "trimmed surface", to_glm(scene.trimmed.vertices), to_faces(scene.trimmed.triangles));
        attach_scalars(trimmed, scene.trimmed_scalars);
        trimmed->setEnabled(false);
    }

    if (!scene.domain_flat.vertices.empty()) {
        auto* domain = polyscope::registerSurfaceMesh("domain triangulation",
                                                      to_glm(scene.domain_flat.vertices),
                                                      to_faces(scene.domain_flat.quads));
        attach_scalars(domain, scene.trimmed_scalars);
        domain->setEnabled(true);
    }

    for (std::size_t i = 0; i < scene.trim_on_surface.size(); ++i) {
        const std::string label = i == 0 ? "trim outer" : fmt::format("trim hole {}", i - 1);
        polyscope::registerCurveNetworkLoop(label + " (surface)",
                                            to_glm(scene.trim_on_surface[i].points));
        polyscope::registerCurveNetworkLoop(label + " (domain)",
                                            to_glm(scene.trim_in_domain[i].points));
    }

    if (scene.has_subdivision) {
        auto* control = polyscope::registerSurfaceMesh("CC control mesh",
                                                       to_glm(scene.control_mesh.vertices),
                                                       to_faces(scene.control_mesh.quads));
        control->setEdgeWidth(1.0);
        control->setEnabled(false);

        auto* limit = polyscope::registerSurfaceMesh(
            "limit surface", to_glm(scene.limit.vertices), to_faces(scene.limit.quads));
        attach_scalars(limit, scene.limit_scalars);
        limit->setEnabled(true);

        if (!scene.extraordinary_vertices.empty()) {
            auto* markers = polyscope::registerPointCloud("extraordinary vertices",
                                                          to_glm(scene.extraordinary_vertices));
            std::vector<double> valences = scene.extraordinary_valences;
            auto* quantity = markers->addScalarQuantity("valence", valences);
            // Valence 3 to 8 covers every layout this project produces; a fixed
            // range keeps the colours meaning the same thing case to case.
            quantity->setMapRange({3.0, 8.0});
            quantity->setColorMap("spectral");
            quantity->setEnabled(true);
        }
    }
}

void rebuild(const n2s::io::Case& test_case) {
    g_viewer.scene = n2s::viewer::build_scene(test_case, g_viewer.options);
    g_viewer.have_scene = true;
    register_scene(g_viewer.scene);

    if (test_case.camera.has_value()) {
        const n2s::io::CameraPose& camera = *test_case.camera;
        polyscope::view::lookAt(glm::vec3{static_cast<float>(camera.eye.x()),
                                          static_cast<float>(camera.eye.y()),
                                          static_cast<float>(camera.eye.z())},
                                glm::vec3{static_cast<float>(camera.target.x()),
                                          static_cast<float>(camera.target.y()),
                                          static_cast<float>(camera.target.z())},
                                glm::vec3{static_cast<float>(camera.up.x()),
                                          static_cast<float>(camera.up.y()),
                                          static_cast<float>(camera.up.z())});
    }
}

void load_case(const std::string& path) {
    try {
        g_viewer.loaded = n2s::io::read_case(path);
        rebuild(*g_viewer.loaded);
        g_viewer.status = fmt::format("Loaded \"{}\": {} limit vertices, {} domain triangles.",
                                      g_viewer.loaded->name,
                                      g_viewer.scene.limit.vertices.size(),
                                      g_viewer.scene.domain.triangles.size());
    } catch (const std::exception& error) {
        g_viewer.loaded.reset();
        g_viewer.have_scene = false;
        g_viewer.status = fmt::format("Failed to load: {}", error.what());
    }
}

void draw_panel() {
    ImGui::PushItemWidth(260);

    ImGui::TextUnformatted("Case");
    ImGui::InputText("path", &g_viewer.case_path[0], g_viewer.case_path.capacity());
    if (ImGui::Button("Load case")) {
        load_case(g_viewer.case_path.c_str());
    }
    ImGui::SameLine();
    if (ImGui::Button("Load built-in saddle case")) {
        n2s::cases::TrimmedCase generated =
            n2s::cases::square_with_circular_hole(n2s::cases::saddle(1.0, 0.6), 0.25);
        n2s::validate_and_repair(generated.region);

        n2s::io::Case built_in{
            .name = generated.name,
            .surface = generated.surface,
            .region = generated.region,
            .layout = n2s::fit::grid_domain_layout(6, 6),
            .control_mesh = std::nullopt,
            .camera = std::nullopt,
            .color_ranges = {},
        };
        g_viewer.loaded = std::move(built_in);
        rebuild(*g_viewer.loaded);
        g_viewer.status = "Loaded the built-in saddle case.";
    }

    ImGui::Separator();
    ImGui::TextWrapped("%s", g_viewer.status.c_str());

    if (g_viewer.loaded.has_value()) {
        ImGui::Separator();
        ImGui::TextUnformatted("Pipeline");

        bool dirty = false;
        dirty |= ImGui::SliderInt("NURBS samples", &g_viewer.options.nurbs_samples, 8, 128);
        dirty |=
            ImGui::SliderInt("limit samples/face", &g_viewer.options.limit_samples_per_face, 1, 24);
        dirty |= ImGui::SliderInt("isophote bands", &g_viewer.options.isophote_bands, 2, 40);
        dirty |= ImGui::InputDouble(
            "trim max sagitta", &g_viewer.options.trim_sampling.max_sagitta, 0.0, 0.0, "%.5f");
        dirty |= ImGui::InputDouble(
            "max triangle area", &g_viewer.options.refinement.max_triangle_area, 0.0, 0.0, "%.5f");
        dirty |=
            ImGui::InputDouble("domain offset", &g_viewer.options.domain_offset, 0.0, 0.0, "%.2f");

        if (ImGui::Button("Rebuild") || dirty) {
            rebuild(*g_viewer.loaded);
        }

        if (ImGui::Button("Compute error to NURBS (slow)")) {
            n2s::viewer::add_error_field(g_viewer.scene, *g_viewer.loaded, g_viewer.options);
            register_scene(g_viewer.scene);
            g_viewer.status = "Error field computed.";
        }
    }

    if (g_viewer.have_scene) {
        ImGui::Separator();
        ImGui::TextUnformatted("Colour ranges");
        ImGui::TextWrapped("Editing a range pins it: the field counts as explicitly ranged "
                           "from then on, because you chose the numbers. Copy them into the "
                           "case file to keep them.");

        bool ranges_changed = false;
        const auto edit = [&](const char* group, std::vector<n2s::viewer::ScalarField>& fields) {
            for (std::size_t i = 0; i < fields.size(); ++i) {
                n2s::viewer::ScalarField& field = fields[i];
                ImGui::PushID(static_cast<int>(i) + (group[0] << 8));

                std::array<double, 2> bounds{field.range.low, field.range.high};
                const std::string label = fmt::format("{} / {}", group, field.name);
                if (ImGui::InputScalarN(label.c_str(),
                                        ImGuiDataType_Double,
                                        bounds.data(),
                                        2,
                                        nullptr,
                                        nullptr,
                                        "%.5g")) {
                    if (bounds[0] < bounds[1]) {
                        field.range.low = bounds[0];
                        field.range.high = bounds[1];
                        field.explicit_range = true;
                        ranges_changed = true;
                    }
                }
                ImGui::SameLine();
                ImGui::TextUnformatted(field.explicit_range ? "pinned" : "AUTO");

                ImGui::PopID();
            }
        };
        edit("nurbs", g_viewer.scene.nurbs_scalars);
        edit("trimmed", g_viewer.scene.trimmed_scalars);
        edit("limit", g_viewer.scene.limit_scalars);

        if (ranges_changed) {
            // Re-registering is cheap: the geometry is already computed and
            // only the map ranges change.
            register_scene(g_viewer.scene);
        }

        ImGui::Separator();
        ImGui::TextUnformatted("Figure export");
        ImGui::InputText("png", &g_viewer.screenshot_path[0], g_viewer.screenshot_path.capacity());

        const bool comparable = g_viewer.scene.all_ranges_explicit();
        if (!comparable) {
            // The enforcement point for the no-auto-ranging rule. Exploring
            // with a computed range is fine; exporting one as a figure is how
            // defect 6 of the thesis happened.
            ImGui::TextWrapped("Screenshot disabled: at least one scalar is auto-ranged. Put "
                               "its range in the case file to make the figure comparable.");
        }
        ImGui::BeginDisabled(!comparable);
        if (ImGui::Button("Screenshot")) {
            polyscope::screenshot(g_viewer.screenshot_path.c_str());
            g_viewer.status = fmt::format("Wrote {}", g_viewer.screenshot_path.c_str());
        }
        ImGui::EndDisabled();
    }

    if (!g_viewer.scene.notes.empty()) {
        ImGui::Separator();
        ImGui::TextUnformatted("Notes");
        for (const std::string& note : g_viewer.scene.notes) {
            ImGui::Bullet();
            ImGui::TextWrapped("%s", note.c_str());
        }
    }

    ImGui::PopItemWidth();
}

/// Builds a scene headlessly and reports what came out. Gives CI something real
/// to check on a component that otherwise only a human can judge.
int selftest() {
    n2s::cases::TrimmedCase generated =
        n2s::cases::square_with_circular_hole(n2s::cases::saddle(1.0, 0.6), 0.25);
    const n2s::TrimReport report = n2s::validate_and_repair(generated.region);
    if (!report.ok()) {
        fmt::print(stderr, "selftest: the built-in case did not validate\n{}", report.to_string());
        return 1;
    }

    const n2s::io::Case test_case{
        .name = generated.name,
        .surface = generated.surface,
        .region = generated.region,
        .layout = n2s::fit::grid_domain_layout(6, 6),
        .control_mesh = std::nullopt,
        .camera = std::nullopt,
        .color_ranges = {{"mean_curvature", {-2.0, 2.0}},
                         {"gaussian_curvature", {-4.0, 1.0}},
                         {"error", {0.0, 0.5}}},
    };

    SceneOptions options;
    options.nurbs_samples = 16;
    options.limit_samples_per_face = 3;

    const Scene scene = n2s::viewer::build_scene(test_case, options);

    if (scene.nurbs.quads.empty() || scene.domain.triangles.empty() ||
        scene.trim_on_surface.size() != 2 || !scene.has_subdivision ||
        scene.limit.vertices.empty()) {
        fmt::print(stderr, "selftest: a layer came back empty\n");
        return 1;
    }
    if (!scene.all_ranges_explicit()) {
        fmt::print(stderr,
                   "selftest: a scalar was auto-ranged despite the case supplying "
                   "every range\n");
        return 1;
    }

    fmt::print("selftest ok: {} nurbs quads, {} domain triangles, {} trim loops, {} limit "
               "vertices, {} scalars\n",
               scene.nurbs.quads.size(),
               scene.domain.triangles.size(),
               scene.trim_on_surface.size(),
               scene.limit.vertices.size(),
               scene.nurbs_scalars.size() + scene.trimmed_scalars.size() +
                   scene.limit_scalars.size());
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    std::string case_argument;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-v" || arg == "--version") {
            fmt::print("{}\n", n2s::build_info_string());
            return 0;
        }
        if (arg == "--selftest") {
            return selftest();
        }
        case_argument = arg;
    }

    g_viewer.case_path.reserve(512);
    g_viewer.case_path = case_argument;
    g_viewer.case_path.resize(512, '\0');
    g_viewer.screenshot_path.reserve(512);
    g_viewer.screenshot_path.resize(512, '\0');
    std::string("figure.png").copy(g_viewer.screenshot_path.data(), 10);

    polyscope::options::programName = "nurbs2subd viewer";
    polyscope::options::verbosity = 0;
    polyscope::init();

    polyscope::state::userCallback = draw_panel;

    if (!case_argument.empty()) {
        load_case(case_argument);
    }

    polyscope::show();
    return 0;
}
