#include "scene.hpp"

#include "n2s/fit/layout.hpp"
#include "n2s/nurbs/differential.hpp"
#include "n2s/nurbs/projection.hpp"
#include "n2s/trim/validate.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace n2s::viewer {

namespace {

/// Looks the range up in the case file. A field with no entry there is still
/// drawn, but marked as not comparable, and the viewer will not screenshot
/// while one is visible.
ScalarField make_field(const std::string& name,
                       std::vector<double> values,
                       const io::Case& test_case,
                       std::vector<std::string>& notes) {
    ScalarField field;
    field.name = name;
    field.values = std::move(values);

    const auto found = test_case.color_ranges.find(name);
    if (found != test_case.color_ranges.end()) {
        field.range = found->second;
        field.explicit_range = true;
        return field;
    }

    double low = std::numeric_limits<double>::infinity();
    double high = -std::numeric_limits<double>::infinity();
    for (const double value : field.values) {
        if (std::isfinite(value)) {
            low = std::min(low, value);
            high = std::max(high, value);
        }
    }
    if (!(low < high)) {
        low = 0.0;
        high = 1.0;
    }

    field.range = io::ScalarRange{low, high};
    field.explicit_range = false;
    notes.push_back(fmt::format(
        "\"{}\" has no range in the case file, so it is shown auto-ranged at [{:.4g}, {:.4g}]. "
        "Not comparable with any other case; add \"color_ranges\": {{\"{}\": [low, high]}} "
        "before using it in a figure.",
        name,
        low,
        high,
        name));
    return field;
}

/// Mean and Gaussian curvature of the NURBS at a list of domain points.
/// Degenerate points report zero, which is wrong but bounded; the alternative
/// is a NaN that poisons the colour range for the whole layer.
void nurbs_curvature_fields(const NurbsSurface& surface,
                            const std::vector<Eigen::Vector2d>& domain_points,
                            std::vector<double>& mean,
                            std::vector<double>& gaussian) {
    mean.clear();
    gaussian.clear();
    mean.reserve(domain_points.size());
    gaussian.reserve(domain_points.size());

    for (const Eigen::Vector2d& p : domain_points) {
        const SurfaceDerivatives derivatives = surface.derivatives(p.x(), p.y(), 2);
        const std::optional<SurfaceCurvature> curvature = surface_curvature(derivatives);
        mean.push_back(curvature.has_value() ? curvature->mean : 0.0);
        gaussian.push_back(curvature.has_value() ? curvature->gaussian : 0.0);
    }
}

double bounding_width(const std::vector<Eigen::Vector3d>& points) {
    if (points.empty()) {
        return 1.0;
    }
    Eigen::Vector3d low = points.front();
    Eigen::Vector3d high = points.front();
    for (const Eigen::Vector3d& p : points) {
        low = low.cwiseMin(p);
        high = high.cwiseMax(p);
    }
    return std::max(1e-6, (high - low).maxCoeff());
}

} // namespace

bool Scene::all_ranges_explicit() const {
    const auto explicit_everywhere = [](const std::vector<ScalarField>& fields) {
        return std::all_of(
            fields.begin(), fields.end(), [](const ScalarField& f) { return f.explicit_range; });
    };
    return explicit_everywhere(nurbs_scalars) && explicit_everywhere(trimmed_scalars) &&
           explicit_everywhere(limit_scalars);
}

Scene build_scene(const io::Case& test_case, const SceneOptions& options) {
    Scene scene;

    // ---- The base NURBS surface, untrimmed ----------------------------------
    const int steps = std::max(2, options.nurbs_samples);
    std::vector<Eigen::Vector2d> nurbs_domain_points;
    nurbs_domain_points.reserve(static_cast<std::size_t>((steps + 1) * (steps + 1)));

    for (int i = 0; i <= steps; ++i) {
        for (int j = 0; j <= steps; ++j) {
            const double u = static_cast<double>(i) / static_cast<double>(steps);
            const double v = static_cast<double>(j) / static_cast<double>(steps);
            nurbs_domain_points.emplace_back(u, v);
            scene.nurbs.vertices.push_back(test_case.surface.evaluate(u, v));
        }
    }
    for (int i = 0; i < steps; ++i) {
        for (int j = 0; j < steps; ++j) {
            const int base = i * (steps + 1) + j;
            scene.nurbs.quads.push_back({base, base + steps + 1, base + steps + 2, base + 1});
        }
    }

    std::vector<double> mean;
    std::vector<double> gaussian;
    nurbs_curvature_fields(test_case.surface, nurbs_domain_points, mean, gaussian);
    scene.nurbs_scalars.push_back(
        make_field("mean_curvature", std::move(mean), test_case, scene.notes));
    scene.nurbs_scalars.push_back(
        make_field("gaussian_curvature", std::move(gaussian), test_case, scene.notes));

    // ---- NURBS control net --------------------------------------------------
    const std::size_t nu = test_case.surface.num_u();
    const std::size_t nv = test_case.surface.num_v();
    scene.control_net_nodes = test_case.surface.control_points();
    for (std::size_t i = 0; i < nu; ++i) {
        for (std::size_t j = 0; j < nv; ++j) {
            const auto here = static_cast<int>(i * nv + j);
            if (i + 1 < nu) {
                scene.control_net_edges.push_back({here, here + static_cast<int>(nv)});
            }
            if (j + 1 < nv) {
                scene.control_net_edges.push_back({here, here + 1});
            }
        }
    }

    // ---- The trimmed region -------------------------------------------------
    TrimRegion region = test_case.region;
    const TrimReport report = validate_and_repair(region);
    for (const std::string& repair : report.repairs) {
        scene.notes.push_back("trim repaired: " + repair);
    }
    for (const std::string& error : report.errors) {
        scene.notes.push_back("trim error: " + error);
    }

    if (report.ok()) {
        scene.domain = triangulate_refined(
            region, test_case.surface, options.trim_sampling, options.refinement);
        scene.trimmed = map_to_surface(scene.domain, test_case.surface);

        std::vector<double> trimmed_mean;
        std::vector<double> trimmed_gaussian;
        nurbs_curvature_fields(
            test_case.surface, scene.domain.vertices, trimmed_mean, trimmed_gaussian);
        scene.trimmed_scalars.push_back(
            make_field("mean_curvature", std::move(trimmed_mean), test_case, scene.notes));
        scene.trimmed_scalars.push_back(
            make_field("gaussian_curvature", std::move(trimmed_gaussian), test_case, scene.notes));

        // The flat domain view, placed beside the model.
        const double width = bounding_width(scene.nurbs.vertices);
        const double shift = options.domain_offset * width;
        scene.domain_flat.quads.clear();
        for (const Eigen::Vector2d& p : scene.domain.vertices) {
            scene.domain_flat.vertices.emplace_back(shift + p.x() * width, p.y() * width, 0.0);
        }
        // DomainMesh is triangles; PolyMesh stores quads, so each triangle is
        // sent as a quad with a repeated last corner. Polyscope renders that as
        // the triangle it is.
        for (const std::array<std::size_t, 3>& t : scene.domain.triangles) {
            scene.domain_flat.quads.push_back({static_cast<int>(t[0]),
                                               static_cast<int>(t[1]),
                                               static_cast<int>(t[2]),
                                               static_cast<int>(t[2])});
        }

        // Trim loops, on the surface and in the flat view.
        const auto add_loop = [&](const TrimLoop& loop) {
            const std::vector<Eigen::Vector2d> points =
                sample_loop(loop, test_case.surface, options.trim_sampling);

            Polyline on_surface;
            Polyline in_domain;
            on_surface.points.reserve(points.size());
            in_domain.points.reserve(points.size());
            for (const Eigen::Vector2d& p : points) {
                on_surface.points.push_back(test_case.surface.evaluate(p.x(), p.y()));
                in_domain.points.emplace_back(shift + p.x() * width, p.y() * width, 0.0);
            }
            scene.trim_on_surface.push_back(std::move(on_surface));
            scene.trim_in_domain.push_back(std::move(in_domain));
        };

        add_loop(region.outer());
        for (const TrimLoop& hole : region.holes()) {
            add_loop(hole);
        }
    } else {
        scene.notes.emplace_back(
            "the trim region did not validate, so the triangulation and every layer built on "
            "it were skipped");
    }

    // ---- Subdivision layers -------------------------------------------------
    const std::optional<fit::ResolvedLayout> resolved =
        fit::resolve_layout(test_case.layout, test_case.control_mesh, test_case.surface);

    if (resolved.has_value()) {
        scene.has_subdivision = true;
        const ControlMesh& control = resolved->mesh;
        const SubdivisionSurface subdivision{control};

        scene.control_mesh = PolyMesh{control.vertices(), control.quads()};

        const TessellatedLimit tessellation =
            tessellate_limit(subdivision, options.limit_samples_per_face);
        scene.limit = tessellation.mesh;

        // Isophotes: bands of constant angle to a fixed light direction. Fixed,
        // not view-dependent, so the same picture comes back every time --
        // which a screen-space zebra does not.
        const Eigen::Vector3d light = options.isophote_light.normalized();
        const auto bands = static_cast<double>(std::max(1, options.isophote_bands));
        std::vector<double> isophote;
        isophote.reserve(tessellation.normals.size());
        for (const Eigen::Vector3d& normal : tessellation.normals) {
            const double shade = 0.5 * (normal.dot(light) + 1.0);
            isophote.push_back(shade * bands - std::floor(shade * bands));
        }
        // An isophote stripe is a phase, so its range is [0, 1] by
        // construction rather than by measurement: no case file entry needed.
        ScalarField stripes;
        stripes.name = "isophote";
        stripes.values = std::move(isophote);
        stripes.range = io::ScalarRange{0.0, 1.0};
        stripes.explicit_range = true;
        scene.limit_scalars.push_back(std::move(stripes));

        // Extraordinary vertices: interior vertices whose valence is not 4.
        for (std::size_t i = 0; i < control.num_vertices(); ++i) {
            const auto vertex = static_cast<int>(i);
            const int faces = control.face_count(vertex);
            const bool boundary = control.is_boundary_vertex(vertex);
            if (!boundary && faces != 4) {
                scene.extraordinary_vertices.push_back(control.vertices()[i]);
                scene.extraordinary_valences.push_back(static_cast<double>(faces));
            }
        }
    }

    return scene;
}

void add_error_field(Scene& scene, const io::Case& test_case, const SceneOptions& options) {
    if (!scene.has_subdivision || scene.limit.vertices.empty()) {
        scene.notes.emplace_back("no limit surface to measure: the case ships no control mesh");
        return;
    }

    // Remove any previous error field so repeated presses do not stack up.
    scene.limit_scalars.erase(
        std::remove_if(scene.limit_scalars.begin(),
                       scene.limit_scalars.end(),
                       [](const ScalarField& f) { return f.name == "error"; }),
        scene.limit_scalars.end());

    ProjectionOptions projection;
    projection.samples_per_span_u = 6;
    projection.samples_per_span_v = 6;

    std::size_t unconverged = 0;
    std::vector<double> distances;
    distances.reserve(scene.limit.vertices.size());

    for (const Eigen::Vector3d& point : scene.limit.vertices) {
        const SurfaceProjection found = project_to_surface(test_case.surface, point, projection);
        if (!found.converged()) {
            ++unconverged;
        }
        distances.push_back(found.distance);
    }

    if (unconverged > 0) {
        // Reported rather than swallowed: an unconverged projection contributes
        // a distance that is an upper bound at best, and a quietly optimistic
        // error map is worse than none.
        scene.notes.push_back(
            fmt::format("{} of {} projections did not converge; those distances are unreliable",
                        unconverged,
                        scene.limit.vertices.size()));
    }

    scene.limit_scalars.push_back(
        make_field("error", std::move(distances), test_case, scene.notes));

    (void)options;
}

} // namespace n2s::viewer
