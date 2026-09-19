#pragma once

#include "n2s/io/case_json.hpp"
#include "n2s/subd/subdivision.hpp"
#include "n2s/trim/sampling.hpp"
#include "n2s/trim/triangulation.hpp"

#include <Eigen/Core>

#include <optional>
#include <string>
#include <vector>

/// Everything the viewer draws, computed from a case with no Polyscope
/// involved. Keeping the two apart means this half can be exercised headlessly
/// -- see `--selftest` -- which matters for a component whose output is
/// otherwise only checkable by eye.
namespace n2s::viewer {

struct SceneOptions {
    /// Grid resolution for the untrimmed NURBS surface layer.
    int nurbs_samples = 48;

    /// Limit surface samples per control mesh face edge.
    int limit_samples_per_face = 8;

    SamplingOptions trim_sampling;

    /// Lighter than the core defaults on purpose. The core defaults are sized
    /// for measurement, and at those settings a single rebuild here produces
    /// tens of thousands of triangles and the viewer stops feeling
    /// interactive. Measurement runs set their own values from a config.
    RefinementOptions refinement{
        .max_triangle_area = 0.02,
        .min_angle_degrees = 25.0,
        .max_vertices = 4000,
        .max_rounds = 12,
        .separation_fraction = 0.5,
    };

    /// Where the flat domain view sits relative to the model, as a multiple of
    /// the model's bounding box width. Polyscope draws one 3D scene, so the
    /// plan's "side-by-side domain and surface views" is done by placing the
    /// domain beside the model rather than in a second window.
    double domain_offset = 1.4;

    int isophote_bands = 12;
    Eigen::Vector3d isophote_light{0.35, 0.25, 1.0};
};

/// A scalar to be drawn, with the colour range it must be drawn at.
///
/// `explicit_range` is the crux. CLAUDE.md forbids auto-ranging in figures
/// meant for comparison, and defect 6 of the thesis was curvature plots with no
/// fixed scale. A field whose range came from the case file is comparable
/// across cases; one whose range was computed from its own values is not, and
/// says so, and the viewer refuses to export a screenshot while any visible
/// field is in that state.
struct ScalarField {
    std::string name;
    std::vector<double> values;
    io::ScalarRange range;
    bool explicit_range = false;
};

/// A polyline, closed unless stated otherwise.
struct Polyline {
    std::vector<Eigen::Vector3d> points;
};

struct Scene {
    /// The base NURBS surface over its whole domain, ignoring the trim.
    PolyMesh nurbs;
    std::vector<ScalarField> nurbs_scalars;

    /// NURBS control net as a wireframe.
    std::vector<Eigen::Vector3d> control_net_nodes;
    std::vector<std::array<int, 2>> control_net_edges;

    /// The trimmed region: the domain triangulation and its image on the
    /// surface. `domain_flat` is the same triangulation placed beside the model
    /// for the side-by-side view.
    DomainMesh domain;
    SurfaceMesh trimmed;
    PolyMesh domain_flat;
    std::vector<ScalarField> trimmed_scalars;

    /// Trim loops, on the surface and in the flat domain view.
    std::vector<Polyline> trim_on_surface;
    std::vector<Polyline> trim_in_domain;

    /// Present only when the case ships a control mesh.
    bool has_subdivision = false;
    PolyMesh control_mesh;
    PolyMesh limit;
    std::vector<ScalarField> limit_scalars;
    std::vector<Eigen::Vector3d> extraordinary_vertices;
    std::vector<double> extraordinary_valences;

    /// Anything the user should know: repairs the trim region needed, scalars
    /// left without an explicit range, stages that could not run.
    std::vector<std::string> notes;

    /// True when every scalar field carries a range from the case file.
    bool all_ranges_explicit() const;
};

/// Builds everything except the error-to-NURBS field, which is left out because
/// it costs a closest-point projection per limit sample and the viewer should
/// come up promptly.
Scene build_scene(const io::Case& test_case, const SceneOptions& options = {});

/// Adds the "error" field to `scene.limit_scalars`: the distance from each
/// limit sample to its closest point on the trimmed NURBS surface. Slow, hence
/// separate and on demand.
void add_error_field(Scene& scene, const io::Case& test_case, const SceneOptions& options = {});

} // namespace n2s::viewer
