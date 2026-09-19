#pragma once

#include "n2s/nurbs/surface.hpp"
#include "n2s/trim/sampling.hpp"
#include "n2s/trim/trim_loop.hpp"

#include <Eigen/Core>

#include <array>
#include <cstddef>
#include <vector>

namespace n2s {

/// A triangulation of the trimmed part of a surface's parametric domain.
///
/// `boundary_vertex_count` is the number of leading vertices that came from the
/// sampled trim loops. They are stored first and never moved, so a caller can
/// tell a boundary vertex from an interior one by index alone -- which matters
/// later, because the boundary discretisation is what two neighbouring patches
/// have to agree on exactly for a watertight join.
struct DomainMesh {
    std::vector<Eigen::Vector2d> vertices;
    std::vector<std::array<std::size_t, 3>> triangles;
    std::size_t boundary_vertex_count = 0;

    /// Total area in the parametric domain.
    double area() const;

    bool is_boundary_vertex(std::size_t index) const { return index < boundary_vertex_count; }
};

/// The same mesh with its vertices lifted onto the surface.
struct SurfaceMesh {
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::array<std::size_t, 3>> triangles;
    std::size_t boundary_vertex_count = 0;

    /// Total model-space area, the sum of the triangle areas. Converges to the
    /// true area of the trimmed surface as the triangulation is refined, always
    /// from below for a convex patch.
    double area() const;
};

/// Constrained Delaunay triangulation of `region`, with the trim loops as
/// constraint edges and the holes removed.
///
/// The boundary is discretised by `sampling`, so the trim curves are respected
/// to the tolerances given there and nowhere else: this function adds no
/// vertices of its own. Interior refinement is a separate pass, because CDT
/// does not provide one.
///
/// Throws `std::runtime_error` if the triangulation comes back empty, which in
/// practice means the loops did not form closed boundaries.
DomainMesh triangulate(const TrimRegion& region,
                       const NurbsSurface& surface,
                       const SamplingOptions& sampling = {});

/// Lifts a domain mesh onto the surface: every vertex `(u, v)` becomes
/// `S(u, v)`. Triangle indices are untouched, so the two meshes stay in
/// correspondence and a domain point can be compared against its 3D image.
SurfaceMesh map_to_surface(const DomainMesh& mesh, const NurbsSurface& surface);

struct RefinementOptions {
    /// Upper bound on the model-space area of one triangle. Measured on the
    /// surface, not in the domain, so a stretched parameterisation gets the
    /// extra vertices it needs. Zero or less disables the area criterion.
    double max_triangle_area = 0.01;

    /// Smallest interior angle a triangle may have, in degrees, measured in the
    /// domain. Zero or less disables the shape criterion.
    double min_angle_degrees = 25.0;

    /// Hard budget on the total vertex count, and on how many refinement rounds
    /// run. Both are safety nets: this is a heuristic refinement, not a
    /// guaranteed-quality one, and a pathological case must stop rather than
    /// grind.
    ///
    /// `max_vertices` does double duty: it also sets the minimum spacing
    /// between inserted vertices, as `separation_fraction * sqrt(domain area /
    /// max_vertices)`. Without an absolute floor on that spacing the
    /// refinement degenerates -- a sliver's circumcentre lands almost on top of
    /// an existing vertex, which creates a thinner sliver, which does it again.
    /// Tying the floor to the budget means the two can never disagree: the
    /// spacing admits about `max_vertices` points and no more.
    std::size_t max_vertices = 20000;
    int max_rounds = 20;

    /// Fraction of the budget-implied spacing that two vertices must be apart.
    /// Lower packs more vertices in and tolerates thinner triangles.
    double separation_fraction = 0.5;
};

/// `triangulate` followed by interior refinement until the area and angle
/// criteria hold.
///
/// **Why this exists.** CDT ships no mesh refinement, so `triangulate` alone
/// returns a triangulation of the boundary polygon with no interior vertices
/// at all. On a curved surface its triangles chord straight across, which
/// underestimates area and leaves error sampling sparse in the middle of a
/// patch exactly where the fitting error is worth measuring.
///
/// **What it does not do.** Boundary vertices are never added, moved or
/// removed. The trim boundary discretisation is fixed by `sampling` and stays
/// that way, because two neighbouring patches have to agree on it exactly for
/// the joins to be watertight; a refinement that split boundary edges on its
/// own schedule would break that agreement. Candidate points that would crowd
/// the boundary are dropped instead. That costs the formal angle guarantee a
/// full Ruppert refinement would give -- hence the budgets above -- and in
/// exchange the boundary is inviolable.
DomainMesh triangulate_refined(const TrimRegion& region,
                               const NurbsSurface& surface,
                               const SamplingOptions& sampling = {},
                               const RefinementOptions& refinement = {});

} // namespace n2s
