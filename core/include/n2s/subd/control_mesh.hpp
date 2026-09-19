#pragma once

#include <Eigen/Core>

#include <array>
#include <cstddef>
#include <filesystem>
#include <vector>

namespace n2s {

/// A quad control mesh for Catmull-Clark subdivision.
///
/// Quads only. Every layout this project produces is a quad layout, and
/// restricting the type means the face-index-to-ptex-face correspondence that
/// the limit evaluation depends on is the identity rather than a mapping that
/// has to be maintained and can silently drift.
///
/// Sharpness values follow the usual convention: 0 is smooth, and values at or
/// above 10 are treated as infinitely sharp by OpenSubdiv.
class ControlMesh {
public:
    /// Throws `std::invalid_argument` on an empty mesh or a face index out of
    /// range.
    ControlMesh(std::vector<Eigen::Vector3d> vertices, std::vector<std::array<int, 4>> quads);

    /// Regular `rows x columns` grid of vertices in the z = 0 plane spanning
    /// [0, 1]^2, with `(rows - 1) * (columns - 1)` quads. The workhorse of the
    /// M3 oracles: its interior vertices are all valence 4, so its limit
    /// surface is an ordinary uniform bicubic B-spline.
    static ControlMesh grid(int rows, int columns);

    /// Cube: eight vertices, six quads, every vertex valence 3. The smallest
    /// closed mesh with extraordinary vertices.
    static ControlMesh cube(double half_size = 0.5);

    /// A fan of `valence` quads around a single interior vertex, which is
    /// vertex 0. Vertices `1 .. valence` are its edge-adjacent 1-ring in
    /// rotational order, and `valence+1 .. 2*valence` are the face diagonals,
    /// diagonal `i` belonging to the quad between edge vertices `i` and `i+1`.
    ///
    /// This is the configuration every extraordinary-vertex oracle is stated
    /// over, and the one R5 measures curvature around. The index ordering is
    /// part of the contract, because the limit weight formulas are stated in
    /// terms of it.
    static ControlMesh vertex_fan(int valence);

    const std::vector<Eigen::Vector3d>& vertices() const noexcept { return vertices_; }

    const std::vector<std::array<int, 4>>& quads() const noexcept { return quads_; }

    std::size_t num_vertices() const noexcept { return vertices_.size(); }

    std::size_t num_quads() const noexcept { return quads_.size(); }

    /// Replaces the vertex positions, keeping the topology. This is what a
    /// fitting solve writes back: the topology is fixed by the layout and only
    /// the positions are solved for.
    void set_vertices(std::vector<Eigen::Vector3d> vertices);

    /// Number of faces meeting at a vertex. For an interior vertex of a quad
    /// mesh this is its valence; on a boundary it is one less.
    int face_count(int vertex) const;

    /// True if the vertex lies on a boundary, i.e. some incident edge is used
    /// by only one face.
    bool is_boundary_vertex(int vertex) const;

    /// Tags an edge as a crease. Order of the two vertices does not matter.
    /// Throws if the vertices do not span an edge of the mesh.
    void set_crease(int v0, int v1, double sharpness);

    /// Tags a vertex as a corner.
    void set_corner(int vertex, double sharpness);

    struct Crease {
        int v0;
        int v1;
        double sharpness;
    };

    struct Corner {
        int vertex;
        double sharpness;
    };

    const std::vector<Crease>& creases() const noexcept { return creases_; }

    const std::vector<Corner>& corners() const noexcept { return corners_; }

    /// Writes a Wavefront OBJ. Text format, 1-based indices.
    void write_obj(const std::filesystem::path& path) const;

private:
    void validate() const;

    std::vector<Eigen::Vector3d> vertices_;
    std::vector<std::array<int, 4>> quads_;
    std::vector<Crease> creases_;
    std::vector<Corner> corners_;
};

/// A triangle or quad mesh with no subdivision semantics, used for display and
/// for OBJ export of refined output.
struct PolyMesh {
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::array<int, 4>> quads;

    void write_obj(const std::filesystem::path& path) const;
};

} // namespace n2s
