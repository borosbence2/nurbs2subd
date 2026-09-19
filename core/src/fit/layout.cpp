#include "n2s/fit/layout.hpp"

#include <fmt/format.h>

#include <stdexcept>
#include <vector>

namespace n2s::fit {

ControlMesh grid_layout(const NurbsSurface& surface, int rows, int columns) {
    if (rows < 2 || columns < 2) {
        throw std::invalid_argument(
            fmt::format("a grid layout needs at least 2 x 2 vertices, got {} x {}", rows, columns));
    }

    // Topology from ControlMesh::grid so that the face ordering, and therefore
    // metrics::grid_domain_map, applies unchanged.
    ControlMesh mesh = ControlMesh::grid(rows, columns);

    std::vector<Eigen::Vector3d> vertices;
    vertices.reserve(static_cast<std::size_t>(rows * columns));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            const double u = static_cast<double>(i) / static_cast<double>(rows - 1);
            const double v = static_cast<double>(j) / static_cast<double>(columns - 1);
            vertices.push_back(surface.evaluate(u, v));
        }
    }

    mesh.set_vertices(std::move(vertices));
    return mesh;
}

} // namespace n2s::fit
