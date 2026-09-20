#include "n2s/fit/layout.hpp"

#include <fmt/format.h>

#include <stdexcept>
#include <vector>

namespace n2s::fit {

namespace {

void validate(const DomainLayout& layout) {
    if (layout.vertices.empty() || layout.quads.empty()) {
        throw std::invalid_argument("a domain layout needs at least one quad");
    }
    const auto count = static_cast<int>(layout.vertices.size());
    for (std::size_t f = 0; f < layout.quads.size(); ++f) {
        for (const int index : layout.quads[f]) {
            if (index < 0 || index >= count) {
                throw std::invalid_argument(fmt::format(
                    "layout quad {} references vertex {}, outside [0, {})", f, index, count));
            }
        }
    }
}

} // namespace

ControlMesh lift(const DomainLayout& layout, const NurbsSurface& surface) {
    validate(layout);

    std::vector<Eigen::Vector3d> points;
    points.reserve(layout.vertices.size());
    for (const Eigen::Vector2d& p : layout.vertices) {
        points.push_back(surface.evaluate(p.x(), p.y()));
    }

    return ControlMesh{std::move(points), layout.quads};
}

metrics::DomainMap bilinear_domain_map(const DomainLayout& layout) {
    validate(layout);

    // Captured by value: the returned map outlives the call, and a layout that
    // changed underneath it would silently start reporting a correspondence
    // that no longer matches the control mesh.
    return [vertices = layout.vertices,
            quads = layout.quads](const LimitLocation& location) -> Eigen::Vector2d {
        const auto face = static_cast<std::size_t>(location.face);
        if (location.face < 0 || face >= quads.size()) {
            throw std::invalid_argument(fmt::format(
                "limit location names face {}, outside [0, {})", location.face, quads.size()));
        }

        const std::array<int, 4>& quad = quads[face];
        const Eigen::Vector2d& a = vertices[static_cast<std::size_t>(quad[0])];
        const Eigen::Vector2d& b = vertices[static_cast<std::size_t>(quad[1])];
        const Eigen::Vector2d& c = vertices[static_cast<std::size_t>(quad[2])];
        const Eigen::Vector2d& d = vertices[static_cast<std::size_t>(quad[3])];

        // OpenSubdiv parameterises a quad face with (0,0) at its first vertex,
        // u running toward the second and v toward the fourth. The blend below
        // has to match that or the correspondence is transposed.
        const double u = location.u;
        const double v = location.v;
        return (1.0 - u) * (1.0 - v) * a + u * (1.0 - v) * b + u * v * c + (1.0 - u) * v * d;
    };
}

std::optional<ResolvedLayout> resolve_layout(const std::optional<DomainLayout>& layout,
                                             const std::optional<ControlMesh>& mesh,
                                             const NurbsSurface& surface) {
    if (layout.has_value()) {
        return ResolvedLayout{lift(*layout, surface), bilinear_domain_map(*layout), true};
    }
    if (mesh.has_value()) {
        return ResolvedLayout{*mesh, metrics::DomainMap{}, false};
    }
    return std::nullopt;
}

DomainLayout grid_domain_layout(int rows, int columns) {
    if (rows < 2 || columns < 2) {
        throw std::invalid_argument(
            fmt::format("a grid layout needs at least 2 x 2 vertices, got {} x {}", rows, columns));
    }

    DomainLayout layout;
    layout.vertices.reserve(static_cast<std::size_t>(rows * columns));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            layout.vertices.emplace_back(static_cast<double>(i) / static_cast<double>(rows - 1),
                                         static_cast<double>(j) / static_cast<double>(columns - 1));
        }
    }

    // Same face order and winding as ControlMesh::grid, so the two agree.
    layout.quads.reserve(static_cast<std::size_t>((rows - 1) * (columns - 1)));
    for (int i = 0; i + 1 < rows; ++i) {
        for (int j = 0; j + 1 < columns; ++j) {
            const int base = i * columns + j;
            layout.quads.push_back({base, base + columns, base + columns + 1, base + 1});
        }
    }

    return layout;
}

ControlMesh grid_layout(const NurbsSurface& surface, int rows, int columns) {
    return lift(grid_domain_layout(rows, columns), surface);
}

} // namespace n2s::fit
