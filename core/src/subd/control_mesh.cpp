#include "n2s/subd/control_mesh.hpp"

#include <fmt/format.h>
#include <fmt/ostream.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <utility>

namespace n2s {

namespace {

std::pair<int, int> ordered_edge(int a, int b) {
    return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
}

void write_obj_file(const std::filesystem::path& path,
                    const std::vector<Eigen::Vector3d>& vertices,
                    const std::vector<std::array<int, 4>>& quads,
                    const char* what) {
    std::ofstream stream(path);
    if (!stream) {
        throw std::runtime_error(fmt::format("cannot open {} for writing", path.string()));
    }

    stream << "# nurbs2subd " << what << "\n";
    for (const Eigen::Vector3d& v : vertices) {
        stream << fmt::format("v {:.17g} {:.17g} {:.17g}\n", v.x(), v.y(), v.z());
    }
    for (const std::array<int, 4>& q : quads) {
        // OBJ indices are 1-based.
        stream << fmt::format("f {} {} {} {}\n", q[0] + 1, q[1] + 1, q[2] + 1, q[3] + 1);
    }
}

} // namespace

ControlMesh::ControlMesh(std::vector<Eigen::Vector3d> vertices,
                         std::vector<std::array<int, 4>> quads)
    : vertices_(std::move(vertices)),
      quads_(std::move(quads)) {
    validate();
}

void ControlMesh::validate() const {
    if (vertices_.empty()) {
        throw std::invalid_argument("a control mesh needs at least one vertex");
    }
    if (quads_.empty()) {
        throw std::invalid_argument("a control mesh needs at least one quad");
    }

    const auto count = static_cast<int>(vertices_.size());
    for (std::size_t f = 0; f < quads_.size(); ++f) {
        const std::array<int, 4>& q = quads_[f];
        for (int corner = 0; corner < 4; ++corner) {
            if (q[static_cast<std::size_t>(corner)] < 0 ||
                q[static_cast<std::size_t>(corner)] >= count) {
                throw std::invalid_argument(
                    fmt::format("quad {} corner {} references vertex {}, outside [0, {})",
                                f,
                                corner,
                                q[static_cast<std::size_t>(corner)],
                                count));
            }
        }
        // A repeated corner collapses an edge, which OpenSubdiv accepts but
        // which always means the layout upstream went wrong.
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                if (q[static_cast<std::size_t>(i)] == q[static_cast<std::size_t>(j)]) {
                    throw std::invalid_argument(
                        fmt::format("quad {} uses vertex {} twice, so it is degenerate",
                                    f,
                                    q[static_cast<std::size_t>(i)]));
                }
            }
        }
    }
}

ControlMesh ControlMesh::grid(int rows, int columns) {
    if (rows < 2 || columns < 2) {
        throw std::invalid_argument(
            fmt::format("a grid needs at least 2 x 2 vertices, got {} x {}", rows, columns));
    }

    std::vector<Eigen::Vector3d> vertices;
    vertices.reserve(static_cast<std::size_t>(rows * columns));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            vertices.emplace_back(static_cast<double>(i) / static_cast<double>(rows - 1),
                                  static_cast<double>(j) / static_cast<double>(columns - 1),
                                  0.0);
        }
    }

    std::vector<std::array<int, 4>> quads;
    quads.reserve(static_cast<std::size_t>((rows - 1) * (columns - 1)));
    for (int i = 0; i + 1 < rows; ++i) {
        for (int j = 0; j + 1 < columns; ++j) {
            const int base = i * columns + j;
            // Counter-clockwise seen from +z.
            quads.push_back({base, base + columns, base + columns + 1, base + 1});
        }
    }

    return ControlMesh{std::move(vertices), std::move(quads)};
}

ControlMesh ControlMesh::cube(double half_size) {
    const double h = half_size;
    std::vector<Eigen::Vector3d> vertices{{-h, -h, -h},
                                          {h, -h, -h},
                                          {h, h, -h},
                                          {-h, h, -h},
                                          {-h, -h, h},
                                          {h, -h, h},
                                          {h, h, h},
                                          {-h, h, h}};

    // Every face wound so its normal points outward.
    std::vector<std::array<int, 4>> quads{
        {0, 3, 2, 1}, // -z
        {4, 5, 6, 7}, // +z
        {0, 1, 5, 4}, // -y
        {2, 3, 7, 6}, // +y
        {1, 2, 6, 5}, // +x
        {0, 4, 7, 3}, // -x
    };

    return ControlMesh{std::move(vertices), std::move(quads)};
}

ControlMesh ControlMesh::vertex_fan(int valence) {
    if (valence < 3) {
        throw std::invalid_argument(
            fmt::format("a vertex fan needs a valence of at least 3, got {}", valence));
    }

    const auto n = static_cast<std::size_t>(valence);
    const double step = 2.0 * std::numbers::pi / static_cast<double>(valence);

    // The diagonal ring sits where the two adjacent edge directions meet, so
    // the quads come out as close to regular as a fan allows. For valence 4
    // this reproduces a unit grid patch exactly.
    const double diagonal_radius = 1.0 / std::cos(std::numbers::pi / static_cast<double>(valence));

    std::vector<Eigen::Vector3d> vertices;
    vertices.reserve(2 * n + 1);
    vertices.emplace_back(0.0, 0.0, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double angle = step * static_cast<double>(i);
        vertices.emplace_back(std::cos(angle), std::sin(angle), 0.0);
    }
    for (std::size_t i = 0; i < n; ++i) {
        const double angle = step * (static_cast<double>(i) + 0.5);
        vertices.emplace_back(
            diagonal_radius * std::cos(angle), diagonal_radius * std::sin(angle), 0.0);
    }

    std::vector<std::array<int, 4>> quads;
    quads.reserve(n);
    for (int i = 0; i < valence; ++i) {
        quads.push_back({0, 1 + i, valence + 1 + i, 1 + (i + 1) % valence});
    }

    return ControlMesh{std::move(vertices), std::move(quads)};
}

void ControlMesh::set_vertices(std::vector<Eigen::Vector3d> vertices) {
    if (vertices.size() != vertices_.size()) {
        throw std::invalid_argument(
            fmt::format("the mesh has {} vertices, got {}", vertices_.size(), vertices.size()));
    }
    vertices_ = std::move(vertices);
}

int ControlMesh::face_count(int vertex) const {
    int count = 0;
    for (const std::array<int, 4>& q : quads_) {
        if (std::find(q.begin(), q.end(), vertex) != q.end()) {
            ++count;
        }
    }
    return count;
}

bool ControlMesh::is_boundary_vertex(int vertex) const {
    // An edge used by exactly one face is a boundary edge; a vertex on one is a
    // boundary vertex.
    std::map<std::pair<int, int>, int> edge_use;
    for (const std::array<int, 4>& q : quads_) {
        for (int i = 0; i < 4; ++i) {
            ++edge_use[ordered_edge(q[static_cast<std::size_t>(i)],
                                    q[static_cast<std::size_t>((i + 1) % 4)])];
        }
    }

    for (const auto& [edge, uses] : edge_use) {
        if (uses == 1 && (edge.first == vertex || edge.second == vertex)) {
            return true;
        }
    }
    return false;
}

void ControlMesh::set_crease(int v0, int v1, double sharpness) {
    const std::pair<int, int> wanted = ordered_edge(v0, v1);

    bool found = false;
    for (const std::array<int, 4>& q : quads_) {
        for (int i = 0; i < 4 && !found; ++i) {
            if (ordered_edge(q[static_cast<std::size_t>(i)],
                             q[static_cast<std::size_t>((i + 1) % 4)]) == wanted) {
                found = true;
            }
        }
    }
    if (!found) {
        throw std::invalid_argument(
            fmt::format("vertices {} and {} do not span an edge of this mesh", v0, v1));
    }

    creases_.push_back(Crease{v0, v1, sharpness});
}

void ControlMesh::set_corner(int vertex, double sharpness) {
    if (vertex < 0 || vertex >= static_cast<int>(vertices_.size())) {
        throw std::invalid_argument(
            fmt::format("vertex {} is outside [0, {})", vertex, vertices_.size()));
    }
    corners_.push_back(Corner{vertex, sharpness});
}

void ControlMesh::write_obj(const std::filesystem::path& path) const {
    write_obj_file(path, vertices_, quads_, "control mesh");
}

void PolyMesh::write_obj(const std::filesystem::path& path) const {
    write_obj_file(path, vertices, quads, "refined mesh");
}

} // namespace n2s
