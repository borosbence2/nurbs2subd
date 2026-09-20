#include "n2s/fit/refine.hpp"

#include <fmt/format.h>

#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace n2s::fit {

namespace {

/// Key for an undirected edge, so that the two faces sharing it agree on which
/// edge it is regardless of the direction each traverses it.
std::pair<int, int> edge_key(int a, int b) {
    return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
}

} // namespace

DomainLayout refine_quads(const DomainLayout& layout, int factor) {
    if (factor < 1) {
        throw std::invalid_argument(
            fmt::format("refinement factor must be at least 1, got {}", factor));
    }
    if (layout.quads.empty() || layout.vertices.empty()) {
        throw std::invalid_argument("cannot refine an empty layout");
    }
    const auto vertex_count = static_cast<int>(layout.vertices.size());
    for (std::size_t f = 0; f < layout.quads.size(); ++f) {
        for (const int index : layout.quads[f]) {
            if (index < 0 || index >= vertex_count) {
                throw std::invalid_argument(
                    fmt::format("layout quad {} references vertex {}, outside [0, {})",
                                f,
                                index,
                                vertex_count));
            }
        }
    }

    if (factor == 1) {
        return layout;
    }

    DomainLayout refined;
    refined.vertices = layout.vertices; // original corners keep their indices

    // New vertices interior to an edge, stored once per undirected edge in the
    // direction of its lower-numbered endpoint.
    std::map<std::pair<int, int>, std::vector<int>> edge_vertices;

    const auto interior_of_edge = [&](int a, int b) -> const std::vector<int>& {
        const std::pair<int, int> key = edge_key(a, b);
        const auto found = edge_vertices.find(key);
        if (found != edge_vertices.end()) {
            return found->second;
        }

        std::vector<int> created;
        created.reserve(static_cast<std::size_t>(factor - 1));
        const Eigen::Vector2d& from = layout.vertices[static_cast<std::size_t>(key.first)];
        const Eigen::Vector2d& to = layout.vertices[static_cast<std::size_t>(key.second)];
        for (int k = 1; k < factor; ++k) {
            const double t = static_cast<double>(k) / static_cast<double>(factor);
            refined.vertices.push_back(from + t * (to - from));
            created.push_back(static_cast<int>(refined.vertices.size()) - 1);
        }
        return edge_vertices.emplace(key, std::move(created)).first->second;
    };

    for (const std::array<int, 4>& quad : layout.quads) {
        const int v0 = quad[0];
        const int v1 = quad[1];
        const int v2 = quad[2];
        const int v3 = quad[3];

        // grid[i][j] is the sample at (i/factor, j/factor) in the face's local
        // (u, v), matching the convention that u runs v0 -> v1 and v runs
        // v0 -> v3.
        std::vector<std::vector<int>> grid(
            static_cast<std::size_t>(factor) + 1,
            std::vector<int>(static_cast<std::size_t>(factor) + 1, -1));

        const auto at = [&](int i, int j) -> int& {
            return grid[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        };

        at(0, 0) = v0;
        at(factor, 0) = v1;
        at(factor, factor) = v2;
        at(0, factor) = v3;

        // The four edges. `interior_of_edge` returns them ordered from the
        // lower-numbered endpoint, so each is walked in the direction this face
        // needs it.
        const auto place_edge = [&](int a, int b, auto&& assign) {
            const std::vector<int>& shared = interior_of_edge(a, b);
            const bool forward = a < b;
            for (int k = 1; k < factor; ++k) {
                const std::size_t index =
                    static_cast<std::size_t>(forward ? k - 1 : factor - 1 - k);
                assign(k, shared[index]);
            }
        };

        place_edge(v0, v1, [&](int k, int index) { at(k, 0) = index; });
        place_edge(v1, v2, [&](int k, int index) { at(factor, k) = index; });
        place_edge(v3, v2, [&](int k, int index) { at(k, factor) = index; });
        place_edge(v0, v3, [&](int k, int index) { at(0, k) = index; });

        // The interior, bilinearly.
        const Eigen::Vector2d& a = layout.vertices[static_cast<std::size_t>(v0)];
        const Eigen::Vector2d& b = layout.vertices[static_cast<std::size_t>(v1)];
        const Eigen::Vector2d& c = layout.vertices[static_cast<std::size_t>(v2)];
        const Eigen::Vector2d& d = layout.vertices[static_cast<std::size_t>(v3)];

        for (int i = 1; i < factor; ++i) {
            for (int j = 1; j < factor; ++j) {
                const double u = static_cast<double>(i) / static_cast<double>(factor);
                const double v = static_cast<double>(j) / static_cast<double>(factor);
                refined.vertices.push_back((1.0 - u) * (1.0 - v) * a + u * (1.0 - v) * b +
                                           u * v * c + (1.0 - u) * v * d);
                at(i, j) = static_cast<int>(refined.vertices.size()) - 1;
            }
        }

        for (int i = 0; i < factor; ++i) {
            for (int j = 0; j < factor; ++j) {
                refined.quads.push_back({at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)});
            }
        }
    }

    return refined;
}

} // namespace n2s::fit
