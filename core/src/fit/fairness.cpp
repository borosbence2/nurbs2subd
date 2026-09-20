#include "n2s/fit/fairness.hpp"

#include <algorithm>
#include <array>
#include <set>
#include <utility>
#include <vector>

namespace n2s::fit {

namespace {

/// Edge-adjacent neighbours of every vertex.
///
/// A quad contributes its four sides and not its two diagonals, which is what
/// makes this the 1-ring of the umbrella rather than the full stencil
/// neighbourhood. `std::set` rather than a vector with a dedup pass: every
/// interior edge is walked twice, once from each of its faces, and the ordering
/// also makes the assembled matrix independent of face order.
std::vector<std::set<int>> edge_neighbours(const ControlMesh& mesh) {
    std::vector<std::set<int>> neighbours(mesh.num_vertices());

    for (const std::array<int, 4>& quad : mesh.quads()) {
        for (int corner = 0; corner < 4; ++corner) {
            const int a = quad[static_cast<std::size_t>(corner)];
            const int b = quad[static_cast<std::size_t>((corner + 1) % 4)];
            neighbours[static_cast<std::size_t>(a)].insert(b);
            neighbours[static_cast<std::size_t>(b)].insert(a);
        }
    }

    return neighbours;
}

} // namespace

FairnessOperator umbrella_operator(const ControlMesh& mesh) {
    const std::vector<std::set<int>> neighbours = edge_neighbours(mesh);

    FairnessOperator fairness;
    std::vector<Eigen::Triplet<double>> triplets;

    for (std::size_t v = 0; v < mesh.num_vertices(); ++v) {
        const auto vertex = static_cast<int>(v);
        if (mesh.is_boundary_vertex(vertex)) {
            continue; // one-sided umbrella; see the header
        }

        const std::set<int>& ring = neighbours[v];
        if (ring.empty()) {
            continue; // isolated vertex: nothing to be fair against
        }

        const auto row = static_cast<int>(fairness.row_vertices.size());
        fairness.row_vertices.push_back(vertex);

        const double weight = -1.0 / static_cast<double>(ring.size());
        triplets.emplace_back(row, vertex, 1.0);
        for (const int neighbour : ring) {
            triplets.emplace_back(row, neighbour, weight);
        }
    }

    fairness.matrix.resize(static_cast<Eigen::Index>(fairness.row_vertices.size()),
                           static_cast<Eigen::Index>(mesh.num_vertices()));
    fairness.matrix.setFromTriplets(triplets.begin(), triplets.end());
    fairness.matrix.makeCompressed();

    return fairness;
}

Eigen::MatrixXd vertex_matrix(const ControlMesh& mesh) {
    Eigen::MatrixXd points(static_cast<Eigen::Index>(mesh.num_vertices()), 3);
    for (std::size_t v = 0; v < mesh.num_vertices(); ++v) {
        points.row(static_cast<Eigen::Index>(v)) = mesh.vertices()[v].transpose();
    }
    return points;
}

double fairness_energy(const FairnessOperator& fairness, const ControlMesh& mesh) {
    if (fairness.rows() == 0) {
        return 0.0;
    }
    return (fairness.matrix * vertex_matrix(mesh)).squaredNorm();
}

} // namespace n2s::fit
