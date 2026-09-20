#include "n2s/fit/fairness.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

using Catch::Approx;
using n2s::ControlMesh;
using n2s::fit::FairnessOperator;

namespace {

/// Moves every vertex by `f(position)`, keeping the topology.
template<typename F>
ControlMesh mapped(const ControlMesh& mesh, F&& f) {
    ControlMesh copy = mesh;
    std::vector<Eigen::Vector3d> moved;
    moved.reserve(mesh.num_vertices());
    for (const Eigen::Vector3d& v : mesh.vertices()) {
        moved.push_back(f(v));
    }
    copy.set_vertices(std::move(moved));
    return copy;
}

} // namespace

TEST_CASE("every umbrella row sums to zero", "[fit][fairness]") {
    // The defining property: the row is a difference between a vertex and an
    // average of its neighbours, so it annihilates constants. An operator whose
    // rows did not sum to zero would penalise translating the whole mesh, and
    // the fit would pull the surface toward the origin.
    const ControlMesh mesh = ControlMesh::grid(5, 6);
    const FairnessOperator fairness = n2s::fit::umbrella_operator(mesh);

    REQUIRE(fairness.rows() > 0);
    for (Eigen::Index r = 0; r < fairness.rows(); ++r) {
        double sum = 0.0;
        for (Eigen::Index c = 0; c < fairness.matrix.cols(); ++c) {
            sum += fairness.matrix.coeff(r, c);
        }
        INFO("row " << r);
        CHECK(sum == Approx(0.0).margin(1e-14));
    }
}

TEST_CASE("the umbrella has a row for every interior vertex and no others", "[fit][fairness]") {
    const ControlMesh mesh = ControlMesh::grid(5, 5);
    const FairnessOperator fairness = n2s::fit::umbrella_operator(mesh);

    std::size_t interior = 0;
    for (std::size_t v = 0; v < mesh.num_vertices(); ++v) {
        if (!mesh.is_boundary_vertex(static_cast<int>(v))) {
            ++interior;
        }
    }

    // A 5x5 grid has a 3x3 interior.
    REQUIRE(interior == 9);
    CHECK(static_cast<std::size_t>(fairness.rows()) == interior);
    CHECK(fairness.row_vertices.size() == interior);

    for (const int v : fairness.row_vertices) {
        INFO("row vertex " << v);
        CHECK_FALSE(mesh.is_boundary_vertex(v));
    }
}

TEST_CASE("a regular grid carries no fairness energy", "[fit][fairness][oracle]") {
    // On a regular grid every interior vertex sits exactly at the centroid of
    // its four neighbours, so the energy is identically zero. This is the
    // analytic oracle: any non-zero answer here is the operator, not the mesh.
    const ControlMesh mesh = ControlMesh::grid(6, 7);
    const FairnessOperator fairness = n2s::fit::umbrella_operator(mesh);

    CHECK(n2s::fit::fairness_energy(fairness, mesh) == Approx(0.0).margin(1e-24));
}

TEST_CASE("an affine image of a regular grid still carries none", "[fit][fairness][oracle]") {
    // The umbrella commutes with affine maps on a symmetric 1-ring, so shearing,
    // scaling, rotating or tilting a regular grid leaves the energy at zero. A
    // fairness term that failed this would charge the fit for a plane that
    // happened not to be axis-aligned.
    const ControlMesh grid = ControlMesh::grid(6, 6);
    const FairnessOperator fairness = n2s::fit::umbrella_operator(grid);

    const ControlMesh sheared = mapped(grid, [](const Eigen::Vector3d& p) {
        return Eigen::Vector3d{2.0 * p.x() + 0.5 * p.y() + 1.0,
                               -0.25 * p.x() + 3.0 * p.y() - 2.0,
                               0.75 * p.x() - 0.5 * p.y() + 4.0};
    });

    CHECK(n2s::fit::fairness_energy(fairness, sheared) == Approx(0.0).margin(1e-20));
}

TEST_CASE("a curved grid carries energy, and more of it when more curved",
          "[fit][fairness][oracle]") {
    // The umbrella of a quadratic is its second difference, which is constant
    // and non-zero, so the energy has to rise with the coefficient. Checked as
    // a ratio rather than a value: doubling the coefficient doubles every
    // umbrella and so quadruples the energy.
    const ControlMesh grid = ControlMesh::grid(6, 6);
    const FairnessOperator fairness = n2s::fit::umbrella_operator(grid);

    const auto paraboloid = [](double k) {
        return [k](const Eigen::Vector3d& p) {
            return Eigen::Vector3d{p.x(), p.y(), k * (p.x() * p.x() + p.y() * p.y())};
        };
    };

    const double one = n2s::fit::fairness_energy(fairness, mapped(grid, paraboloid(1.0)));
    const double two = n2s::fit::fairness_energy(fairness, mapped(grid, paraboloid(2.0)));

    CHECK(one > 0.0);
    CHECK(two == Approx(4.0 * one).epsilon(1e-12));
}

TEST_CASE("the umbrella weights are the centroid of the 1-ring", "[fit][fairness][oracle]") {
    // Stated directly on the smallest mesh where it can be checked by hand: the
    // centre of a 3x3 grid has four edge-adjacent neighbours, so its row is 1 on
    // itself and -1/4 on each of them, and nothing on the four diagonals.
    const ControlMesh mesh = ControlMesh::grid(3, 3);
    const FairnessOperator fairness = n2s::fit::umbrella_operator(mesh);

    REQUIRE(fairness.rows() == 1);
    const int centre = fairness.row_vertices.front();

    // Grid vertex (i, j) is index i * columns + j, so the centre of a 3x3 is 4
    // and its edge neighbours are 1, 3, 5 and 7. The diagonals are 0, 2, 6, 8.
    REQUIRE(centre == 4);
    CHECK(fairness.matrix.coeff(0, 4) == Approx(1.0));
    for (const int neighbour : {1, 3, 5, 7}) {
        INFO("edge neighbour " << neighbour);
        CHECK(fairness.matrix.coeff(0, neighbour) == Approx(-0.25));
    }
    for (const int diagonal : {0, 2, 6, 8}) {
        INFO("diagonal " << diagonal);
        CHECK(fairness.matrix.coeff(0, diagonal) == Approx(0.0).margin(1e-15));
    }
}

TEST_CASE("an extraordinary vertex gets a row of its own valence", "[fit][fairness]") {
    // The fan puts one interior vertex of valence n at index 0, with its
    // edge-adjacent ring at 1..n and the face diagonals after it. The umbrella
    // must use the ring and ignore the diagonals whatever n is.
    for (const int valence : {3, 5, 6}) {
        const ControlMesh fan = ControlMesh::vertex_fan(valence);
        const FairnessOperator fairness = n2s::fit::umbrella_operator(fan);

        INFO("valence " << valence);
        REQUIRE(fairness.rows() == 1);
        CHECK(fairness.row_vertices.front() == 0);
        CHECK(fairness.matrix.coeff(0, 0) == Approx(1.0));

        for (int ring = 1; ring <= valence; ++ring) {
            INFO("ring vertex " << ring);
            CHECK(fairness.matrix.coeff(0, ring) == Approx(-1.0 / valence));
        }
        for (int diagonal = valence + 1; diagonal <= 2 * valence; ++diagonal) {
            INFO("diagonal vertex " << diagonal);
            CHECK(fairness.matrix.coeff(0, diagonal) == Approx(0.0).margin(1e-15));
        }
    }
}

TEST_CASE("the operator depends on topology alone", "[fit][fairness]") {
    // It is assembled once and reused across a solve, so it must not depend on
    // the positions it will be applied to.
    const ControlMesh grid = ControlMesh::grid(4, 5);
    const FairnessOperator flat = n2s::fit::umbrella_operator(grid);
    const FairnessOperator moved =
        n2s::fit::umbrella_operator(mapped(grid, [](const Eigen::Vector3d& p) {
            return Eigen::Vector3d{p.x(), p.y(), std::sin(5.0 * p.x()) * p.y()};
        }));

    REQUIRE(flat.rows() == moved.rows());
    CHECK(flat.row_vertices == moved.row_vertices);
    CHECK((Eigen::MatrixXd{flat.matrix} - Eigen::MatrixXd{moved.matrix}).norm() ==
          Approx(0.0).margin(1e-15));
}

TEST_CASE("a mesh with no interior vertex has an empty operator", "[fit][fairness]") {
    // A single quad is all boundary. The energy is then identically zero and
    // the fit is pure least squares, which is the right answer rather than an
    // error: there is nothing to fair.
    const ControlMesh single = ControlMesh::grid(2, 2);
    const FairnessOperator fairness = n2s::fit::umbrella_operator(single);

    CHECK(fairness.rows() == 0);
    CHECK(fairness.row_vertices.empty());
    CHECK(n2s::fit::fairness_energy(fairness, single) == Approx(0.0));
}
