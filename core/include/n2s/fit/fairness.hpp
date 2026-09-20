#pragma once

#include "n2s/subd/control_mesh.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <vector>

namespace n2s::fit {

/// The discrete fairness term of the R2 objective, as a sparse matrix.
///
/// `E_fair(V) = ||L V||^2`, with `L` the uniform umbrella operator: the row for
/// vertex `i` holds `1` on the diagonal and `-1/d` on each of its `d`
/// edge-adjacent neighbours, so `(L V)_i` is the offset of vertex `i` from the
/// centroid of its 1-ring. Minimising the squared umbrella is the standard
/// discrete thin-plate energy on a mesh.
///
/// **Three choices worth stating, because each one could have gone otherwise.**
///
/// *Uniform weights, not cotangent.* Cotangent weights depend on the current
/// vertex positions, which are the unknowns; using them would make `E_fair`
/// non-quadratic and the fit a nonlinear solve rather than the one sparse
/// system the plan asks for. The uniform umbrella is purely topological, so `L`
/// is assembled once and stays constant through the solve.
///
/// *Rows for interior vertices only.* The umbrella at a boundary vertex is
/// one-sided: its 1-ring lies to one side, so penalising the offset from that
/// centroid pulls the boundary inward, which is a shrinkage term disguised as a
/// fairness term. It would also fight the boundary constraint, since R2 holds
/// the boundary control points fixed.
///
/// **The caveat this inherits.** The uniform umbrella annihilates affine data
/// only where the 1-ring is symmetric about the vertex, which on a quad mesh
/// means a regular interior vertex with evenly spaced neighbours. At an
/// extraordinary vertex, or anywhere the layout is stretched, a perfectly flat
/// configuration still carries a non-zero umbrella, so `lambda > 0` biases
/// those regions even when there is nothing to fair. That is a real effect and
/// not a small one near the thesis layout's valence-5 vertices; R2's lambda
/// sweep measures it rather than assuming it away.
struct FairnessOperator {
    /// `rows() x mesh.num_vertices()`. Applied to the `n x 3` matrix of control
    /// point coordinates.
    Eigen::SparseMatrix<double> matrix;

    /// The mesh vertex each row belongs to, in row order. Boundary vertices do
    /// not appear.
    std::vector<int> row_vertices;

    Eigen::Index rows() const { return matrix.rows(); }
};

/// Builds the umbrella operator for a mesh's topology.
///
/// Depends only on the connectivity, so two meshes with the same quads produce
/// the same operator whatever their vertex positions.
FairnessOperator umbrella_operator(const ControlMesh& mesh);

/// `||L V||^2` for the mesh's current positions: zero when every interior
/// vertex sits at the centroid of its 1-ring.
double fairness_energy(const FairnessOperator& fairness, const ControlMesh& mesh);

/// The `n x 3` matrix of a mesh's vertex positions, which is what both the
/// fairness operator and the limit stencil matrix are applied to.
Eigen::MatrixXd vertex_matrix(const ControlMesh& mesh);

} // namespace n2s::fit
