#pragma once

#include "n2s/subd/control_mesh.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <memory>
#include <vector>

namespace n2s {

/// How the limit surface behaves at a mesh boundary.
///
/// `EdgeAndCorner` is the default and the one the whole watertightness argument
/// rests on: with it, the limit curve along a boundary depends only on the
/// boundary control points, so two patches that share boundary control points
/// share a boundary curve exactly. M3 has a test pinning that property down.
enum class BoundaryInterpolation {
    None,          ///< Boundaries float; the limit surface shrinks away from them.
    EdgeOnly,      ///< Boundary edges are creased, corners are not.
    EdgeAndCorner, ///< Boundary edges creased and corners pinned.
};

struct SubdivisionOptions {
    BoundaryInterpolation boundary = BoundaryInterpolation::EdgeAndCorner;

    /// Depth limit for the adaptive refinement backing the limit evaluation.
    /// Higher costs memory and buys accuracy only near extraordinary vertices.
    int max_isolation_level = 6;
};

/// A point at which to evaluate the limit surface: a face of the control mesh
/// and a parametric coordinate within it.
struct LimitLocation {
    int face; ///< Index into `ControlMesh::quads()`.
    double u;
    double v;
};

/// How many orders of derivative stencils to generate. Second derivatives cost
/// extra to build and are only needed for curvature, so they are opt-in.
enum class DerivativeOrder {
    None,
    First,
    Second,
};

/// Limit surface quantities at one location. Derivatives are with respect to
/// the face-local `(u, v)`.
struct LimitSample {
    Eigen::Vector3d position;
    Eigen::Vector3d du;
    Eigen::Vector3d dv;
    /// Unit normal, `du x dv` normalised. Zero where the two are parallel,
    /// which on a limit surface happens only at a degenerate configuration.
    Eigen::Vector3d normal;

    /// Second derivatives. Zero unless the sample was requested with
    /// `DerivativeOrder::Second`.
    Eigen::Vector3d duu = Eigen::Vector3d::Zero();
    Eigen::Vector3d duv = Eigen::Vector3d::Zero();
    Eigen::Vector3d dvv = Eigen::Vector3d::Zero();
};

/// Sparse matrices mapping control points to limit quantities:
/// `limit = A * P`, where `P` is the `num_control_vertices x 3` matrix of
/// control point positions and each row of `A` corresponds to one requested
/// location.
///
/// This is the backbone of every fitting method in Part B. Interpolation solves
/// `A V = target` for the interior control points; least squares minimises
/// `|A V - P|^2`. Both need `A` and nothing else, which is why it is built once
/// here rather than rediscovered per method.
struct LimitMatrices {
    Eigen::SparseMatrix<double, Eigen::RowMajor> position;
    Eigen::SparseMatrix<double, Eigen::RowMajor> du;
    Eigen::SparseMatrix<double, Eigen::RowMajor> dv;
    Eigen::SparseMatrix<double, Eigen::RowMajor> duu;
    Eigen::SparseMatrix<double, Eigen::RowMajor> duv;
    Eigen::SparseMatrix<double, Eigen::RowMajor> dvv;

    /// The order actually generated. Matrices above that order are left empty
    /// rather than filled with zeros, so using one by mistake fails loudly on a
    /// dimension mismatch instead of quietly producing zero curvature.
    DerivativeOrder order = DerivativeOrder::None;
};

/// A Catmull-Clark subdivision surface over a quad control mesh.
///
/// Every stencil and limit weight comes from OpenSubdiv. CLAUDE.md forbids
/// hand-derived subdivision weights outright, because getting the
/// extraordinary-vertex limit weights wrong by hand is exactly what invalidated
/// the thesis this project supersedes. The hand-written formulas in the tests
/// are oracles checking OpenSubdiv, never a substitute for it.
class SubdivisionSurface {
public:
    explicit SubdivisionSurface(ControlMesh mesh, SubdivisionOptions options = {});
    ~SubdivisionSurface();

    SubdivisionSurface(SubdivisionSurface&&) noexcept;
    SubdivisionSurface& operator=(SubdivisionSurface&&) noexcept;
    SubdivisionSurface(const SubdivisionSurface&) = delete;
    SubdivisionSurface& operator=(const SubdivisionSurface&) = delete;

    const ControlMesh& control_mesh() const noexcept;
    const SubdivisionOptions& options() const noexcept;

    /// Limit position, derivatives and normal at each location.
    std::vector<LimitSample> evaluate_limit(const std::vector<LimitLocation>& locations,
                                            DerivativeOrder order = DerivativeOrder::First) const;

    /// Single-location convenience. Builds a table for one point, so prefer the
    /// batch form in any loop.
    LimitSample evaluate_limit(const LimitLocation& location,
                               DerivativeOrder order = DerivativeOrder::First) const;

    /// Builds `A` and the requested derivative matrices for the given
    /// locations. Columns correspond to control mesh vertices in order.
    LimitMatrices build_limit_matrices(const std::vector<LimitLocation>& locations,
                                       DerivativeOrder order = DerivativeOrder::First) const;

    /// Uniform Catmull-Clark refinement to the given level, for display and for
    /// OBJ export. Level 0 returns the control mesh itself.
    PolyMesh refine_uniform(int level) const;

    /// The local subdivision matrix about an interior vertex: the linear map
    /// taking that vertex and its 1-ring to the corresponding vertex and 1-ring
    /// after one Catmull-Clark step.
    ///
    /// Rows and columns share one ordering: the centre, then the
    /// edge-adjacent 1-ring in rotational order, then the face diagonals in the
    /// matching order. Row `i` is the image of the entry that column `i` names,
    /// so the matrix is square and its eigenvalues are meaningful.
    ///
    /// Its eigenstructure governs how the limit surface behaves near an
    /// extraordinary vertex -- the subdominant eigenvalue sets the rate at which
    /// curvature converges, or fails to -- which is what R5 is about. The M3
    /// oracle checks the valence-4 spectrum against the textbook values.
    ///
    /// Throws `std::invalid_argument` if `vertex` is on a boundary, where the
    /// 1-ring is not closed and the matrix is not defined.
    Eigen::MatrixXd local_subdivision_matrix(int vertex) const;

    /// Limit positions of the control mesh vertices themselves. The natural
    /// "where does this vertex end up" query, and the one the regular and
    /// valence-n limit weight oracles check.
    std::vector<Eigen::Vector3d> limit_positions_of_vertices() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

/// Control point positions as a `num_vertices x 3` matrix, so that
/// `matrices.position * control_point_matrix(mesh)` gives the limit points.
Eigen::MatrixXd control_point_matrix(const ControlMesh& mesh);

/// Inverse of the above.
std::vector<Eigen::Vector3d> to_points(const Eigen::MatrixXd& rows);

} // namespace n2s
