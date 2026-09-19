#include "n2s/nurbs/curve.hpp"
#include "n2s/nurbs/surface.hpp"
#include "n2s/subd/subdivision.hpp"

#include <Eigen/Eigenvalues>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Catch::Approx;
using n2s::ControlMesh;
using n2s::KnotVector;
using n2s::LimitLocation;
using n2s::LimitMatrices;
using n2s::SubdivisionSurface;

namespace {

/// Uniform (unclamped) cubic knot vector for `count` control points, knots at
/// the integers. The interior of the resulting B-spline is exactly what
/// Catmull-Clark produces over a regular region, which is what makes it usable
/// as an oracle.
KnotVector uniform_cubic(std::size_t count) {
    std::vector<double> knots;
    knots.reserve(count + 4);
    for (std::size_t i = 0; i < count + 4; ++i) {
        knots.push_back(static_cast<double>(i));
    }
    return KnotVector{3, std::move(knots)};
}

/// Heights making a grid non-planar and asymmetric, so that a transposed index
/// or a mirrored parameterisation cannot pass unnoticed.
double bumpy(int i, int j) {
    return 0.31 * std::sin(0.9 * i) + 0.22 * std::cos(1.4 * j) + 0.07 * i * j;
}

ControlMesh bumpy_grid(int rows, int columns) {
    ControlMesh mesh = ControlMesh::grid(rows, columns);
    std::vector<Eigen::Vector3d> vertices = mesh.vertices();
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            vertices[static_cast<std::size_t>(i * columns + j)].z() = bumpy(i, j);
        }
    }
    mesh.set_vertices(std::move(vertices));
    return mesh;
}

} // namespace

// ---------------------------------------------------------------------------
// The limit weight oracles. These hand-written formulas are the one place
// CLAUDE.md permits them: they check OpenSubdiv, they never stand in for it.
// ---------------------------------------------------------------------------

TEST_CASE("regular vertex limit weights are 4/9, 1/9 and 1/36", "[subd][oracle]") {
    // The valence-4 case, stated separately because it is the one every text
    // quotes and the one a reader will check first.
    const SubdivisionSurface surface{ControlMesh::vertex_fan(4)};
    const LimitMatrices matrices =
        surface.build_limit_matrices({LimitLocation{0, 0.0, 0.0}}, n2s::DerivativeOrder::None);

    const Eigen::MatrixXd dense = Eigen::MatrixXd(matrices.position);
    REQUIRE(dense.rows() == 1);
    REQUIRE(dense.cols() == 9);

    CHECK(dense(0, 0) == Approx(4.0 / 9.0).epsilon(1e-12));
    for (int i = 1; i <= 4; ++i) {
        INFO("edge neighbour " << i);
        CHECK(dense(0, i) == Approx(1.0 / 9.0).epsilon(1e-12));
    }
    for (int i = 5; i <= 8; ++i) {
        INFO("face diagonal " << i);
        CHECK(dense(0, i) == Approx(1.0 / 36.0).epsilon(1e-12));
    }
}

TEST_CASE("valence-n interior limit weights match the analytic mask", "[subd][oracle]") {
    // Centre n/(n+5), each edge-adjacent 1-ring vertex 4/(n(n+5)), each face
    // diagonal 1/(n(n+5)). Getting this wrong by hand for n != 4 is precisely
    // the defect that invalidated the thesis, so it is checked for every
    // valence the plan names.
    for (const int valence : {3, 4, 5, 6}) {
        const SubdivisionSurface surface{ControlMesh::vertex_fan(valence)};
        const LimitMatrices matrices =
            surface.build_limit_matrices({LimitLocation{0, 0.0, 0.0}}, n2s::DerivativeOrder::None);
        const Eigen::MatrixXd dense = Eigen::MatrixXd(matrices.position);

        const auto n = static_cast<double>(valence);
        const double centre = n / (n + 5.0);
        const double edge = 4.0 / (n * (n + 5.0));
        const double diagonal = 1.0 / (n * (n + 5.0));

        INFO("valence " << valence);
        REQUIRE(dense.cols() == 2 * valence + 1);

        CHECK(dense(0, 0) == Approx(centre).epsilon(1e-12));
        for (int i = 1; i <= valence; ++i) {
            INFO("edge neighbour " << i);
            CHECK(dense(0, i) == Approx(edge).epsilon(1e-12));
        }
        for (int i = valence + 1; i <= 2 * valence; ++i) {
            INFO("face diagonal " << i);
            CHECK(dense(0, i) == Approx(diagonal).epsilon(1e-12));
        }

        // The mask is affine, as every subdivision mask must be.
        CHECK(dense.row(0).sum() == Approx(1.0).epsilon(1e-12));
    }
}

TEST_CASE("the regular local subdivision matrix has the textbook spectrum", "[subd][oracle]") {
    // {1, 1/2, 1/2, 1/4, 1/4, 1/4, 1/8, 1/8, 1/16}. The leading 1 is the affine
    // invariant; the subdominant 1/2 governs how fast the surface flattens
    // toward the limit point, which is the quantity R5 cares about at
    // extraordinary vertices.
    const SubdivisionSurface surface{ControlMesh::vertex_fan(4)};
    const Eigen::MatrixXd matrix = surface.local_subdivision_matrix(0);

    REQUIRE(matrix.rows() == 9);
    REQUIRE(matrix.cols() == 9);

    const Eigen::EigenSolver<Eigen::MatrixXd> solver(matrix);
    std::vector<double> eigenvalues;
    for (Eigen::Index i = 0; i < solver.eigenvalues().size(); ++i) {
        const std::complex<double> value = solver.eigenvalues()(i);
        INFO("eigenvalue " << i << " = " << value.real() << " + " << value.imag() << "i");
        CHECK(std::abs(value.imag()) < 1e-10); // the regular case is diagonalisable over the reals
        eigenvalues.push_back(value.real());
    }

    std::sort(eigenvalues.begin(), eigenvalues.end(), std::greater<>());

    const std::vector<double> expected{1.0, 0.5, 0.5, 0.25, 0.25, 0.25, 0.125, 0.125, 0.0625};
    REQUIRE(eigenvalues.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        INFO("sorted eigenvalue " << i);
        CHECK(eigenvalues[i] == Approx(expected[i]).margin(1e-10));
    }
}

// ---------------------------------------------------------------------------
// Reproduction of the regular B-spline.
// ---------------------------------------------------------------------------

TEST_CASE("a regular grid reproduces the uniform bicubic B-spline exactly", "[subd][oracle]") {
    // Away from extraordinary vertices the Catmull-Clark limit surface *is* the
    // uniform bicubic B-spline over the same control net. Exactly, not
    // approximately, so the tolerance is machine precision.
    const int size = 6;
    const ControlMesh mesh = bumpy_grid(size, size);
    const SubdivisionSurface surface{mesh};

    const KnotVector knots = uniform_cubic(static_cast<std::size_t>(size));
    const n2s::NurbsSurface bspline = n2s::NurbsSurface::bspline(knots, knots, mesh.vertices());

    double worst = 0.0;
    // Quad (i, j) is regular when its whole 4x4 neighbourhood exists, i.e. for
    // i, j in [1, size-3]. Its 4x4 net is control points (i-1..i+2), which is
    // the B-spline patch over the knot interval [i+2, i+3].
    for (int i = 1; i <= size - 3; ++i) {
        for (int j = 1; j <= size - 3; ++j) {
            const int face = i * (size - 1) + j;

            for (int a = 0; a <= 4; ++a) {
                for (int b = 0; b <= 4; ++b) {
                    const double u = static_cast<double>(a) / 4.0;
                    const double v = static_cast<double>(b) / 4.0;

                    const Eigen::Vector3d limit =
                        surface.evaluate_limit(LimitLocation{face, u, v}).position;
                    const Eigen::Vector3d reference = bspline.evaluate(
                        static_cast<double>(i) + 2.0 + u, static_cast<double>(j) + 2.0 + v);

                    INFO("face (" << i << ", " << j << ") at (" << u << ", " << v << ")");
                    REQUIRE((limit - reference).norm() < 1e-12);
                    worst = std::max(worst, (limit - reference).norm());
                }
            }
        }
    }
    INFO("worst deviation " << worst);
    CHECK(worst < 1e-12);
}

// ---------------------------------------------------------------------------
// The boundary property that watertightness rests on.
// ---------------------------------------------------------------------------

TEST_CASE("the boundary limit curve depends only on the boundary control points",
          "[subd][oracle][watertight]") {
    // This is the whole basis of R3. If two patches share their boundary
    // control points, and the limit boundary curve is a function of those
    // points alone, then the two limit surfaces meet with no gap at all --
    // not a small gap, none. The test asserts the independence directly, on
    // the stencil weights, rather than inferring it from a measured gap.
    const int size = 6;
    const ControlMesh mesh = bumpy_grid(size, size);
    const SubdivisionSurface surface{mesh};

    // The j = 0 boundary: vertices with index divisible by `size`.
    std::vector<LimitLocation> locations;
    for (int i = 0; i + 1 < size; ++i) {
        const int face = i * (size - 1) + 0;
        for (int a = 0; a <= 4; ++a) {
            locations.push_back(LimitLocation{face, static_cast<double>(a) / 4.0, 0.0});
        }
    }

    const LimitMatrices matrices =
        surface.build_limit_matrices(locations, n2s::DerivativeOrder::None);
    const Eigen::MatrixXd dense = Eigen::MatrixXd(matrices.position);

    for (Eigen::Index row = 0; row < dense.rows(); ++row) {
        for (Eigen::Index column = 0; column < dense.cols(); ++column) {
            if (column % size == 0) {
                continue; // a boundary control point, allowed to contribute
            }
            INFO("row " << row << " draws on interior control point " << column);
            REQUIRE(std::abs(dense(row, column)) < 1e-14);
        }
    }
}

TEST_CASE("the boundary limit curve is the cubic B-spline of the boundary points",
          "[subd][oracle][watertight]") {
    const int size = 6;
    const ControlMesh mesh = bumpy_grid(size, size);
    const SubdivisionSurface surface{mesh};

    std::vector<Eigen::Vector3d> boundary;
    for (int i = 0; i < size; ++i) {
        boundary.push_back(mesh.vertices()[static_cast<std::size_t>(i * size)]);
    }
    const n2s::NurbsCurve curve =
        n2s::NurbsCurve::bspline(uniform_cubic(boundary.size()), boundary);

    // Only the spans with a full four-point neighbourhood are ordinary cubic
    // B-spline spans. The two end spans are governed by the corner rule
    // instead, which pins the endpoint rather than following the uniform
    // B-spline, and are checked by the interpolation assertion below.
    for (int i = 1; i <= size - 3; ++i) {
        const int face = i * (size - 1) + 0;
        for (int a = 0; a <= 8; ++a) {
            const double u = static_cast<double>(a) / 8.0;
            const Eigen::Vector3d limit =
                surface.evaluate_limit(LimitLocation{face, u, 0.0}).position;
            const Eigen::Vector3d reference = curve.evaluate(static_cast<double>(i) + 2.0 + u);

            INFO("face " << face << " at u = " << u);
            REQUIRE((limit - reference).norm() < 1e-12);
        }
    }
}

TEST_CASE("EDGE_AND_CORNER interpolates the corner control points", "[subd][oracle]") {
    const int size = 5;
    const ControlMesh mesh = bumpy_grid(size, size);
    const SubdivisionSurface surface{mesh};

    const std::vector<Eigen::Vector3d> limits = surface.limit_positions_of_vertices();
    REQUIRE(limits.size() == mesh.num_vertices());

    const std::vector<std::size_t> corners{0,
                                           static_cast<std::size_t>(size - 1),
                                           static_cast<std::size_t>((size - 1) * size),
                                           static_cast<std::size_t>(size * size - 1)};

    for (const std::size_t corner : corners) {
        INFO("corner vertex " << corner);
        CHECK((limits[corner] - mesh.vertices()[corner]).norm() < 1e-12);
    }
}

TEST_CASE("BoundaryInterpolation::None leaves the boundary with no limit surface",
          "[subd][subdivision]") {
    // Not a quirk worth working around, but the definition of the mode: with no
    // boundary rule the limit surface pulls inward and simply does not exist
    // over the boundary region, so OpenSubdiv produces no stencil there. That
    // is the concrete reason this project defaults to EdgeAndCorner -- R3 needs
    // a limit curve *on* the boundary to make two patches agree along it.
    const ControlMesh mesh = bumpy_grid(5, 5);

    const SubdivisionSurface pinned{
        mesh, n2s::SubdivisionOptions{n2s::BoundaryInterpolation::EdgeAndCorner, 6}};
    CHECK((pinned.limit_positions_of_vertices()[0] - mesh.vertices()[0]).norm() < 1e-12);

    const SubdivisionSurface floating{mesh,
                                      n2s::SubdivisionOptions{n2s::BoundaryInterpolation::None, 6}};
    CHECK_THROWS_WITH(floating.limit_positions_of_vertices(),
                      Catch::Matchers::ContainsSubstring("no limit surface"));

    // The interior still evaluates perfectly well under either mode.
    const int centre_face = 1 * 4 + 1;
    CHECK_NOTHROW(floating.evaluate_limit(LimitLocation{centre_face, 0.5, 0.5}));
}

// ---------------------------------------------------------------------------
// Properties of the stencil matrix itself.
// ---------------------------------------------------------------------------

TEST_CASE("rows of the limit matrix sum to one", "[subd][oracle]") {
    // Affine invariance: translating every control point must translate the
    // limit surface by the same amount and nothing else.
    const std::vector<ControlMesh> meshes{
        ControlMesh::grid(5, 5), ControlMesh::cube(), ControlMesh::vertex_fan(5)};

    for (const ControlMesh& mesh : meshes) {
        const SubdivisionSurface surface{mesh};

        std::vector<LimitLocation> locations;
        for (std::size_t f = 0; f < mesh.num_quads(); ++f) {
            for (const double u : {0.0, 0.25, 0.5, 1.0}) {
                for (const double v : {0.0, 0.5, 0.75, 1.0}) {
                    locations.push_back(LimitLocation{static_cast<int>(f), u, v});
                }
            }
        }

        const LimitMatrices matrices =
            surface.build_limit_matrices(locations, n2s::DerivativeOrder::First);
        const Eigen::VectorXd sums =
            Eigen::MatrixXd(matrices.position) * Eigen::VectorXd::Ones(matrices.position.cols());

        for (Eigen::Index row = 0; row < sums.size(); ++row) {
            INFO("row " << row << " of a " << mesh.num_quads() << "-quad mesh");
            REQUIRE(sums(row) == Approx(1.0).epsilon(1e-12));
        }

        // Derivative stencils annihilate constants: a translated surface has
        // the same tangents.
        const Eigen::VectorXd du_sums =
            Eigen::MatrixXd(matrices.du) * Eigen::VectorXd::Ones(matrices.du.cols());
        for (Eigen::Index row = 0; row < du_sums.size(); ++row) {
            INFO("du row " << row);
            REQUIRE(std::abs(du_sums(row)) < 1e-12);
        }
    }
}

TEST_CASE("the limit matrix reproduces direct evaluation", "[subd][subdivision]") {
    const ControlMesh mesh = bumpy_grid(5, 5);
    const SubdivisionSurface surface{mesh};

    const std::vector<LimitLocation> locations{
        {0, 0.25, 0.75}, {3, 0.5, 0.5}, {7, 1.0, 0.0}, {11, 0.125, 0.875}};

    const LimitMatrices matrices =
        surface.build_limit_matrices(locations, n2s::DerivativeOrder::First);
    const Eigen::MatrixXd points = n2s::control_point_matrix(mesh);
    const Eigen::MatrixXd positions = matrices.position * points;

    const std::vector<n2s::LimitSample> samples = surface.evaluate_limit(locations);
    REQUIRE(samples.size() == locations.size());

    for (std::size_t i = 0; i < samples.size(); ++i) {
        INFO("location " << i);
        CHECK(
            (samples[i].position - positions.row(static_cast<Eigen::Index>(i)).transpose()).norm() <
            1e-14);
        CHECK(samples[i].normal.norm() == Approx(1.0).epsilon(1e-12));
    }
}

TEST_CASE("derivative stencils match central differences of the limit surface",
          "[subd][subdivision]") {
    const ControlMesh mesh = bumpy_grid(6, 6);
    const SubdivisionSurface surface{mesh};

    const int face = 2 * 5 + 2; // an interior, fully regular quad
    const double h = 1e-6;

    for (const double u : {0.3, 0.5, 0.7}) {
        for (const double v : {0.25, 0.6}) {
            const n2s::LimitSample sample = surface.evaluate_limit(LimitLocation{face, u, v});

            const Eigen::Vector3d fd_u =
                (surface.evaluate_limit(LimitLocation{face, u + h, v}).position -
                 surface.evaluate_limit(LimitLocation{face, u - h, v}).position) /
                (2.0 * h);
            const Eigen::Vector3d fd_v =
                (surface.evaluate_limit(LimitLocation{face, u, v + h}).position -
                 surface.evaluate_limit(LimitLocation{face, u, v - h}).position) /
                (2.0 * h);

            INFO("u = " << u << ", v = " << v);
            CHECK((sample.du - fd_u).norm() < 1e-6);
            CHECK((sample.dv - fd_v).norm() < 1e-6);
        }
    }
}

TEST_CASE("second derivative stencils match central differences of the first",
          "[subd][subdivision]") {
    // Differentiating the analytic first derivatives rather than the positions:
    // a second central difference of position loses half the available digits,
    // while differencing the first derivative once does not.
    const ControlMesh mesh = bumpy_grid(6, 6);
    const SubdivisionSurface surface{mesh};

    const int face = 2 * 5 + 2; // interior, fully regular
    const double h = 1e-6;

    for (const double u : {0.3, 0.55}) {
        for (const double v : {0.35, 0.7}) {
            const n2s::LimitSample sample =
                surface.evaluate_limit(LimitLocation{face, u, v}, n2s::DerivativeOrder::Second);

            const auto first_at = [&](double a, double b) {
                return surface.evaluate_limit(LimitLocation{face, a, b},
                                              n2s::DerivativeOrder::First);
            };

            const Eigen::Vector3d fd_uu =
                (first_at(u + h, v).du - first_at(u - h, v).du) / (2.0 * h);
            const Eigen::Vector3d fd_vv =
                (first_at(u, v + h).dv - first_at(u, v - h).dv) / (2.0 * h);
            const Eigen::Vector3d fd_uv =
                (first_at(u, v + h).du - first_at(u, v - h).du) / (2.0 * h);

            INFO("u = " << u << ", v = " << v);
            CHECK((sample.duu - fd_uu).norm() < 1e-5);
            CHECK((sample.dvv - fd_vv).norm() < 1e-5);
            CHECK((sample.duv - fd_uv).norm() < 1e-5);

            // The patch is genuinely curved here, so these must not be the
            // zeros a silently-unpopulated weight array would give.
            CHECK(sample.duu.norm() > 1e-6);
        }
    }
}

TEST_CASE("derivative matrices above the requested order are left empty", "[subd][subdivision]") {
    // Empty rather than zero-filled, so that using one by mistake fails on a
    // dimension mismatch instead of quietly reporting zero curvature.
    const SubdivisionSurface surface{bumpy_grid(5, 5)};
    const std::vector<LimitLocation> locations{{0, 0.5, 0.5}, {1, 0.25, 0.75}};

    const LimitMatrices none = surface.build_limit_matrices(locations, n2s::DerivativeOrder::None);
    CHECK(none.order == n2s::DerivativeOrder::None);
    CHECK(none.position.rows() == 2);
    CHECK(none.du.size() == 0);

    const LimitMatrices first =
        surface.build_limit_matrices(locations, n2s::DerivativeOrder::First);
    CHECK(first.du.rows() == 2);
    CHECK(first.duu.size() == 0);

    const LimitMatrices second =
        surface.build_limit_matrices(locations, n2s::DerivativeOrder::Second);
    CHECK(second.duu.rows() == 2);
    CHECK(second.dvv.rows() == 2);

    // Second derivative stencils annihilate linear functions, so their rows sum
    // to zero just as the first derivative rows do.
    const Eigen::VectorXd sums =
        Eigen::MatrixXd(second.duu) * Eigen::VectorXd::Ones(second.duu.cols());
    for (Eigen::Index row = 0; row < sums.size(); ++row) {
        CHECK(std::abs(sums(row)) < 1e-12);
    }
}

// ---------------------------------------------------------------------------
// Refinement and export.
// ---------------------------------------------------------------------------

TEST_CASE("uniform refinement quadruples the face count per level", "[subd][subdivision]") {
    const ControlMesh mesh = ControlMesh::cube();
    const SubdivisionSurface surface{mesh};

    CHECK(surface.refine_uniform(0).quads.size() == mesh.num_quads());
    CHECK(surface.refine_uniform(1).quads.size() == mesh.num_quads() * 4);
    CHECK(surface.refine_uniform(2).quads.size() == mesh.num_quads() * 16);

    // The refined mesh converges onto the limit surface, so it must stay within
    // the convex hull of the control mesh: a cube of half size 0.5.
    for (const Eigen::Vector3d& v : surface.refine_uniform(3).vertices) {
        CHECK(v.cwiseAbs().maxCoeff() <= 0.5 + 1e-12);
    }
}

TEST_CASE("refined vertices approach the limit surface", "[subd][subdivision]") {
    const ControlMesh mesh = bumpy_grid(5, 5);
    const SubdivisionSurface surface{mesh};

    // The centre of face 0 is a fixed point of refinement in parameter space,
    // so the refined mesh vertex nearest it must converge onto the limit point.
    const Eigen::Vector3d limit = surface.evaluate_limit(LimitLocation{0, 0.5, 0.5}).position;

    double previous = 1e9;
    for (const int level : {1, 2, 3, 4}) {
        const n2s::PolyMesh refined = surface.refine_uniform(level);
        double closest = 1e9;
        for (const Eigen::Vector3d& v : refined.vertices) {
            closest = std::min(closest, (v - limit).norm());
        }
        INFO("level " << level << ", closest refined vertex " << closest);
        CHECK(closest < previous);
        previous = closest;
    }
    CHECK(previous < 1e-3);
}

TEST_CASE("OBJ export round-trips through a file", "[subd][io]") {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "n2s_test_control_mesh.obj";
    std::filesystem::remove(path);

    const ControlMesh mesh = bumpy_grid(4, 4);
    mesh.write_obj(path);
    REQUIRE(std::filesystem::exists(path));

    std::ifstream stream(path);
    REQUIRE(stream);

    std::size_t vertex_lines = 0;
    std::size_t face_lines = 0;
    std::vector<Eigen::Vector3d> read_back;
    std::string line;
    while (std::getline(stream, line)) {
        if (line.rfind("v ", 0) == 0) {
            ++vertex_lines;
            std::istringstream values(line.substr(2));
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            values >> x >> y >> z;
            read_back.emplace_back(x, y, z);
        } else if (line.rfind("f ", 0) == 0) {
            ++face_lines;
            std::istringstream indices(line.substr(2));
            int index = 0;
            int corners = 0;
            while (indices >> index) {
                // OBJ is 1-based, so every index must be in [1, numVertices].
                CHECK(index >= 1);
                CHECK(index <= static_cast<int>(mesh.num_vertices()));
                ++corners;
            }
            CHECK(corners == 4);
        }
    }

    // Close before removing: on Windows an open handle makes the file
    // undeletable, and the failure surfaces as a filesystem error from the
    // cleanup rather than from anything the test is actually checking.
    stream.close();

    CHECK(vertex_lines == mesh.num_vertices());
    CHECK(face_lines == mesh.num_quads());

    // Written with enough digits to survive the trip: an OBJ at default
    // precision would quietly lose the geometry it is meant to preserve.
    REQUIRE(read_back.size() == mesh.num_vertices());
    for (std::size_t i = 0; i < read_back.size(); ++i) {
        INFO("vertex " << i);
        CHECK((read_back[i] - mesh.vertices()[i]).norm() < 1e-15);
    }

    std::filesystem::remove(path);
}

TEST_CASE("refined meshes export too", "[subd][io]") {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "n2s_test_refined.obj";
    std::filesystem::remove(path);

    const SubdivisionSurface surface{ControlMesh::cube()};
    const n2s::PolyMesh refined = surface.refine_uniform(2);
    refined.write_obj(path);

    REQUIRE(std::filesystem::exists(path));
    CHECK(std::filesystem::file_size(path) > 0);

    std::filesystem::remove(path);
}

TEST_CASE("limit tessellation is consistent with direct evaluation", "[subd][subdivision]") {
    const ControlMesh mesh = bumpy_grid(4, 4);
    const SubdivisionSurface surface{mesh};

    const int per_edge = 4;
    const n2s::TessellatedLimit tessellation = n2s::tessellate_limit(surface, per_edge);

    const std::size_t per_face = static_cast<std::size_t>((per_edge + 1) * (per_edge + 1));
    CHECK(tessellation.mesh.vertices.size() == mesh.num_quads() * per_face);
    CHECK(tessellation.locations.size() == tessellation.mesh.vertices.size());
    CHECK(tessellation.normals.size() == tessellation.mesh.vertices.size());
    CHECK(tessellation.mesh.quads.size() ==
          mesh.num_quads() * static_cast<std::size_t>(per_edge * per_edge));

    // Every vertex really is the limit point of the location recorded for it,
    // which is what lets a scalar computed per location be drawn per vertex.
    for (std::size_t i = 0; i < tessellation.locations.size(); i += 7) {
        const Eigen::Vector3d direct = surface.evaluate_limit(tessellation.locations[i]).position;
        INFO("vertex " << i);
        REQUIRE((tessellation.mesh.vertices[i] - direct).norm() < 1e-14);
    }

    for (const Eigen::Vector3d& n : tessellation.normals) {
        CHECK(n.norm() == Approx(1.0).epsilon(1e-9));
    }
}

TEST_CASE("neighbouring tessellated faces agree along their shared edge", "[subd][subdivision]") {
    // The mesh is deliberately not welded, so the duplicated vertices along a
    // shared edge have to agree numerically instead. If they ever stop
    // agreeing, the limit evaluation is face-dependent and the whole
    // watertightness argument is in trouble.
    const ControlMesh mesh = bumpy_grid(5, 5);
    const SubdivisionSurface surface{mesh};

    // Faces 0 and 1 of a 4x4 quad grid share the edge u in [0,1] at v = 1 of
    // face 0, which is v = 0 of face 1.
    for (int a = 0; a <= 8; ++a) {
        const double u = static_cast<double>(a) / 8.0;
        const Eigen::Vector3d from_first =
            surface.evaluate_limit(LimitLocation{0, u, 1.0}).position;
        const Eigen::Vector3d from_second =
            surface.evaluate_limit(LimitLocation{1, u, 0.0}).position;

        INFO("u = " << u);
        CHECK((from_first - from_second).norm() < 1e-13);
    }
}

TEST_CASE("invalid subdivision inputs are rejected", "[subd][subdivision]") {
    SECTION("a degenerate quad") {
        const std::vector<Eigen::Vector3d> vertices(4, Eigen::Vector3d::Zero());
        CHECK_THROWS_AS(ControlMesh(vertices, {{0, 1, 1, 2}}), std::invalid_argument);
    }

    SECTION("a face index out of range") {
        CHECK_THROWS_AS(
            ControlMesh(std::vector<Eigen::Vector3d>(4, Eigen::Vector3d::Zero()), {{0, 1, 2, 9}}),
            std::invalid_argument);
    }

    SECTION("a limit location naming a face that does not exist") {
        const SubdivisionSurface surface{ControlMesh::grid(3, 3)};
        CHECK_THROWS_AS(surface.evaluate_limit(LimitLocation{99, 0.5, 0.5}), std::invalid_argument);
    }

    SECTION("the subdivision matrix of a boundary vertex") {
        const SubdivisionSurface surface{ControlMesh::grid(4, 4)};
        CHECK_THROWS_AS(surface.local_subdivision_matrix(0), std::invalid_argument);
    }

    SECTION("a fan with too low a valence") {
        CHECK_THROWS_AS(ControlMesh::vertex_fan(2), std::invalid_argument);
    }
}
