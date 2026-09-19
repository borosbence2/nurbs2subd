#include "n2s/subd/subdivision.hpp"

#include <Eigen/Geometry>
#include <fmt/format.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/stencilTableFactory.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/topologyRefiner.h>
#include <opensubdiv/far/topologyRefinerFactory.h>

#include <memory>
#include <stdexcept>
#include <utility>

namespace n2s {

namespace {

namespace Far = OpenSubdiv::Far;
namespace Sdc = OpenSubdiv::Sdc;

/// Double-precision limit stencils. OpenSubdiv's convenience typedefs are all
/// float; the `Real` templates are explicitly instantiated for double as well,
/// and CLAUDE.md requires double throughout the geometry code, so the Real
/// forms are used directly here.
using LimitStencilTableD = Far::LimitStencilTableReal<double>;
using LimitStencilFactoryD = Far::LimitStencilTableFactoryReal<double>;

Sdc::Options::VtxBoundaryInterpolation to_sdc(BoundaryInterpolation boundary) {
    switch (boundary) {
    case BoundaryInterpolation::None:
        return Sdc::Options::VTX_BOUNDARY_NONE;
    case BoundaryInterpolation::EdgeOnly:
        return Sdc::Options::VTX_BOUNDARY_EDGE_ONLY;
    case BoundaryInterpolation::EdgeAndCorner:
        return Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER;
    }
    return Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER;
}

/// Vertex type for PrimvarRefiner, which interpolates through Clear and
/// AddWithWeight rather than through a matrix.
struct RefinablePoint {
    void Clear() { position.setZero(); }

    void AddWithWeight(const RefinablePoint& other, double weight) {
        position += weight * other.position;
    }

    Eigen::Vector3d position = Eigen::Vector3d::Zero();
};

/// The same interface carrying a whole row of weights instead of a position, so
/// that refining the identity yields the subdivision matrix itself.
struct WeightRow {
    void Clear() { weights.setZero(); }

    void AddWithWeight(const WeightRow& other, double weight) { weights += weight * other.weights; }

    Eigen::VectorXd weights;
};

} // namespace

/// Holds the topology arrays alive, because Far::TopologyDescriptor stores bare
/// pointers into them, and caches the adaptively refined refiner used for all
/// limit work.
struct SubdivisionSurface::Impl {
    ControlMesh mesh;
    SubdivisionOptions options;

    std::vector<int> verts_per_face;
    std::vector<Far::Index> face_verts;
    std::vector<Far::Index> crease_pairs;
    std::vector<float> crease_weights;
    std::vector<Far::Index> corner_indices;
    std::vector<float> corner_weights;

    std::unique_ptr<Far::TopologyRefiner> adaptive;

    Impl(ControlMesh m, SubdivisionOptions o)
        : mesh(std::move(m)),
          options(o) {
        verts_per_face.assign(mesh.num_quads(), 4);

        face_verts.reserve(mesh.num_quads() * 4);
        for (const std::array<int, 4>& q : mesh.quads()) {
            for (const int index : q) {
                face_verts.push_back(index);
            }
        }

        for (const ControlMesh::Crease& crease : mesh.creases()) {
            crease_pairs.push_back(crease.v0);
            crease_pairs.push_back(crease.v1);
            // Sharpness is float in OpenSubdiv's descriptor whatever the
            // precision of the stencils; it is a tag, not a coordinate.
            crease_weights.push_back(static_cast<float>(crease.sharpness));
        }
        for (const ControlMesh::Corner& corner : mesh.corners()) {
            corner_indices.push_back(corner.vertex);
            corner_weights.push_back(static_cast<float>(corner.sharpness));
        }

        adaptive = make_refiner();
        adaptive->RefineAdaptive(
            Far::TopologyRefiner::AdaptiveOptions(options.max_isolation_level));

        // The limit locations below address faces by ptex index. For an
        // all-quad mesh each face is exactly one ptex face, in order, so the
        // two indices coincide -- which is why ControlMesh is quads only.
        const Far::PtexIndices ptex(*adaptive);
        if (ptex.GetNumFaces() != static_cast<int>(mesh.num_quads())) {
            throw std::runtime_error(
                fmt::format("expected one ptex face per quad, got {} ptex faces for {} quads",
                            ptex.GetNumFaces(),
                            mesh.num_quads()));
        }
    }

    /// A refiner can only be refined once, and adaptive and uniform refinement
    /// are mutually exclusive on one instance, so uniform refinement builds its
    /// own.
    std::unique_ptr<Far::TopologyRefiner> make_refiner() const {
        Far::TopologyDescriptor descriptor;
        descriptor.numVertices = static_cast<int>(mesh.num_vertices());
        descriptor.numFaces = static_cast<int>(mesh.num_quads());
        descriptor.numVertsPerFace = verts_per_face.data();
        descriptor.vertIndicesPerFace = face_verts.data();

        if (!crease_weights.empty()) {
            descriptor.numCreases = static_cast<int>(crease_weights.size());
            descriptor.creaseVertexIndexPairs = crease_pairs.data();
            descriptor.creaseWeights = crease_weights.data();
        }
        if (!corner_weights.empty()) {
            descriptor.numCorners = static_cast<int>(corner_weights.size());
            descriptor.cornerVertexIndices = corner_indices.data();
            descriptor.cornerWeights = corner_weights.data();
        }

        Sdc::Options sdc_options;
        sdc_options.SetVtxBoundaryInterpolation(to_sdc(options.boundary));

        using Factory = Far::TopologyRefinerFactory<Far::TopologyDescriptor>;
        Factory::Options factory_options(Sdc::SCHEME_CATMARK, sdc_options);

        std::unique_ptr<Far::TopologyRefiner> refiner{Factory::Create(descriptor, factory_options)};
        if (!refiner) {
            throw std::runtime_error(
                "OpenSubdiv rejected the control mesh topology. The usual cause is a "
                "non-manifold configuration: an edge shared by more than two faces, or "
                "faces meeting at a vertex in more than one fan.");
        }
        return refiner;
    }

    /// Builds the limit stencil table for the given locations. One
    /// LocationArray per face, because a LocationArray carries a single ptex
    /// index.
    std::unique_ptr<const LimitStencilTableD>
    make_limit_table(const std::vector<LimitLocation>& locations, DerivativeOrder order) const {
        if (locations.empty()) {
            throw std::invalid_argument("no limit locations requested");
        }

        const auto face_count = static_cast<int>(mesh.num_quads());

        // The coordinates must outlive the factory call: LocationArray holds
        // bare pointers into them.
        std::vector<double> s;
        std::vector<double> t;
        s.reserve(locations.size());
        t.reserve(locations.size());

        std::vector<LimitStencilFactoryD::LocationArray> arrays;
        for (std::size_t i = 0; i < locations.size(); ++i) {
            const LimitLocation& location = locations[i];
            if (location.face < 0 || location.face >= face_count) {
                throw std::invalid_argument(
                    fmt::format("limit location {} names face {}, outside [0, {})",
                                i,
                                location.face,
                                face_count));
            }
            s.push_back(location.u);
            t.push_back(location.v);
        }

        // Group consecutive locations sharing a face into one array; callers
        // that sample face by face then pay for one array per face rather than
        // one per point.
        std::size_t run_start = 0;
        for (std::size_t i = 1; i <= locations.size(); ++i) {
            if (i == locations.size() || locations[i].face != locations[run_start].face) {
                LimitStencilFactoryD::LocationArray array;
                array.ptexIdx = locations[run_start].face;
                array.numLocations = static_cast<int>(i - run_start);
                array.s = s.data() + run_start;
                array.t = t.data() + run_start;
                arrays.push_back(array);
                run_start = i;
            }
        }

        LimitStencilFactoryD::Options factory_options;
        factory_options.generate1stDerivatives = order == DerivativeOrder::None ? 0u : 1u;
        factory_options.generate2ndDerivatives = order == DerivativeOrder::Second ? 1u : 0u;

        std::unique_ptr<const LimitStencilTableD> table{
            LimitStencilFactoryD::Create(*adaptive, arrays, nullptr, nullptr, factory_options)};
        if (!table) {
            throw std::runtime_error("OpenSubdiv could not build a limit stencil table for "
                                     "these locations");
        }
        return table;
    }
};

SubdivisionSurface::SubdivisionSurface(ControlMesh mesh, SubdivisionOptions options)
    : impl_(std::make_unique<Impl>(std::move(mesh), options)) {}

SubdivisionSurface::~SubdivisionSurface() = default;
SubdivisionSurface::SubdivisionSurface(SubdivisionSurface&&) noexcept = default;
SubdivisionSurface& SubdivisionSurface::operator=(SubdivisionSurface&&) noexcept = default;

const ControlMesh& SubdivisionSurface::control_mesh() const noexcept {
    return impl_->mesh;
}

const SubdivisionOptions& SubdivisionSurface::options() const noexcept {
    return impl_->options;
}

LimitMatrices SubdivisionSurface::build_limit_matrices(const std::vector<LimitLocation>& locations,
                                                       DerivativeOrder order) const {
    const std::unique_ptr<const LimitStencilTableD> table =
        impl_->make_limit_table(locations, order);

    const auto rows = static_cast<Eigen::Index>(locations.size());
    const auto columns = static_cast<Eigen::Index>(impl_->mesh.num_vertices());

    const std::vector<int>& sizes = table->GetSizes();
    const std::vector<Far::Index>& offsets = table->GetOffsets();
    const std::vector<Far::Index>& indices = table->GetControlIndices();

    if (static_cast<std::size_t>(table->GetNumStencils()) != locations.size()) {
        // OpenSubdiv silently drops locations that have no limit surface under
        // them rather than failing, so the count has to be checked. The usual
        // cause is BoundaryInterpolation::None, under which the limit surface
        // does not reach the boundary at all and every location on or near one
        // is discarded.
        throw std::runtime_error(fmt::format(
            "asked for {} limit points but OpenSubdiv produced {} stencils. Locations with no "
            "limit surface beneath them are dropped; with BoundaryInterpolation::None that is "
            "true of the whole boundary region, which is one reason this project defaults to "
            "EdgeAndCorner.",
            locations.size(),
            table->GetNumStencils()));
    }

    using Triplet = Eigen::Triplet<double>;

    // One pass per weight array, all sharing the same sparsity pattern.
    const auto gather = [&](const std::vector<double>& weights) {
        std::vector<Triplet> triplets;
        triplets.reserve(weights.size());
        for (Eigen::Index row = 0; row < rows; ++row) {
            const auto index = static_cast<std::size_t>(row);
            const int size = sizes[index];
            const int offset = offsets[index];
            for (int k = 0; k < size; ++k) {
                const auto at = static_cast<std::size_t>(offset + k);
                triplets.emplace_back(row, static_cast<Eigen::Index>(indices[at]), weights[at]);
            }
        }

        Eigen::SparseMatrix<double, Eigen::RowMajor> matrix(rows, columns);
        matrix.setFromTriplets(triplets.begin(), triplets.end());
        return matrix;
    };

    LimitMatrices matrices;
    matrices.order = order;
    matrices.position = gather(table->GetWeights());

    if (order != DerivativeOrder::None) {
        matrices.du = gather(table->GetDuWeights());
        matrices.dv = gather(table->GetDvWeights());
    }
    if (order == DerivativeOrder::Second) {
        matrices.duu = gather(table->GetDuuWeights());
        matrices.duv = gather(table->GetDuvWeights());
        matrices.dvv = gather(table->GetDvvWeights());
    }

    return matrices;
}

std::vector<LimitSample>
SubdivisionSurface::evaluate_limit(const std::vector<LimitLocation>& locations,
                                   DerivativeOrder order) const {
    // Evaluated through the same matrices the fitting uses, rather than through
    // a separate PatchTable path. One code path means the matrix and the
    // evaluator cannot disagree, and a disagreement between them would be
    // invisible until it showed up as an unexplained fitting residual.
    const DerivativeOrder effective =
        order == DerivativeOrder::None ? DerivativeOrder::First : order;
    const LimitMatrices matrices = build_limit_matrices(locations, effective);
    const Eigen::MatrixXd points = control_point_matrix(impl_->mesh);

    const Eigen::MatrixXd positions = matrices.position * points;
    const Eigen::MatrixXd du = matrices.du * points;
    const Eigen::MatrixXd dv = matrices.dv * points;

    Eigen::MatrixXd duu;
    Eigen::MatrixXd duv;
    Eigen::MatrixXd dvv;
    if (effective == DerivativeOrder::Second) {
        duu = matrices.duu * points;
        duv = matrices.duv * points;
        dvv = matrices.dvv * points;
    }

    std::vector<LimitSample> samples;
    samples.reserve(locations.size());
    for (Eigen::Index row = 0; row < positions.rows(); ++row) {
        LimitSample sample;
        sample.position = positions.row(row);
        sample.du = du.row(row);
        sample.dv = dv.row(row);

        const Eigen::Vector3d cross = sample.du.cross(sample.dv);
        const double length = cross.norm();
        sample.normal = length > 0.0 ? Eigen::Vector3d{cross / length} : Eigen::Vector3d::Zero();

        if (effective == DerivativeOrder::Second) {
            sample.duu = duu.row(row);
            sample.duv = duv.row(row);
            sample.dvv = dvv.row(row);
        }

        samples.push_back(sample);
    }
    return samples;
}

LimitSample SubdivisionSurface::evaluate_limit(const LimitLocation& location,
                                               DerivativeOrder order) const {
    return evaluate_limit(std::vector<LimitLocation>{location}, order).front();
}

std::vector<Eigen::Vector3d> SubdivisionSurface::limit_positions_of_vertices() const {
    // Each vertex is addressed through one of its incident faces, at that
    // face's corner. Which face is immaterial: the limit point is a property of
    // the vertex, and the regular and valence-n oracles check exactly that.
    std::vector<LimitLocation> locations(impl_->mesh.num_vertices(), LimitLocation{-1, 0.0, 0.0});

    const std::vector<std::array<int, 4>>& quads = impl_->mesh.quads();
    for (std::size_t f = 0; f < quads.size(); ++f) {
        static constexpr double kCorners[4][2] = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
        for (std::size_t corner = 0; corner < 4; ++corner) {
            const auto vertex = static_cast<std::size_t>(quads[f][corner]);
            if (locations[vertex].face < 0) {
                locations[vertex] =
                    LimitLocation{static_cast<int>(f), kCorners[corner][0], kCorners[corner][1]};
            }
        }
    }

    for (std::size_t i = 0; i < locations.size(); ++i) {
        if (locations[i].face < 0) {
            throw std::runtime_error(
                fmt::format("vertex {} belongs to no quad, so it has no limit point", i));
        }
    }

    const std::vector<LimitSample> samples = evaluate_limit(locations);
    std::vector<Eigen::Vector3d> positions;
    positions.reserve(samples.size());
    for (const LimitSample& sample : samples) {
        positions.push_back(sample.position);
    }
    return positions;
}

Eigen::MatrixXd SubdivisionSurface::local_subdivision_matrix(int vertex) const {
    const ControlMesh& mesh = impl_->mesh;
    if (vertex < 0 || vertex >= static_cast<int>(mesh.num_vertices())) {
        throw std::invalid_argument(
            fmt::format("vertex {} is outside [0, {})", vertex, mesh.num_vertices()));
    }
    if (mesh.is_boundary_vertex(vertex)) {
        throw std::invalid_argument(
            fmt::format("vertex {} is on a boundary, where the 1-ring is not closed and the local "
                        "subdivision matrix is undefined",
                        vertex));
    }

    const std::unique_ptr<Far::TopologyRefiner> refiner = impl_->make_refiner();
    refiner->RefineUniform(Far::TopologyRefiner::UniformOptions(1));

    const Far::TopologyLevel& base = refiner->GetLevel(0);
    const OpenSubdiv::Far::ConstIndexArray edges = base.GetVertexEdges(vertex);
    const OpenSubdiv::Far::ConstIndexArray faces = base.GetVertexFaces(vertex);

    if (edges.size() != faces.size()) {
        throw std::runtime_error(
            fmt::format("vertex {} has {} edges but {} faces, so its 1-ring is not closed",
                        vertex,
                        edges.size(),
                        faces.size()));
    }
    const int valence = edges.size();
    const auto ring_size = static_cast<Eigen::Index>(2 * valence + 1);

    // Source indices at level 0: the centre, the far end of each incident edge,
    // then the vertex diagonally opposite across each incident face.
    std::vector<int> source;
    source.reserve(static_cast<std::size_t>(ring_size));
    source.push_back(vertex);
    for (int i = 0; i < valence; ++i) {
        const OpenSubdiv::Far::ConstIndexArray ends = base.GetEdgeVertices(edges[i]);
        source.push_back(ends[0] == vertex ? ends[1] : ends[0]);
    }
    for (int i = 0; i < valence; ++i) {
        const OpenSubdiv::Far::ConstIndexArray corners = base.GetFaceVertices(faces[i]);
        // On a quad the diagonal is two steps around from the centre.
        int at = 0;
        for (int c = 0; c < corners.size(); ++c) {
            if (corners[c] == vertex) {
                at = c;
            }
        }
        source.push_back(corners[(at + 2) % corners.size()]);
    }

    // Target indices at level 1, under the natural correspondence: the centre
    // maps to its child, each incident edge to its edge child (the new
    // 1-ring), each incident face to its face child (the new diagonals).
    std::vector<int> target;
    target.reserve(static_cast<std::size_t>(ring_size));
    target.push_back(base.GetVertexChildVertex(vertex));
    for (int i = 0; i < valence; ++i) {
        target.push_back(base.GetEdgeChildVertex(edges[i]));
    }
    for (int i = 0; i < valence; ++i) {
        target.push_back(base.GetFaceChildVertex(faces[i]));
    }

    // Refine the identity: every base vertex carries its own unit vector, so
    // after one step each level-1 vertex holds its row of the global
    // subdivision matrix. The local matrix is then a submatrix of that, which
    // means it comes from OpenSubdiv's stencils rather than from any formula
    // written here.
    const auto base_count = static_cast<Eigen::Index>(mesh.num_vertices());
    std::vector<WeightRow> buffer(static_cast<std::size_t>(refiner->GetNumVerticesTotal()),
                                  WeightRow{Eigen::VectorXd::Zero(base_count)});
    for (Eigen::Index i = 0; i < base_count; ++i) {
        buffer[static_cast<std::size_t>(i)].weights(i) = 1.0;
    }

    Far::PrimvarRefinerReal<double> primvar(*refiner);
    WeightRow* level_0 = buffer.data();
    WeightRow* level_1 = level_0 + base.GetNumVertices();
    primvar.Interpolate(1, level_0, level_1);

    Eigen::MatrixXd matrix(ring_size, ring_size);
    for (Eigen::Index row = 0; row < ring_size; ++row) {
        const WeightRow& refined = level_1[target[static_cast<std::size_t>(row)]];
        for (Eigen::Index column = 0; column < ring_size; ++column) {
            matrix(row, column) = refined.weights(source[static_cast<std::size_t>(column)]);
        }
    }
    return matrix;
}

PolyMesh SubdivisionSurface::refine_uniform(int level) const {
    if (level < 0) {
        throw std::invalid_argument(
            fmt::format("refinement level must not be negative, got {}", level));
    }
    if (level == 0) {
        return PolyMesh{impl_->mesh.vertices(), impl_->mesh.quads()};
    }

    const std::unique_ptr<Far::TopologyRefiner> refiner = impl_->make_refiner();
    refiner->RefineUniform(Far::TopologyRefiner::UniformOptions(level));

    // Positions for every level, laid out one level after another.
    std::vector<RefinablePoint> buffer(static_cast<std::size_t>(refiner->GetNumVerticesTotal()));
    for (std::size_t i = 0; i < impl_->mesh.num_vertices(); ++i) {
        buffer[i].position = impl_->mesh.vertices()[i];
    }

    Far::PrimvarRefinerReal<double> primvar(*refiner);
    RefinablePoint* source = buffer.data();
    for (int l = 1; l <= level; ++l) {
        RefinablePoint* destination = source + refiner->GetLevel(l - 1).GetNumVertices();
        primvar.Interpolate(l, source, destination);
        source = destination;
    }

    const Far::TopologyLevel& last = refiner->GetLevel(level);

    PolyMesh mesh;
    mesh.vertices.reserve(static_cast<std::size_t>(last.GetNumVertices()));
    for (int i = 0; i < last.GetNumVertices(); ++i) {
        mesh.vertices.push_back(source[i].position);
    }

    mesh.quads.reserve(static_cast<std::size_t>(last.GetNumFaces()));
    for (int f = 0; f < last.GetNumFaces(); ++f) {
        const OpenSubdiv::Far::ConstIndexArray face = last.GetFaceVertices(f);
        // Catmull-Clark produces only quads from level 1 onwards.
        if (face.size() != 4) {
            throw std::runtime_error(fmt::format(
                "expected quads after Catmull-Clark refinement, face {} has {} vertices",
                f,
                face.size()));
        }
        mesh.quads.push_back({face[0], face[1], face[2], face[3]});
    }

    return mesh;
}

TessellatedLimit tessellate_limit(const SubdivisionSurface& surface, int samples_per_edge) {
    if (samples_per_edge < 1) {
        throw std::invalid_argument(
            fmt::format("samples_per_edge must be at least 1, got {}", samples_per_edge));
    }

    const std::size_t faces = surface.control_mesh().num_quads();
    const int side = samples_per_edge + 1;
    const auto per_face = static_cast<std::size_t>(side * side);

    TessellatedLimit result;
    result.locations.reserve(faces * per_face);

    for (std::size_t f = 0; f < faces; ++f) {
        for (int a = 0; a < side; ++a) {
            for (int b = 0; b < side; ++b) {
                result.locations.push_back(
                    LimitLocation{static_cast<int>(f),
                                  static_cast<double>(a) / static_cast<double>(samples_per_edge),
                                  static_cast<double>(b) / static_cast<double>(samples_per_edge)});
            }
        }
    }

    const std::vector<LimitSample> samples = surface.evaluate_limit(result.locations);

    result.mesh.vertices.reserve(samples.size());
    result.normals.reserve(samples.size());
    for (const LimitSample& sample : samples) {
        result.mesh.vertices.push_back(sample.position);
        result.normals.push_back(sample.normal);
    }

    result.mesh.quads.reserve(faces *
                              static_cast<std::size_t>(samples_per_edge * samples_per_edge));
    for (std::size_t f = 0; f < faces; ++f) {
        const auto base = static_cast<int>(f * per_face);
        for (int a = 0; a < samples_per_edge; ++a) {
            for (int b = 0; b < samples_per_edge; ++b) {
                const int corner = base + a * side + b;
                result.mesh.quads.push_back({corner, corner + side, corner + side + 1, corner + 1});
            }
        }
    }

    return result;
}

Eigen::MatrixXd control_point_matrix(const ControlMesh& mesh) {
    Eigen::MatrixXd rows(static_cast<Eigen::Index>(mesh.num_vertices()), 3);
    for (std::size_t i = 0; i < mesh.num_vertices(); ++i) {
        rows.row(static_cast<Eigen::Index>(i)) = mesh.vertices()[i];
    }
    return rows;
}

std::vector<Eigen::Vector3d> to_points(const Eigen::MatrixXd& rows) {
    if (rows.cols() != 3) {
        throw std::invalid_argument(
            fmt::format("expected a matrix with 3 columns, got {}", rows.cols()));
    }
    std::vector<Eigen::Vector3d> points;
    points.reserve(static_cast<std::size_t>(rows.rows()));
    for (Eigen::Index i = 0; i < rows.rows(); ++i) {
        points.emplace_back(rows(i, 0), rows(i, 1), rows(i, 2));
    }
    return points;
}

} // namespace n2s
