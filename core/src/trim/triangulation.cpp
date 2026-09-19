#include "n2s/trim/triangulation.hpp"

#include "n2s/tolerances.hpp"

#include <CDT.h>
#include <Eigen/Geometry> // Vector3d::cross
#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace n2s {

namespace {

double triangle_area(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c) {
    const Eigen::Vector2d ab = b - a;
    const Eigen::Vector2d ac = c - a;
    return 0.5 * std::abs(ab.x() * ac.y() - ab.y() * ac.x());
}

/// Appends one sampled loop as vertices plus a closed ring of constraint edges.
void append_loop(const TrimLoop& loop,
                 const NurbsSurface& surface,
                 const SamplingOptions& sampling,
                 std::vector<CDT::V2d<double>>& vertices,
                 std::vector<CDT::Edge>& edges) {
    const std::vector<Eigen::Vector2d> points = sample_loop(loop, surface, sampling);
    if (points.size() < 3) {
        throw std::runtime_error(
            fmt::format("a trim loop sampled to only {} points, too few to bound a region; "
                        "loosen nothing and check the loop is not degenerate",
                        points.size()));
    }

    const auto first = static_cast<CDT::VertInd>(vertices.size());
    for (const Eigen::Vector2d& p : points) {
        vertices.emplace_back(p.x(), p.y());
    }

    // The ring closes back onto `first`; that last edge is what makes the loop
    // a boundary rather than an open polyline, and without it CDT cannot tell
    // inside from outside.
    const auto count = static_cast<CDT::VertInd>(points.size());
    for (CDT::VertInd i = 0; i < count; ++i) {
        edges.emplace_back(first + i, first + (i + 1) % count);
    }
}

} // namespace

double DomainMesh::area() const {
    double total = 0.0;
    for (const std::array<std::size_t, 3>& t : triangles) {
        total += triangle_area(vertices[t[0]], vertices[t[1]], vertices[t[2]]);
    }
    return total;
}

double SurfaceMesh::area() const {
    double total = 0.0;
    for (const std::array<std::size_t, 3>& t : triangles) {
        const Eigen::Vector3d ab = vertices[t[1]] - vertices[t[0]];
        const Eigen::Vector3d ac = vertices[t[2]] - vertices[t[0]];
        total += 0.5 * ab.cross(ac).norm();
    }
    return total;
}

namespace {

/// Shared by `triangulate` and the refinement rounds: constrain to the sampled
/// loops, add any accumulated interior points, triangulate, drop the outside
/// and the holes.
DomainMesh triangulate_with(const TrimRegion& region,
                            const NurbsSurface& surface,
                            const SamplingOptions& sampling,
                            const std::vector<Eigen::Vector2d>& interior) {
    std::vector<CDT::V2d<double>> vertices;
    std::vector<CDT::Edge> edges;

    append_loop(region.outer(), surface, sampling, vertices, edges);
    for (const TrimLoop& hole : region.holes()) {
        append_loop(hole, surface, sampling, vertices, edges);
    }

    const std::size_t boundary_count = vertices.size();

    // Interior points carry no constraint edges: they are free to be moved
    // around by the Delaunay criterion, unlike the boundary.
    for (const Eigen::Vector2d& p : interior) {
        vertices.emplace_back(p.x(), p.y());
    }

    // Two sampled loops can land a vertex on the same point -- a hole touching
    // the outer boundary, say. CDT rejects duplicates outright, so they are
    // merged and the edges remapped before insertion rather than after it
    // fails.
    CDT::RemoveDuplicatesAndRemapEdges(vertices, edges);

    CDT::Triangulation<double> cdt;
    cdt.insertVertices(vertices);
    cdt.insertEdges(edges);
    cdt.eraseOuterTrianglesAndHoles();

    if (cdt.triangles.empty()) {
        throw std::runtime_error(
            "the constrained triangulation came back empty. The usual cause is trim loops that "
            "do not close, so CDT cannot separate inside from outside; run validate_and_repair "
            "and check its report first.");
    }

    DomainMesh mesh;
    mesh.boundary_vertex_count = std::min(boundary_count, cdt.vertices.size());

    mesh.vertices.reserve(cdt.vertices.size());
    for (const CDT::V2d<double>& v : cdt.vertices) {
        mesh.vertices.emplace_back(v.x, v.y);
    }

    mesh.triangles.reserve(cdt.triangles.size());
    for (const CDT::Triangle& t : cdt.triangles) {
        mesh.triangles.push_back({static_cast<std::size_t>(t.vertices[0]),
                                  static_cast<std::size_t>(t.vertices[1]),
                                  static_cast<std::size_t>(t.vertices[2])});
    }

    return mesh;
}

} // namespace

void SurfaceMesh::write_obj(const std::filesystem::path& path) const {
    std::ofstream stream(path);
    if (!stream) {
        throw std::runtime_error(fmt::format("cannot open {} for writing", path.string()));
    }

    stream << "# nurbs2subd trimmed surface\n";
    for (const Eigen::Vector3d& v : vertices) {
        stream << fmt::format("v {:.17g} {:.17g} {:.17g}\n", v.x(), v.y(), v.z());
    }
    for (const std::array<std::size_t, 3>& t : triangles) {
        // OBJ indices are 1-based.
        stream << fmt::format("f {} {} {}\n", t[0] + 1, t[1] + 1, t[2] + 1);
    }
}

DomainMesh triangulate(const TrimRegion& region,
                       const NurbsSurface& surface,
                       const SamplingOptions& sampling) {
    return triangulate_with(region, surface, sampling, {});
}

SurfaceMesh map_to_surface(const DomainMesh& mesh, const NurbsSurface& surface) {
    SurfaceMesh lifted;
    lifted.triangles = mesh.triangles;
    lifted.boundary_vertex_count = mesh.boundary_vertex_count;

    lifted.vertices.reserve(mesh.vertices.size());
    for (const Eigen::Vector2d& p : mesh.vertices) {
        lifted.vertices.push_back(surface.evaluate(p.x(), p.y()));
    }

    return lifted;
}

namespace {

/// Smallest interior angle of a domain triangle, in degrees.
double smallest_angle_degrees(const Eigen::Vector2d& a,
                              const Eigen::Vector2d& b,
                              const Eigen::Vector2d& c) {
    const double ab = (b - a).norm();
    const double bc = (c - b).norm();
    const double ca = (a - c).norm();
    if (ab < tol::kMinParametricSpan || bc < tol::kMinParametricSpan ||
        ca < tol::kMinParametricSpan) {
        return 0.0;
    }

    // Law of cosines on each vertex; the smallest angle faces the shortest side,
    // so only that one needs computing.
    const double shortest = std::min({ab, bc, ca});
    double cosine = 0.0;
    if (shortest == ab) {
        cosine = (bc * bc + ca * ca - ab * ab) / (2.0 * bc * ca);
    } else if (shortest == bc) {
        cosine = (ab * ab + ca * ca - bc * bc) / (2.0 * ab * ca);
    } else {
        cosine = (ab * ab + bc * bc - ca * ca) / (2.0 * ab * bc);
    }

    return std::acos(std::clamp(cosine, -1.0, 1.0)) * 180.0 / std::numbers::pi;
}

/// Circumcentre of a domain triangle, or nullopt if it is degenerate.
std::optional<Eigen::Vector2d>
circumcentre(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c) {
    const double d =
        2.0 * (a.x() * (b.y() - c.y()) + b.x() * (c.y() - a.y()) + c.x() * (a.y() - b.y()));
    if (std::abs(d) < tol::kMinParametricSpan) {
        return std::nullopt;
    }

    const double a2 = a.squaredNorm();
    const double b2 = b.squaredNorm();
    const double c2 = c.squaredNorm();

    return Eigen::Vector2d{(a2 * (b.y() - c.y()) + b2 * (c.y() - a.y()) + c2 * (a.y() - b.y())) / d,
                           (a2 * (c.x() - b.x()) + b2 * (a.x() - c.x()) + c2 * (b.x() - a.x())) /
                               d};
}

double surface_triangle_area(const NurbsSurface& surface,
                             const Eigen::Vector2d& a,
                             const Eigen::Vector2d& b,
                             const Eigen::Vector2d& c) {
    const Eigen::Vector3d pa = surface.evaluate(a.x(), a.y());
    const Eigen::Vector3d pb = surface.evaluate(b.x(), b.y());
    const Eigen::Vector3d pc = surface.evaluate(c.x(), c.y());
    return 0.5 * (pb - pa).cross(pc - pa).norm();
}

} // namespace

DomainMesh triangulate_refined(const TrimRegion& region,
                               const NurbsSurface& surface,
                               const SamplingOptions& sampling,
                               const RefinementOptions& refinement) {
    std::vector<Eigen::Vector2d> interior;
    DomainMesh mesh = triangulate_with(region, surface, sampling, interior);

    const bool check_area = refinement.max_triangle_area > 0.0;
    const bool check_angle = refinement.min_angle_degrees > 0.0;
    if (!check_area && !check_angle) {
        return mesh;
    }

    // Candidates closer than this to an existing vertex are dropped. Without
    // it, a triangle that cannot be improved -- one pinned against the
    // boundary, say -- would be offered the same circumcentre every round and
    // the loop would never settle.
    const double crowding = 0.25;

    // Polygonise the loops once. TrimRegion::contains re-evaluates every trim
    // curve on every call, which is fine for a handful of queries and ruinous
    // here: this runs per candidate per round, and doing it the lazy way made
    // the test suite thirty times slower.
    const std::vector<Eigen::Vector2d> outer_polygon =
        region.outer().polygonise(kLoopTestSamplesPerCurve);
    std::vector<std::vector<Eigen::Vector2d>> hole_polygons;
    hole_polygons.reserve(region.holes().size());
    for (const TrimLoop& hole : region.holes()) {
        hole_polygons.push_back(hole.polygonise(kLoopTestSamplesPerCurve));
    }

    // Minimum spacing between vertices, derived from the budget so that the two
    // cannot contradict each other. This is the floor that stops the sliver
    // feedback loop described on RefinementOptions.
    const double min_separation =
        refinement.separation_fraction *
        std::sqrt(std::max(mesh.area(), tol::kMinLoopArea) /
                  static_cast<double>(std::max<std::size_t>(refinement.max_vertices, 1)));

    // Uniform grid over the domain at one cell per separation distance, so a
    // candidate only has to be compared against the nine cells around it rather
    // than against every vertex placed so far. The linear scan that preceded
    // this was quadratic in the vertex count and dominated the whole test suite.
    const auto cell_of = [&](const Eigen::Vector2d& p) {
        const auto x = static_cast<std::int64_t>(std::floor(p.x() / min_separation));
        const auto y = static_cast<std::int64_t>(std::floor(p.y() / min_separation));
        return std::make_pair(x, y);
    };
    const auto hash_cell = [](std::pair<std::int64_t, std::int64_t> c) {
        return c.first * 73856093LL ^ c.second * 19349663LL;
    };

    std::unordered_map<std::int64_t, std::vector<Eigen::Vector2d>> grid;
    const auto grid_insert = [&](const Eigen::Vector2d& p) {
        grid[hash_cell(cell_of(p))].push_back(p);
    };
    const auto too_close = [&](const Eigen::Vector2d& p) {
        const auto [cx, cy] = cell_of(p);
        for (std::int64_t dx = -1; dx <= 1; ++dx) {
            for (std::int64_t dy = -1; dy <= 1; ++dy) {
                const auto it = grid.find(hash_cell({cx + dx, cy + dy}));
                if (it == grid.end()) {
                    continue;
                }
                for (const Eigen::Vector2d& other : it->second) {
                    if ((other - p).norm() < min_separation) {
                        return true;
                    }
                }
            }
        }
        return false;
    };

    const auto inside_region = [&](const Eigen::Vector2d& p) {
        if (!point_in_polygon(outer_polygon, p)) {
            return false;
        }
        for (const std::vector<Eigen::Vector2d>& hole : hole_polygons) {
            if (point_in_polygon(hole, p)) {
                return false;
            }
        }
        return true;
    };

    for (int round = 0; round < refinement.max_rounds; ++round) {
        if (mesh.vertices.size() >= refinement.max_vertices) {
            break;
        }

        std::vector<Eigen::Vector2d> candidates;

        for (const std::array<std::size_t, 3>& t : mesh.triangles) {
            const Eigen::Vector2d& a = mesh.vertices[t[0]];
            const Eigen::Vector2d& b = mesh.vertices[t[1]];
            const Eigen::Vector2d& c = mesh.vertices[t[2]];

            const bool too_big = check_area && surface_triangle_area(surface, a, b, c) >
                                                   refinement.max_triangle_area;
            const bool too_thin =
                check_angle && smallest_angle_degrees(a, b, c) < refinement.min_angle_degrees;
            if (!too_big && !too_thin) {
                continue;
            }

            // The circumcentre is what gives a refinement its shape quality,
            // but it can land outside the triangle and so outside the region.
            // The centroid always lies inside, so it is the fallback.
            const std::optional<Eigen::Vector2d> centre = circumcentre(a, b, c);
            const Eigen::Vector2d centroid = (a + b + c) / 3.0;
            Eigen::Vector2d candidate = centre.value_or(centroid);
            if (!inside_region(candidate)) {
                candidate = centroid;
            }
            if (!inside_region(candidate)) {
                continue;
            }

            // Keep the candidate clear of the triangle's own corners, which is
            // what stops it crowding the fixed boundary discretisation.
            const double shortest_edge = std::min({(b - a).norm(), (c - b).norm(), (a - c).norm()});
            const double keep_away = crowding * shortest_edge;
            if ((candidate - a).norm() < keep_away || (candidate - b).norm() < keep_away ||
                (candidate - c).norm() < keep_away) {
                continue;
            }

            candidates.push_back(candidate);
        }

        if (candidates.empty()) {
            break;
        }

        // Rebuild the occupancy grid from the current mesh, so that candidates
        // are kept clear of the boundary vertices as well as of each other.
        grid.clear();
        for (const Eigen::Vector2d& v : mesh.vertices) {
            grid_insert(v);
        }

        std::size_t accepted = 0;
        for (const Eigen::Vector2d& candidate : candidates) {
            if (interior.size() + mesh.boundary_vertex_count >= refinement.max_vertices) {
                break;
            }
            if (too_close(candidate)) {
                continue;
            }
            grid_insert(candidate);
            interior.push_back(candidate);
            ++accepted;
        }

        // Every candidate was rejected for crowding: the mesh is as fine as the
        // separation floor allows and further rounds would change nothing.
        if (accepted == 0) {
            break;
        }

        mesh = triangulate_with(region, surface, sampling, interior);
    }

    return mesh;
}

} // namespace n2s
