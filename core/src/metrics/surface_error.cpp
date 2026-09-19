#include "n2s/metrics/surface_error.hpp"

#include "n2s/nurbs/differential.hpp"
#include "n2s/tolerances.hpp"
#include "n2s/trim/sampling.hpp"

#include <Eigen/Geometry>
#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace n2s::metrics {

namespace {

/// Angle between two unit vectors, in degrees, computed from the cross product
/// rather than the dot product.
///
/// `acos(dot)` loses precision exactly where this metric matters most: for
/// nearly parallel normals the dot product is 1 - O(theta^2), so half the
/// digits of a small angle are gone before the arccos is taken. `atan2` of the
/// cross and dot magnitudes is accurate across the whole range.
double angle_degrees(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    const double cross = a.cross(b).norm();
    const double dot = a.dot(b);
    return std::atan2(cross, dot) * 180.0 / std::numbers::pi;
}

/// Closest distance from a point to a polyline, measured against the segments
/// rather than only the vertices.
double distance_to_polyline(const Eigen::Vector3d& point,
                            const std::vector<Eigen::Vector3d>& polyline,
                            bool closed) {
    if (polyline.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (polyline.size() == 1) {
        return (point - polyline.front()).norm();
    }

    double best = std::numeric_limits<double>::infinity();
    const std::size_t last = closed ? polyline.size() : polyline.size() - 1;
    for (std::size_t i = 0; i < last; ++i) {
        const Eigen::Vector3d& a = polyline[i];
        const Eigen::Vector3d& b = polyline[(i + 1) % polyline.size()];

        const Eigen::Vector3d along = b - a;
        const double length_squared = along.squaredNorm();
        double distance = 0.0;
        if (length_squared < tol::kDegenerateDerivative * tol::kDegenerateDerivative) {
            distance = (point - a).norm();
        } else {
            const double s = std::clamp((point - a).dot(along) / length_squared, 0.0, 1.0);
            distance = (point - (a + s * along)).norm();
        }
        best = std::min(best, distance);
    }
    return best;
}

/// Closest point on a triangle to `point` (Ericson, Real-Time Collision
/// Detection, section 5.1.5): checks the three vertex regions, the three edge
/// regions, then the interior.
Eigen::Vector3d closest_point_on_triangle(const Eigen::Vector3d& point,
                                          const Eigen::Vector3d& a,
                                          const Eigen::Vector3d& b,
                                          const Eigen::Vector3d& c) {
    const Eigen::Vector3d ab = b - a;
    const Eigen::Vector3d ac = c - a;
    const Eigen::Vector3d ap = point - a;

    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0) {
        return a;
    }

    const Eigen::Vector3d bp = point - b;
    const double d3 = ab.dot(bp);
    const double d4 = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3) {
        return b;
    }

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        return a + (d1 / (d1 - d3)) * ab;
    }

    const Eigen::Vector3d cp = point - c;
    const double d5 = ab.dot(cp);
    const double d6 = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6) {
        return c;
    }

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        return a + (d2 / (d2 - d6)) * ac;
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        return b + ((d4 - d3) / ((d4 - d3) + (d5 - d6))) * (c - b);
    }

    const double denom = 1.0 / (va + vb + vc);
    return a + ab * (vb * denom) + ac * (vc * denom);
}

/// Distance from a point to a quad mesh, measured against the faces rather
/// than the vertices.
///
/// Measuring to the nearest vertex instead puts a floor of half the sample
/// spacing under every result, so the reverse Hausdorff distance would report
/// the tessellation density rather than the coverage error it is supposed to
/// measure. That floor was 0.07 on the first version of this code, against a
/// true distance of zero.
double distance_to_mesh(const Eigen::Vector3d& point, const PolyMesh& mesh) {
    double best = std::numeric_limits<double>::infinity();
    for (const std::array<int, 4>& quad : mesh.quads) {
        const Eigen::Vector3d& a = mesh.vertices[static_cast<std::size_t>(quad[0])];
        const Eigen::Vector3d& b = mesh.vertices[static_cast<std::size_t>(quad[1])];
        const Eigen::Vector3d& c = mesh.vertices[static_cast<std::size_t>(quad[2])];
        const Eigen::Vector3d& d = mesh.vertices[static_cast<std::size_t>(quad[3])];

        best = std::min(best, (point - closest_point_on_triangle(point, a, b, c)).squaredNorm());
        best = std::min(best, (point - closest_point_on_triangle(point, a, c, d)).squaredNorm());
    }
    return std::sqrt(best);
}

/// The four boundary edges of a quad face in `(u, v)`, as a parameter walk.
Eigen::Vector2d edge_point(int edge, double t) {
    switch (edge) {
    case 0:
        return {t, 0.0};
    case 1:
        return {1.0, t};
    case 2:
        return {1.0 - t, 1.0};
    default:
        return {0.0, 1.0 - t};
    }
}

} // namespace

DomainMap grid_domain_map(int rows, int columns) {
    if (rows < 2 || columns < 2) {
        throw std::invalid_argument(fmt::format(
            "a grid domain map needs at least a 2 x 2 grid, got {} x {}", rows, columns));
    }

    const int faces_u = rows - 1;
    const int faces_v = columns - 1;

    return [faces_u, faces_v](const LimitLocation& location) {
        // ControlMesh::grid emits faces in row-major order with u as the slow
        // index, so face f sits at (f / faces_v, f % faces_v).
        const int i = location.face / faces_v;
        const int j = location.face % faces_v;

        return Eigen::Vector2d{(static_cast<double>(i) + location.u) / static_cast<double>(faces_u),
                               (static_cast<double>(j) + location.v) /
                                   static_cast<double>(faces_v)};
    };
}

std::vector<ErrorSample> sample_surface_error(const SubdivisionSurface& limit,
                                              const NurbsSurface& nurbs,
                                              const TrimRegion& region,
                                              const DomainMap& domain_map,
                                              const SurfaceErrorOptions& options) {
    if (options.samples_per_face < 1) {
        throw std::invalid_argument(
            fmt::format("samples_per_face must be at least 1, got {}", options.samples_per_face));
    }

    const TessellatedLimit tessellation = tessellate_limit(limit, options.samples_per_face);
    const std::vector<LimitSample> samples =
        limit.evaluate_limit(tessellation.locations, DerivativeOrder::First);

    // Polygonise the trim loops once for the containment test; TrimRegion
    // re-evaluates every curve per call, which is ruinous in a loop this long.
    const std::vector<Eigen::Vector2d> outer = region.outer().polygonise();
    std::vector<std::vector<Eigen::Vector2d>> holes;
    holes.reserve(region.holes().size());
    for (const TrimLoop& hole : region.holes()) {
        holes.push_back(hole.polygonise());
    }
    const auto inside = [&](const Eigen::Vector2d& p) {
        if (!point_in_polygon(outer, p)) {
            return false;
        }
        return std::none_of(
            holes.begin(), holes.end(), [&](const std::vector<Eigen::Vector2d>& hole) {
                return point_in_polygon(hole, p);
            });
    };

    std::vector<ErrorSample> result;
    result.reserve(samples.size());

    for (std::size_t i = 0; i < samples.size(); ++i) {
        ErrorSample record;
        record.location = tessellation.locations[i];
        record.limit_point = samples[i].position;

        if (domain_map) {
            record.domain_point = domain_map(record.location);

            if (options.restrict_to_trimmed_region && !inside(record.domain_point)) {
                record.measured = false;
                result.push_back(record);
                continue;
            }

            const Eigen::Vector3d corresponding =
                nurbs.evaluate(record.domain_point.x(), record.domain_point.y());
            record.parametric = (record.limit_point - corresponding).norm();

            const SurfaceDerivatives derivatives =
                nurbs.derivatives(record.domain_point.x(), record.domain_point.y(), 1);
            const std::optional<Eigen::Vector3d> nurbs_normal =
                normal_from_derivatives(derivatives.du(), derivatives.dv());
            const std::optional<Eigen::Vector3d> limit_normal =
                normal_from_derivatives(samples[i].du, samples[i].dv);

            if (nurbs_normal.has_value() && limit_normal.has_value()) {
                // The two surfaces may be parameterised with opposite
                // orientation, which would report every normal as 180 degrees
                // out. The deviation wanted is between the tangent planes, so
                // the smaller of the angle and its supplement is the answer.
                const double raw = angle_degrees(*limit_normal, *nurbs_normal);
                record.normal_degrees = std::min(raw, 180.0 - raw);
            } else {
                record.measured = false;
            }
        }

        const SurfaceProjection projection =
            project_to_surface(nurbs, record.limit_point, options.projection);
        if (projection.converged()) {
            record.geometric = projection.distance;
        } else {
            record.measured = false;
        }

        result.push_back(record);
    }

    return result;
}

SurfaceError measure_surface_error(const SubdivisionSurface& limit,
                                   const NurbsSurface& nurbs,
                                   const TrimRegion& region,
                                   const DomainMap& domain_map,
                                   const SurfaceErrorOptions& options) {
    SurfaceError error;

    Accumulator parametric;
    Accumulator geometric;
    Accumulator normals;
    Accumulator mean_curvature;
    Accumulator gaussian_curvature;

    const std::vector<ErrorSample> samples =
        sample_surface_error(limit, nurbs, region, domain_map, options);

    for (const ErrorSample& sample : samples) {
        if (!sample.measured) {
            geometric.add_unmeasured();
            if (domain_map) {
                parametric.add_unmeasured();
                normals.add_unmeasured();
            }
            continue;
        }
        geometric.add(sample.geometric);
        if (domain_map) {
            parametric.add(sample.parametric);
            normals.add(sample.normal_degrees);
        }
    }

    // Curvature deviation needs second derivatives on both surfaces, so it is
    // a second pass rather than folded into the first.
    if (domain_map) {
        const TessellatedLimit tessellation = tessellate_limit(limit, options.samples_per_face);
        const std::vector<LimitSample> second =
            limit.evaluate_limit(tessellation.locations, DerivativeOrder::Second);

        for (std::size_t i = 0; i < second.size(); ++i) {
            if (!samples[i].measured) {
                mean_curvature.add_unmeasured();
                gaussian_curvature.add_unmeasured();
                continue;
            }

            const Eigen::Vector2d& domain_point = samples[i].domain_point;
            const std::optional<SurfaceCurvature> nurbs_curvature =
                surface_curvature(nurbs.derivatives(domain_point.x(), domain_point.y(), 2));
            const std::optional<SurfaceCurvature> limit_curvature = curvature_from_derivatives(
                second[i].du, second[i].dv, second[i].duu, second[i].duv, second[i].dvv);

            if (!nurbs_curvature.has_value() || !limit_curvature.has_value()) {
                mean_curvature.add_unmeasured();
                gaussian_curvature.add_unmeasured();
                continue;
            }

            // Magnitudes: the two surfaces may carry opposite normal
            // orientation, which flips the sign of mean curvature without any
            // geometric difference. Gaussian curvature is orientation-free and
            // is compared signed.
            mean_curvature.add(
                std::abs(std::abs(limit_curvature->mean) - std::abs(nurbs_curvature->mean)));
            gaussian_curvature.add(std::abs(limit_curvature->gaussian - nurbs_curvature->gaussian));
        }
    }

    // The other direction: NURBS samples to the limit surface, measured
    // against the tessellated limit rather than by projecting onto the true
    // limit surface. Still an upper bound, but a tight one: measuring to the
    // faces makes the discretisation error second order in the sample spacing
    // instead of first.
    //
    // Cost is the product of the sample count and the face count, with no
    // spatial index. Fine at the densities the experiments use; if the reverse
    // direction ever dominates a run, that is the thing to fix.
    const TessellatedLimit tessellation = tessellate_limit(limit, options.samples_per_face);
    Accumulator reverse;

    const int side = std::max(2, options.nurbs_samples_per_side);
    const std::vector<Eigen::Vector2d> outer = region.outer().polygonise();
    std::vector<std::vector<Eigen::Vector2d>> holes;
    for (const TrimLoop& hole : region.holes()) {
        holes.push_back(hole.polygonise());
    }

    for (int i = 0; i <= side; ++i) {
        for (int j = 0; j <= side; ++j) {
            const Eigen::Vector2d domain_point{static_cast<double>(i) / static_cast<double>(side),
                                               static_cast<double>(j) / static_cast<double>(side)};
            if (!point_in_polygon(outer, domain_point)) {
                continue;
            }
            const bool in_hole = std::any_of(
                holes.begin(), holes.end(), [&](const std::vector<Eigen::Vector2d>& hole) {
                    return point_in_polygon(hole, domain_point);
                });
            if (in_hole) {
                continue;
            }

            const Eigen::Vector3d target = nurbs.evaluate(domain_point.x(), domain_point.y());
            reverse.add(distance_to_mesh(target, tessellation.mesh));
        }
    }

    error.parametric = parametric.result();
    error.geometric = geometric.result();
    error.reverse_geometric = reverse.result();
    error.normal_degrees = normals.result();
    error.mean_curvature = mean_curvature.result();
    error.gaussian_curvature = gaussian_curvature.result();
    error.hausdorff = std::max(error.geometric.max, error.reverse_geometric.max);

    return error;
}

Stats measure_boundary_error(const SubdivisionSurface& limit,
                             const NurbsSurface& nurbs,
                             const TrimRegion& region,
                             const BoundaryErrorOptions& options) {
    if (options.samples_per_edge < 1) {
        throw std::invalid_argument(
            fmt::format("samples_per_edge must be at least 1, got {}", options.samples_per_edge));
    }

    // The trim curves' image on the surface, as dense polylines.
    SamplingOptions trim_sampling;
    trim_sampling.mode = SamplingMode::Uniform;
    trim_sampling.samples_per_span = std::max(1, options.trim_samples_per_curve);

    std::vector<std::vector<Eigen::Vector3d>> trim_images;
    const auto add_loop = [&](const TrimLoop& loop) {
        const std::vector<Eigen::Vector2d> points = sample_loop(loop, nurbs, trim_sampling);
        std::vector<Eigen::Vector3d> image;
        image.reserve(points.size());
        for (const Eigen::Vector2d& p : points) {
            image.push_back(nurbs.evaluate(p.x(), p.y()));
        }
        trim_images.push_back(std::move(image));
    };
    add_loop(region.outer());
    for (const TrimLoop& hole : region.holes()) {
        add_loop(hole);
    }

    // Every boundary edge of the control mesh, i.e. every face edge used by
    // exactly one face.
    const ControlMesh& mesh = limit.control_mesh();
    std::vector<LimitLocation> boundary_samples;

    for (std::size_t f = 0; f < mesh.num_quads(); ++f) {
        const std::array<int, 4>& quad = mesh.quads()[f];
        for (int edge = 0; edge < 4; ++edge) {
            const int a = quad[static_cast<std::size_t>(edge)];
            const int b = quad[static_cast<std::size_t>((edge + 1) % 4)];

            int uses = 0;
            for (const std::array<int, 4>& other : mesh.quads()) {
                for (int k = 0; k < 4; ++k) {
                    const int p = other[static_cast<std::size_t>(k)];
                    const int q = other[static_cast<std::size_t>((k + 1) % 4)];
                    if ((p == a && q == b) || (p == b && q == a)) {
                        ++uses;
                    }
                }
            }
            if (uses != 1) {
                continue; // interior edge
            }

            for (int s = 0; s <= options.samples_per_edge; ++s) {
                const double t =
                    static_cast<double>(s) / static_cast<double>(options.samples_per_edge);
                const Eigen::Vector2d local = edge_point(edge, t);
                boundary_samples.push_back(
                    LimitLocation{static_cast<int>(f), local.x(), local.y()});
            }
        }
    }

    Accumulator accumulator;
    if (boundary_samples.empty()) {
        // A closed control mesh has no boundary, so there is nothing to
        // measure. Reported as an empty statistic rather than as zero, which
        // would read as a perfect result.
        return accumulator.result();
    }

    const std::vector<LimitSample> samples =
        limit.evaluate_limit(boundary_samples, DerivativeOrder::None);

    for (const LimitSample& sample : samples) {
        double best = std::numeric_limits<double>::infinity();
        for (const std::vector<Eigen::Vector3d>& image : trim_images) {
            best = std::min(best, distance_to_polyline(sample.position, image, true));
        }
        accumulator.add(best);
    }

    return accumulator.result();
}

} // namespace n2s::metrics
