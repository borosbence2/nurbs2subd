#include "n2s/metrics/seam_error.hpp"

#include "n2s/nurbs/differential.hpp"

#include <Eigen/Geometry>
#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace n2s::metrics {

namespace {

double angle_degrees(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return std::atan2(a.cross(b).norm(), a.dot(b)) * 180.0 / std::numbers::pi;
}

} // namespace

Eigen::Vector2d seam_parameter(const SeamEdge& edge, double t) {
    const double s = edge.reversed ? 1.0 - t : t;
    switch (edge.edge) {
    case FaceEdge::VMin:
        return {s, 0.0};
    case FaceEdge::UMax:
        return {1.0, s};
    case FaceEdge::VMax:
        return {s, 1.0};
    case FaceEdge::UMin:
        return {0.0, s};
    }
    return {s, 0.0};
}

std::vector<SeamSample> sample_seam(const SubdivisionSurface& a,
                                    const SeamEdge& edge_a,
                                    const SubdivisionSurface& b,
                                    const SeamEdge& edge_b,
                                    int samples) {
    if (samples < 1) {
        throw std::invalid_argument(
            fmt::format("a seam needs at least one interval, got {}", samples));
    }

    std::vector<LimitLocation> locations_a;
    std::vector<LimitLocation> locations_b;
    std::vector<double> parameters;
    locations_a.reserve(static_cast<std::size_t>(samples) + 1);
    locations_b.reserve(static_cast<std::size_t>(samples) + 1);
    parameters.reserve(static_cast<std::size_t>(samples) + 1);

    for (int i = 0; i <= samples; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(samples);
        parameters.push_back(t);

        const Eigen::Vector2d local_a = seam_parameter(edge_a, t);
        const Eigen::Vector2d local_b = seam_parameter(edge_b, t);
        locations_a.push_back(LimitLocation{edge_a.face, local_a.x(), local_a.y()});
        locations_b.push_back(LimitLocation{edge_b.face, local_b.x(), local_b.y()});
    }

    const std::vector<LimitSample> on_a = a.evaluate_limit(locations_a, DerivativeOrder::First);
    const std::vector<LimitSample> on_b = b.evaluate_limit(locations_b, DerivativeOrder::First);

    std::vector<SeamSample> result;
    result.reserve(parameters.size());

    for (std::size_t i = 0; i < parameters.size(); ++i) {
        SeamSample sample;
        sample.t = parameters[i];
        sample.point_a = on_a[i].position;
        sample.point_b = on_b[i].position;
        sample.gap = (sample.point_a - sample.point_b).norm();

        const std::optional<Eigen::Vector3d> normal_a =
            normal_from_derivatives(on_a[i].du, on_a[i].dv);
        const std::optional<Eigen::Vector3d> normal_b =
            normal_from_derivatives(on_b[i].du, on_b[i].dv);

        if (normal_a.has_value() && normal_b.has_value()) {
            // Two patches meeting along a seam may be parameterised with
            // opposite orientation, which would report a perfectly tangent
            // join as 180 degrees out. The quantity wanted is the angle
            // between the tangent planes.
            const double raw = angle_degrees(*normal_a, *normal_b);
            sample.normal_degrees = std::min(raw, 180.0 - raw);
        } else {
            sample.normal_degrees = std::numeric_limits<double>::quiet_NaN();
        }

        result.push_back(sample);
    }

    return result;
}

SeamError measure_seam(const SubdivisionSurface& a,
                       const SeamEdge& edge_a,
                       const SubdivisionSurface& b,
                       const SeamEdge& edge_b,
                       int samples) {
    Accumulator gap;
    Accumulator normals;

    for (const SeamSample& sample : sample_seam(a, edge_a, b, edge_b, samples)) {
        gap.add(sample.gap);
        // A NaN here is a point where a normal was undefined; the accumulator
        // counts it as unmeasured rather than poisoning the statistics.
        normals.add(sample.normal_degrees);
    }

    return SeamError{gap.result(), normals.result()};
}

} // namespace n2s::metrics
