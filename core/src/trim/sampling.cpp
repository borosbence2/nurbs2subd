#include "n2s/trim/sampling.hpp"

#include "n2s/tolerances.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace n2s {

namespace {

/// The trim curve lifted onto the surface: t -> S(C(t)).
Eigen::Vector3d on_surface(const NurbsCurve2& curve, const NurbsSurface& surface, double t) {
    const Eigen::Vector2d domain_point = curve.evaluate(t);
    return surface.evaluate(domain_point.x(), domain_point.y());
}

/// Distance from `p` to the segment ab. Used as the sagitta: how far the true
/// curve strays from the chord that is about to stand in for it.
double
distance_to_segment(const Eigen::Vector3d& p, const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    const Eigen::Vector3d along = b - a;
    const double length_squared = along.squaredNorm();
    if (length_squared < tol::kDegenerateDerivative * tol::kDegenerateDerivative) {
        return (p - a).norm();
    }
    const double s = std::clamp((p - a).dot(along) / length_squared, 0.0, 1.0);
    return (p - (a + s * along)).norm();
}

/// The distinct knots of the curve inside its domain, endpoints included.
/// These are the parameters a sampler may never step over.
std::vector<double> distinct_knots(const NurbsCurve2& curve) {
    const KnotVector& kv = curve.knots();
    const double a = kv.domain_start();
    const double b = kv.domain_end();

    std::vector<double> values{a};
    for (std::size_t i = 0; i < kv.size(); ++i) {
        const double knot = kv[i];
        if (knot > a + tol::kKnot && knot < b - tol::kKnot && knot > values.back() + tol::kKnot) {
            values.push_back(knot);
        }
    }
    values.push_back(b);
    return values;
}

/// Recursively bisects [a, b] until both tolerances hold, appending the
/// interior parameters and `b` to `out`. `a` is assumed already present.
void subdivide(const NurbsCurve2& curve,
               const NurbsSurface& surface,
               const SamplingOptions& options,
               double a,
               double b,
               const Eigen::Vector3d& pa,
               const Eigen::Vector3d& pb,
               int depth,
               std::vector<double>& out) {
    const double mid = 0.5 * (a + b);
    const Eigen::Vector3d pm = on_surface(curve, surface, mid);

    const double chord = (pb - pa).norm();
    const double sagitta = distance_to_segment(pm, pa, pb);

    const bool too_long = chord > options.max_segment_length;
    const bool too_bent = sagitta > options.max_sagitta;

    if ((too_long || too_bent) && depth < options.max_subdivisions) {
        subdivide(curve, surface, options, a, mid, pa, pm, depth + 1, out);
        subdivide(curve, surface, options, mid, b, pm, pb, depth + 1, out);
        return;
    }

    out.push_back(b);
}

std::vector<double> uniform_parameters(const NurbsCurve2& curve, const SamplingOptions& options) {
    const std::vector<double> knots = distinct_knots(curve);
    const int per_span = std::max(1, options.samples_per_span);

    std::vector<double> out{knots.front()};
    for (std::size_t i = 0; i + 1 < knots.size(); ++i) {
        const double a = knots[i];
        const double b = knots[i + 1];
        for (int k = 1; k <= per_span; ++k) {
            out.push_back(a + (b - a) * static_cast<double>(k) / static_cast<double>(per_span));
        }
    }
    return out;
}

} // namespace

std::vector<double> sample_parameters(const NurbsCurve2& curve,
                                      const NurbsSurface& surface,
                                      const SamplingOptions& options) {
    if (options.max_segment_length <= 0.0 || options.max_sagitta <= 0.0) {
        throw std::invalid_argument(
            fmt::format("sampling tolerances must be positive, got max_segment_length = {} and "
                        "max_sagitta = {}",
                        options.max_segment_length,
                        options.max_sagitta));
    }
    if (options.max_subdivisions < 0) {
        throw std::invalid_argument(
            fmt::format("max_subdivisions must not be negative, got {}", options.max_subdivisions));
    }

    if (options.mode == SamplingMode::Uniform) {
        return uniform_parameters(curve, options);
    }

    // Start from the knots and refine each span independently. Refining spans
    // rather than the whole domain is what guarantees the knots survive into
    // the output.
    const std::vector<double> knots = distinct_knots(curve);

    std::vector<double> out{knots.front()};
    for (std::size_t i = 0; i + 1 < knots.size(); ++i) {
        const double a = knots[i];
        const double b = knots[i + 1];
        subdivide(curve,
                  surface,
                  options,
                  a,
                  b,
                  on_surface(curve, surface, a),
                  on_surface(curve, surface, b),
                  0,
                  out);
    }
    return out;
}

std::vector<Eigen::Vector2d> sample_curve(const NurbsCurve2& curve,
                                          const NurbsSurface& surface,
                                          const SamplingOptions& options) {
    const std::vector<double> parameters = sample_parameters(curve, surface, options);

    std::vector<Eigen::Vector2d> points;
    points.reserve(parameters.size());
    for (const double t : parameters) {
        points.push_back(curve.evaluate(t));
    }
    return points;
}

std::vector<Eigen::Vector2d>
sample_loop(const TrimLoop& loop, const NurbsSurface& surface, const SamplingOptions& options) {
    std::vector<Eigen::Vector2d> points;

    for (const NurbsCurve2& curve : loop.curves()) {
        std::vector<Eigen::Vector2d> one = sample_curve(curve, surface, options);
        // Drop the closing point of each curve: it is the next curve's opening
        // point, and a duplicated vertex is a zero-length constraint edge that
        // the triangulator will refuse.
        one.pop_back();
        points.insert(points.end(), one.begin(), one.end());
    }

    return points;
}

double sampled_loop_length(const TrimLoop& loop,
                           const NurbsSurface& surface,
                           const SamplingOptions& options) {
    const std::vector<Eigen::Vector2d> points = sample_loop(loop, surface, options);

    double length = 0.0;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const Eigen::Vector2d& a = points[i];
        const Eigen::Vector2d& b = points[(i + 1) % points.size()];
        length += (surface.evaluate(b.x(), b.y()) - surface.evaluate(a.x(), a.y())).norm();
    }
    return length;
}

} // namespace n2s
