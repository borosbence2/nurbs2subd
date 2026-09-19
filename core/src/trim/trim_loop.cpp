#include "n2s/trim/trim_loop.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace n2s {

TrimLoop::TrimLoop(std::vector<NurbsCurve2> curves)
    : curves_(std::move(curves)) {
    if (curves_.empty()) {
        throw std::invalid_argument("a trim loop needs at least one curve");
    }
}

double TrimLoop::worst_closure_gap() const {
    double worst = 0.0;
    for (std::size_t i = 0; i < curves_.size(); ++i) {
        const NurbsCurve2& current = curves_[i];
        const NurbsCurve2& next = curves_[(i + 1) % curves_.size()];

        const Eigen::Vector2d end = current.evaluate(current.knots().domain_end());
        const Eigen::Vector2d start = next.evaluate(next.knots().domain_start());
        worst = std::max(worst, (end - start).norm());
    }
    return worst;
}

std::vector<Eigen::Vector2d> TrimLoop::polygonise(int samples_per_curve) const {
    const int steps = std::max(1, samples_per_curve);

    std::vector<Eigen::Vector2d> points;
    points.reserve(curves_.size() * static_cast<std::size_t>(steps) + 1);

    for (const NurbsCurve2& curve : curves_) {
        const double a = curve.knots().domain_start();
        const double b = curve.knots().domain_end();
        // Each curve contributes its start but not its end: the next curve's
        // start is the same point, and duplicating it would create zero-length
        // segments that upset the intersection tests.
        for (int i = 0; i < steps; ++i) {
            const double t = a + (b - a) * static_cast<double>(i) / static_cast<double>(steps);
            points.push_back(curve.evaluate(t));
        }
    }

    return points;
}

double polygon_signed_area(const std::vector<Eigen::Vector2d>& polygon) {
    if (polygon.size() < 3) {
        return 0.0;
    }

    // Shoelace formula. Summing the cross products about the first vertex
    // rather than about the origin keeps the terms the size of the polygon
    // instead of the size of its offset from the origin, which matters for a
    // small trim loop far from (0, 0).
    const Eigen::Vector2d& origin = polygon.front();
    double twice_area = 0.0;
    for (std::size_t i = 1; i + 1 < polygon.size(); ++i) {
        const Eigen::Vector2d a = polygon[i] - origin;
        const Eigen::Vector2d b = polygon[i + 1] - origin;
        twice_area += a.x() * b.y() - a.y() * b.x();
    }
    return 0.5 * twice_area;
}

double TrimLoop::signed_area(int samples_per_curve) const {
    return polygon_signed_area(polygonise(samples_per_curve));
}

Orientation TrimLoop::orientation(int samples_per_curve) const {
    return signed_area(samples_per_curve) >= 0.0 ? Orientation::CounterClockwise
                                                 : Orientation::Clockwise;
}

void TrimLoop::reverse() {
    // Both halves are needed: reversing the order of the curves alone would
    // leave each individual curve still running head-to-tail the old way, so
    // the joins would no longer meet.
    std::reverse(curves_.begin(), curves_.end());
    for (NurbsCurve2& curve : curves_) {
        curve = reversed(curve);
    }
}

bool point_in_polygon(const std::vector<Eigen::Vector2d>& polygon, const Eigen::Vector2d& point) {
    if (polygon.size() < 3) {
        return false;
    }

    // Even-odd crossing count of a ray cast in +x. The `(yi > y) != (yj > y)`
    // form counts a vertex exactly once even when the ray passes through it,
    // which is what keeps the test consistent on shared boundaries.
    bool inside = false;
    for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
        const Eigen::Vector2d& a = polygon[i];
        const Eigen::Vector2d& b = polygon[j];

        const bool straddles = (a.y() > point.y()) != (b.y() > point.y());
        if (!straddles) {
            continue;
        }

        const double crossing_x = a.x() + (point.y() - a.y()) / (b.y() - a.y()) * (b.x() - a.x());
        if (point.x() < crossing_x) {
            inside = !inside;
        }
    }
    return inside;
}

TrimRegion::TrimRegion(TrimLoop outer, std::vector<TrimLoop> holes)
    : outer_(std::move(outer)),
      holes_(std::move(holes)) {}

bool TrimRegion::contains(const Eigen::Vector2d& point, int samples_per_curve) const {
    if (!point_in_polygon(outer_.polygonise(samples_per_curve), point)) {
        return false;
    }
    for (const TrimLoop& hole : holes_) {
        if (point_in_polygon(hole.polygonise(samples_per_curve), point)) {
            return false;
        }
    }
    return true;
}

double TrimRegion::area(int samples_per_curve) const {
    // Magnitudes, so the answer does not depend on the loops having been
    // through the orientation repair yet.
    double total = std::abs(outer_.signed_area(samples_per_curve));
    for (const TrimLoop& hole : holes_) {
        total -= std::abs(hole.signed_area(samples_per_curve));
    }
    return total;
}

} // namespace n2s
