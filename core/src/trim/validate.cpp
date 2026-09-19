#include "n2s/trim/validate.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>

namespace n2s {

namespace {

/// True if the open segments ab and cd cross properly, i.e. at a single point
/// interior to both. Touching endpoints do not count: consecutive edges of a
/// polygon share one by construction.
bool segments_properly_cross(const Eigen::Vector2d& a,
                             const Eigen::Vector2d& b,
                             const Eigen::Vector2d& c,
                             const Eigen::Vector2d& d) {
    const auto cross = [](const Eigen::Vector2d& p, const Eigen::Vector2d& q) {
        return p.x() * q.y() - p.y() * q.x();
    };

    const double d1 = cross(b - a, c - a);
    const double d2 = cross(b - a, d - a);
    const double d3 = cross(d - c, a - c);
    const double d4 = cross(d - c, b - c);

    // Strict sign opposition on both sides. Collinear and touching cases fall
    // through as "no crossing", which is the right answer for a polygon whose
    // vertices are sampled from smooth curves: a genuine self-intersection
    // there is transversal.
    return ((d1 > 0.0) != (d2 > 0.0)) && ((d3 > 0.0) != (d4 > 0.0));
}

/// Index of the first properly crossing pair of non-adjacent edges, or nullopt.
std::optional<std::pair<std::size_t, std::size_t>>
find_self_intersection(const std::vector<Eigen::Vector2d>& polygon) {
    const std::size_t n = polygon.size();
    if (n < 4) {
        return std::nullopt;
    }

    for (std::size_t i = 0; i < n; ++i) {
        const Eigen::Vector2d& a = polygon[i];
        const Eigen::Vector2d& b = polygon[(i + 1) % n];

        for (std::size_t j = i + 2; j < n; ++j) {
            // Skip the wrap-around pair, which is adjacent to edge 0.
            if (i == 0 && j + 1 == n) {
                continue;
            }
            const Eigen::Vector2d& c = polygon[j];
            const Eigen::Vector2d& d = polygon[(j + 1) % n];

            if (segments_properly_cross(a, b, c, d)) {
                return std::make_pair(i, j);
            }
        }
    }
    return std::nullopt;
}

/// Moves the shared endpoint of two consecutive curves to their midpoint. Both
/// curves are clamped, so the endpoint *is* the first/last control point and
/// moving it moves nothing else.
void snap_join(NurbsCurve2& before, NurbsCurve2& after) {
    std::vector<Eigen::Vector2d> before_points = before.control_points();
    std::vector<Eigen::Vector2d> after_points = after.control_points();

    const Eigen::Vector2d meeting = 0.5 * (before_points.back() + after_points.front());
    before_points.back() = meeting;
    after_points.front() = meeting;

    before = NurbsCurve2{before.knots(), std::move(before_points), before.weights()};
    after = NurbsCurve2{after.knots(), std::move(after_points), after.weights()};
}

/// Checks one loop, repairing what it can. `expected` is the orientation the
/// convention demands of this loop.
void check_loop(TrimLoop& loop,
                const std::string& label,
                Orientation expected,
                const TrimValidationOptions& options,
                TrimReport& report) {
    for (std::size_t i = 0; i < loop.size(); ++i) {
        if (!loop.curves()[i].knots().is_clamped()) {
            report.errors.push_back(
                fmt::format("{}: curve {} is not clamped, so its endpoints are not its end control "
                            "points and the loop cannot be closed by snapping",
                            label,
                            i));
        }
    }
    if (!report.errors.empty()) {
        return;
    }

    // Closure gaps, join by join.
    std::vector<NurbsCurve2>& curves = loop.mutable_curves();
    for (std::size_t i = 0; i < curves.size(); ++i) {
        const std::size_t next = (i + 1) % curves.size();

        const Eigen::Vector2d end = curves[i].evaluate(curves[i].knots().domain_end());
        const Eigen::Vector2d start = curves[next].evaluate(curves[next].knots().domain_start());
        const double gap = (end - start).norm();

        if (gap <= 0.0) {
            continue;
        }
        if (gap <= tol::kNegligibleClosureGap) {
            // Rounding noise, not a defect. Close it without saying so.
            if (options.snap_closure_gaps) {
                snap_join(curves[i], curves[next]);
            }
            continue;
        }
        if (gap > options.closure_tolerance) {
            report.errors.push_back(
                fmt::format("{}: gap of {:.3e} between curve {} and curve {}, above the closure "
                            "tolerance of {:.3e}",
                            label,
                            gap,
                            i,
                            next,
                            options.closure_tolerance));
            continue;
        }
        if (options.snap_closure_gaps) {
            snap_join(curves[i], curves[next]);
            report.repairs.push_back(fmt::format(
                "{}: snapped a gap of {:.3e} between curve {} and curve {}", label, gap, i, next));
        }
    }

    // Self-intersection is checked before orientation on purpose. A loop that
    // crosses itself has no meaningful enclosed area -- a symmetric bowtie has
    // a signed area of exactly zero -- so testing orientation first would
    // report a degenerate area and bury the real defect.
    if (options.check_self_intersection) {
        const std::vector<Eigen::Vector2d> polygon = loop.polygonise(options.samples_per_curve);
        if (const auto crossing = find_self_intersection(polygon)) {
            report.errors.push_back(fmt::format(
                "{}: self-intersects, between polygonised edges {} and {}. Sampled at {} points "
                "per curve, so a near-tangency may be a sampling artefact, but a genuine "
                "crossing makes the enclosed region undefined",
                label,
                crossing->first,
                crossing->second,
                options.samples_per_curve));
            return;
        }
    }

    // Orientation.
    const double area = loop.signed_area(options.samples_per_curve);
    if (std::abs(area) < tol::kMinLoopArea) {
        report.errors.push_back(fmt::format(
            "{}: encloses an area of only {:.3e}, too small to orient reliably; the loop is "
            "degenerate",
            label,
            std::abs(area)));
        return;
    }

    const Orientation actual = area >= 0.0 ? Orientation::CounterClockwise : Orientation::Clockwise;
    if (actual != expected) {
        if (options.fix_orientation) {
            loop.reverse();
            report.repairs.push_back(fmt::format(
                "{}: reversed, it was wound {}",
                label,
                actual == Orientation::CounterClockwise ? "counter-clockwise" : "clockwise"));
        } else {
            report.errors.push_back(
                fmt::format("{}: wound the wrong way and orientation repair is disabled", label));
        }
    }
}

} // namespace

std::string TrimReport::to_string() const {
    std::string out;
    for (const std::string& repair : repairs) {
        out += fmt::format("repaired: {}\n", repair);
    }
    for (const std::string& error : errors) {
        out += fmt::format("error:    {}\n", error);
    }
    if (out.empty()) {
        out = "trim region is valid, nothing to repair\n";
    }
    return out;
}

TrimReport validate_and_repair(TrimRegion& region, const TrimValidationOptions& options) {
    TrimReport report;

    check_loop(
        region.mutable_outer(), "outer loop", Orientation::CounterClockwise, options, report);

    std::vector<TrimLoop>& holes = region.mutable_holes();
    for (std::size_t i = 0; i < holes.size(); ++i) {
        check_loop(holes[i], fmt::format("hole {}", i), Orientation::Clockwise, options, report);
    }

    if (!report.ok()) {
        // Containment between loops is meaningless once a loop is known bad.
        return report;
    }

    const std::vector<Eigen::Vector2d> outer_polygon =
        region.outer().polygonise(options.samples_per_curve);

    for (std::size_t i = 0; i < holes.size(); ++i) {
        const std::vector<Eigen::Vector2d> hole_polygon =
            holes[i].polygonise(options.samples_per_curve);

        // Every vertex, not just one: a hole that pokes out through the
        // boundary has some vertices inside and some outside, and testing a
        // single representative point would miss exactly that case.
        const bool all_inside =
            std::all_of(hole_polygon.begin(), hole_polygon.end(), [&](const Eigen::Vector2d& p) {
                return point_in_polygon(outer_polygon, p);
            });
        if (!all_inside) {
            report.errors.push_back(fmt::format("hole {} is not contained in the outer loop", i));
        }

        for (std::size_t j = i + 1; j < holes.size(); ++j) {
            const std::vector<Eigen::Vector2d> other =
                holes[j].polygonise(options.samples_per_curve);
            const bool overlaps =
                std::any_of(hole_polygon.begin(),
                            hole_polygon.end(),
                            [&](const Eigen::Vector2d& p) { return point_in_polygon(other, p); }) ||
                std::any_of(other.begin(), other.end(), [&](const Eigen::Vector2d& p) {
                    return point_in_polygon(hole_polygon, p);
                });
            if (overlaps) {
                report.errors.push_back(fmt::format("holes {} and {} overlap", i, j));
            }
        }
    }

    return report;
}

} // namespace n2s
