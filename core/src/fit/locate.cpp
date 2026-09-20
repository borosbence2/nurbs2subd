#include "n2s/fit/locate.hpp"

#include "n2s/tolerances.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace n2s::fit {

namespace {

/// The 2D cross product, a signed area. The whole inversion below is written in
/// terms of it.
double cross(const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
    return a.x() * b.y() - a.y() * b.x();
}

/// How far `(u, v)` lies outside the unit square: zero inside, growing with the
/// distance out. The ordering key that decides which quad claims a point lying
/// on a shared edge.
double outsideness(const Eigen::Vector2d& uv) {
    const double du = std::max({0.0, -uv.x(), uv.x() - 1.0});
    const double dv = std::max({0.0, -uv.y(), uv.y() - 1.0});
    return std::max(du, dv);
}

void validate(const DomainLayout& layout) {
    if (layout.vertices.empty() || layout.quads.empty()) {
        throw std::invalid_argument("a domain layout needs at least one quad");
    }
    const auto count = static_cast<int>(layout.vertices.size());
    for (std::size_t f = 0; f < layout.quads.size(); ++f) {
        for (const int index : layout.quads[f]) {
            if (index < 0 || index >= count) {
                throw std::invalid_argument(fmt::format(
                    "layout quad {} references vertex {}, outside [0, {})", f, index, count));
            }
        }
    }
}

} // namespace

std::optional<Eigen::Vector2d> invert_bilinear(const Eigen::Vector2d& corner00,
                                               const Eigen::Vector2d& corner10,
                                               const Eigen::Vector2d& corner11,
                                               const Eigen::Vector2d& corner01,
                                               const Eigen::Vector2d& point) {
    // Writing the forward map as
    //
    //     P(u,v) = a + u B + v D + u v C
    //
    // with B = b - a, D = d - a and C = a - b + c - d, the problem is to solve
    // u B + v D + u v C = q for q = point - a. Grouping it as
    // u (B + v C) + v D and crossing both sides with (B + v C) kills the u
    // term, leaving a quadratic in v alone. That is the whole method: it is
    // closed form, so there is no iteration to half-converge, and the two
    // roots are exactly the two branches of a bilinear patch's inverse.
    const Eigen::Vector2d B = corner10 - corner00;
    const Eigen::Vector2d D = corner01 - corner00;
    const Eigen::Vector2d C = corner00 - corner10 + corner11 - corner01;
    const Eigen::Vector2d q = point - corner00;

    const double scale = std::max({B.norm(), D.norm(), C.norm()});
    if (scale <= 0.0) {
        return std::nullopt; // all four corners coincide
    }

    const double a2 = cross(C, D);
    const double a1 = cross(B, D) - cross(C, q);
    const double a0 = -cross(B, q);

    // Given v, recover u from u (B + v C) = q - v D. Solved as a projection
    // rather than by dividing one chosen component, so there is no component to
    // choose badly when B + v C happens to be axis-aligned.
    const auto solve_u = [&](double v) -> std::optional<double> {
        const Eigen::Vector2d w = B + v * C;
        const double ww = w.squaredNorm();
        if (ww <= tol::kRelativeBilinearDegeneracy * scale * scale) {
            return std::nullopt; // the quad pinches to a point along this row
        }
        return w.dot(q - v * D) / ww;
    };

    const auto finish = [&](double v) -> std::optional<Eigen::Vector2d> {
        const std::optional<double> u = solve_u(v);
        if (!u.has_value()) {
            return std::nullopt;
        }
        return Eigen::Vector2d{*u, v};
    };

    // A parallelogram has C = 0 exactly and the quadratic degenerates to a
    // linear equation. Relative threshold: an absolute one would classify the
    // same quad differently at two model scales.
    if (std::abs(a2) <= tol::kRelativeBilinearDegeneracy * scale * scale) {
        if (std::abs(a1) <= tol::kRelativeBilinearDegeneracy * scale * scale) {
            return std::nullopt; // degenerate: no isolated solution
        }
        return finish(-a0 / a1);
    }

    const double discriminant = a1 * a1 - 4.0 * a2 * a0;
    if (discriminant < 0.0) {
        return std::nullopt; // the point is off the bilinear patch entirely
    }

    // The numerically stable pair of roots: forming both as -a1 +/- sqrt(...)
    // over 2 a2 loses the small root to cancellation whenever a0 a2 << a1^2,
    // which is exactly the near-parallelogram case this has to handle well.
    const double root = std::sqrt(discriminant);
    const double qq = -0.5 * (a1 + std::copysign(root, a1));
    const double v1 = qq / a2;
    const double v2 = (std::abs(qq) > 0.0) ? a0 / qq : v1;

    // Both roots are genuine inverses of the bilinear map; only one of them is
    // inside the quad for a point that is inside the quad. Pick by how far
    // outside the unit square each lands, so that a point marginally outside
    // still resolves to the branch that belongs to this quad.
    std::optional<Eigen::Vector2d> best;
    double best_outside = std::numeric_limits<double>::infinity();
    for (const double v : {v1, v2}) {
        const std::optional<Eigen::Vector2d> uv = finish(v);
        if (!uv.has_value()) {
            continue;
        }
        const double out = outsideness(*uv);
        if (out < best_outside) {
            best_outside = out;
            best = uv;
        }
    }
    return best;
}

LayoutLocator::LayoutLocator(const DomainLayout& layout)
    : layout_(layout) {
    validate(layout_);

    Eigen::Vector2d low = layout_.vertices.front();
    Eigen::Vector2d high = layout_.vertices.front();
    for (const Eigen::Vector2d& v : layout_.vertices) {
        low = low.cwiseMin(v);
        high = high.cwiseMax(v);
    }
    origin_ = low;

    // One cell per quad on average: enough to make the scan short without
    // spending more time filling buckets than searching them.
    const Eigen::Vector2d extent = high - low;
    const double area = std::max(extent.x() * extent.y(), tol::kMinLoopArea);
    cell_size_ = std::sqrt(area / static_cast<double>(layout_.num_quads()));
    if (!(cell_size_ > 0.0)) {
        cell_size_ = 1.0;
    }

    // Padded by the containment slack, so a point that a quad accepts is never
    // filed in a cell that quad was not registered in. The search would
    // otherwise be tighter than the test it feeds.
    const double pad = tol::kLayoutContainment * std::max(extent.maxCoeff(), 1.0);

    for (std::size_t f = 0; f < layout_.quads.size(); ++f) {
        Eigen::Vector2d quad_low = layout_.vertices[static_cast<std::size_t>(layout_.quads[f][0])];
        Eigen::Vector2d quad_high = quad_low;
        for (const int index : layout_.quads[f]) {
            const Eigen::Vector2d& v = layout_.vertices[static_cast<std::size_t>(index)];
            quad_low = quad_low.cwiseMin(v);
            quad_high = quad_high.cwiseMax(v);
        }
        quad_low.array() -= pad;
        quad_high.array() += pad;

        const auto x0 =
            static_cast<std::int64_t>(std::floor((quad_low.x() - origin_.x()) / cell_size_));
        const auto x1 =
            static_cast<std::int64_t>(std::floor((quad_high.x() - origin_.x()) / cell_size_));
        const auto y0 =
            static_cast<std::int64_t>(std::floor((quad_low.y() - origin_.y()) / cell_size_));
        const auto y1 =
            static_cast<std::int64_t>(std::floor((quad_high.y() - origin_.y()) / cell_size_));

        for (std::int64_t x = x0; x <= x1; ++x) {
            for (std::int64_t y = y0; y <= y1; ++y) {
                buckets_[x * 73856093LL ^ y * 19349663LL].push_back(static_cast<int>(f));
            }
        }
    }
}

std::int64_t LayoutLocator::cell_key(const Eigen::Vector2d& point) const {
    const auto x = static_cast<std::int64_t>(std::floor((point.x() - origin_.x()) / cell_size_));
    const auto y = static_cast<std::int64_t>(std::floor((point.y() - origin_.y()) / cell_size_));
    return x * 73856093LL ^ y * 19349663LL;
}

const std::vector<int>& LayoutLocator::candidates(const Eigen::Vector2d& point) const {
    const auto found = buckets_.find(cell_key(point));
    return found == buckets_.end() ? empty_ : found->second;
}

std::optional<LimitLocation> LayoutLocator::locate(const Eigen::Vector2d& point) const {
    std::optional<LimitLocation> best;
    double best_outside = std::numeric_limits<double>::infinity();

    for (const int face : candidates(point)) {
        const std::array<int, 4>& quad = layout_.quads[static_cast<std::size_t>(face)];
        const std::optional<Eigen::Vector2d> uv =
            invert_bilinear(layout_.vertices[static_cast<std::size_t>(quad[0])],
                            layout_.vertices[static_cast<std::size_t>(quad[1])],
                            layout_.vertices[static_cast<std::size_t>(quad[2])],
                            layout_.vertices[static_cast<std::size_t>(quad[3])],
                            point);
        if (!uv.has_value()) {
            continue;
        }

        const double out = outsideness(*uv);
        if (out > tol::kLayoutContainment) {
            continue;
        }

        // Strictly less, and the candidates arrive in ascending face order, so
        // a tie on a shared edge always resolves to the lower face index.
        if (out < best_outside) {
            best_outside = out;
            best =
                LimitLocation{face, std::clamp(uv->x(), 0.0, 1.0), std::clamp(uv->y(), 0.0, 1.0)};
        }
    }

    return best;
}

} // namespace n2s::fit
