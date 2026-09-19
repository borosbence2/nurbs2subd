#include "n2s/nurbs/knot_insertion.hpp"

#include "n2s/tolerances.hpp"

#include <fmt/format.h>

#include <stdexcept>
#include <vector>

namespace n2s {

namespace detail {

InsertionPlan plan_insertion(const KnotVector& kv, double value, int times) {
    if (times < 1) {
        throw std::invalid_argument(
            fmt::format("knot insertion count must be at least 1, got {}", times));
    }
    if (value < kv.domain_start() - tol::kKnot || value > kv.domain_end() + tol::kKnot) {
        throw std::invalid_argument(fmt::format("cannot insert knot {} outside the domain [{}, {}]",
                                                value,
                                                kv.domain_start(),
                                                kv.domain_end()));
    }

    const int multiplicity = kv.multiplicity(value);
    const int degree = kv.degree();

    // Multiplicity above the degree makes the curve discontinuous there, which
    // KnotVector rejects on construction anyway; failing here gives the caller
    // a message that names the real problem.
    if (multiplicity + times > degree) {
        throw std::invalid_argument(fmt::format(
            "inserting knot {} {} times would raise its multiplicity to {}, above the degree {}",
            value,
            times,
            multiplicity + times,
            degree));
    }

    const std::size_t span = kv.find_span(value);

    std::vector<double> knots;
    knots.reserve(kv.size() + static_cast<std::size_t>(times));
    for (std::size_t i = 0; i <= span; ++i) {
        knots.push_back(kv[i]);
    }
    for (int i = 0; i < times; ++i) {
        knots.push_back(value);
    }
    for (std::size_t i = span + 1; i < kv.size(); ++i) {
        knots.push_back(kv[i]);
    }

    return InsertionPlan{span, multiplicity, times, std::move(knots)};
}

} // namespace detail

NurbsSurface
insert_knot(const NurbsSurface& surface, Direction direction, double value, int times) {
    using Homogeneous = Eigen::Vector4d;

    const KnotVector& target = direction == Direction::U ? surface.knots_u() : surface.knots_v();
    const detail::InsertionPlan plan = detail::plan_insertion(target, value, times);

    const std::size_t num_u = surface.num_u();
    const std::size_t num_v = surface.num_v();
    const std::size_t new_count =
        (direction == Direction::U ? num_u : num_v) + static_cast<std::size_t>(plan.times);

    const std::size_t rows = direction == Direction::U ? num_v : num_u;
    const std::size_t new_num_u = direction == Direction::U ? new_count : num_u;
    const std::size_t new_num_v = direction == Direction::U ? num_v : new_count;

    std::vector<Eigen::Vector3d> points(new_num_u * new_num_v, Eigen::Vector3d::Zero());
    std::vector<double> weights(new_num_u * new_num_v, 1.0);

    // Insert into every line of the control net that runs along `direction`:
    // a tensor product surface is a curve in each direction, so the curve
    // algorithm applied line by line is exactly the surface algorithm.
    for (std::size_t line = 0; line < rows; ++line) {
        std::vector<Homogeneous> homogeneous;
        homogeneous.reserve(direction == Direction::U ? num_u : num_v);

        if (direction == Direction::U) {
            for (std::size_t i = 0; i < num_u; ++i) {
                homogeneous.push_back(detail::to_homogeneous<3>(surface.control_point(i, line),
                                                                surface.weight(i, line)));
            }
        } else {
            for (std::size_t j = 0; j < num_v; ++j) {
                homogeneous.push_back(detail::to_homogeneous<3>(surface.control_point(line, j),
                                                                surface.weight(line, j)));
            }
        }

        const std::vector<Homogeneous> inserted =
            detail::insert_into<3>(target, homogeneous, value, plan);

        for (std::size_t n = 0; n < inserted.size(); ++n) {
            const Homogeneous& h = inserted[n];
            const std::size_t flat =
                direction == Direction::U ? n * new_num_v + line : line * new_num_v + n;
            points[flat] = h.head<3>() / h.w();
            weights[flat] = h.w();
        }
    }

    KnotVector knots_u =
        direction == Direction::U ? KnotVector{surface.degree_u(), plan.knots} : surface.knots_u();
    KnotVector knots_v =
        direction == Direction::V ? KnotVector{surface.degree_v(), plan.knots} : surface.knots_v();

    return NurbsSurface{
        std::move(knots_u), std::move(knots_v), std::move(points), std::move(weights)};
}

} // namespace n2s
