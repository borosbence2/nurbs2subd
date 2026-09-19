#include "n2s/nurbs/differential.hpp"

#include "n2s/tolerances.hpp"

#include <Eigen/Geometry>
#include <fmt/format.h>

#include <cmath>
#include <stdexcept>

namespace n2s {

namespace {

void require_second_order(const SurfaceDerivatives& d) {
    if (d.max_order() < 2) {
        throw std::invalid_argument(
            fmt::format("curvature needs derivatives up to order 2, got {}", d.max_order()));
    }
}

} // namespace

std::optional<Eigen::Vector3d> surface_normal(const SurfaceDerivatives& d) {
    const Eigen::Vector3d cross = d.du().cross(d.dv());
    const double length = cross.norm();

    // Catches both a vanishing partial and two parallel partials; either way
    // there is no tangent plane here.
    if (length < tol::kDegenerateDerivative) {
        return std::nullopt;
    }
    return cross / length;
}

FirstFundamentalForm first_fundamental_form(const SurfaceDerivatives& d) {
    return FirstFundamentalForm{
        .e = d.du().squaredNorm(),
        .f = d.du().dot(d.dv()),
        .g = d.dv().squaredNorm(),
    };
}

std::optional<SecondFundamentalForm> second_fundamental_form(const SurfaceDerivatives& d) {
    require_second_order(d);

    const std::optional<Eigen::Vector3d> normal = surface_normal(d);
    if (!normal.has_value()) {
        return std::nullopt;
    }

    return SecondFundamentalForm{
        .l = normal->dot(d.duu()),
        .m = normal->dot(d.duv()),
        .n = normal->dot(d.dvv()),
    };
}

std::optional<SurfaceCurvature> surface_curvature(const SurfaceDerivatives& d) {
    require_second_order(d);

    const std::optional<SecondFundamentalForm> second = second_fundamental_form(d);
    if (!second.has_value()) {
        return std::nullopt;
    }

    const FirstFundamentalForm first = first_fundamental_form(d);
    const double det = first.determinant();
    if (std::abs(det) < tol::kDegenerateDerivative) {
        return std::nullopt;
    }

    const double gaussian = (second->l * second->n - second->m * second->m) / det;
    const double mean =
        (first.e * second->n - 2.0 * first.f * second->m + first.g * second->l) / (2.0 * det);

    // k = H +- sqrt(H^2 - K). The discriminant is non-negative in exact
    // arithmetic (it is the squared half-difference of the principal
    // curvatures) but rounds below zero at umbilic points, where the two
    // curvatures coincide; clamping there is the correct answer, not a fudge.
    const double discriminant = mean * mean - gaussian;
    const double root = discriminant > 0.0 ? std::sqrt(discriminant) : 0.0;

    return SurfaceCurvature{
        .mean = mean,
        .gaussian = gaussian,
        .k1 = mean + root,
        .k2 = mean - root,
    };
}

} // namespace n2s
