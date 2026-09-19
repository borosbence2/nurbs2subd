#pragma once

#include "n2s/nurbs/surface.hpp"

#include <Eigen/Core>

#include <optional>

namespace n2s {

/// Coefficients of the first fundamental form:
/// `E = Su.Su`, `F = Su.Sv`, `G = Sv.Sv`.
struct FirstFundamentalForm {
    double e;
    double f;
    double g;

    /// `EG - F^2`, the squared area of the parallelogram spanned by the
    /// partials. Zero exactly where the surface is degenerate.
    double determinant() const { return e * g - f * f; }
};

/// Coefficients of the second fundamental form:
/// `L = n.Suu`, `M = n.Suv`, `N = n.Svv`, with `n` the unit normal below.
struct SecondFundamentalForm {
    double l;
    double m;
    double n;
};

/// Curvature of a surface at a point.
///
/// The sign of `mean`, `k1` and `k2` follows the orientation of the unit
/// normal, which is fixed by the parametric direction of the surface: swapping
/// u and v flips all three. `gaussian` is orientation-independent. Code that
/// compares curvature across surfaces with unrelated parameterisations should
/// compare magnitudes, or fix the orientation first.
/// Accuracy note: `k1` and `k2` come from `H +- sqrt(H^2 - K)`, and near an
/// umbilic point (`k1 == k2`, e.g. anywhere on a sphere) that discriminant
/// cancels down to rounding noise. The square root then costs half the
/// available digits, so the principal curvatures there carry a relative error
/// around `sqrt(machine epsilon)`, roughly 1e-8, even though `mean` and
/// `gaussian` are accurate to near machine precision. Prefer `mean` and
/// `gaussian` for error metrics and figures.
struct SurfaceCurvature {
    double mean;     ///< H = (EN - 2FM + GL) / (2 (EG - F^2))
    double gaussian; ///< K = (LN - M^2) / (EG - F^2)
    double k1;       ///< Larger principal curvature, H + sqrt(H^2 - K).
    double k2;       ///< Smaller principal curvature, H - sqrt(H^2 - K).
};

/// Unit normal `(Su x Sv) / |Su x Sv|`.
///
/// Returns `nullopt` where the partials vanish or are parallel, which is where
/// the tangent plane genuinely does not exist: the pole of a sphere of
/// revolution, a collapsed patch edge, a degenerate control net. Callers have
/// to decide what that means for them, so it is not reported as a zero vector.
std::optional<Eigen::Vector3d> surface_normal(const SurfaceDerivatives& d);

/// Always defined: the first fundamental form needs only the first partials.
FirstFundamentalForm first_fundamental_form(const SurfaceDerivatives& d);

/// Needs the unit normal, so it inherits its degenerate case.
/// `d` must have been computed with `max_order >= 2`.
std::optional<SecondFundamentalForm> second_fundamental_form(const SurfaceDerivatives& d);

/// Mean, Gaussian and principal curvatures.
/// `d` must have been computed with `max_order >= 2`.
std::optional<SurfaceCurvature> surface_curvature(const SurfaceDerivatives& d);

} // namespace n2s
