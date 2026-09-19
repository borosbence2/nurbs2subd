#pragma once

#include "n2s/nurbs/surface.hpp"
#include "n2s/trim/trim_loop.hpp"

#include <Eigen/Core>

#include <string>
#include <vector>

/// Synthetic test cases: a base surface plus a trimmed domain, built in code so
/// that every experiment can be reproduced without a data file.
///
/// The surfaces are *exact*. The paraboloid and saddle are quadratic in each
/// direction and are therefore represented exactly by a biquadratic Bezier
/// patch, and the cylinder and sphere exactly by their rational forms. That
/// matters: an approximate base surface would put a floor under every
/// approximation error this project sets out to measure.
namespace n2s::cases {

// ---------------------------------------------------------------------------
// Base surfaces. All are parameterised over (u, v) in [0, 1]^2.
// ---------------------------------------------------------------------------

/// Planar rectangle in the z = 0 plane, spanning [0, width] x [0, height].
NurbsSurface plane(double width = 1.0, double height = 1.0);

/// Elliptic paraboloid `z = curvature * (x^2 + y^2)` over
/// [-half_width, half_width]^2. Gaussian curvature is positive everywhere.
NurbsSurface paraboloid(double half_width = 1.0, double curvature = 0.5);

/// Hyperbolic paraboloid `z = curvature * (x^2 - y^2)` over the same square.
/// Gaussian curvature is negative everywhere, which is the interesting case for
/// subdivision fitting.
NurbsSurface saddle(double half_width = 1.0, double curvature = 0.5);

/// Circular cylinder, axis along z, base at z = 0. Rational quadratic around,
/// linear along.
NurbsSurface cylinder(double radius = 1.0, double height = 1.0);

/// Sphere centred at the origin, as a surface of revolution. Note that v = 0
/// and v = 1 are degenerate poles.
NurbsSurface sphere(double radius = 1.0);

// ---------------------------------------------------------------------------
// Trim loops, in the (u, v) domain.
//
// Each is built with the orientation the convention wants, but nothing depends
// on that: `validate_and_repair` fixes a reversed loop, and the tests rely on
// it doing so.
// ---------------------------------------------------------------------------

/// Axis-aligned rectangle, counter-clockwise, as four degree-1 curves.
TrimLoop rectangle_loop(const Eigen::Vector2d& low, const Eigen::Vector2d& high);

/// Exact circle as four rational quadratic segments, counter-clockwise.
/// Reverse it (or let the validator do it) to use as a hole.
TrimLoop circle_loop(const Eigen::Vector2d& centre, double radius);

/// L-shape: `outer` with the corner block of size `notch` removed from its
/// high-u, high-v corner. Six degree-1 curves, counter-clockwise.
TrimLoop
l_shape_loop(const Eigen::Vector2d& low, const Eigen::Vector2d& high, const Eigen::Vector2d& notch);

/// A rectangle with one corner replaced by a concave circular arc, the domain
/// of Shen et al. (2014) figure 1. The arc is centred on the low corner with
/// the given radius, so the trimmed region is the rectangle *minus* a quarter
/// disc. Counter-clockwise.
TrimLoop
quarter_arc_corner_loop(const Eigen::Vector2d& low, const Eigen::Vector2d& high, double radius);

// ---------------------------------------------------------------------------
// Whole cases.
// ---------------------------------------------------------------------------

/// A named surface and the trimmed part of its domain.
struct TrimmedCase {
    std::string name;
    NurbsSurface surface;
    TrimRegion region;

    /// Analytic area of the trimmed domain in `(u, v)`, where one is known in
    /// closed form. Negative means "no closed form for this case". This is the
    /// oracle the triangulation area test converges against.
    double analytic_domain_area = -1.0;
};

/// The unit domain with a circular hole in the middle. Domain area is
/// `1 - pi r^2` exactly, which makes this the convergence oracle for the
/// triangulator.
TrimmedCase square_with_circular_hole(const NurbsSurface& surface, double radius = 0.25);

/// The Shen et al. figure 1 domain on the given surface.
TrimmedCase quarter_arc_corner(const NurbsSurface& surface, double radius = 0.5);

/// The L-shaped domain on the given surface.
TrimmedCase l_shape(const NurbsSurface& surface, const Eigen::Vector2d& notch = {0.5, 0.5});

/// Every synthetic case, each base surface crossed with each trim, for sweeps
/// and for the "all synthetic cases triangulate cleanly" exit criterion.
std::vector<TrimmedCase> all_synthetic_cases();

} // namespace n2s::cases
