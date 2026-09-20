#pragma once

#include "n2s/fit/domain_layout.hpp"

namespace n2s::fit {

/// Splits every quad of a layout into a `factor x factor` grid of quads,
/// placing the new vertices bilinearly within each original quad.
///
/// `factor == 3` is the *edge thirding* the thesis used, and the baseline
/// refinement R1 measures against.
///
/// Why refinement is needed at all: a hand-authored layout is mostly boundary.
/// The thesis DoubleVB layout has 26 vertices of which 22 lie on the boundary,
/// so a fit over it has four free interior control points and its error says
/// far more about the boundary treatment than about the fitting method.
/// Thirding once turns those 14 quads into 126 and gives the interior real
/// freedom.
///
/// Vertices on a shared edge are created once and referenced by both faces.
/// That matters beyond tidiness: a duplicated vertex along a shared edge is a
/// crack in the control mesh, and Catmull-Clark would treat the two sides as
/// unrelated boundaries. Because the new positions on an edge depend only on
/// that edge's two endpoints, both faces compute them identically, so the
/// sharing is exact rather than within a tolerance.
///
/// Throws `std::invalid_argument` for a factor below 1 or an invalid layout.
DomainLayout refine_quads(const DomainLayout& layout, int factor);

/// `refine_quads(layout, 3)`, named for what the thesis called it.
inline DomainLayout third_edges(const DomainLayout& layout) {
    return refine_quads(layout, 3);
}

} // namespace n2s::fit
