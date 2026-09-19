#pragma once

#include "n2s/nurbs/surface.hpp"
#include "n2s/subd/control_mesh.hpp"

namespace n2s::fit {

/// A quad layout placed by sampling the surface on a regular grid of its
/// parametric domain: control vertex `(i, j)` sits at `S(i/(rows-1),
/// j/(columns-1))`.
///
/// This is the *naive* layout, and deliberately so. Its limit surface does not
/// interpolate those points -- Catmull-Clark pulls the limit surface inside its
/// control net -- so the approximation error it produces is dominated by that
/// shrinkage rather than by anything about the trim or the layout. It exists to
/// give the harness a baseline to measure before R1 supplies a fitted layout,
/// and to make the improvement R1 brings quantifiable rather than asserted.
///
/// Throws `std::invalid_argument` for a grid smaller than 2 x 2.
ControlMesh grid_layout(const NurbsSurface& surface, int rows, int columns);

} // namespace n2s::fit
