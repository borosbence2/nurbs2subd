#pragma once

#include "n2s/fit/domain_layout.hpp"
#include "n2s/metrics/surface_error.hpp"
#include "n2s/nurbs/surface.hpp"
#include "n2s/subd/control_mesh.hpp"

#include <optional>

namespace n2s::fit {

/// Lifts a layout onto a surface: control vertex `(u, v)` becomes `S(u, v)`.
///
/// Note what this is *not*: the limit surface of the result does not pass
/// through those points, because Catmull-Clark pulls its limit surface inside
/// the control net. This is the starting point a fit improves on, not a fit.
ControlMesh lift(const DomainLayout& layout, const NurbsSurface& surface);

/// The correspondence induced by a layout: within each quad, `(u, v)` maps
/// bilinearly onto the quad's four domain corners.
///
/// Bilinear is a choice, and the plan asks for it to be documented. It is the
/// only map determined by the four corners alone, it is exact when the quad is
/// an axis-aligned rectangle in the domain (which is the case for every grid
/// layout), and it agrees with its neighbours along shared edges, so the
/// correspondence is continuous across the whole layout. What it is not is
/// area-preserving: a strongly non-rectangular quad is sampled unevenly in the
/// domain, and that shows up as parametric error even where the surfaces
/// coincide.
metrics::DomainMap bilinear_domain_map(const DomainLayout& layout);

/// A layout resolved into the two things every downstream stage needs: the
/// control mesh to subdivide, and the domain correspondence to measure
/// against.
struct ResolvedLayout {
    ControlMesh mesh;

    /// Empty when the case supplied bare control points, in which case the
    /// correspondence-dependent metrics are unavailable and must be reported
    /// as absent rather than guessed.
    metrics::DomainMap domain_map;

    /// True when the layout came from the domain and therefore carries a
    /// correspondence.
    bool has_correspondence = false;
};

/// Chooses between a domain layout and a bare control mesh.
///
/// The domain layout wins whenever both are present: it carries a
/// correspondence and the control mesh does not, and lifting the layout
/// reproduces the control mesh anyway for any case where the two agree.
/// Returns `nullopt` when neither is supplied.
std::optional<ResolvedLayout> resolve_layout(const std::optional<DomainLayout>& layout,
                                             const std::optional<ControlMesh>& mesh,
                                             const NurbsSurface& surface);

/// The layout implied by `ControlMesh::grid`, in the unit domain. Equivalent to
/// `metrics::grid_domain_map` but as an explicit layout.
DomainLayout grid_domain_layout(int rows, int columns);

/// A quad layout placed by sampling the surface on a regular grid of its
/// parametric domain: control vertex `(i, j)` sits at `S(i/(rows-1),
/// j/(columns-1))`.
///
/// This is the *naive* layout, and deliberately so. Its limit surface does not
/// interpolate those points, so the approximation error it produces is
/// dominated by that shrinkage rather than by anything about the trim or the
/// layout. It exists to give the harness a baseline to measure before R1
/// supplies a fitted layout, and to make the improvement R1 brings quantifiable
/// rather than asserted.
///
/// Throws `std::invalid_argument` for a grid smaller than 2 x 2.
ControlMesh grid_layout(const NurbsSurface& surface, int rows, int columns);

} // namespace n2s::fit
