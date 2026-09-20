#pragma once

#include <Eigen/Core>

#include <array>
#include <cstddef>
#include <vector>

namespace n2s::fit {

/// A quad layout authored in the parametric domain of a surface.
///
/// This, rather than a mesh of 3D points, is what a layout actually is. Two
/// things follow from keeping it in the domain and only lifting it when
/// needed:
///
///  - the domain correspondence is known by construction, so the parametric,
///    normal and curvature metrics are available rather than absent (a control
///    mesh arriving as bare 3D points carries no correspondence, and there is
///    no way to recover one);
///  - the layout can be checked against the trim region it is supposed to
///    cover, which is a domain-space question.
///
/// Kept in a header of its own, free of the surface and metrics includes, so
/// that the case file format can carry a layout without pulling in the whole
/// fitting and measurement stack.
struct DomainLayout {
    std::vector<Eigen::Vector2d> vertices;
    std::vector<std::array<int, 4>> quads;

    std::size_t num_vertices() const { return vertices.size(); }

    std::size_t num_quads() const { return quads.size(); }
};

} // namespace n2s::fit
