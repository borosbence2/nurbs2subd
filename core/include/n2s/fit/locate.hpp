#pragma once

#include "n2s/fit/domain_layout.hpp"
#include "n2s/subd/subdivision.hpp"

#include <Eigen/Core>

#include <cstdint>
#include <optional>
#include <unordered_map>
#include <vector>

namespace n2s::fit {

/// Maps a point of the parametric domain back to the `(face, u, v)` the
/// layout's bilinear map sends there: the inverse of `bilinear_domain_map`.
///
/// This is what a least-squares fit needs and an interpolating fit does not.
/// R1 had one target per layout vertex, and a vertex already knows which faces
/// it belongs to. R2 samples the domain densely instead -- far more samples
/// than control points -- and every one of those samples has to be turned into
/// a location at which the limit surface can be evaluated before it can appear
/// in a row of the system.
///
/// **The choice the plan asks to be documented.** The inverse is taken of the
/// layout's *bilinear* map, the same one `bilinear_domain_map` defines, and not
/// of some other parameterisation of the quad. Three reasons:
///
///  - it is the map the metrics already use, so a sample's location means the
///    same thing to the fit and to the measurement that judges it. Fitting
///    against one correspondence and measuring against another would produce a
///    parametric error that is partly a disagreement between two conventions,
///    and nothing in the numbers would separate the two;
///  - it agrees with its neighbours along shared edges, so a sample on an edge
///    gets the same domain point from either side and the fit sees no seam
///    where the layout has none;
///  - it is invertible in closed form (one quadratic), so locating a sample
///    needs no iteration and cannot half-converge.
///
/// What it inherits from the forward map is the same caveat: bilinear is not
/// area-preserving, so a strongly non-rectangular quad is sampled unevenly in
/// the domain. A uniform grid of domain samples therefore does *not* produce a
/// uniform density of samples in `(u, v)`, and a fit weighting every sample
/// equally weights those regions unequally. R2 measures the effect rather than
/// assuming it away.
class LayoutLocator {
public:
    /// Builds the index. Throws `std::invalid_argument` for a layout with no
    /// quads or with out-of-range vertex indices.
    explicit LayoutLocator(const DomainLayout& layout);

    /// The location whose forward map is `point`, or `nullopt` when the point
    /// lies outside every quad.
    ///
    /// A point on a shared edge belongs to both neighbouring quads and both
    /// answers are correct. The one returned is chosen deterministically --
    /// the quad it sits furthest inside, ties broken by the lower face index --
    /// so that a fit built from the same samples twice produces the same
    /// system both times.
    std::optional<LimitLocation> locate(const Eigen::Vector2d& point) const;

    const DomainLayout& layout() const { return layout_; }

private:
    /// Quads whose padded bounding box overlaps the cell a point falls in.
    const std::vector<int>& candidates(const Eigen::Vector2d& point) const;

    std::int64_t cell_key(const Eigen::Vector2d& point) const;

    DomainLayout layout_;

    /// Uniform bucket grid over the domain, one bucket per cell, holding the
    /// quads whose bounding box reaches it. Without it, locating `n` samples in
    /// a layout of `f` quads is `n * f`, and R2's sweeps push both up at once.
    std::unordered_map<std::int64_t, std::vector<int>> buckets_;
    std::vector<int> empty_;

    Eigen::Vector2d origin_{0.0, 0.0};
    double cell_size_ = 1.0;
};

/// Inverts the bilinear map of a single quad, given its four corners in the
/// order `bilinear_domain_map` uses: `(0,0)`, `(1,0)`, `(1,1)`, `(0,1)`.
///
/// Returns the `(u, v)` solving the bilinear system, whether or not it lands
/// inside `[0,1]^2`; `nullopt` only when the quad is too degenerate to invert.
/// Exposed because it is the part worth testing against an analytic oracle on
/// its own, separately from the search that finds the quad.
std::optional<Eigen::Vector2d> invert_bilinear(const Eigen::Vector2d& corner00,
                                               const Eigen::Vector2d& corner10,
                                               const Eigen::Vector2d& corner11,
                                               const Eigen::Vector2d& corner01,
                                               const Eigen::Vector2d& point);

} // namespace n2s::fit
