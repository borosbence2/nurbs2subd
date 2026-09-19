#pragma once

#include "n2s/tolerances.hpp"
#include "n2s/trim/trim_loop.hpp"

#include <string>
#include <vector>

namespace n2s {

struct TrimValidationOptions {
    /// Joins closer than this are snapped shut; anything wider is an error.
    double closure_tolerance = tol::kTrimClosure;

    /// Density of the polygonisation used for the orientation, containment and
    /// self-intersection tests.
    int samples_per_curve = kLoopTestSamplesPerCurve;

    /// Snap closure gaps that are within tolerance. Turning this off makes the
    /// pass report-only.
    bool snap_closure_gaps = true;

    /// Reverse loops whose orientation disagrees with the convention (outer
    /// counter-clockwise, holes clockwise).
    bool fix_orientation = true;

    /// Self-intersection is O(n^2) in the polygonised vertex count. It is the
    /// expensive check and the one most worth skipping on input already known
    /// to be clean.
    bool check_self_intersection = true;
};

/// What a validation pass found. Repairs are things that were silently wrong
/// and have been put right; errors are things that make the region unusable and
/// that this code will not guess its way past.
struct TrimReport {
    std::vector<std::string> repairs;
    std::vector<std::string> errors;

    bool ok() const { return errors.empty(); }

    /// Multi-line summary for logs and for the `meta.json` of an experiment.
    std::string to_string() const;
};

/// Validates `region` and repairs what it can, in place.
///
/// Checks, in order:
///  - every curve is clamped, so that its endpoints are its end control points
///    (the snap below relies on this, and unclamped trim curves do not occur in
///    CAD output);
///  - closure gaps at every join, snapped shut when within tolerance and
///    reported as an error when not;
///  - orientation, reversed to match the convention when it disagrees;
///  - self-intersection within each loop;
///  - holes lying inside the outer loop and not overlapping each other.
///
/// Every repair is reported. The defect this guards against is the silent one:
/// a hole loop wound the wrong way triangulates into a region with the hole
/// filled in and the outside hollow, which looks plausible until an error
/// metric is computed over it.
TrimReport validate_and_repair(TrimRegion& region, const TrimValidationOptions& options = {});

} // namespace n2s
