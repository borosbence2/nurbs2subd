#pragma once

#include "n2s/fit/domain_layout.hpp"
#include "n2s/nurbs/surface.hpp"
#include "n2s/subd/control_mesh.hpp"
#include "n2s/trim/trim_loop.hpp"

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <array>
#include <filesystem>
#include <map>
#include <optional>
#include <string>

namespace n2s::io {

/// A camera, stored in the case file so that a figure can be regenerated with
/// the same viewpoint months later. A screenshot whose camera lives only in
/// whoever took it is not reproducible, and the plan requires that every figure
/// be regenerable from committed inputs.
struct CameraPose {
    Eigen::Vector3d eye{3.0, 3.0, 3.0};
    Eigen::Vector3d target{0.0, 0.0, 0.0};
    Eigen::Vector3d up{0.0, 0.0, 1.0};
};

/// Lower and upper bound of a scalar colour map.
///
/// CLAUDE.md forbids auto-ranging in any figure meant for comparison, and
/// defect 6 of the thesis was curvature plots with no fixed scale. Ranges
/// therefore live in the case file rather than being chosen by whatever the
/// viewer happened to compute, so two cases plotted side by side are actually
/// comparable.
struct ScalarRange {
    double low = 0.0;
    double high = 1.0;

    bool contains(double value) const { return value >= low && value <= high; }
};

/// Everything needed to run and to draw one test case.
struct Case {
    std::string name;
    NurbsSurface surface;
    TrimRegion region;

    /// The quad layout, authored in the parametric domain.
    ///
    /// Preferred over `control_mesh` whenever both could be given: a layout
    /// carries its own domain correspondence, so the parametric, normal and
    /// curvature metrics are available for it. A bare control mesh does not,
    /// and no correspondence can be recovered from one after the fact.
    std::optional<fit::DomainLayout> layout;

    /// The Catmull-Clark control mesh, for a case that supplies control points
    /// directly rather than a domain layout. Absent for cases that only
    /// exercise the trimming pipeline.
    std::optional<ControlMesh> control_mesh;

    std::optional<CameraPose> camera;

    /// Colour ranges by scalar name, e.g. "error", "mean_curvature",
    /// "gaussian_curvature". A viewer must refuse to draw a scalar with no
    /// range rather than invent one.
    std::map<std::string, ScalarRange> color_ranges;
};

/// Document tag and schema version for case files.
inline constexpr int kCaseJsonVersion = 1;

nlohmann::json to_json(const Case& test_case);

/// Throws `std::runtime_error` naming the offending field if the document is
/// malformed. Does *not* validate or repair the trim region: that is
/// `validate_and_repair`'s job, and doing it here would hide from the caller
/// that the file on disk needed fixing.
Case case_from_json(const nlohmann::json& document);

void write_case(const std::filesystem::path& path, const Case& test_case);
Case read_case(const std::filesystem::path& path);

} // namespace n2s::io
