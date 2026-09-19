#pragma once

#include "n2s/nurbs/curve.hpp"
#include "n2s/nurbs/surface.hpp"

#include <nlohmann/json.hpp>

#include <filesystem>

/// JSON serialisation for NURBS geometry.
///
/// The document shapes are deliberately explicit rather than compact: every
/// redundant field (`num_u`, `num_v`) is validated against the others on read,
/// because a silently transposed control net produces a plausible-looking
/// surface and a nonsense error metric three milestones later.
///
/// Curve document:
/// ```json
/// {
///   "format": "n2s-curve", "version": 1,
///   "degree": 3,
///   "knots": [0, 0, 0, 0, 0.5, 1, 1, 1, 1],
///   "control_points": [[0, 0, 0], ...],
///   "weights": [1, 1, ...]
/// }
/// ```
///
/// Surface document:
/// ```json
/// {
///   "format": "n2s-surface", "version": 1,
///   "degree_u": 3, "degree_v": 3,
///   "knots_u": [...], "knots_v": [...],
///   "num_u": 4, "num_v": 4,
///   "control_points": [[x, y, z], ...],
///   "weights": [...]
/// }
/// ```
///
/// `control_points` is row-major with u as the slow index, matching
/// `NurbsSurface`. `weights` is optional; omitting it means a non-rational
/// curve or surface.
namespace n2s::io {

/// Schema version written by this build. Readers reject anything else rather
/// than guessing at a format they do not know.
inline constexpr int kNurbsJsonVersion = 1;

nlohmann::json to_json(const NurbsCurve& curve);

/// Domain curves carry the tag "n2s-curve2" and two-component control points.
/// A separate tag rather than a dimension field, so that handing a trim curve
/// to a reader expecting model geometry fails immediately instead of producing
/// a curve in the wrong space.
nlohmann::json to_json(const NurbsCurve2& curve);

nlohmann::json to_json(const NurbsSurface& surface);

/// Throws `std::runtime_error` naming the offending field if the document is
/// malformed, and `std::invalid_argument` from the geometry constructors if the
/// document is well-formed but describes an invalid curve or surface.
NurbsCurve curve_from_json(const nlohmann::json& document);
NurbsCurve2 curve2_from_json(const nlohmann::json& document);
NurbsSurface surface_from_json(const nlohmann::json& document);

/// File-level convenience. Writes pretty-printed JSON, since these files are
/// read and hand-edited as often as they are parsed.
void write_curve(const std::filesystem::path& path, const NurbsCurve& curve);
void write_surface(const std::filesystem::path& path, const NurbsSurface& surface);

NurbsCurve read_curve(const std::filesystem::path& path);
NurbsSurface read_surface(const std::filesystem::path& path);

} // namespace n2s::io
