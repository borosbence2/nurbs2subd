#pragma once

#include <string_view>

namespace n2s {

/// Provenance of the running binary. Every experiment result records these
/// fields so that a metrics file can be traced back to the exact source state
/// that produced it.
struct BuildInfo {
    std::string_view version;      ///< Project version, e.g. "0.1.0".
    std::string_view git_sha;      ///< Full commit hash at configure time.
    std::string_view git_describe; ///< `git describe --always --dirty --tags`.
    std::string_view git_dirty;    ///< "clean", "dirty" or "unknown".
    std::string_view build_type;   ///< CMAKE_BUILD_TYPE.
    std::string_view compiler;     ///< Compiler id and version.
    std::string_view sanitizers;   ///< Active sanitizer set, or "OFF".
};

/// Returns the build provenance captured when CMake last configured the tree.
const BuildInfo& build_info();

/// One-line human-readable summary, as printed by `nurbs2subd --version`.
std::string_view build_info_string();

} // namespace n2s
