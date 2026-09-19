#include "n2s/build_info.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <string>
#include <string_view>

namespace {

bool is_full_sha(std::string_view s) {
    return s.size() == 40 && std::all_of(s.begin(), s.end(), [](char c) {
               return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
           });
}

} // namespace

TEST_CASE("build provenance is populated at configure time", "[build_info]") {
    const n2s::BuildInfo& info = n2s::build_info();

    REQUIRE_FALSE(info.version.empty());
    REQUIRE_FALSE(info.git_sha.empty());
    REQUIRE_FALSE(info.git_describe.empty());
    REQUIRE_FALSE(info.compiler.empty());

    // Either a real commit, or one of the two documented placeholders. Anything
    // else means GitVersion.cmake produced a value that would silently mislabel
    // an experiment result.
    const bool sha_is_valid =
        is_full_sha(info.git_sha) || info.git_sha == "unknown" || info.git_sha == "uncommitted";
    CHECK(sha_is_valid);

    const bool dirty_is_valid =
        info.git_dirty == "clean" || info.git_dirty == "dirty" || info.git_dirty == "unknown";
    CHECK(dirty_is_valid);
}

TEST_CASE("the version banner names the version and the commit", "[build_info]") {
    const std::string banner{n2s::build_info_string()};
    const n2s::BuildInfo& info = n2s::build_info();

    CHECK(banner.find(std::string{info.version}) != std::string::npos);
    CHECK(banner.find(std::string{info.git_describe}) != std::string::npos);
}
