# All third-party code is fetched and built from source at pinned tags, per
# CLAUDE.md. Nothing here may be replaced by a system package without updating
# the stack section of CLAUDE.md first.

include(FetchContent)

# Several upstream projects still declare `cmake_minimum_required(VERSION <3.5)`,
# which CMake 4.x rejects outright. This makes those projects configure with 3.5
# semantics rather than failing. It only affects the fetched subprojects.
if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0)
    set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
endif()

# FetchContent gained the SYSTEM keyword in 3.25; without it dependency headers
# are compiled with our -Werror set, which no upstream project is expected to
# survive.
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.25)
    set(_n2s_fc_system SYSTEM)
else()
    set(_n2s_fc_system "")
endif()

set(FETCHCONTENT_QUIET OFF)

# Dependencies are static; the project ships as a set of executables, not a
# distributable library.
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

# --------------------------------------------------------------------------
# Eigen 3.4 - header only.
#
# Eigen's own CMakeLists pulls in test/doc/install machinery we do not want, so
# the archive is only downloaded (SOURCE_SUBDIR names a directory that holds no
# CMakeLists.txt, which makes FetchContent_MakeAvailable skip add_subdirectory)
# and wrapped in an interface target here.
# --------------------------------------------------------------------------
FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR do-not-configure-eigen)
FetchContent_MakeAvailable(eigen)

add_library(n2s_eigen INTERFACE)
target_include_directories(n2s_eigen SYSTEM INTERFACE "${eigen_SOURCE_DIR}")
add_library(Eigen3::Eigen ALIAS n2s_eigen)

# --------------------------------------------------------------------------
# OpenSubdiv - Catmull-Clark refinement and exact limit evaluation.
# Only the CPU evaluation backend is needed; every GPU/threading backend and
# all samples are disabled to keep configure and CI times sane.
# --------------------------------------------------------------------------
set(NO_EXAMPLES ON CACHE BOOL "" FORCE)
set(NO_TUTORIALS ON CACHE BOOL "" FORCE)
set(NO_REGRESSION ON CACHE BOOL "" FORCE)
set(NO_DOC ON CACHE BOOL "" FORCE)
set(NO_TESTS ON CACHE BOOL "" FORCE)
set(NO_GLTESTS ON CACHE BOOL "" FORCE)
set(NO_OMP ON CACHE BOOL "" FORCE)
set(NO_TBB ON CACHE BOOL "" FORCE)
set(NO_CUDA ON CACHE BOOL "" FORCE)
set(NO_OPENCL ON CACHE BOOL "" FORCE)
set(NO_CLEW ON CACHE BOOL "" FORCE)
set(NO_OPENGL ON CACHE BOOL "" FORCE)
set(NO_METAL ON CACHE BOOL "" FORCE)
set(NO_DX ON CACHE BOOL "" FORCE)
set(NO_PTEX ON CACHE BOOL "" FORCE)
set(NO_GLEW ON CACHE BOOL "" FORCE)
set(NO_GLFW ON CACHE BOOL "" FORCE)
set(NO_GLFW_X11 ON CACHE BOOL "" FORCE)
set(NO_MACOS_FRAMEWORK ON CACHE BOOL "" FORCE)

FetchContent_Declare(opensubdiv
    ${_n2s_fc_system}
    GIT_REPOSITORY https://github.com/PixarAnimationStudios/OpenSubdiv.git
    GIT_TAG v3_6_1
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(opensubdiv)

# The CPU library target has been spelled differently across OpenSubdiv
# releases; resolve it instead of hard-coding one name.
set(N2S_OSD_TARGET "")
foreach(_candidate osdCPU_static osdCPU osd_static_cpu osd_dynamic_cpu)
    if(TARGET ${_candidate})
        set(N2S_OSD_TARGET ${_candidate})
        break()
    endif()
endforeach()
if(N2S_OSD_TARGET STREQUAL "")
    message(FATAL_ERROR
        "Could not find the OpenSubdiv CPU library target. Inspect "
        "${opensubdiv_SOURCE_DIR}/opensubdiv/CMakeLists.txt and extend the "
        "candidate list in cmake/Dependencies.cmake.")
endif()
message(STATUS "OpenSubdiv CPU target: ${N2S_OSD_TARGET}")

# OpenSubdiv headers are included as <opensubdiv/far/...>, which lives one level
# above the library's own include path.
add_library(n2s_opensubdiv INTERFACE)
target_link_libraries(n2s_opensubdiv INTERFACE ${N2S_OSD_TARGET})
target_include_directories(n2s_opensubdiv SYSTEM INTERFACE "${opensubdiv_SOURCE_DIR}")
add_library(n2s::opensubdiv ALIAS n2s_opensubdiv)

# --------------------------------------------------------------------------
# CDT - constrained Delaunay triangulation, header-only.
# Shewchuk's Triangle is deliberately not used (non-free licence).
# --------------------------------------------------------------------------
set(CDT_USE_AS_COMPILED_LIBRARY OFF CACHE BOOL "" FORCE)
FetchContent_Declare(cdt
    ${_n2s_fc_system}
    GIT_REPOSITORY https://github.com/artem-ogre/CDT.git
    GIT_TAG 1.4.5
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR CDT)
FetchContent_MakeAvailable(cdt)

# --------------------------------------------------------------------------
# Small utility libraries.
# --------------------------------------------------------------------------
set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
FetchContent_Declare(nlohmann_json
    ${_n2s_fc_system}
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.12.0
    GIT_SHALLOW TRUE)

set(FMT_INSTALL OFF CACHE BOOL "" FORCE)
set(FMT_TEST OFF CACHE BOOL "" FORCE)
set(FMT_DOC OFF CACHE BOOL "" FORCE)
FetchContent_Declare(fmt
    ${_n2s_fc_system}
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 11.2.0
    GIT_SHALLOW TRUE)

set(CLI11_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(CLI11_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(CLI11_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(CLI11_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_Declare(cli11
    ${_n2s_fc_system}
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.7.2
    GIT_SHALLOW TRUE)

FetchContent_MakeAvailable(nlohmann_json fmt cli11)

# --------------------------------------------------------------------------
# Catch2 v3 - tests only.
# --------------------------------------------------------------------------
if(N2S_BUILD_TESTS)
    set(CATCH_INSTALL_DOCS OFF CACHE BOOL "" FORCE)
    set(CATCH_INSTALL_EXTRAS ON CACHE BOOL "" FORCE)  # needed for catch_discover_tests
    FetchContent_Declare(Catch2
        ${_n2s_fc_system}
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.9.1
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
endif()

# --------------------------------------------------------------------------
# Polyscope (+ bundled Dear ImGui) - viewer only. Never linked into core/.
# --------------------------------------------------------------------------
if(N2S_BUILD_VIEWER)
    FetchContent_Declare(polyscope
        ${_n2s_fc_system}
        GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
        GIT_TAG v2.6.1
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(polyscope)
endif()
