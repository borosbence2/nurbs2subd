cmake_minimum_required(VERSION 3.24)

project(nurbs2subd
    VERSION 0.1.0
    DESCRIPTION "Trimmed NURBS to Catmull-Clark subdivision surface conversion with measured error"
    LANGUAGES CXX)

if(PROJECT_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed. Use `cmake --preset dev`.")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(MINGW)
    # Git for Windows ships its own, older libstdc++-6.dll in <Git>/mingw64/bin,
    # which precedes the real toolchain on PATH inside a Git Bash shell. Freshly
    # linked executables then fail to start with STATUS_ENTRYPOINT_NOT_FOUND
    # (0xC0000139) before main() is reached. Linking the runtime statically
    # removes the DLL search order from the picture and makes the binaries
    # self-contained. Applied before the dependencies are added so that it
    # covers every target in the build, not just ours.
    add_link_options(-static-libgcc -static-libstdc++ -static)
endif()

# Single output directory keeps `./build/<preset>/apps/cli/nurbs2subd` predictable
# across generators (MSVC otherwise appends a per-config subdirectory).
foreach(_cfg "" _DEBUG _RELEASE _RELWITHDEBINFO _MINSIZEREL)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY${_cfg} "${CMAKE_BINARY_DIR}/bin")
endforeach()

option(N2S_BUILD_VIEWER "Build the Polyscope viewer application" ON)
option(N2S_BUILD_TESTS "Build the Catch2 test suite" ON)
option(N2S_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" ON)
option(N2S_ENABLE_SANITIZERS "Enable ASan + UBSan (Linux/macOS, GCC/Clang only)" OFF)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include(CompilerWarnings)
include(Sanitizers)
include(GitVersion)
include(Dependencies)

if(N2S_BUILD_TESTS)
    enable_testing()
endif()

add_subdirectory(core)
add_subdirectory(apps/cli)

if(N2S_BUILD_VIEWER)
    add_subdirectory(apps/viewer)
endif()

if(N2S_BUILD_TESTS)
    add_subdirectory(tests)
endif()

message(STATUS "")
message(STATUS "nurbs2subd ${PROJECT_VERSION} (${N2S_GIT_DESCRIBE})")
message(STATUS "  build type ......... ${CMAKE_BUILD_TYPE}")
message(STATUS "  compiler ........... ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  viewer ............. ${N2S_BUILD_VIEWER}")
message(STATUS "  tests .............. ${N2S_BUILD_TESTS}")
message(STATUS "  warnings as errors . ${N2S_WARNINGS_AS_ERRORS}")
message(STATUS "  sanitizers ......... ${N2S_SANITIZERS_ACTIVE}")
message(STATUS "")
