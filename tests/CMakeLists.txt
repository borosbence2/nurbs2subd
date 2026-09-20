include(Catch)

add_executable(n2s_tests
    test_build_info.cpp
    io/test_case_json.cpp
    io/test_data_cases.cpp
    io/test_nurbs_json.cpp
    metrics/test_seam_error.cpp
    metrics/test_runner.cpp
    metrics/test_layout.cpp
    metrics/test_surface_error.cpp
    nurbs/test_basis.cpp
    nurbs/test_curve.cpp
    nurbs/test_differential.cpp
    nurbs/test_knot_insertion.cpp
    nurbs/test_knot_vector.cpp
    nurbs/test_projection.cpp
    nurbs/test_surface.cpp
    subd/test_subdivision.cpp
    trim/test_sampling.cpp
    trim/test_trim_loop.cpp
    trim/test_triangulation.cpp
    trim/test_validate.cpp
    support/analytic_nurbs.cpp)

# Lets tests include "support/analytic_nurbs.hpp" regardless of their own depth.
target_include_directories(n2s_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

# The committed case files are test input, so the tests have to be able to find
# them wherever the build tree happens to live.
target_compile_definitions(n2s_tests PRIVATE
    N2S_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../data")

target_link_libraries(n2s_tests PRIVATE
    n2s::core
    Catch2::Catch2WithMain
    n2s::warnings
    n2s::sanitizers)

catch_discover_tests(n2s_tests)
