include(Catch)

add_executable(n2s_tests
    test_build_info.cpp
    io/test_nurbs_json.cpp
    nurbs/test_basis.cpp
    nurbs/test_curve.cpp
    nurbs/test_differential.cpp
    nurbs/test_knot_insertion.cpp
    nurbs/test_knot_vector.cpp
    nurbs/test_projection.cpp
    nurbs/test_surface.cpp
    trim/test_sampling.cpp
    trim/test_trim_loop.cpp
    trim/test_triangulation.cpp
    trim/test_validate.cpp
    support/analytic_nurbs.cpp)

# Lets tests include "support/analytic_nurbs.hpp" regardless of their own depth.
target_include_directories(n2s_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(n2s_tests PRIVATE
    n2s::core
    Catch2::Catch2WithMain
    n2s::warnings
    n2s::sanitizers)

catch_discover_tests(n2s_tests)
