# core/ is the geometry library. Hard rule: it must build and test with no
# graphics dependency of any kind, so nothing here may link Polyscope, ImGui,
# GLFW or OpenGL.

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/build_info.cpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/build_info.cpp"
    @ONLY)

add_library(n2s_core STATIC
    include/n2s/build_info.hpp
    include/n2s/io/nurbs_json.hpp
    include/n2s/nurbs/basis.hpp
    include/n2s/nurbs/curve.hpp
    include/n2s/nurbs/detail/binomial.hpp
    include/n2s/nurbs/differential.hpp
    include/n2s/nurbs/knot_insertion.hpp
    include/n2s/nurbs/knot_vector.hpp
    include/n2s/nurbs/projection.hpp
    include/n2s/nurbs/surface.hpp
    include/n2s/tolerances.hpp
    include/n2s/trim/cases.hpp
    include/n2s/trim/sampling.hpp
    include/n2s/trim/trim_loop.hpp
    include/n2s/trim/triangulation.hpp
    include/n2s/trim/validate.hpp
    src/io/nurbs_json.cpp
    src/nurbs/basis.cpp
    src/nurbs/differential.cpp
    src/nurbs/knot_insertion.cpp
    src/nurbs/knot_vector.cpp
    src/nurbs/projection.cpp
    src/nurbs/surface.cpp
    src/trim/cases.cpp
    src/trim/sampling.cpp
    src/trim/trim_loop.cpp
    src/trim/triangulation.cpp
    src/trim/validate.cpp
    "${CMAKE_CURRENT_BINARY_DIR}/generated/build_info.cpp")

add_library(n2s::core ALIAS n2s_core)

target_include_directories(n2s_core
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
    PRIVATE
        # Implementation-only headers under src/, e.g. nurbs/detail/binomial.hpp.
        "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_link_libraries(n2s_core
    PUBLIC
        Eigen3::Eigen
        n2s::opensubdiv
        CDT::CDT
        nlohmann_json::nlohmann_json
        fmt::fmt
        n2s::sanitizers
    PRIVATE
        n2s::warnings)
