# core/ is the geometry library. Hard rule: it must build and test with no
# graphics dependency of any kind, so nothing here may link Polyscope, ImGui,
# GLFW or OpenGL.

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/build_info.cpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/build_info.cpp"
    @ONLY)

add_library(n2s_core STATIC
    include/n2s/build_info.hpp
    "${CMAKE_CURRENT_BINARY_DIR}/generated/build_info.cpp")

add_library(n2s::core ALIAS n2s_core)

target_include_directories(n2s_core PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include")

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
