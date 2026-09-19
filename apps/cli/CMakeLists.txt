add_executable(nurbs2subd main.cpp)

target_link_libraries(nurbs2subd PRIVATE
    n2s::core
    CLI11::CLI11
    fmt::fmt
    n2s::warnings
    n2s::sanitizers)
