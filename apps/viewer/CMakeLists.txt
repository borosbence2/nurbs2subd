add_executable(nurbs2subd_viewer
    main.cpp
    scene.cpp
    scene.hpp)

target_link_libraries(nurbs2subd_viewer PRIVATE
    n2s::core
    polyscope
    fmt::fmt
    n2s::warnings
    n2s::sanitizers)
