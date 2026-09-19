include(Catch)

add_executable(n2s_tests
    test_build_info.cpp)

target_link_libraries(n2s_tests PRIVATE
    n2s::core
    Catch2::Catch2WithMain
    n2s::warnings
    n2s::sanitizers)

catch_discover_tests(n2s_tests)
