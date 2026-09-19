# Warning configuration. Applied to first-party targets only: dependencies are
# pulled in with SYSTEM include directories so their headers stay quiet.

add_library(n2s_warnings INTERFACE)
add_library(n2s::warnings ALIAS n2s_warnings)

if(MSVC)
    target_compile_options(n2s_warnings INTERFACE
        /W4
        /permissive-        # strict standard conformance
        /Zc:__cplusplus     # report the real __cplusplus value
        /Zc:preprocessor
        /utf-8
        /bigobj             # Eigen/OpenSubdiv templates overflow the default section limit
        /wd4127             # conditional expression is constant: fires inside Eigen
    )
    if(N2S_WARNINGS_AS_ERRORS)
        target_compile_options(n2s_warnings INTERFACE /WX)
    endif()
else()
    target_compile_options(n2s_warnings INTERFACE
        -Wall
        -Wextra
        -Wpedantic
        -Wshadow
        -Wnon-virtual-dtor
        -Wold-style-cast
        -Wcast-align
        -Wunused
        -Woverloaded-virtual
        -Wconversion
        -Wsign-conversion
        -Wdouble-promotion
        -Wformat=2
        -Wimplicit-fallthrough
    )
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(n2s_warnings INTERFACE
            -Wduplicated-cond
            -Wduplicated-branches
            -Wlogical-op
            -Wuseless-cast
        )
    endif()
    if(N2S_WARNINGS_AS_ERRORS)
        target_compile_options(n2s_warnings INTERFACE -Werror)
    endif()
endif()
