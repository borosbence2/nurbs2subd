# ASan + UBSan for Debug builds. The research plan mandates these on Linux; the
# MinGW toolchain ships no sanitizer runtime and MSVC has no UBSan, so on those
# platforms the option degrades to a no-op instead of breaking the `dev` preset.

add_library(n2s_sanitizers INTERFACE)
add_library(n2s::sanitizers ALIAS n2s_sanitizers)

set(N2S_SANITIZERS_ACTIVE "OFF" CACHE INTERNAL "" FORCE)

if(N2S_ENABLE_SANITIZERS)
    if(WIN32)
        message(STATUS "Sanitizers requested but unsupported on this platform - skipping")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set(_n2s_san_flags
            -fsanitize=address,undefined
            -fno-omit-frame-pointer
            -fno-sanitize-recover=all)
        target_compile_options(n2s_sanitizers INTERFACE ${_n2s_san_flags})
        target_link_options(n2s_sanitizers INTERFACE ${_n2s_san_flags})
        set(N2S_SANITIZERS_ACTIVE "address,undefined" CACHE INTERNAL "" FORCE)
    else()
        message(STATUS "Sanitizers requested but ${CMAKE_CXX_COMPILER_ID} is unsupported - skipping")
    endif()
endif()
