# Capture the git revision at configure time so that every experiment result can
# record the exact commit it came from (hard rule in CLAUDE.md). Re-run CMake
# after committing to refresh it; `meta.json` writing in M5 re-reads it at
# runtime as well, so a stale value here never silently mislabels a result.

find_package(Git QUIET)

set(N2S_GIT_SHA "unknown")
set(N2S_GIT_DESCRIBE "unknown")
set(N2S_GIT_DIRTY "unknown")

if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        OUTPUT_VARIABLE N2S_GIT_SHA
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)

    execute_process(
        COMMAND "${GIT_EXECUTABLE}" describe --always --dirty --tags
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        OUTPUT_VARIABLE N2S_GIT_DESCRIBE
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)

    execute_process(
        COMMAND "${GIT_EXECUTABLE}" status --porcelain
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        OUTPUT_VARIABLE _n2s_git_status
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)

    if(_n2s_git_status STREQUAL "")
        set(N2S_GIT_DIRTY "clean")
    else()
        set(N2S_GIT_DIRTY "dirty")
    endif()

    # On an unborn branch `git rev-parse HEAD` fails but still echoes the
    # literal string "HEAD" on stdout, so the result has to be validated rather
    # than merely tested for emptiness.
    if(NOT N2S_GIT_SHA MATCHES "^[0-9a-f]+$")
        set(N2S_GIT_SHA "uncommitted")
        set(N2S_GIT_DESCRIBE "uncommitted")
    endif()
endif()
