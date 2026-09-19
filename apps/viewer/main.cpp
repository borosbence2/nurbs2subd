// Polyscope viewer. The layer set and the ImGui control panel arrive with M4;
// for now this opens an empty window so that the graphics dependency chain is
// verified by the build.
//
// `--version` deliberately returns before polyscope::init so that CI, which has
// no display, can still smoke-test that the binary links and starts.

#include "n2s/build_info.hpp"

#include <fmt/core.h>
#include <polyscope/polyscope.h>

#include <string>

int main(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-v" || arg == "--version") {
            fmt::print("{}\n", n2s::build_info_string());
            return 0;
        }
    }

    polyscope::options::programName = "nurbs2subd viewer";
    polyscope::options::verbosity = 0;
    polyscope::init();

    polyscope::warning("No case loaded. Loading and the pipeline layers arrive with plan "
                       "milestone M4.");

    polyscope::show();
    return 0;
}
