# nurbs2subd

Conversion of **trimmed NURBS surfaces** to **Catmull–Clark subdivision surfaces**,
with approximation error that is actually *measured* and patch joins that are
provably watertight.

This is the successor to an MSc thesis (BME VIK, 2023) whose implementation had
no quantitative validation. The research question, the baseline literature, the
known defects being corrected, and the milestone plan all live in
[`docs/RESEARCH_PLAN.md`](docs/RESEARCH_PLAN.md) — read that before touching
anything.

**Status: M0 (build infrastructure).** The geometry pipeline is not implemented
yet; `nurbs2subd run` and the viewer are deliberate stubs.

## Requirements

| | |
|---|---|
| Compiler | C++20: GCC 12+, Clang 15+, or MSVC 19.3x (VS 2022) |
| CMake | 3.24 or newer |
| Build tool | Ninja |
| Python | 3.11+ — analysis scripts only, not needed to build |

Every C++ dependency is fetched and built from source by CMake at a pinned tag;
nothing needs to be installed by hand. The first configure clones roughly 400 MB
of sources and takes a few minutes.

On Linux the viewer additionally needs the X11 and GL development headers:

```sh
sudo apt-get install -y xorg-dev libglu1-mesa-dev
```

Building with `-DN2S_BUILD_VIEWER=OFF` (or the `ci-nogfx` preset) removes that
requirement entirely — `core/` never links a graphics library.

## Build

```sh
cmake --preset dev            # Debug, ASan + UBSan where supported
cmake --build --preset dev
ctest --preset dev --output-on-failure
```

| Preset | Build type | Sanitizers | Viewer | Purpose |
|---|---|---|---|---|
| `dev` | Debug | ASan + UBSan | on | day-to-day work |
| `release` | Release | off | on | running experiments; the only preset whose timings mean anything |
| `ci` | RelWithDebInfo | off | on | what GitHub Actions builds |
| `ci-nogfx` | RelWithDebInfo | off | off | guards the "core builds without graphics" rule |

Binaries land in `build/<preset>/bin/`.

Sanitizers are enabled on Linux and macOS with GCC or Clang. On Windows the
`dev` preset still configures and builds, with sanitizers reported as `OFF` — the
MinGW toolchain ships no sanitizer runtime and MSVC has no UBSan. The CI matrix
therefore runs the sanitized Debug build on Linux.

## Run

```sh
./build/release/bin/nurbs2subd --version
./build/release/bin/nurbs2subd run experiments/configs/<name>.json
./build/release/bin/nurbs2subd_viewer
python experiments/scripts/plot_<name>.py results/<run-id>
```

Experiment output goes to `results/<run-id>/`, which is gitignored. Each run
records its config, the git commit and a timestamp, so any number in the write-up
can be traced back to the source state that produced it.

## Layout

```
core/         geometry library, no UI dependencies
              nurbs/ trim/ subd/ fit/ metrics/ io/
apps/cli/     nurbs2subd — headless experiment runner
apps/viewer/  Polyscope + ImGui viewer
tests/        Catch2 v3 tests, mirroring core/
data/         test cases (JSON, legacy .bsc); small files only
experiments/  configs/*.json, scripts/*.py; results/ is gitignored
docs/         research plan, notes, paper draft
cmake/        warning, sanitizer, git-provenance and dependency modules
```

## Dependencies

All pinned in [`cmake/Dependencies.cmake`](cmake/Dependencies.cmake).

| Library | Version | Licence | Used for |
|---|---|---|---|
| [Eigen](https://eigen.tuxfamily.org) | 3.4.0 | MPL-2.0 | linear algebra, sparse solvers |
| [OpenSubdiv](https://github.com/PixarAnimationStudios/OpenSubdiv) | v3_6_1 | modified Apache-2.0 | Catmull–Clark refinement, exact limit evaluation |
| [CDT](https://github.com/artem-ogre/CDT) | 1.4.5 | MPL-2.0 | constrained Delaunay triangulation |
| [Polyscope](https://polyscope.run) (+ Dear ImGui) | v2.6.1 | MIT | viewer |
| [nlohmann/json](https://github.com/nlohmann/json) | v3.12.0 | MIT | case and config IO |
| [fmt](https://github.com/fmtlib/fmt) | 11.2.0 | MIT | formatting |
| [CLI11](https://github.com/CLIUtils/CLI11) | v2.7.2 | BSD-3-Clause | CLI parsing |
| [Catch2](https://github.com/catchorg/Catch2) | v3.9.1 | BSL-1.0 | tests |

Subdivision stencils and limit weights come from OpenSubdiv only. Hand-derived
formulas were a defect of the thesis implementation and are permitted in this
repository solely as test oracles.

Shewchuk's *Triangle* is deliberately **not** used — its licence is not free.

## Licence

MIT, see [LICENSE](LICENSE). All dependency licences above are compatible with it.
