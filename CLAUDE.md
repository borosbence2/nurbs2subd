# CLAUDE.md — nurbs2subd

Research codebase: convert trimmed NURBS surfaces to Catmull–Clark subdivision
surfaces with **measured** approximation error and **watertight** joins between
neighbouring patches. Successor to an MSc thesis (BME VIK, 2023) that lacked
quantitative validation. The full plan is in `docs/RESEARCH_PLAN.md` — read it
before starting any task.

## How to work
- Work milestone by milestone, in plan order. Do not start a research milestone
  (R*) before every tooling milestone (M*) meets its exit criteria.
- One task = one focused commit (Conventional Commits: `feat:`, `fix:`, `test:`…).
- For math code: write the test (analytic oracle) first, then the implementation.
- After each task, update the `## Progress` section of `docs/RESEARCH_PLAN.md`
  (checkbox + one line of notes, including any deviation from the plan).
- If a task is ambiguous or the plan looks wrong, stop and ask. Do not silently
  change scope.
- Ask before adding any dependency not listed below.

## Stack (fixed)
- C++20, CMake ≥ 3.24, dependencies via `FetchContent` (vcpkg only if a package
  fails to build that way).
- Eigen 3.4 (linear algebra, sparse solvers)
- OpenSubdiv 3.6+ (Catmull–Clark refinement, **exact limit evaluation**)
- CDT (artem-ogre, header-only) — constrained Delaunay. **Never use Shewchuk's
  Triangle** (non-free license).
- Polyscope + Dear ImGui (viewer). **No Qt, no libQGLViewer, no legacy OpenGL.**
- nlohmann/json, fmt, CLI11, Catch2 v3
- Python 3.11+ (numpy, pandas, matplotlib) for analysis scripts only.
- Optional, later and behind a CMake option: OpenMesh or geometry-central,
  OpenCASCADE (STEP import).

## Layout
```
core/        geometry library (no UI deps): nurbs/, trim/, subd/, fit/, metrics/, io/
apps/viewer/ Polyscope app
apps/cli/    headless experiment runner (nurbs2subd)
tests/       Catch2 tests, mirror core/ structure
data/        test cases (JSON), small files only
experiments/ configs/*.json, scripts/*.py; results/ is gitignored
docs/        plan, notes, paper draft
```

## Hard rules
- `core/` must build and test without any graphics library.
- **Never hand-derive subdivision stencils or limit weights.** Get them from
  OpenSubdiv (`Far::LimitStencilTable` / `Far::PatchTable`). Hand-written
  formulas may exist only as test oracles (see plan M3).
- Every scalar visualization uses an explicit, displayed color range. No
  auto-ranging in figures meant for comparison.
- Every experiment is reproducible from a config file: results record the config,
  git hash, and timestamp.
- Warnings are errors (`-Wall -Wextra -Wpedantic` / `/W4`). Debug builds run with
  ASan + UBSan on Linux.
- `double` everywhere in geometry code. Tolerances are named constants, never
  magic numbers.
- No files > 1 MB committed to `data/`.

## Commands
```
cmake --preset dev && cmake --build --preset dev
ctest --preset dev --output-on-failure
./build/dev/apps/cli/nurbs2subd run experiments/configs/<name>.json
python experiments/scripts/plot_<name>.py results/<run-id>
```
