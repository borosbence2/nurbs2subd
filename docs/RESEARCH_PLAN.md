# Research plan — nurbs2subd

**Research question.** Can trimmed NURBS patches be converted to Catmull–Clark
subdivision surfaces such that (a) neighbouring patches join with *exactly zero*
gap, (b) the approximation error to the original NURBS is bounded and measured,
and (c) surface quality near extraordinary vertices (EVs) is acceptable?

**Baseline literature.**
- Shen, Kosinka, Sabin, Dodgson (2014), *Conversion of trimmed NURBS surfaces to
  Catmull–Clark subdivision surfaces*, CAGD 31(7–8).
- Chen et al. (2008), *Progressive interpolation based on Catmull–Clark
  subdivision surfaces*, CGF 27(7).
- Halstead, Kass, DeRose (1993), *Efficient, fair interpolation using
  Catmull–Clark surfaces*, SIGGRAPH.
- Sederberg et al. (2008), *Watertight trimmed NURBS*, ACM TOG 27(3).
- Stam (1998), *Exact evaluation of Catmull–Clark subdivision surfaces at
  arbitrary parameter values*, SIGGRAPH.

**Known defects of the thesis implementation (do not repeat).**
1. Extraordinary-vertex limit weights were derived by hand and are wrong.
2. Trim curves were sampled uniformly per curve, ignoring length and curvature.
3. No distance or error metric against the original NURBS.
4. Only single patches were tested; watertightness was never demonstrated.
5. The progressive-iteration and direct-solve results were compared as if they
   were different methods. They solve the same problem; differences only
   measured non-convergence.
6. Curvature plots had no fixed color scale.

---

## Part A — Tooling

### M0 — Repository and build infrastructure
Tasks:
- [x] CMake project with presets `dev` (Debug + ASan/UBSan), `release`, `ci`.
- [x] FetchContent for all fixed dependencies, with pinned versions/tags.
- [x] `.clang-format`, `.clang-tidy`, `.editorconfig`, `.gitignore`.
- [x] GitHub Actions: Linux (GCC + Clang) and Windows (MSVC), build + ctest.
- [x] Empty `core` library, `viewer` and `cli` executables, one passing test.
- [x] README with build instructions. License: MIT (dependencies are MIT, MPL-2.0,
      or modified Apache 2.0, all compatible).

Exit: CI green on all three configurations.

### M1 — NURBS core
Tasks:
- [x] Knot vector type with validation (non-decreasing, correct length, clamped check).
- [x] B-spline basis functions and derivatives (Piegl & Tiller A2.2 / A2.3).
- [x] Rational curve and surface evaluation with 1st and 2nd derivatives.
- [x] Surface normal, first/second fundamental forms, mean and Gaussian curvature.
- [x] Closest-point projection onto a surface: grid initialization + Newton with
      domain clamping; report convergence status.
- [x] Knot insertion (needed later for sampling and tests).
- [x] IO: JSON format (degree, knots, control points, weights). **Scope change
      2026-09-19:** the legacy thesis `.bsc` reader is dropped; JSON is the only
      format `core/` reads. See Progress for the reasoning and for how the
      thesis surfaces get in.

Tests (analytic oracles):
- [x] Partition of unity; derivatives against central finite differences.
- [x] Exact circle/cylinder via rational quadratic NURBS: radius error < 1e-12.
- [x] Bilinear and bicubic Bézier patches against closed-form Bernstein evaluation.
- [x] Curvature of a sphere patch (1/r) and a cylinder (1/(2r) mean).
- [x] Projection round trip: project `S(u,v)` recovers `(u,v)`.

Exit: all tests pass; a rational surface survives a JSON file round trip
unchanged.

### M2 — Trimmed domain
Tasks:
- [ ] Trim loop = ordered list of 2D B-spline curves in the (u,v) domain.
- [ ] Validation: closure gaps (report, snap below tolerance), self-intersection
      check, orientation (outer CCW, holes CW; fix and warn if reversed).
- [ ] Adaptive sampling of trim curves by 3D arc length *on the surface* and by
      curvature, with a max-segment-length parameter. Uniform sampling only as an
      option for comparison. Curve knots are always included as samples.
- [ ] Constrained Delaunay triangulation of the trimmed domain (CDT), with holes,
      plus a refinement pass (max area / min angle).
- [ ] Map domain triangulation to 3D.
- [ ] Synthetic test-case generator: plane, paraboloid, saddle, cylinder
      patches × trims (circular hole, quarter-arc corner as in Shen et al. Fig. 1,
      L-shape, the two thesis domains).

Tests:
- Triangulated area equals the analytic trimmed area (e.g. square minus disc)
  within sampling tolerance, and converges as tolerance shrinks.
- Every domain triangle lies inside the trimmed region.
- Reversed-orientation input yields the same triangulation after auto-fix.

Exit: all synthetic cases and both thesis domains triangulate cleanly.

### M3 — Subdivision core (OpenSubdiv wrapper)
Tasks:
- [ ] `ControlMesh` (vertices, quad faces, boundary tags, crease/corner tags).
- [ ] Wrapper: `ControlMesh` → `Far::TopologyRefiner` (SCHEME_CATMARK, boundary
      interpolation EDGE_AND_CORNER, configurable).
- [ ] Limit evaluation at `(face, u, v)` returning position, first derivatives, normal,
      and second derivatives where available (`Far::PatchTable` + `PatchMap`).
- [ ] **Limit stencil export as a sparse matrix**: for a set of `(face,u,v)` sample
      locations, build the matrix `A` such that `limit_points = A * control_points`.
      This is the backbone of all fitting.
- [ ] Uniform refinement to level k for display.
- [ ] Export OBJ (control mesh and refined mesh).

Tests (oracles; these are the only allowed hand-written formulas):
- Regular vertex limit weights: centre 4/9, edge neighbours 1/9, diagonal 1/36.
- Valence-n interior vertex limit position: centre n/(n+5), each edge-connected
  1-ring vertex 4/(n(n+5)), each face-diagonal vertex 1/(n(n+5)). Check n = 3, 5, 6
  against OpenSubdiv's limit stencils.
- Eigenvalues of the regular (valence-4) local subdivision matrix:
  {1, 1/2, 1/2, 1/4, 1/4, 1/4, 1/8, 1/8, 1/16}.
- A regular grid of control points taken from a uniform bicubic B-spline
  reproduces that B-spline surface exactly (< 1e-12).
- Boundary limit curve with EDGE_AND_CORNER equals the cubic B-spline of the
  boundary control points (the property watertightness relies on).
- Rows of `A` sum to 1.

Exit: all oracle tests pass for valences 3–6.

### M4 — Viewer
Tasks:
- [ ] Polyscope app with an ImGui panel: load a case (JSON), run pipeline stages
      individually, toggle layers.
- [ ] Layers: NURBS surface, control net, trim curves (3D and domain view), domain
      triangulation, quad layout, Catmull–Clark control mesh, limit surface.
- [ ] Scalar maps: error-to-NURBS, mean curvature, Gaussian curvature. Fixed,
      user-editable range, colorbar always shown.
- [ ] Isophote / reflection-line stripes (static light direction, stripe colormap)
      as a view-independent zebra substitute.
- [ ] EV markers colored by valence.
- [ ] Side-by-side domain (2D) and surface (3D) views.
- [ ] Screenshot export with a fixed camera stored in the case file.

Exit: every figure from the thesis can be reproduced with correct color scales.

### M5 — Metrics and experiment harness
Tasks:
- [ ] Metrics module:
  - Parametric error: `|L(x) − S(φ(x))|` using the known domain correspondence.
  - Geometric error: limit sample → closest point on the trimmed NURBS.
  - Two-sided Hausdorff and RMS distance (dense sampling both ways).
  - Boundary deviation: limit boundary vs the trim curve image on the surface.
  - Normal deviation (degrees), curvature deviation.
  - Seam metrics for two patches: positional gap and normal angle along the seam.
- [ ] CLI `nurbs2subd run <config.json>`: runs the pipeline, writes
      `results/<run-id>/{config.json, metrics.json, samples.csv, *.obj, meta.json}`;
      `meta.json` holds the git hash, timestamp, and timings.
- [ ] `nurbs2subd sweep <config.json>` over parameter grids.
- [ ] Python: common loader, table generator (Markdown + LaTeX), plot helpers.
- [ ] Regression test: a small reference run whose metrics must stay within tolerance.

Exit: one command regenerates every table and plot from configs.

**GATE A.** All M0–M5 exit criteria met, CI green. Only then start Part B.

---

## Part B — Research

### R1 — Reproduce and correct the thesis baseline
Setup: both thesis test cases, with hand-authored quad layouts stored in JSON
(layout generation is deliberately out of scope here). Include the edge-thirding
refinement from the thesis as the baseline layout refinement.

**Prerequisite (from the M1 scope change):** the two thesis surfaces must exist
as `data/*.json` before this milestone starts. They are converted once, by a
throwaway script, from whatever the thesis wrote; the converter is not committed
and never enters `core/`. Blocked until the thesis files are to hand.

Tasks:
- [ ] Interpolation (square system): interior control points solved so that the
      limit surface passes through the layout vertices mapped to the NURBS; boundary
      handled as in Shen et al.
- [ ] Progressive iterative approximation (PIA), with a stopping criterion on the
      max update.
- [ ] Convergence study: PIA error vs iteration count, overlaid with the direct solve.
- [ ] Record: surface error, boundary error, curvature maps (fixed scale).

Expected finding: PIA converges to the direct solution; the thesis differences
were non-convergence. The report states this explicitly.

Exit: corrected baseline numbers for both thesis cases.

### R2 — Least-squares fitting with fairness
Formulation: minimise `‖A·V − P‖² + λ·E_fair(V)` over interior control points `V`,
with boundary control points constrained (elimination or KKT). `P` = NURBS points
at dense domain samples mapped to `(face,u,v)`; `E_fair` = discrete thin-plate
energy on the control mesh (optionally a limit-surface bending energy).

Tasks:
- [ ] Sample-location mapping from domain points to `(face,u,v)` via the layout's
      bilinear/parametric map; document the choice.
- [ ] Solver: sparse normal equations (`SimplicialLDLT`); fall back to `SparseQR`
      when ill-conditioned; report the condition estimate.
- [ ] Parameter sweeps: sample density, λ, refinement level of the layout.
- [ ] Plots: error vs control-point count (log–log), error vs λ, curvature quality vs λ.
- [ ] Compare against R1 on identical layouts.

Exit: a clear recommendation for the fitting method, with numbers.

### R3 — Watertight joins
Key property: with boundary rules enabled, the limit boundary curve depends only on
the boundary control points. If two patches share identical boundary control
points along the seam, the gap is exactly zero.

Tasks:
- [ ] Case 1 (split domain): one NURBS surface, domain cut by a curve into two
      trimmed patches. Fit the shared seam first as a cubic B-spline (the limit
      boundary curve), then fit both patches with the seam control points fixed.
- [ ] Case 2 (true intersection): two NURBS surfaces intersecting (e.g. cylinder and
      bicubic patch). Compute the intersection curve (marching + Newton refinement),
      preimages in both domains, fit a shared 3D seam curve, then fit both patches.
- [ ] Seam refinement compatibility: both sides must use the same number of seam
      control points; resolve mismatched layouts by refining the coarser side.
- [ ] Metrics: seam gap (must be ≤ 1e-12 × model size), normal angle along the
      seam, deviation of the seam from the true intersection curve.
- [ ] Tangent-continuity improvement (optional): add a soft constraint on normals
      across the seam; measure the G1 defect before and after.

Exit: zero gap demonstrated on both cases, with a bounded surface error.

### DECISION GATE B
Write `docs/notes/gate-b.md`: results so far, where the error concentrates
(boundary with high curvature? EVs? layout distortion?), and a choice of R4, R5, or
both. Ask the user before continuing.

### R4 — Automatic, error-driven quad layout (primary contribution candidate)
Tasks:
- [ ] Baseline automatic layout for simple domains (boundary-aligned, few EVs),
      or integrate an external layout if licensing allows (ask first).
- [ ] Replace edge thirding with an **error-driven adaptive refinement**: refine
      layout regions whose fitting error or boundary deviation exceeds a tolerance,
      while maintaining the constraint that EVs stay ≥ 2 rings from the boundary.
- [ ] Refinement placement driven by trim-curve curvature instead of knot positions.
- [ ] Compare: thirding vs uniform vs adaptive, at an equal control-point budget.

Exit: adaptive refinement beats thirding at an equal budget on all test cases,
or a documented negative result.

### R5 — Quality near extraordinary vertices
Tasks:
- [ ] Measure curvature fluctuation in rings around EVs, as a function of valence
      and refinement.
- [ ] Evaluate mitigation options: EV placement away from high-curvature regions,
      local fairness weighting, and (if available in OpenSubdiv) alternative EV
      rules or crease tags.
- [ ] Optional comparison against multi-sided/transfinite patches for the same regions.

Exit: quantified EV quality with at least one effective mitigation, or a
documented limitation.

### R6 — Write-up
Tasks:
- [ ] `docs/paper/` LaTeX skeleton (4–6 pages): problem, method, experiments,
      limitations.
- [ ] Every figure and table generated by scripts from committed configs.
- [ ] Test data and reference results packaged for release.
- [ ] Limitations section: exact conics are lost, EV quality, layout robustness,
      no STEP input yet.

Exit: a draft ready to send to a co-author / reviewer.

---

## Out of scope (until explicitly requested)
- STEP/IGES import, Blender add-on, GUI polish, GPU evaluation, performance work
  beyond "experiments finish in minutes".

## Progress
<!-- Claude Code: tick tasks above and add one dated line here per completed task. -->

### M0 — Repository and build infrastructure

- 2026-09-19 — Repository initialised (branch `master`). `RESEARCH_PLAN.md` moved
  to `docs/RESEARCH_PLAN.md`, the location CLAUDE.md already referenced.
- 2026-09-19 — CMake 3.24 project, C++20, Ninja. Presets `dev` (Debug,
  ASan+UBSan), `release`, `ci` (RelWithDebInfo). **Deviation:** added a fourth
  preset `ci-nogfx` (viewer off) so that the "core builds without graphics" hard
  rule is enforced by a CI job rather than by convention.
- 2026-09-19 — `cmake/Dependencies.cmake` fetches every fixed dependency at a
  pinned tag: Eigen 3.4.0, OpenSubdiv v3_6_1, CDT 1.4.5, Polyscope v2.6.1,
  nlohmann/json v3.12.0, fmt 11.2.0, CLI11 v2.7.2, Catch2 v3.9.1. Polyscope and
  Catch2 are only fetched when the viewer / tests are enabled. Dependency headers
  are included as SYSTEM so `-Werror` applies to first-party code only.
- 2026-09-19 — Warnings: `/W4 /WX` on MSVC, `-Wall -Wextra -Wpedantic` plus
  conversion, shadow and old-style-cast warnings as errors on GCC/Clang.
  Sanitizers degrade to a no-op on Windows (no MinGW runtime, no MSVC UBSan)
  instead of breaking `--preset dev`; CI runs the sanitized Debug build on Linux.
- 2026-09-19 — `cmake/GitVersion.cmake` captures commit, describe and dirty state
  at configure time and exposes them through `n2s::build_info()`. **Deviation
  (small, forward-looking):** not literally required by M0, but the M5 rule that
  every result records its git hash needs this plumbing, and it gives the one
  required passing test something real to assert.
- 2026-09-19 — `core` library (build provenance only), `nurbs2subd` CLI with
  `run`/`sweep` stubs, Polyscope viewer stub, 2 passing Catch2 tests.
  `.clang-format`, `.clang-tidy`, `.editorconfig`, `.gitignore`, MIT `LICENSE`,
  README with build instructions and a dependency/licence table.
- 2026-09-19 — GitHub Actions matrix: linux-gcc, linux-clang, windows-msvc
  (the three required configurations), plus linux-core-nogfx and a sanitized
  linux-asan job, plus a clang-format check.
- 2026-09-19 — Gotcha worth recording: on Windows, Git for Windows' bundled
  `mingw64/bin` precedes the real MinGW toolchain on PATH in a Git Bash shell, so
  freshly linked executables die with `STATUS_ENTRYPOINT_NOT_FOUND` before
  `main()`. The top-level `CMakeLists.txt` links the MinGW runtime statically to
  remove the DLL search order from the picture.
- 2026-09-19 — Verified locally (Windows, GCC 14.2 / MinGW, Ninja): `dev`,
  `release` and `ci-nogfx` all configure, build and pass `ctest`.
  **Exit criterion met:** on the first push, linux-gcc, linux-clang and
  windows-msvc all passed, along with linux-core-nogfx and the sanitized
  linux-asan job. The clang-format job failed on one line in
  `build_info.cpp.in`; the file is now formatted and the job pins clang-format
  to 19.1.7 so that the check cannot fail on version drift instead of on real
  formatting.

### M1 - NURBS core

- 2026-09-19 - `KnotVector`: validation as a class invariant (degree, length,
  monotonicity, degenerate domain, interior multiplicity above the degree),
  `find_span` per Piegl & Tiller A2.1 including the `u == U[m-p]` special case,
  multiplicity, clamped check, `uniform_clamped` factory. Tested against the
  worked example in The NURBS Book section 2.5.
- 2026-09-19 - B-spline basis functions (A2.2) and derivatives (A2.3). Oracles:
  partition of unity, non-negativity, reduction to Bernstein polynomials on a
  Bezier knot vector, derivative orders summing to zero, orders above the degree
  vanishing, and central finite differences.
- 2026-09-19 - `NurbsCurve` and `NurbsSurface`: rational evaluation and mixed
  partials via the homogeneous derivatives (A3.2 / A3.6) and the rational
  quotient rule (A4.2 / A4.4). Exactness oracles: circle, cylinder and a sphere
  of revolution all reproduced to better than 1e-12; bilinear and bicubic Bezier
  patches to 1e-14 against closed-form Bernstein evaluation.
- 2026-09-19 - Differential geometry: unit normal, first and second fundamental
  forms, mean, Gaussian and principal curvatures. Sphere gives 1/r, cylinder
  1/(2r) mean and zero Gaussian, plane zero. Degenerate points (the poles of the
  sphere of revolution) return `nullopt` rather than a silently wrong normal.
  Recorded in the header: near an umbilic the principal curvatures lose half
  their digits to cancellation under the square root, so error metrics and
  figures should use mean and Gaussian curvature.
- 2026-09-19 - Closest-point projection: grid search seeded per knot span, then
  Newton on the two perpendicularity conditions with clamping to the domain.
  Reports a status (`PointCoincident`, `Perpendicular`, `ClampedToBoundary`,
  `Stalled`, `IterationLimit`, `DegenerateJacobian`) rather than a bare bool, so
  that an unconverged projection cannot quietly flatter an error metric.
  Round-trip oracle from the plan passes to 1e-9 in parameter space.
- 2026-09-19 - Knot insertion (A5.1) for curves, and for surfaces by applying it
  line by line to the control net. Carried out on the weighted control points,
  which the exact-circle and exact-cylinder tests after insertion are there to
  enforce. Also tested: geometry unchanged, and Bezier extraction putting a
  control point exactly on the curve at full multiplicity.
- 2026-09-19 - JSON IO for curves and surfaces, with `num_u`/`num_v` stated
  redundantly so that a transposed control net is rejected instead of silently
  producing a plausible surface. Malformed documents name the offending field.
- 2026-09-19 - **Scope change, agreed with the user: `.bsc` support dropped.**
  The plan described the format only as "degree; knots; count; control points",
  which does not say whether a file holds a curve or a surface, whether a
  surface repeats those fields per direction, whether weights appear at all, or
  whether it is text or binary. Rather than guess, JSON becomes the only format
  `core/` reads, because it is strictly better for what this project needs: it
  is versioned and self-describing, it carries weights so rational surfaces
  survive at all, it states the control net shape redundantly so a transposed
  net is rejected, it is hand-editable and diffable in review, and it is already
  covered by round-trip tests. A legacy binary parser in `core/` would have to
  be carried, tested and trusted forever for two input files.
  **How the thesis surfaces get in:** once the files are to hand they are
  converted *once* by a throwaway script under `experiments/scripts/`, and the
  resulting JSON is committed to `data/`. The converter is not part of `core/`
  and does not need to be robust, only correct once, with the round trip through
  `n2s::io` checking the result. Recorded as a prerequisite on R1.
- 2026-09-19 - **M1 complete.** 53 tests pass on both the `dev` and `ci-nogfx`
  presets. M2's synthetic case generator supplies plane, paraboloid, saddle and
  cylinder test cases, so no milestone before R1 depends on the thesis data.
