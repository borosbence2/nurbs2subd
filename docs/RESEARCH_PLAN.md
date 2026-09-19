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
- [x] Trim loop = ordered list of 2D B-spline curves in the (u,v) domain.
- [x] Validation: closure gaps (report, snap below tolerance), self-intersection
      check, orientation (outer CCW, holes CW; fix and warn if reversed).
- [x] Adaptive sampling of trim curves by 3D arc length *on the surface* and by
      curvature, with a max-segment-length parameter. Uniform sampling only as an
      option for comparison. Curve knots are always included as samples.
- [x] Constrained Delaunay triangulation of the trimmed domain (CDT), with holes,
      plus a refinement pass (max area / min angle). **CDT ships no refinement;
      ours is hand-written. See Progress.**
- [x] Map domain triangulation to 3D.
- [~] Synthetic test-case generator: plane, paraboloid, saddle, cylinder
      patches × trims (circular hole, quarter-arc corner as in Shen et al. Fig. 1,
      L-shape, the two thesis domains). All but the thesis domains, which need
      the thesis data (same prerequisite as R1).

Tests:
- [x] Triangulated area equals the analytic trimmed area (e.g. square minus disc)
      within sampling tolerance, and converges as tolerance shrinks.
- [x] Every domain triangle lies inside the trimmed region.
- [x] Reversed-orientation input yields the same triangulation after auto-fix.

Exit: all synthetic cases triangulate cleanly (**met**: 12 of 12). Both thesis
domains: blocked on the thesis data, as for R1.

### M3 — Subdivision core (OpenSubdiv wrapper)
Tasks:
- [x] `ControlMesh` (vertices, quad faces, boundary tags, crease/corner tags).
- [x] Wrapper: `ControlMesh` → `Far::TopologyRefiner` (SCHEME_CATMARK, boundary
      interpolation EDGE_AND_CORNER, configurable).
- [x] Limit evaluation at `(face, u, v)` returning position, first derivatives, normal,
      and second derivatives where available. **Via the limit stencil matrices
      rather than `PatchTable` + `PatchMap`; see Progress.**
- [x] **Limit stencil export as a sparse matrix**: for a set of `(face,u,v)` sample
      locations, build the matrix `A` such that `limit_points = A * control_points`.
      This is the backbone of all fitting.
- [x] Uniform refinement to level k for display.
- [x] Export OBJ (control mesh and refined mesh).

Tests (oracles; these are the only allowed hand-written formulas):
- [x] Regular vertex limit weights: centre 4/9, edge neighbours 1/9, diagonal 1/36.
- [x] Valence-n interior vertex limit position: centre n/(n+5), each edge-connected
      1-ring vertex 4/(n(n+5)), each face-diagonal vertex 1/(n(n+5)). Checked for
      n = 3, 4, 5, 6 against OpenSubdiv's limit stencils.
- [x] Eigenvalues of the regular (valence-4) local subdivision matrix:
      {1, 1/2, 1/2, 1/4, 1/4, 1/4, 1/8, 1/8, 1/16}.
- [x] A regular grid of control points taken from a uniform bicubic B-spline
      reproduces that B-spline surface exactly (< 1e-12).
- [x] Boundary limit curve with EDGE_AND_CORNER equals the cubic B-spline of the
      boundary control points (the property watertightness relies on).
- [x] Rows of `A` sum to 1.

Exit: all oracle tests pass for valences 3-6. **Met.**

### M4 — Viewer
Tasks:
- [x] Polyscope app with an ImGui panel: load a case (JSON), run pipeline stages
      individually, toggle layers.
- [x] Layers: NURBS surface, control net, trim curves (3D and domain view), domain
      triangulation, quad layout, Catmull–Clark control mesh, limit surface.
      (The quad layout *is* the Catmull–Clark control mesh in this pipeline, so
      they are one layer, not two.)
- [x] Scalar maps: error-to-NURBS, mean curvature, Gaussian curvature. Fixed,
      user-editable range, colorbar always shown.
- [x] Isophote / reflection-line stripes (static light direction, stripe colormap)
      as a view-independent zebra substitute.
- [x] EV markers colored by valence.
- [x] Side-by-side domain (2D) and surface (3D) views. **Deviation:** one 3D
      scene with the domain placed beside the model, not two windows; see Progress.
- [x] Screenshot export with a fixed camera stored in the case file.

Exit: every figure from the thesis can be reproduced with correct color scales.
**Partially blocked**: the machinery is in place and colour scales are enforced,
but the thesis figures themselves need the thesis data, as R1 does.

### M5 — Metrics and experiment harness
Tasks:
- [x] Metrics module:
  - [x] Parametric error: `|L(x) − S(φ(x))|` using the known domain correspondence.
  - [x] Geometric error: limit sample → closest point on the trimmed NURBS.
  - [x] Two-sided Hausdorff and RMS distance (dense sampling both ways).
  - [x] Boundary deviation: limit boundary vs the trim curve image on the surface.
  - [x] Normal deviation (degrees), curvature deviation.
  - [x] Seam metrics for two patches: positional gap and normal angle along the seam.
- [x] CLI `nurbs2subd run <config.json>`: runs the pipeline, writes
      `results/<run-id>/{config.json, metrics.json, samples.csv, *.obj, meta.json}`;
      `meta.json` holds the git hash, timestamp, and timings.
- [x] `nurbs2subd sweep <config.json>` over parameter grids.
- [x] Python: common loader, table generator (Markdown + LaTeX), plot helpers.
- [x] Regression test: a small reference run whose metrics must stay within tolerance.

Exit: one command regenerates every table and plot from configs. **Met.**

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
- 2026-09-19 - M1 CI caught a real portability defect that neither GCC nor Clang
  reports: MSVC treats the `return` after Catch2's `FAIL` as unreachable
  (C4702), and `/WX` turns that into an error. Restructured so no branch is
  unreachable. The Windows job is earning its place.

### M2 - Trimmed domain

- 2026-09-19 - `NurbsCurveT<Dim>` replaces the 3D-only `NurbsCurve`, with
  `NurbsCurve = NurbsCurveT<3>` and `NurbsCurve2 = NurbsCurveT<2>`. Trim curves
  live in the `(u, v)` domain, and duplicating the rational evaluation and
  derivative code per dimension was the alternative. Knot insertion was
  templated with it. All 53 M1 tests passed unchanged through the refactor.
- 2026-09-19 - `TrimLoop` and `TrimRegion`, plus `reversed()` for curves (knot
  vector mirrored within its own domain, control points and weights reversed).
- 2026-09-19 - `validate_and_repair`: clamped-curve check, closure gaps (snapped
  below tolerance, reported above), self-intersection, orientation, hole
  containment and hole-hole overlap. Two deliberate ordering decisions:
  self-intersection is tested *before* orientation, because a symmetric bowtie
  has a signed area of exactly zero and would otherwise be reported as
  "degenerate area" instead of as the self-intersection it is; and closure gaps
  at or below 1e-15 are snapped silently, because a repair log full of 1e-17
  entries is a log nobody reads.
- 2026-09-19 - Adaptive trim sampling by recursive bisection, bounding both the
  model-space chord length and the sagitta, with every knot guaranteed to
  survive into the output. Uniform sampling kept as an option purely so the
  experiments can quantify what adaptive sampling buys, that being defect 2 of
  the thesis implementation. Measured on a cylinder: the curved edge collects
  more than four times the samples of the straight ruling, where uniform
  sampling gives both the same.
- 2026-09-19 - CDT triangulation with holes, and `map_to_surface`. Boundary
  vertices are stored first and the count is exposed, so boundary and interior
  can be told apart by index; R3 needs exactly that to make two patches agree on
  a seam.
- 2026-09-19 - **Deviation: CDT 1.4.5 provides no mesh refinement.** The plan
  assumed the library covered "max area / min angle"; it has no such entry
  point. Written by hand as `triangulate_refined`: circumcentre insertion with a
  centroid fallback, driven by model-space triangle area and domain-space
  minimum angle. Deliberately *not* a full Ruppert refinement -- boundary
  vertices are never added, moved or split, because the boundary discretisation
  is what two patches must agree on for a watertight join, and a refinement
  free to split boundary edges on its own schedule would break that agreement.
  The cost is the formal angle guarantee, covered instead by a vertex budget.
- 2026-09-19 - Bug worth recording, found by a failing area test. The first
  refinement kept candidate points clear of a triangle's corners by a *relative*
  margin, a fraction of that triangle's own shortest edge. For a sliver that
  margin is tiny, so the circumcentre landed almost on an existing vertex and
  produced a thinner sliver, which did it again. The mesh stayed a valid
  triangulation with exactly the right domain area throughout, so nothing
  complained -- but the worst 3D aspect ratio reached 1e13 and the cylinder area
  came out 0.45% high. Fixed with an absolute minimum separation derived from
  the vertex budget, `0.5 * sqrt(domain area / max_vertices)`, so the two cannot
  contradict each other. Worst aspect ratio fell to 13, area error to 0.086%.
  There is now a test asserting the aspect ratio directly, because the area
  test only caught this by luck.
- 2026-09-19 - **Finding that M5 needs to know about: a triangulation inscribed
  in a smooth surface can have more area than the surface, not less.** That is
  Schwarz's lantern, and anisotropic triangles on a cylinder are precisely the
  construction. Every refined cylinder mesh here overshoots. Worse, refining is
  not monotonically an improvement: the unrefined mesh is about 0.4% low, while
  the refined mesh at the smallest vertex budget is nearly 3% high. Any
  area-weighted error metric in M5 must bound the magnitude of its
  discretisation error and never assume its sign.
- 2026-09-19 - Performance: the refinement pass first ran the candidate-thinning
  test as a linear scan over every accumulated vertex, which is quadratic in the
  vertex count and took the test suite from 13 s to 436 s. Replaced with a
  uniform grid at one cell per separation distance; the suite is back to 34 s.
- 2026-09-19 - **M2 exit criterion met for the synthetic cases**: all 12 (4
  surfaces x 3 trims) validate and triangulate cleanly, and the triangulated
  area converges to the analytic `1 - pi r^2`. 90 tests pass on both the `dev`
  and `ci-nogfx` presets. The two thesis domains remain blocked on the thesis
  data, exactly as R1 is.
- 2026-09-19 - M2 CI fully green on all six jobs, windows-msvc included: the
  C4702 fix held.

### M3 - Subdivision core (OpenSubdiv wrapper)

- 2026-09-19 - `ControlMesh`: quad faces, crease and corner tags, OBJ export,
  and named constructors `grid`, `cube` and `vertex_fan(valence)`. Quads only,
  deliberately: for an all-quad mesh each face is exactly one ptex face, so the
  face index the limit locations use *is* the ptex index, rather than a mapping
  that has to be maintained and can silently drift. The constructor verifies the
  correspondence against `Far::PtexIndices` rather than assuming it.
- 2026-09-19 - `SubdivisionSurface`: topology descriptor to
  `Far::TopologyRefiner` (SCHEME_CATMARK, configurable boundary interpolation,
  default EDGE_AND_CORNER), adaptive refinement, limit stencil matrices, limit
  evaluation, uniform refinement, OBJ export.
- 2026-09-19 - **Double precision preserved.** Every `Far::...Table` convenience
  typedef in OpenSubdiv is `float`. The underlying `...Real<REAL>` templates are
  explicitly instantiated for `double` as well, so
  `LimitStencilTableFactoryReal<double>` is used directly and the "double
  everywhere in geometry code" rule holds through the subdivision layer. Worth
  knowing: using the obvious `Far::LimitStencilTableFactory` would have quietly
  capped every limit position, derivative and fitting matrix at single
  precision, and the resulting ~1e-7 error floor would have been indistinguishable
  from a genuine approximation error in M5.
- 2026-09-19 - **Deviation: limit evaluation goes through the limit stencil
  matrices, not `PatchTable` + `PatchMap`.** The plan named the latter. Both are
  correct, but routing evaluation through the same matrices the fitting uses
  means the two cannot disagree -- and a disagreement between the evaluator and
  the fitting matrix would not show up as a crash, only as an unexplained
  residual somewhere in R1 or R2.
- 2026-09-19 - Added `local_subdivision_matrix(vertex)`, the map from a vertex
  and its 1-ring to the same after one Catmull-Clark step. Needed for the
  eigenvalue oracle, and R5 needs it again to study behaviour at extraordinary
  vertices. Built by refining the identity through `PrimvarRefinerReal<double>`,
  so it comes from OpenSubdiv's stencils rather than from any formula written
  here.
- 2026-09-19 - **Every oracle passed on the first run.** In particular the
  valence-3, 5 and 6 limit masks agree with `n/(n+5)`, `4/(n(n+5))`,
  `1/(n(n+5))` to 1e-12, the regular spectrum is exactly
  {1, 1/2, 1/2, 1/4, 1/4, 1/4, 1/8, 1/8, 1/16}, and a regular grid reproduces
  the uniform bicubic B-spline to better than 1e-12. Defect 1 of the thesis --
  hand-derived extraordinary-vertex limit weights that were wrong -- is now
  covered by a test that would catch it.
- 2026-09-19 - The watertightness precondition is asserted directly rather than
  inferred: for locations on a boundary, every stencil weight on a non-boundary
  control point is zero to 1e-14. That is the property R3 depends on, stated as
  a fact about the matrix rather than as a measured gap that happened to be
  small.
- 2026-09-19 - Build gotcha: OpenSubdiv's `Sdc` scheme headers use `M_PI`, which
  is not standard C++ and is absent under the strict conformance this project
  builds with. Fixed by putting `_USE_MATH_DEFINES` and `_DEFAULT_SOURCE` on the
  dependency's interface target, covering MSVC, MinGW and glibc.
- 2026-09-19 - Behaviour worth knowing: OpenSubdiv *silently drops* limit
  locations that have no limit surface beneath them rather than failing. With
  `BoundaryInterpolation::None` that is the entire boundary region. The wrapper
  compares the stencil count against the request and raises an error naming the
  cause; there is a test pinning the behaviour down, because a silently short
  result would otherwise misalign every downstream sample array.
- 2026-09-19 - **M3 exit criterion met.** 108 tests pass on both the `dev` and
  `ci-nogfx` presets.

### M4 - Viewer

- 2026-09-19 - **Case file format** (`n2s-case`), which M4 needed and M5 needs
  again: surface, trim region, optional control mesh, optional camera, and
  colour ranges. The camera lives in the file because a screenshot whose
  viewpoint exists only in whoever took it is not reproducible, and the plan
  requires every figure to be regenerable from committed inputs.
- 2026-09-19 - Curve JSON now covers domain curves as well as model curves,
  under a separate `n2s-curve2` tag rather than a dimension field. A trim curve
  handed to a reader expecting model geometry then fails immediately instead of
  producing a curve in the wrong space. There is a test for both directions.
- 2026-09-19 - `nurbs2subd export-cases` writes all 12 synthetic cases to
  `data/` (76 KB total, well under the 1 MB rule), so `data/` is generated from
  the same code the tests use rather than hand-maintained. A test loads every
  committed file, validates, triangulates and checks the size cap, which is what
  catches a change to the case definitions invalidating what is on disk.
- 2026-09-19 - Polyscope viewer with every layer the plan lists, an ImGui panel
  for the pipeline parameters, editable colour ranges, isophote stripes about a
  fixed light direction, and EV markers coloured by valence on a fixed 3-to-8
  scale so the colours mean the same thing from case to case.
- 2026-09-19 - **The no-auto-ranging rule is enforced, not just documented.**
  Defect 6 of the thesis was curvature plots with no fixed colour scale. A
  scalar with no range in the case file is still drawn -- exploring is useful --
  but it is labelled `AUTO - not comparable`, and **the screenshot button is
  disabled while any such field exists**. Exporting a figure is the moment the
  rule has to bite, so that is where it bites. Editing a range in the panel
  counts as pinning it, because the user chose the numbers.
- 2026-09-19 - **Deviation: side-by-side views.** Polyscope draws a single 3D
  scene, so the domain triangulation and the domain-space trim curves are placed
  beside the model rather than shown in a second window. Same information,
  one viewport; a second window would have meant a second renderer.
- 2026-09-19 - **Deviation: error-to-NURBS is computed on demand**, behind a
  button, because it costs a closest-point projection per limit sample and the
  viewer should come up promptly. Unconverged projections are counted and
  reported in the panel rather than folded silently into the map, since an
  unconverged projection gives an upper bound at best and a quietly optimistic
  error picture is worse than none.
- 2026-09-19 - Verification note: **the drawing itself cannot be checked
  automatically.** The geometry behind it can, so the viewer is split in two --
  `scene.cpp` knows nothing about Polyscope, and `nurbs2subd_viewer --selftest`
  builds a whole scene headlessly, asserts every layer came back non-empty and
  that no scalar was auto-ranged, then exits. CI runs it. Whether the picture
  *looks* right still needs a human, and nothing here claims otherwise.
- 2026-09-19 - CI gained three smoke steps: the case exporter, the viewer
  selftest, and the existing CLI version check.
- 2026-09-19 - **M4 exit criterion partially met.** 121 tests pass on both
  presets and the viewer selftest passes. Reproducing the thesis figures
  specifically still needs the thesis data, exactly as R1 does.

### M5 - Metrics and experiment harness

- 2026-09-19 - `metrics/statistics`: max, mean and RMS, with failed
  measurements counted separately and **excluded** rather than folded in.
  Folding a failed projection in as zero flatters the result; folding it in as
  something large invents data. A non-finite sample is treated as a failed
  measurement for the same reason -- one NaN in the sums turns every statistic
  into NaN and hides which sample was at fault. `max` is reported before `rms`
  throughout: a fitting method is easy to make look good on RMS while leaving
  one region badly wrong.
- 2026-09-19 - `metrics/surface_error`: parametric, geometric, two-sided
  Hausdorff and RMS, normal deviation and curvature deviation, plus boundary
  deviation against the trim curve's *image on the surface*. The domain
  correspondence is a `std::function` supplied by the caller; when none is
  given, the correspondence-dependent statistics come back absent rather than
  invented from a guessed mapping.
- 2026-09-19 - `curvature_from_derivatives` was factored out of the NURBS
  differential code so the limit surface's curvature is computed by *the same
  formulas* as the NURBS it is compared against. Two implementations of the same
  mathematics would make the curvature-deviation metric measure the difference
  between two formulas rather than between two surfaces.
- 2026-09-19 - Bug found by an exact oracle. The reverse (NURBS → limit)
  distance first measured to the nearest tessellation *vertex*, which puts a
  floor of half the sample spacing under every result: it read 0.0707 where the
  true distance was 0, so the two-sided Hausdorff reported the tessellation
  density rather than the coverage error. Fixed to measure against the faces.
- 2026-09-19 - `metrics/seam_error`. **Two patches whose shared boundary control
  points are identical join with a measured gap below 1e-14** -- zero, not
  small. That is R3's central claim, and it is now measurable before R3 starts.
  Recorded alongside it: the normal deviation across such a seam is *not* zero.
  Shared control points give G0 and nothing more, and a zero gap is easy to
  mistake for a smooth join.
- 2026-09-19 - `fit::grid_layout`, the naive baseline layout, and
  `experiment::run`/`sweep`. A run writes config.json, metrics.json,
  samples.csv, three OBJs and meta.json with the git hash, the UTC timestamp and
  per-stage timings. The provenance is asserted by a test rather than trusted.
- 2026-09-19 - **Deviation: a new `core/experiment/` directory.** The run and
  sweep harness is neither IO nor geometry, and filing it under `io/` would have
  mislabelled it. CLAUDE.md's layout section and the README were updated to
  match. `core/fit/` was created at the same time and *is* in the documented
  layout.
- 2026-09-19 - Bug the sweep tests now pin down: a sweep serialises its base
  config, which already holds a resolved case path, then re-reads it. Resolving
  against the base directory a second time produced a doubled-up path and the
  sweep could not open its own case.
- 2026-09-19 - **Finding worth carrying into R2 and R4.** A sweep over layout
  density on the saddle case shows that an *anisotropic* layout can be worse
  than a coarser isotropic one: 7x5 (35 control points) has a higher error than
  5x5 (25). The mechanism is that a uniform cubic B-spline smooths a quadratic
  by about `h^2/6 * f''`, so on `z = c(x^2 - y^2)` the errors from the two
  directions **cancel when the spacings match** and stop cancelling when they do
  not. Predicted `(c/3)|h_x^2 - h_y^2|` against measured RMS: 0.0081 vs 0.0081
  at 9x7, 0.0231 vs 0.0197 at 7x5. R4's error-driven refinement must therefore
  not refine one direction in isolation, and this is a concrete case where
  "more control points" is not "less error".
- 2026-09-19 - Python: `n2s_results.py` (loader, Markdown and LaTeX tables from
  one set of formatted values so the two cannot disagree), `plot_sweep.py`,
  `make_tables.py`. Tables carry their provenance underneath -- run id, git
  describe, timestamp -- and say so explicitly when a run had unmeasured
  samples.
- 2026-09-19 - A regression run is pinned in `tests/metrics/test_runner.cpp`:
  a small fixed config whose metrics must stay within a relative tolerance. It
  deliberately turns the interior refinement off, that being the part of the
  pipeline most sensitive to floating-point differences between platforms.
- 2026-09-19 - **M5 exit criterion met.** 150 tests pass on both the `dev` and
  `ci-nogfx` presets, and one command regenerates the tables and the figure from
  a result directory.

---

## GATE A status

M0-M5 are complete except for the parts that need the MSc thesis data, which
are the same parts R1 is blocked on:

- M1: the legacy `.bsc` reader was dropped by agreement; the thesis surfaces
  will be converted once to JSON when the files are to hand.
- M2: the two thesis domains are not among the synthetic cases.
- M4: the thesis figures cannot be reproduced without the thesis cases.

Everything that does not depend on that data is done and tested. **Before Part B
starts, the user has to supply the two thesis test cases**; without them R1
cannot reproduce or correct the baseline it exists to correct.
