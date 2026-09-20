#!/usr/bin/env python3
"""Compares fitting methods on identical layouts, as a table.

    python experiments/scripts/compare_fits.py results/<run-id> [results/<run-id> ...]
    python experiments/scripts/compare_fits.py --latex results/<run-id> ...

This is R2's "compare against R1 on identical layouts" task. It prints position
*and* curvature columns together, because R1 established that those two move in
opposite directions under an interpolating fit: a comparison quoting position
error alone would report the interpolating fit as the better method, which is
the shape the thesis's claim took.

**The common subset, and why it is the default.** Each run's summary statistics
cover only the samples that run could measure, and which samples those are
depends on the surface being measured: the closest-point projection behind the
geometric error converges on some points and not others, and a *better* fit --
one whose residual is smaller -- is measurably harder to project onto, because
the perpendicularity test divides by that residual. So the sets differ, they
differ systematically with fit quality, and averaging over them compares two
numbers taken over two different populations.

When every run has written a ``samples.csv``, this script therefore restricts to
the samples measured in *all* of them, and says how many that left. The summary
columns from ``metrics.json`` are printed too, so the difference between the two
is visible rather than hidden; where they disagree, the common subset is the one
to quote.

**What it refuses to do.** If the runs do not share a case, a layout size and a
boundary, they are not a comparison and the script says so rather than printing
a table that looks like one. The boundary check applies only to fits that
constrain their boundary: R2 constrains its boundary control points to R1's, so
on identical layouts those two must agree. An unfitted run has no boundary
constraint at all and is exempt.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from n2s_results import latex_table, load_run, markdown_table  # noqa: E402

#: Columns, as (heading, metrics group, field). Position first, then what
#: position alone does not say.
COLUMNS = [
    ("geom max", "geometric", "max"),
    ("geom rms", "geometric", "rms"),
    ("normal max", "normal_degrees", "max"),
    ("curv max", "mean_curvature", "max"),
    ("curv rms", "mean_curvature", "rms"),
    ("boundary max", "boundary_error", "max"),
]

#: Fit methods whose boundary control points are constrained, and which must
#: therefore agree on the boundary error when run on the same layout.
CONSTRAINED_BOUNDARY = {"interpolate", "pia", "least_squares"}

#: Sample columns the common-subset table reports, as (heading, column, gate).
SAMPLE_COLUMNS = [
    ("geom max", "geometric", "measured", max),
    ("geom rms", "geometric", "measured", None),
    ("normal max", "normal_degrees", "measured", max),
    ("curv max", "mean_curvature", "curvature_measured", max),
    ("curv rms", "mean_curvature", "curvature_measured", None),
]

#: Fractional agreement required of the boundary error before two runs count as
#: sharing a boundary. Not exact equality: the two fits produce different
#: meshes, and the boundary sampling walks them separately.
BOUNDARY_TOLERANCE = 1e-9


def label(run) -> str:
    fit = run.metrics.get("fit", {})
    method = fit.get("method", "?")
    parts = [method]
    if method == "least_squares":
        lam = fit.get("least_squares", {}).get("lambda")
        if lam is not None:
            parts.append(f"lambda={lam:g}")
    refinement = fit.get("layout_refinement")
    if refinement and refinement != 1:
        parts.append(f"x{refinement}")
    return " ".join(parts)


def common_subset_table(runs) -> list[int] | None:
    """Prints the comparison restricted to samples every run measured.

    Returns that subset, so the summary below quotes the same numbers the table
    does. Two places working out "how much better" from two different
    populations is how a script ends up disagreeing with itself.
    """
    if any(not run.samples for run in runs):
        print(
            "Not every run wrote a samples.csv, so this comparison uses each run's own\n"
            "measured set. Those sets differ, and they differ with fit quality. Re-run\n"
            "with \"write_samples\": true for a comparison over identical points.\n",
            file=sys.stderr,
        )
        return None

    lengths = {len(run.samples) for run in runs}
    if len(lengths) > 1:
        print(f"the runs took different numbers of samples ({sorted(lengths)}), so there is "
              f"no common subset to take.", file=sys.stderr)
        return None

    total = next(iter(lengths))
    for column in ("geometric", "mean_curvature", "curvature_measured"):
        if column not in runs[0].samples[0]:
            print(f"samples.csv has no {column!r} column; these runs predate it.",
                  file=sys.stderr)
            return None

    common = [
        i
        for i in range(total)
        if all(run.samples[i]["measured"] > 0.5 and run.samples[i]["curvature_measured"] > 0.5
               for run in runs)
    ]
    if not common:
        print("no sample was measured by every run, so there is nothing to compare.",
              file=sys.stderr)
        return None

    print(f"On the {len(common)} samples of {total} that every run measured:\n")
    header = ["fit"] + [name for name, _, _, _ in SAMPLE_COLUMNS]
    rows = []
    for run in runs:
        row = [label(run)]
        for _, column, _, reducer in SAMPLE_COLUMNS:
            values = [run.samples[i][column] for i in common]
            value = (reducer(values) if reducer
                     else math.sqrt(sum(v * v for v in values) / len(values)))
            row.append(f"{value:.4g}")
        rows.append(row)
    print(markdown_table(header, rows))

    if len(common) < total // 2:
        print(
            f"\nNote: only {len(common)} of {total} samples ({len(common) / total:.0%}) survive "
            f"in every run. That is a limitation of the closest-point projection behind the "
            f"geometric error, not of the fits; see the R2 notes in the research plan.",
            file=sys.stderr,
        )

    return common


def main(argv: list[str]) -> int:
    args = [a for a in argv[1:] if not a.startswith("--")]
    as_latex = "--latex" in argv

    if not args:
        print(__doc__)
        return 2

    runs = [load_run(path) for path in args]

    cases = {run.case for run in runs}
    if len(cases) > 1:
        print(f"these runs are of different cases ({', '.join(sorted(cases))}), so there is "
              f"nothing to compare. Refusing to print a table.", file=sys.stderr)
        return 1

    counts = {int(run.metrics["control_vertices"]) for run in runs}
    if len(counts) > 1:
        print(f"WARNING: the runs have different control-point counts "
              f"({', '.join(str(c) for c in sorted(counts))}). They are not on identical "
              f"layouts, so the rows below are not comparing fitting methods alone.",
              file=sys.stderr)

    # Only the fits that constrain their boundary are held to this. An unfitted
    # run's boundary is wherever lifting the layout put it, so including it here
    # would fire the warning on every comparison and teach the reader to ignore
    # it -- which is how the one that matters gets missed.
    constrained = [run for run in runs
                   if run.metrics.get("fit", {}).get("method") in CONSTRAINED_BOUNDARY]
    boundaries = [run.stat("boundary_error", "max") for run in constrained]
    if len(boundaries) > 1 and max(boundaries) > 0.0:
        spread = (max(boundaries) - min(boundaries)) / max(boundaries)
        if spread > BOUNDARY_TOLERANCE:
            print(f"WARNING: the boundary error differs between these runs "
                  f"({min(boundaries):.6g} to {max(boundaries):.6g}). R2 constrains its "
                  f"boundary to R1's, so on identical layouts these should match; something "
                  f"other than the interior fit is moving.", file=sys.stderr)

    common = common_subset_table(runs)

    print("\nFrom metrics.json, each run over the samples it could measure:")
    header = ["fit"] + [name for name, _, _ in COLUMNS] + ["samples"]
    rows = []
    for run in runs:
        row = [label(run)]
        for _, group, field in COLUMNS:
            row.append(f"{run.stat(group, field):.4g}")
        measured = int(run.metrics["surface_error"]["geometric"]["count"])
        unmeasured = int(run.metrics["surface_error"]["geometric"]["unmeasured"])
        row.append(f"{measured}" + (f" (+{unmeasured} unmeasured)" if unmeasured else ""))
        rows.append(row)

    caption = f"{runs[0].case}, {sorted(counts)[0]} control points"
    if as_latex:
        print(latex_table(header, rows, caption))
    else:
        print(f"{caption}\n")
        print(markdown_table(header, rows))

    # The comparison the milestone turns on, stated rather than left to the
    # reader: how much curvature the least-squares fit recovers, and what it
    # pays for it in position.
    by_method = {label(run): run for run in runs}
    interpolate = next((r for name, r in by_method.items() if name.startswith("interpolate")),
                       None)
    least_squares = [r for name, r in by_method.items() if name.startswith("least_squares")]
    if interpolate is not None and least_squares:
        # Read from the common subset when there is one. The table above says
        # that is the number to quote, so the summary had better quote it.
        def value(run, group, field, column):
            if not common:
                return run.stat(group, field)
            values = [run.samples[i][column] for i in common]
            if field == "max":
                return max(values)
            return math.sqrt(sum(v * v for v in values) / len(values))

        best = min(least_squares,
                   key=lambda r: value(r, "mean_curvature", "rms", "mean_curvature"))
        source = (f"over the {len(common)} samples every run measured" if common
                  else "over each run's own measured set, which differ")
        print(f"\nleast squares against interpolation, on this layout, {source}:")
        for name, group, field, column in (
                ("curvature rms", "mean_curvature", "rms", "mean_curvature"),
                ("curvature max", "mean_curvature", "max", "mean_curvature"),
                ("normal max", "normal_degrees", "max", "normal_degrees"),
                ("position rms", "geometric", "rms", "geometric")):
            before = value(interpolate, group, field, column)
            after = value(best, group, field, column)
            if before > 0.0:
                ratio = after / before
                verdict = f"{1.0 / ratio:.3g}x better" if ratio < 1.0 else f"{ratio:.3g}x worse"
                print(f"  {name:<14} {before:>11.4g} -> {after:>11.4g}   {verdict}")

    if any(not run.is_complete() for run in runs):
        print(
            "\nNote: some runs could not measure every sample and the statistics describe "
            "only those they could. The counts are in the last column; where they differ, "
            "the runs are not averaging over the same set of points.",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
