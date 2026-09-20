#!/usr/bin/env python3
"""Plots the fairness trade-off: what raising lambda costs and what it buys.

    python experiments/scripts/plot_lambda.py results/<sweep-id>

Reads a sweep whose axis is ``fit.least_squares.lambda``. Two panels sharing an
x axis, because the point is the trade-off and not either curve alone:

  * position error against the NURBS (what lambda costs)
  * curvature deviation and normal deviation (what lambda buys)

**Why this figure and not a single "best lambda".** There is no lambda that
minimises everything; there is a knee, past which curvature stops improving and
position error starts paying for it. A script that printed one number would be
hiding the choice rather than making it. The knee is reported as the largest
lambda whose position RMS is still within a stated fraction of the unfaired
fit -- that fraction is a judgement, so it is printed alongside the answer.

Writes ``lambda_tradeoff.png`` into the sweep directory.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from n2s_results import load_sweep  # noqa: E402

#: How much position RMS the knee is allowed to give up against lambda = 0.
#: A judgement, printed with the answer rather than buried in it.
POSITION_BUDGET = 1.05

AXIS = "fit.least_squares.lambda"


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2

    directory = Path(argv[1])
    sweep = load_sweep(directory)

    rows = []
    for point, run in sweep:
        if AXIS not in point:
            print(f"{directory} is not a lambda sweep; its axes are {sorted(point)}",
                  file=sys.stderr)
            return 1
        rows.append(
            {
                "lambda": float(point[AXIS]),
                "geometric_rms": run.stat("geometric", "rms"),
                "geometric_max": run.stat("geometric", "max"),
                "mean_curvature_rms": run.stat("mean_curvature", "rms"),
                "mean_curvature_max": run.stat("mean_curvature", "max"),
                "normal_max": run.stat("normal_degrees", "max"),
                "boundary_max": run.stat("boundary_error", "max"),
                "energy": float(run.metrics["fit"]["least_squares"]["fairness_energy"]),
                "complete": run.is_complete(),
            }
        )

    rows.sort(key=lambda r: r["lambda"])

    # The boundary is constrained to R1's values at every point, so it must not
    # move across the sweep. If it does, the constraint is not being applied and
    # every comparison in this figure is between two things at once.
    boundaries = {round(r["boundary_max"], 12) for r in rows}
    if len(boundaries) > 1:
        print(
            "WARNING: the boundary error moves across this sweep "
            f"({min(boundaries):.6g} to {max(boundaries):.6g}). It is supposed to be "
            "constrained to R1's solve and identical at every lambda, so something is "
            "varying that should not be.",
            file=sys.stderr,
        )

    baseline = rows[0]["geometric_rms"]
    knee = None
    for row in rows:
        if row["geometric_rms"] <= baseline * POSITION_BUDGET:
            knee = row

    print(f"{sweep.index['sweep_id']}: {len(rows)} points")
    print(f"{'lambda':>10} {'geom rms':>11} {'curv rms':>11} {'curv max':>11} "
          f"{'normal max':>11} {'energy':>11}")
    for row in rows:
        flag = "" if row["complete"] else "  (incomplete)"
        print(f"{row['lambda']:>10.3g} {row['geometric_rms']:>11.4g} "
              f"{row['mean_curvature_rms']:>11.4g} {row['mean_curvature_max']:>11.4g} "
              f"{row['normal_max']:>11.4g} {row['energy']:>11.4g}{flag}")

    if knee is not None:
        print(
            f"\nknee: lambda = {knee['lambda']:.3g}, the largest whose position RMS stays "
            f"within {POSITION_BUDGET:.2f}x of the unfaired fit."
        )
        print(f"  position RMS   {baseline:.4g} -> {knee['geometric_rms']:.4g}")
        print(f"  curvature RMS  {rows[0]['mean_curvature_rms']:.4g} -> "
              f"{knee['mean_curvature_rms']:.4g}")
        print(f"  normal max     {rows[0]['normal_max']:.4g} -> {knee['normal_max']:.4g}")

    if any(not row["complete"] for row in rows):
        print(
            "\nNote: some runs could not measure every sample, and the statistics above "
            "describe only those they could. The sample counts differ between points, so "
            "small differences here are not all signal.",
            file=sys.stderr,
        )

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib is not installed; the table above is all there is.",
              file=sys.stderr)
        return 0

    # Lambda spans zero, which a log axis cannot show. The zero point is drawn
    # at the left edge as a separate marker and labelled, rather than dropped:
    # it is the baseline every other point is judged against.
    positive = [r for r in rows if r["lambda"] > 0.0]
    zero = [r for r in rows if r["lambda"] == 0.0]

    figure, (top, bottom) = plt.subplots(2, 1, figsize=(7.0, 7.0), sharex=True)

    top.semilogx([r["lambda"] for r in positive], [r["geometric_rms"] for r in positive],
                 marker="o", label="geometric rms")
    top.semilogx([r["lambda"] for r in positive], [r["geometric_max"] for r in positive],
                 marker="s", linestyle="--", label="geometric max")
    for r in zero:
        top.axhline(r["geometric_rms"], color="0.6", linewidth=0.8,
                    label=f"lambda = 0 (rms {r['geometric_rms']:.3g})")
    top.set_ylabel("position error")
    top.set_yscale("log")
    top.grid(True, which="both", linewidth=0.3)
    top.legend(fontsize=8)
    top.set_title(f"{sweep.runs[0].case}: what the fairness weight costs and buys")

    bottom.loglog([r["lambda"] for r in positive],
                  [r["mean_curvature_rms"] for r in positive],
                  marker="o", label="mean curvature deviation, rms")
    bottom.loglog([r["lambda"] for r in positive],
                  [r["normal_max"] for r in positive],
                  marker="^", linestyle="--", label="normal deviation, max (deg)")
    for r in zero:
        bottom.axhline(r["mean_curvature_rms"], color="0.6", linewidth=0.8,
                       label=f"lambda = 0 (curv rms {r['mean_curvature_rms']:.3g})")
    if knee is not None and knee["lambda"] > 0.0:
        bottom.axvline(knee["lambda"], color="0.3", linewidth=0.8, linestyle=":")
        top.axvline(knee["lambda"], color="0.3", linewidth=0.8, linestyle=":")

    bottom.set_xlabel("lambda")
    bottom.set_ylabel("quality")
    bottom.grid(True, which="both", linewidth=0.3)
    bottom.legend(fontsize=8)

    figure.tight_layout()
    output = directory / "lambda_tradeoff.png"
    figure.savefig(output, dpi=150)
    print(f"\nwrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
