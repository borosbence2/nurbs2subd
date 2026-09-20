#!/usr/bin/env python3
"""Error against control-point count, log-log, one curve per fit method.

    python experiments/scripts/plot_error_vs_controlpoints.py results/<sweep-id>

Reads a sweep over ``fit.layout_refinement`` and, optionally, ``fit.method``.
The x axis is the control-point count each run actually produced, not the
refinement factor that produced it: the factor is an input, the count is what a
reader comparing two methods at equal cost needs.

The slope is what the figure is for. On a log-log plot a method converging at
order ``p`` in the control-point spacing appears as a line of slope ``-p/2``,
since the count grows as the square of the density. The fitted slope is printed
per method, over the last three points, where the asymptotic rate is visible and
the coarse-layout transient is not.

The curvature column is printed but not plotted, and should not be compared
down the column: each row is a different layout, so its samples sit at different
places and a different subset of them survives the projection. Use
``compare_fits.py`` for curvature, which restricts to the samples every run
measured -- that is only possible between runs that share a layout.

Writes ``error_vs_controlpoints.png`` into the sweep directory.
"""

from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from n2s_results import load_sweep  # noqa: E402

#: How many of the finest points the slope is fitted over.
SLOPE_POINTS = 3

#: Statistics drawn. Maximum as well as RMS on purpose: a method is easy to make
#: look good on RMS while leaving one region badly wrong, and the maximum is
#: what a tolerance claim rests on.
SERIES = [
    ("geometric", "rms", "-", "o"),
    ("geometric", "max", "--", "s"),
]


def fitted_slope(xs: list[float], ys: list[float]) -> float | None:
    """Least-squares slope of log(y) against log(x), or None if undefined."""
    pairs = [(x, y) for x, y in zip(xs, ys) if x > 0.0 and y > 0.0]
    if len(pairs) < 2:
        return None

    lx = [math.log(x) for x, _ in pairs]
    ly = [math.log(y) for _, y in pairs]
    n = float(len(pairs))
    mean_x = sum(lx) / n
    mean_y = sum(ly) / n
    denominator = sum((x - mean_x) ** 2 for x in lx)
    if denominator <= 0.0:
        return None
    return sum((x - mean_x) * (y - mean_y) for x, y in zip(lx, ly)) / denominator


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2

    directory = Path(argv[1])
    sweep = load_sweep(directory)

    by_method: dict[str, list[dict[str, float]]] = defaultdict(list)
    for point, run in sweep:
        method = str(point.get("fit.method", run.metrics["fit"]["method"]))
        by_method[method].append(
            {
                "control_vertices": float(run.metrics["control_vertices"]),
                "geometric_rms": run.stat("geometric", "rms"),
                "geometric_max": run.stat("geometric", "max"),
                "mean_curvature_rms": run.stat("mean_curvature", "rms"),
                "complete": run.is_complete(),
            }
        )

    for rows in by_method.values():
        rows.sort(key=lambda r: r["control_vertices"])

    print(f"{sweep.index['sweep_id']}: {len(sweep.runs)} runs, "
          f"{len(by_method)} method(s)")

    slopes: dict[str, dict[str, float | None]] = {}
    for method, rows in sorted(by_method.items()):
        print(f"\n{method}")
        print(f"{'control pts':>12} {'geom rms':>11} {'geom max':>11} {'curv rms':>11}")
        for row in rows:
            flag = "" if row["complete"] else "  (incomplete)"
            print(f"{row['control_vertices']:>12.0f} {row['geometric_rms']:>11.4g} "
                  f"{row['geometric_max']:>11.4g} {row['mean_curvature_rms']:>11.4g}{flag}")

        # Curvature is printed for context but deliberately not plotted, and it
        # must not be read down the column: each row is a different layout, so
        # its error samples sit at different places, and the projection behind
        # the geometric error keeps a different subset of them in each run. Two
        # rows are therefore averages over two unrelated populations. The
        # position columns carry the same caveat but are dominated by the
        # convergence, which is what this figure is for.
        tail = rows[-SLOPE_POINTS:]
        slopes[method] = {
            "rms": fitted_slope([r["control_vertices"] for r in tail],
                                [r["geometric_rms"] for r in tail]),
            "max": fitted_slope([r["control_vertices"] for r in tail],
                                [r["geometric_max"] for r in tail]),
        }
        for name, slope in slopes[method].items():
            if slope is None:
                print(f"  slope ({name}): not enough usable points")
            else:
                # -p/2 in the count is -p in the spacing, which is the number
                # anyone comparing against a published convergence order wants.
                print(f"  slope ({name}): {slope:.3f} in the count, "
                      f"so order {-2.0 * slope:.2f} in the spacing")

    if len(by_method) > 1:
        print("\nAt equal control-point count:")
        methods = sorted(by_method)
        common = sorted(
            {r["control_vertices"] for r in by_method[methods[0]]}.intersection(
                *({r["control_vertices"] for r in by_method[m]} for m in methods[1:])
            )
        )
        for count in common:
            parts = []
            for method in methods:
                row = next(r for r in by_method[method] if r["control_vertices"] == count)
                parts.append(f"{method} {row['geometric_rms']:.4g}")
            print(f"  {count:.0f} control points: " + ", ".join(parts))

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib is not installed; the table above is all there is.",
              file=sys.stderr)
        return 0

    figure, axes = plt.subplots(figsize=(7.0, 5.0))
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for index, (method, rows) in enumerate(sorted(by_method.items())):
        colour = colours[index % len(colours)]
        for group, field, style, marker in SERIES:
            key = f"{group}_{field}"
            axes.loglog([r["control_vertices"] for r in rows],
                        [r[key] for r in rows],
                        linestyle=style,
                        marker=marker,
                        color=colour,
                        label=f"{method}, {group} {field}")

    axes.set_xlabel("control points")
    axes.set_ylabel("geometric error")
    axes.grid(True, which="both", linewidth=0.3)
    axes.legend(fontsize=8)

    subtitle = "; ".join(
        f"{method} rms slope {s['rms']:.2f}" for method, s in sorted(slopes.items())
        if s["rms"] is not None
    )
    axes.set_title(f"{sweep.runs[0].case}: error against control-point count\n{subtitle}",
                   fontsize=10)

    figure.tight_layout()
    output = directory / "error_vs_controlpoints.png"
    figure.savefig(output, dpi=150)
    print(f"\nwrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
