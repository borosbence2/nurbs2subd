#!/usr/bin/env python3
"""Plots error against layout size for a sweep.

    python experiments/scripts/plot_sweep.py results/sweep_layout_density-20260919-195519

Writes ``error_vs_control_points.png`` into the sweep directory and prints the
same numbers as a table, so the figure and the text cannot disagree.

Colour and axis ranges are fixed by the data rather than chosen per figure, and
both axes are logarithmic: the quantity of interest is the *rate* at which
error falls with control point count, which a linear axis hides.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from n2s_results import load_sweep, markdown_table, metrics_table  # noqa: E402


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2

    sweep_directory = Path(argv[1])
    sweep = load_sweep(sweep_directory)

    header, rows = metrics_table(sweep.runs)
    header = ["point"] + header[1:]
    for point, row in zip(sweep.points, rows):
        row[0] = " ".join(f"{key}={value}" for key, value in sorted(point.items()))
    print(markdown_table(header, rows))

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "\nmatplotlib is not installed, so the table above is all there is.\n"
            "Install it to get the figure: pip install matplotlib",
            file=sys.stderr,
        )
        return 0

    control_points = [run.metrics["control_vertices"] for run in sweep.runs]
    geometric_max = [run.stat("geometric", "max") for run in sweep.runs]
    geometric_rms = [run.stat("geometric", "rms") for run in sweep.runs]

    # Square layouts are separated out. Mixing them with anisotropic ones on a
    # single curve is misleading: on a saddle the two directions' smoothing
    # errors cancel when the spacings match, so an anisotropic layout can be
    # worse than a coarser square one despite having more control points.
    def is_square(point: dict[str, object]) -> bool:
        rows_value = point.get("layout_rows")
        columns_value = point.get("layout_columns")
        return rows_value is not None and rows_value == columns_value

    square = [i for i, point in enumerate(sweep.points) if is_square(point)]
    other = [i for i in range(len(sweep.points)) if i not in square]

    figure, axes = plt.subplots(figsize=(6.5, 4.5))
    if other:
        axes.scatter(
            [control_points[i] for i in other],
            [geometric_max[i] for i in other],
            marker="x",
            label="max, anisotropic layout",
        )
    if square:
        axes.plot(
            [control_points[i] for i in square],
            [geometric_max[i] for i in square],
            marker="o",
            label="max, square layout",
        )
        axes.plot(
            [control_points[i] for i in square],
            [geometric_rms[i] for i in square],
            marker="s",
            linestyle="--",
            label="rms, square layout",
        )

    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlabel("control points")
    axes.set_ylabel("geometric error to the NURBS")
    axes.set_title(sweep.index["sweep_id"])
    axes.grid(True, which="both", linewidth=0.3)
    axes.legend()
    figure.tight_layout()

    output = sweep_directory / "error_vs_control_points.png"
    figure.savefig(output, dpi=150)
    print(f"\nwrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
