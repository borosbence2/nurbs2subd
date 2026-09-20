#!/usr/bin/env python3
"""Plots PIA convergence against the direct solve.

    python experiments/scripts/plot_convergence.py results/<pia-run-id>

This is the figure for defect 5 of the thesis: progressive iteration and the
direct solve were compared as though they were different methods, when they are
the same square system. The curve to read is `distance to direct solve` -- if
the two were genuinely different methods it would flatten out at some non-zero
value instead of falling to machine precision.

Writes `pia_convergence.png` into the run directory.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from n2s_results import load_run  # noqa: E402


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2

    directory = Path(argv[1])
    run = load_run(directory)

    history_path = directory / "pia_convergence.csv"
    if not history_path.exists():
        print(
            f"{directory} has no pia_convergence.csv. Only a run whose config sets\n"
            f'  "fit": {{ "method": "pia" }}\n'
            f"produces one.",
            file=sys.stderr,
        )
        return 1

    with history_path.open(newline="") as handle:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(handle)]

    iterations = [int(r["iteration"]) for r in rows]
    error = [r["max_error"] for r in rows]
    distance = [r["distance_to_direct"] for r in rows]

    fit = run.metrics.get("fit", {})
    print(f"{run.run_id}: {fit.get('method')} , {len(rows)} iterations")
    print(f"  final interpolation error   {error[-1]:.4e}")
    print(f"  final distance to direct    {distance[-1]:.4e}")
    print(f"  distance at iteration 1     {distance[0]:.4e}")
    if distance[0] > 0:
        print(f"  reduction factor            {distance[0] / max(distance[-1], 1e-300):.3g}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib is not installed; the numbers above are all there is.",
              file=sys.stderr)
        return 0

    figure, axes = plt.subplots(figsize=(6.5, 4.5))
    axes.semilogy(iterations, error, label="interpolation error")
    axes.semilogy(iterations, distance, label="distance to direct solve", linestyle="--")

    # The direct solve is a single solve, so its own residual is a horizontal
    # line. Drawing it makes the point of the figure visible without reading
    # the caption: PIA is walking toward it, not somewhere else.
    direct_residual = float(fit.get("max_interpolation_error", 0.0))
    if fit.get("method") == "pia" and direct_residual > 0:
        axes.axhline(
            direct_residual,
            color="0.4",
            linewidth=0.8,
            label=f"converged residual ({direct_residual:.1e})",
        )

    axes.set_xlabel("PIA iteration")
    axes.set_ylabel("error")
    axes.set_title(f"{run.case}: PIA against the direct solve")
    axes.grid(True, which="both", linewidth=0.3)
    axes.legend()
    figure.tight_layout()

    output = directory / "pia_convergence.png"
    figure.savefig(output, dpi=150)
    print(f"\nwrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
