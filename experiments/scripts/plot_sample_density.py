#!/usr/bin/env python3
"""Error and fairness energy against sample density.

    python experiments/scripts/plot_sample_density.py results/<sweep-id>

Reads a sweep over ``fit.least_squares.samples_per_side``. Two questions, and
the figure exists to answer the second as much as the first:

  1. **Has the fit converged in the sample count?** The error should flatten.
     A fit still moving at the densest point was built from too few samples,
     and every number taken from it is a number about the sampling.

  2. **Is ``normalise_lambda`` doing its job?** It divides the data term by the
     sample count and the fairness term by its row count so that lambda means
     the same thing at every density. If it works, the fairness energy is flat
     across this sweep. If it drifts, then sweeping sample density silently
     sweeps the fairness weight too, and the lambda sweep beside it is
     measuring a mixture of the two.

Writes ``sample_density.png`` into the sweep directory.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from n2s_results import load_sweep  # noqa: E402

AXIS = "fit.least_squares.samples_per_side"

#: Fractional spread in the fairness energy above which normalisation is
#: reported as not holding. Generous: the sample set changes with density, so
#: some movement is expected; an order of magnitude is not.
ENERGY_DRIFT_LIMIT = 0.25


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2

    directory = Path(argv[1])
    sweep = load_sweep(directory)

    rows = []
    for point, run in sweep:
        if AXIS not in point:
            print(f"{directory} is not a sample-density sweep; its axes are {sorted(point)}",
                  file=sys.stderr)
            return 1
        least_squares = run.metrics["fit"]["least_squares"]
        rows.append(
            {
                "per_side": int(point[AXIS]),
                "used": int(least_squares["samples_used"]),
                "free": int(least_squares["free_vertices"]),
                "geometric_rms": run.stat("geometric", "rms"),
                "geometric_max": run.stat("geometric", "max"),
                "mean_curvature_rms": run.stat("mean_curvature", "rms"),
                "energy": float(least_squares["fairness_energy"]),
                "condition": float(least_squares["condition_estimate"]),
            }
        )

    rows.sort(key=lambda r: r["per_side"])

    print(f"{sweep.index['sweep_id']}: {len(rows)} points")
    print(f"{'per side':>9} {'used':>7} {'per free':>9} {'geom rms':>11} "
          f"{'curv rms':>11} {'energy':>11} {'condition':>11}")
    for row in rows:
        per_free = row["used"] / row["free"] if row["free"] else 0.0
        print(f"{row['per_side']:>9} {row['used']:>7} {per_free:>9.1f} "
              f"{row['geometric_rms']:>11.4g} {row['mean_curvature_rms']:>11.4g} "
              f"{row['energy']:>11.4g} {row['condition']:>11.4g}")

    # Under-determined runs are named rather than merely plotted. Fewer samples
    # than free control points means the fairness term alone is holding the
    # system down, and the fit is not the one the config asked for.
    starved = [r for r in rows if r["used"] < r["free"]]
    if starved:
        print(
            "\nWARNING: "
            + ", ".join(f"{r['per_side']} per side ({r['used']} samples < {r['free']} free "
                        f"control points)" for r in starved)
            + " -- under-determined by the data term alone; only the fairness term makes "
              "these solvable, so they are not measuring what the others measure.",
            file=sys.stderr,
        )

    energies = [r["energy"] for r in rows]
    if energies and min(energies) > 0.0:
        drift = (max(energies) - min(energies)) / min(energies)
        verdict = "holds" if drift <= ENERGY_DRIFT_LIMIT else "does NOT hold"
        print(f"\nfairness energy spread across the sweep: {drift:.1%} -- "
              f"normalise_lambda {verdict}")
        if drift > ENERGY_DRIFT_LIMIT:
            print(
                "  Sweeping sample density is therefore also sweeping the effective "
                "fairness weight, and the lambda sweep beside this one is measuring a "
                "mixture of the two.",
                file=sys.stderr,
            )

    finest = rows[-1]
    previous = rows[-2] if len(rows) > 1 else None
    if previous is not None and finest["geometric_rms"] > 0.0:
        change = abs(finest["geometric_rms"] - previous["geometric_rms"]) / finest["geometric_rms"]
        print(f"position RMS change over the last doubling: {change:.2%}"
              + ("  (converged)" if change < 0.02 else "  (still moving)"))

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib is not installed; the table above is all there is.",
              file=sys.stderr)
        return 0

    figure, (top, bottom) = plt.subplots(2, 1, figsize=(7.0, 7.0), sharex=True)

    top.loglog([r["used"] for r in rows], [r["geometric_rms"] for r in rows],
               marker="o", label="geometric rms")
    top.loglog([r["used"] for r in rows], [r["geometric_max"] for r in rows],
               marker="s", linestyle="--", label="geometric max")
    top.loglog([r["used"] for r in rows], [r["mean_curvature_rms"] for r in rows],
               marker="^", linestyle=":", label="mean curvature deviation, rms")
    top.set_ylabel("error")
    top.grid(True, which="both", linewidth=0.3)
    top.legend(fontsize=8)
    top.set_title(f"{sweep.runs[0].case}: convergence in the sample count", fontsize=10)

    bottom.semilogx([r["used"] for r in rows], [r["energy"] for r in rows],
                    marker="o", label="fairness energy")
    # Drawn against a flat reference, because flat is the claim being checked
    # and a curve without one just looks like data.
    bottom.axhline(rows[0]["energy"], color="0.6", linewidth=0.8,
                   label=f"energy at the coarsest ({rows[0]['energy']:.3g})")
    bottom.set_xlabel("samples used")
    bottom.set_ylabel("fairness energy")
    bottom.grid(True, which="both", linewidth=0.3)
    bottom.legend(fontsize=8)

    figure.tight_layout()
    output = directory / "sample_density.png"
    figure.savefig(output, dpi=150)
    print(f"\nwrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
