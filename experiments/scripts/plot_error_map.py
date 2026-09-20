#!/usr/bin/env python3
"""Draws a scalar field over the trimmed domain, at a fixed colour scale.

    python experiments/scripts/plot_error_map.py results/<run-id> geometric
    python experiments/scripts/plot_error_map.py results/<run-id> geometric --range 0 0.01

The range comes from the case file's ``color_ranges`` or from ``--range``. If
neither supplies one the script **refuses to draw**.

That is deliberate, and it is the same rule the viewer enforces at its
screenshot button. Defect 6 of the thesis was curvature plots with no fixed
colour scale: every figure auto-ranged to its own data, so two of them placed
side by side compared nothing, while looking exactly as though they did. A
figure that cannot be compared is worse than no figure, because it is
persuasive.

For the same reason the curvature fields read ``mean_curvature_deviation``
rather than the viewer's ``mean_curvature``. They are different quantities on
different scales: the viewer colours the curvature of a surface, which is
signed, while these columns are ``|limit - nurbs|``, which is not. Sharing one
key would put two incomparable figures on what looks like a common scale.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from n2s_results import load_run  # noqa: E402

#: Columns in samples.csv that can be drawn. For each: the ``color_ranges`` key
#: in the case file, and the column saying whether that sample holds a real
#: measurement. Curvature has its own gate because it needs second derivatives
#: on both surfaces and can fail where position and normal are perfectly fine.
FIELDS = {
    "geometric": ("error", "measured"),
    "parametric": ("error", "measured"),
    "normal_degrees": ("normal_degrees", "measured"),
    "mean_curvature": ("mean_curvature_deviation", "curvature_measured"),
    "gaussian_curvature": ("gaussian_curvature_deviation", "curvature_measured"),
}


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print(__doc__)
        return 2

    directory = Path(argv[1])
    field = argv[2]

    explicit_range = None
    if "--range" in argv:
        index = argv.index("--range")
        explicit_range = (float(argv[index + 1]), float(argv[index + 2]))

    if field not in FIELDS:
        print(f"unknown field {field!r}; expected one of {', '.join(FIELDS)}", file=sys.stderr)
        return 2

    range_key, gate = FIELDS[field]

    run = load_run(directory)
    samples_path = directory / "samples.csv"
    if not samples_path.exists():
        print(f"{directory} has no samples.csv; set \"write_samples\": true in the config",
              file=sys.stderr)
        return 1

    with samples_path.open(newline="") as handle:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(handle)]

    if rows and gate not in rows[0]:
        # A samples.csv written before the curvature columns existed. Saying so
        # beats drawing nothing and beats a KeyError.
        print(f"{samples_path} has no {gate!r} column, so it predates this field; "
              f"re-run the experiment to get it.", file=sys.stderr)
        return 1

    # Unmeasured samples are dropped rather than drawn as zero. A failed
    # projection plotted as zero error is a picture of a perfect fit in exactly
    # the region where the measurement failed.
    measured = [r for r in rows if r[gate] > 0.5]
    dropped = len(rows) - len(measured)
    if not measured:
        print("every sample in this run was unmeasured; there is nothing to draw",
              file=sys.stderr)
        return 1

    # Resolve the colour range.
    value_range = explicit_range
    source = "--range"
    if value_range is None:
        # The config records the case path as it was resolved when the run
        # started, which is relative to the working directory of that run --
        # not to the result directory. Try it as given first, then beside the
        # result, so the script works from either.
        candidates = [Path(run.config["case"])]
        candidates.append(directory / run.config["case"])
        case_path = next((c for c in candidates if c.exists()), candidates[0])
        if case_path.exists():
            import json

            with case_path.open() as handle:
                ranges = json.load(handle).get("color_ranges", {})
            if range_key in ranges:
                value_range = (float(ranges[range_key][0]), float(ranges[range_key][1]))
                source = f"case file ({range_key})"

    if value_range is None:
        values = [r[field] for r in measured]
        print(
            f"No colour range for {field!r}.\n"
            f"  The case file has no \"color_ranges\" entry for {range_key!r} and no\n"
            f"  --range was given, so this figure would auto-range to its own data and\n"
            f"  could not be compared with any other. Refusing to draw it.\n\n"
            f"  This run's values span [{min(values):.4g}, {max(values):.4g}]; pick a range\n"
            f"  covering every run you mean to compare, then either add it to the case file\n"
            f"  or pass --range LOW HIGH.",
            file=sys.stderr,
        )
        return 1

    print(f"{run.run_id}: {field} over {len(measured)} samples")
    print(f"  colour range [{value_range[0]:.4g}, {value_range[1]:.4g}] from {source}")
    if dropped:
        print(f"  {dropped} unmeasured samples dropped")

    over = sum(1 for r in measured if r[field] > value_range[1])
    if over:
        # Clipping is normal, but silent clipping hides the worst of the error
        # behind the top colour, which is the half of the picture that matters.
        print(f"  {over} samples exceed the top of the range and are clipped")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib is not installed; the summary above is all there is.",
              file=sys.stderr)
        return 0

    figure, axes = plt.subplots(figsize=(6.0, 5.0))
    scatter = axes.scatter(
        [r["domain_u"] for r in measured],
        [r["domain_v"] for r in measured],
        c=[r[field] for r in measured],
        vmin=value_range[0],
        vmax=value_range[1],
        s=8,
        cmap="viridis",
    )

    # The colourbar is not optional: a scalar figure without its scale shown is
    # the defect this whole script exists to prevent.
    # Labelled by the range key, not the column name, so that a curvature
    # figure says on its face that it shows a deviation rather than a curvature.
    bar = figure.colorbar(scatter, ax=axes)
    bar.set_label(f"{range_key}  [{value_range[0]:.4g}, {value_range[1]:.4g}]")

    axes.set_xlabel("domain u")
    axes.set_ylabel("domain v")
    axes.set_aspect("equal")
    axes.set_title(f"{run.case}\n{field}, {run.metrics.get('fit', {}).get('method', '?')} fit")
    figure.tight_layout()

    output = directory / f"map_{field}.png"
    figure.savefig(output, dpi=150)
    print(f"\nwrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
