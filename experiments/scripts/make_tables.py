#!/usr/bin/env python3
"""Generates the Markdown and LaTeX tables for a set of runs.

    python experiments/scripts/make_tables.py results/run-a results/run-b ...
    python experiments/scripts/make_tables.py --sweep results/sweep-xyz

Writes ``tables.md`` and ``tables.tex`` beside the first input, and prints the
Markdown. Both formats come from one set of formatted values, so the paper and
the notes cannot end up quoting different numbers for the same run.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from n2s_results import (  # noqa: E402
    latex_table,
    load_run,
    load_sweep,
    markdown_table,
    metrics_table,
)


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2

    if argv[1] == "--sweep":
        if len(argv) != 3:
            print(__doc__)
            return 2
        sweep = load_sweep(argv[2])
        runs = sweep.runs
        output_directory = Path(argv[2])
        caption = f"Sweep {sweep.index['sweep_id']}"
    else:
        runs = [load_run(path) for path in argv[1:]]
        output_directory = Path(argv[1])
        caption = "Run comparison"

    header, rows = metrics_table(runs)

    markdown = markdown_table(header, rows)
    print(markdown)

    # Provenance under the table, not in a separate file nobody opens. A table
    # that cannot be traced to the commit that produced it is not usable in a
    # write-up, and the commit is the thing most easily lost.
    provenance = ["", "Provenance:", ""]
    incomplete = []
    for run in runs:
        provenance.append(
            f"- `{run.run_id}` from `{run.meta.get('git_describe', 'unknown')}` "
            f"at {run.meta.get('timestamp_utc', 'unknown')}"
        )
        if not run.is_complete():
            incomplete.append(run.run_id)

    if incomplete:
        provenance += [
            "",
            "> Some samples could not be measured in: "
            + ", ".join(f"`{name}`" for name in incomplete)
            + ". The statistics above describe only the samples that were measured.",
        ]

    (output_directory / "tables.md").write_text(
        markdown + "\n" + "\n".join(provenance) + "\n", encoding="utf-8"
    )
    (output_directory / "tables.tex").write_text(
        latex_table(header, rows, caption) + "\n", encoding="utf-8"
    )

    print("\n".join(provenance))
    print(f"\nwrote {output_directory / 'tables.md'} and {output_directory / 'tables.tex'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
