"""Loading and tabulating nurbs2subd experiment results.

Every table and figure in the write-up is generated from a result directory
through this module, so that a number in the paper can be traced to the run
that produced it and, through that run's ``meta.json``, to the commit.

Usage::

    from n2s_results import load_run, load_sweep, metrics_table

    run = load_run("results/baseline_saddle-20260919-195348")
    print(run.summary_line())
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator

__all__ = [
    "Run",
    "Sweep",
    "load_run",
    "load_sweep",
    "metrics_table",
    "markdown_table",
    "latex_table",
]


@dataclass
class Run:
    """One result directory: its config, metrics and provenance."""

    directory: Path
    config: dict[str, Any]
    metrics: dict[str, Any]
    meta: dict[str, Any]
    _samples: list[dict[str, float]] | None = field(default=None, repr=False)

    @property
    def run_id(self) -> str:
        return self.metrics["run_id"]

    @property
    def case(self) -> str:
        return self.metrics["case"]

    @property
    def notes(self) -> list[str]:
        return self.metrics.get("notes", [])

    def stat(self, group: str, field_name: str = "max") -> float:
        """A single statistic, e.g. ``stat("geometric", "rms")``.

        ``group`` is a key under ``surface_error``, or ``"boundary_error"``;
        ``field_name`` is one of ``max``, ``mean``, ``rms``, ``count``,
        ``unmeasured``.
        """
        if group == "boundary_error":
            return float(self.metrics["boundary_error"][field_name])
        return float(self.metrics["surface_error"][group][field_name])

    def is_complete(self, group: str = "geometric") -> bool:
        """False when some samples could not be measured.

        Worth checking before quoting a number: the statistics describe only
        the samples that *were* measured, and a run with exclusions is
        describing a smaller set than it was asked to.
        """
        if group == "boundary_error":
            return int(self.metrics["boundary_error"]["unmeasured"]) == 0
        return int(self.metrics["surface_error"][group]["unmeasured"]) == 0

    @property
    def samples(self) -> list[dict[str, float]]:
        """Per-sample rows from ``samples.csv``, loaded on first access."""
        if self._samples is None:
            path = self.directory / "samples.csv"
            if not path.exists():
                self._samples = []
            else:
                with path.open(newline="") as handle:
                    self._samples = [
                        {key: float(value) for key, value in row.items()}
                        for row in csv.DictReader(handle)
                    ]
        return self._samples

    def summary_line(self) -> str:
        provenance = self.meta.get("git_describe", "unknown")
        flag = "" if self.is_complete() else "  (incomplete)"
        return (
            f"{self.run_id}  [{provenance}]  "
            f"geometric max {self.stat('geometric'):.4g}  "
            f"hausdorff {self.metrics['surface_error']['hausdorff']:.4g}"
            f"{flag}"
        )


@dataclass
class Sweep:
    """A sweep directory: the index plus the runs it names."""

    directory: Path
    index: dict[str, Any]
    runs: list[Run]
    points: list[dict[str, Any]]

    def __iter__(self) -> Iterator[tuple[dict[str, Any], Run]]:
        return zip(self.points, self.runs)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def load_run(directory: str | Path) -> Run:
    """Loads one result directory."""
    directory = Path(directory)
    if not directory.is_dir():
        raise FileNotFoundError(f"no such result directory: {directory}")

    missing = [
        name
        for name in ("config.json", "metrics.json", "meta.json")
        if not (directory / name).exists()
    ]
    if missing:
        # A result without its provenance is not a result this project accepts.
        raise FileNotFoundError(
            f"{directory} is missing {', '.join(missing)}; it was not written by "
            f"`nurbs2subd run`, or the run did not finish"
        )

    return Run(
        directory=directory,
        config=_read_json(directory / "config.json"),
        metrics=_read_json(directory / "metrics.json"),
        meta=_read_json(directory / "meta.json"),
    )


def load_sweep(directory: str | Path) -> Sweep:
    """Loads a sweep directory and every run it indexes."""
    directory = Path(directory)
    index = _read_json(directory / "sweep.json")

    results_root = directory.parent
    runs: list[Run] = []
    points: list[dict[str, Any]] = []
    for entry in index["runs"]:
        runs.append(load_run(results_root / entry["run_id"]))
        points.append(entry["point"])

    return Sweep(directory=directory, index=index, runs=runs, points=points)


#: The columns quoted by default. Maximum before RMS on purpose: a fitting
#: method is easy to make look good on RMS while leaving one region badly
#: wrong, and the maximum is what a tolerance claim rests on.
DEFAULT_COLUMNS: list[tuple[str, str, str]] = [
    ("geometric max", "geometric", "max"),
    ("geometric rms", "geometric", "rms"),
    ("hausdorff", "hausdorff", ""),
    ("normal max (deg)", "normal_degrees", "max"),
    ("boundary max", "boundary_error", "max"),
]


def metrics_table(
    runs: list[Run],
    columns: list[tuple[str, str, str]] | None = None,
    label: str = "run",
) -> tuple[list[str], list[list[str]]]:
    """Builds a header and rows for the given runs.

    Returns plain strings so that the Markdown and LaTeX writers share exactly
    the same numbers; formatting a table twice is how two versions of a figure
    end up disagreeing.
    """
    columns = columns or DEFAULT_COLUMNS
    header = [label] + [name for name, _, _ in columns]

    rows: list[list[str]] = []
    for run in runs:
        row = [run.run_id]
        for _, group, field_name in columns:
            if group == "hausdorff":
                value = float(run.metrics["surface_error"]["hausdorff"])
            else:
                value = run.stat(group, field_name)
            row.append(f"{value:.4g}")
        rows.append(row)

    return header, rows


def markdown_table(header: list[str], rows: list[list[str]]) -> str:
    widths = [
        max(len(header[i]), *(len(row[i]) for row in rows)) if rows else len(header[i])
        for i in range(len(header))
    ]
    lines = [
        "| " + " | ".join(h.ljust(w) for h, w in zip(header, widths)) + " |",
        "|" + "|".join("-" * (w + 2) for w in widths) + "|",
    ]
    lines += [
        "| " + " | ".join(c.ljust(w) for c, w in zip(row, widths)) + " |" for row in rows
    ]
    return "\n".join(lines)


def latex_table(header: list[str], rows: list[list[str]], caption: str = "") -> str:
    def escape(text: str) -> str:
        return text.replace("_", r"\_")

    spec = "l" + "r" * (len(header) - 1)
    lines = [
        r"\begin{tabular}{" + spec + "}",
        r"\hline",
        " & ".join(escape(h) for h in header) + r" \\",
        r"\hline",
    ]
    lines += [" & ".join(escape(c) for c in row) + r" \\" for row in rows]
    lines += [r"\hline", r"\end{tabular}"]

    table = "\n".join(lines)
    if caption:
        table = (
            "\\begin{table}[h]\n\\centering\n"
            + table
            + f"\n\\caption{{{escape(caption)}}}\n\\end{{table}}"
        )
    return table
