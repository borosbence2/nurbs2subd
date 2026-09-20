#!/usr/bin/env python3
"""Converts the MSc thesis test data into nurbs2subd case files.

    python experiments/scripts/import_thesis_cases.py "<assets dir>" data

Run once; the JSON it writes is committed and is what the project actually
uses. Keeping this as a script rather than as C++ was a deliberate decision:
the legacy formats are read exactly once, and a parser for them in ``core/``
would have to be carried, tested and trusted forever for four input files.

WHAT IS RECONSTRUCTED, AND WHY
------------------------------
The thesis files are incomplete and were partly hand-written, so some of this
is reconstruction rather than transcription. Every such decision is listed
here, because a baseline built on undocumented guesses is exactly what this
project exists to replace.

1. *Surfaces.* Each ``surface_*.txt`` gives four boundary curves of five
   control points each -- not a control net. The interior is filled by the
   discrete Coons construction (a bilinear blend of the boundary control
   points, minus the bilinear blend of the corners), giving a 5x5 net read as
   a bicubic B-spline with clamped knots. Coons is the natural reading of
   "here are the four edges", it reproduces the boundaries exactly, and it is
   the only interior determined by the boundaries alone.

2. *The DoubleVB domain is rescaled.* Its trim loop and layout span roughly
   450 x 307 units, while a surface's parametric domain is [0,1]^2. The
   mapping is a single uniform scale plus a translation, so the aspect ratio
   is preserved and the domain sits strictly inside [0,1]^2 with a small
   margin. A per-axis fit to [0,1]^2 would have squashed it.

3. *DoubleVB is paired with the folder-3 surface.* The original pairing is not
   recorded anywhere in the assets, and no surface accompanies the DoubleVB
   files. The two live in the same folder, which is the only evidence there
   is, and the pairing is stated in the case name.

4. *Four trim joins are snapped.* Eighteen of the twenty-two joins in the
   DoubleVB loop are bitwise exact; four are off by 0.4 to 0.9 units, which is
   0.1 to 0.2 percent of the model. Those four endpoints are snapped to the
   layout's boundary vertices, the layout being the authority: the loop is the
   layout's boundary, and every other endpoint already coincides with one.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

# --------------------------------------------------------------------------
# Readers for the legacy formats.
# --------------------------------------------------------------------------


def read_bsc(path: Path) -> list[dict[str, Any]]:
    """Reads the thesis curve format.

    Per curve: a degree, then a count followed by that many knots, then a
    count followed by that many ``x y z`` control points.
    """
    lines = path.read_text().split("\n")
    curves: list[dict[str, Any]] = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        degree = int(lines[i].strip())
        tokens = lines[i + 1].split()
        knot_count = int(tokens[0])
        knots = [float(t) for t in tokens[1 : 1 + knot_count]]
        point_count = int(lines[i + 2].strip())
        points = [
            tuple(float(v) for v in lines[i + 3 + k].split()) for k in range(point_count)
        ]

        implied = knot_count - degree - 1
        if implied != point_count:
            raise ValueError(
                f"{path.name} curve {len(curves)}: {knot_count} knots at degree {degree} "
                f"imply {implied} control points but the file says {point_count}"
            )

        curves.append({"degree": degree, "knots": knots, "points": points})
        i += 3 + point_count
    return curves


def read_bezier_trim(path: Path) -> list[dict[str, Any]]:
    """Reads the ``x0,y0,x1,y1,x2,y2,x3,y3`` per line form: cubic Bezier segments."""
    curves = []
    for line in path.read_text().split("\n"):
        values = [float(v) for v in line.split(",") if v.strip()]
        if not values:
            continue
        if len(values) != 8:
            raise ValueError(f"{path.name}: expected 8 numbers per line, got {len(values)}")
        points = [(values[2 * k], values[2 * k + 1], 0.0) for k in range(4)]
        curves.append({"degree": 3, "knots": [0, 0, 0, 0, 1, 1, 1, 1], "points": points})
    return curves


def read_obj(path: Path) -> tuple[list[tuple[float, float, float]], list[list[int]]]:
    vertices, faces = [], []
    for line in path.read_text().split("\n"):
        if line.startswith("v "):
            vertices.append(tuple(float(t) for t in line.split()[1:4]))
        elif line.startswith("f "):
            faces.append([int(t.split("/")[0]) - 1 for t in line.split()[1:]])
    return vertices, faces


def read_boundary_surface(path: Path) -> list[list[tuple[float, float, float]]]:
    """Reads four boundary curves, five control points each.

    The four lines are *not* in a consistent order or direction between files,
    so each is identified by which coordinate it holds constant and then
    oriented by the one that varies. Trusting the file order silently
    transposes two of the four surfaces.
    """
    rows = []
    for line in path.read_text().split("\n"):
        values = [float(v) for v in line.split(",") if v.strip()]
        if values:
            rows.append([tuple(values[3 * k : 3 * k + 3]) for k in range(len(values) // 3)])
    if len(rows) != 4:
        raise ValueError(f"{path.name}: expected 4 boundary curves, got {len(rows)}")

    x_values = [v[0] for row in rows for v in row]
    y_values = [v[1] for row in rows for v in row]
    x_lo, x_hi = min(x_values), max(x_values)
    y_lo, y_hi = min(y_values), max(y_values)

    found: dict[str, list[tuple[float, float, float]]] = {}
    for row in rows:
        xs = {round(v[0], 9) for v in row}
        ys = {round(v[1], 9) for v in row}
        if len(xs) == 1:
            key = "u0" if abs(row[0][0] - x_lo) < 1e-9 else "u1"
            found[key] = sorted(row, key=lambda v: v[1])
        elif len(ys) == 1:
            key = "v0" if abs(row[0][1] - y_lo) < 1e-9 else "v1"
            found[key] = sorted(row, key=lambda v: v[0])
        else:
            raise ValueError(f"{path.name}: a boundary curve holds neither x nor y constant")

    missing = {"u0", "u1", "v0", "v1"} - set(found)
    if missing:
        raise ValueError(f"{path.name}: missing boundaries {sorted(missing)}")

    _ = (x_hi, y_hi)
    return [found["u0"], found["u1"], found["v0"], found["v1"]]


# --------------------------------------------------------------------------
# Reconstruction.
# --------------------------------------------------------------------------


def coons_control_net(boundaries: list[list[tuple[float, float, float]]]) -> list[list[Any]]:
    """Fills a 5x5 control net from its four boundary control polygons.

    The discrete Coons formula: each interior point is the bilinear blend of
    the two opposing boundary points, minus the bilinear blend of the four
    corners that the first term counts twice.
    """
    u0, u1, v0, v1 = boundaries
    n = len(u0)
    if any(len(b) != n for b in boundaries):
        raise ValueError("the four boundaries must have the same control point count")

    net = [[None] * n for _ in range(n)]
    for j in range(n):
        net[0][j] = u0[j]
        net[n - 1][j] = u1[j]
    for i in range(n):
        net[i][0] = v0[i]
        net[i][n - 1] = v1[i]

    for i in range(n):
        for j in range(n):
            if net[i][j] is not None:
                continue
            s = i / (n - 1)
            t = j / (n - 1)
            point = []
            for k in range(3):
                edges = (
                    (1 - s) * net[0][j][k]
                    + s * net[n - 1][j][k]
                    + (1 - t) * net[i][0][k]
                    + t * net[i][n - 1][k]
                )
                corners = (
                    (1 - s) * (1 - t) * net[0][0][k]
                    + (1 - s) * t * net[0][n - 1][k]
                    + s * (1 - t) * net[n - 1][0][k]
                    + s * t * net[n - 1][n - 1][k]
                )
                point.append(edges - corners)
            net[i][j] = tuple(point)

    return net


def clamped_cubic_knots(count: int) -> list[float]:
    interior = count - 4
    if interior < 0:
        raise ValueError(f"a cubic B-spline needs at least 4 control points, got {count}")
    return [0.0] * 4 + [(k + 1) / (interior + 1) for k in range(interior)] + [1.0] * 4


def surface_json(net: list[list[Any]]) -> dict[str, Any]:
    rows, columns = len(net), len(net[0])
    return {
        "format": "n2s-surface",
        "version": 1,
        "degree_u": 3,
        "degree_v": 3,
        "knots_u": clamped_cubic_knots(rows),
        "knots_v": clamped_cubic_knots(columns),
        "num_u": rows,
        "num_v": columns,
        # Row major with u as the slow index, matching NurbsSurface.
        "control_points": [list(net[i][j]) for i in range(rows) for j in range(columns)],
        "weights": [1.0] * (rows * columns),
    }


def curve2_json(curve: dict[str, Any]) -> dict[str, Any]:
    return {
        "format": "n2s-curve2",
        "version": 1,
        "degree": curve["degree"],
        "knots": curve["knots"],
        "control_points": [[p[0], p[1]] for p in curve["points"]],
        "weights": [1.0] * len(curve["points"]),
    }


def uniform_fit(points: list[tuple[float, ...]], margin: float = 0.02):
    """A uniform scale and translation carrying `points` into [0,1]^2.

    Uniform, not per-axis: squashing the domain to fill the unit square would
    change the shape of every trim curve and every layout quad in it.
    """
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    span = max(max(xs) - min(xs), max(ys) - min(ys))
    if span <= 0:
        raise ValueError("degenerate domain")

    scale = (1.0 - 2.0 * margin) / span
    cx = 0.5 * (min(xs) + max(xs))
    cy = 0.5 * (min(ys) + max(ys))

    def transform(p):
        return (0.5 + (p[0] - cx) * scale, 0.5 + (p[1] - cy) * scale)

    return transform


def snap_to_layout(curves, layout_vertices, tolerance, report):
    """Snaps stray trim endpoints onto the layout vertex they belong to.

    The loop *is* the layout's boundary, and all but a handful of endpoints
    already coincide with a layout vertex exactly, so the layout is the
    authority on where the strays were meant to be.
    """
    for index, curve in enumerate(curves):
        for which in (0, -1):
            p = curve["points"][which]
            best, best_distance = None, float("inf")
            for v in layout_vertices:
                d = ((p[0] - v[0]) ** 2 + (p[1] - v[1]) ** 2) ** 0.5
                if d < best_distance:
                    best, best_distance = v, d
            if best is not None and 0.0 < best_distance <= tolerance:
                points = list(curve["points"])
                points[which] = (best[0], best[1], 0.0)
                curve["points"] = points
                report.append(
                    f"    snapped curve {index} "
                    f"{'start' if which == 0 else 'end'} by {best_distance:.4g}"
                )


def worst_join_gap(curves) -> float:
    worst = 0.0
    for i, curve in enumerate(curves):
        end = curve["points"][-1]
        start = curves[(i + 1) % len(curves)]["points"][0]
        worst = max(worst, ((end[0] - start[0]) ** 2 + (end[1] - start[1]) ** 2) ** 0.5)
    return worst


def orient_ccw(curves):
    """Reverses the loop if it runs clockwise.

    An exported case should load without needing repair. Leaving a clockwise
    outer loop in the file would mean every reader silently reverses it, and a
    silent repair on load is how an orientation defect survives unnoticed --
    which is defect territory this project exists to stay out of.
    """
    points = [p for c in curves for p in c["points"]]
    area = sum(
        points[i][0] * points[(i + 1) % len(points)][1]
        - points[(i + 1) % len(points)][0] * points[i][1]
        for i in range(len(points))
    ) / 2.0
    if area >= 0.0:
        return curves, False

    reversed_curves = []
    for curve in reversed(curves):
        knots = curve["knots"]
        a, b = knots[0], knots[-1]
        # Mirror the knot vector within its own domain so the reversed curve is
        # parameterised over the same interval.
        reversed_curves.append(
            {
                "degree": curve["degree"],
                "knots": [a + b - k for k in reversed(knots)],
                "points": list(reversed(curve["points"])),
            }
        )
    return reversed_curves, True


DEFAULT_RANGES = {
    "mean_curvature": [-2.0, 2.0],
    "gaussian_curvature": [-4.0, 4.0],
    "error": [0.0, 0.05],
}


def write_case(out: Path, name: str, surface, trim_curves, layout=None) -> None:
    document = {
        "format": "n2s-case",
        "version": 1,
        "name": name,
        "surface": surface,
        "trim": {"outer": [curve2_json(c) for c in trim_curves], "holes": []},
        "color_ranges": DEFAULT_RANGES,
    }
    if layout is not None:
        document["layout"] = layout

    path = out / f"{name}.json"
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {path} ({path.stat().st_size / 1024:.1f} KB)")


# --------------------------------------------------------------------------


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__)
        return 2

    assets = Path(argv[1])
    out = Path(argv[2])
    out.mkdir(parents=True, exist_ok=True)

    # ---- Case 1: the small Bezier-trimmed patch --------------------------
    print("thesis_slot (folder 1):")
    surface = surface_json(coons_control_net(read_boundary_surface(assets / "1" / "surface_1x1.txt")))
    trim = read_bezier_trim(assets / "1" / "bspline_teszt_2.txt")
    trim, flipped = orient_ccw(trim)
    print(f"    {len(trim)} cubic Bezier segments, worst join gap {worst_join_gap(trim):.3g}"
          f"{', reversed to counter-clockwise' if flipped else ''}")
    write_case(out, "thesis_slot", surface, trim)

    # ---- Folder 2 (pelda1) is deliberately not imported -------------------
    #
    # Neither of its trim files is an ordered loop. Splitting them wherever
    # consecutive curves fail to meet gives 31 and 9 fragments respectively,
    # with gaps of 6 to 100 units between them and no consistent orientation.
    # A closed loop can only be recovered from that by deciding which fragment
    # follows which and bridging the gaps -- which is guessing, and a trimmed
    # region built on a guessed boundary would make every number measured
    # against it meaningless. The cage (`cage_1.obj`) is 3D control points with
    # no domain correspondence, so it cannot stand in for the missing loop
    # either.
    #
    # If the intended boundary can be recovered -- from the thesis text, or by
    # re-exporting it -- this is a three-line addition. Until then the two
    # cases below are the usable ones.

    # ---- Case 2: DoubleVB, the one with a hand-authored layout -----------
    print("thesis_doublevb (folder 3):")
    trim = read_bsc(assets / "3" / "DoubleVB_trim_ordered.txt")
    vertices, faces = read_obj(assets / "3" / "DoubleVB.obj")
    if any(len(f) != 4 for f in faces):
        raise ValueError("the DoubleVB layout is not all quads")

    print(f"    {len(trim)} trim curves, {len(vertices)} layout vertices, {len(faces)} quads")
    print(f"    worst join gap before repair: {worst_join_gap(trim):.4g}")

    report: list[str] = []
    span = max(
        max(v[0] for v in vertices) - min(v[0] for v in vertices),
        max(v[1] for v in vertices) - min(v[1] for v in vertices),
    )
    snap_to_layout(trim, vertices, 0.01 * span, report)
    for line in report:
        print(line)
    print(f"    worst join gap after repair:  {worst_join_gap(trim):.4g}")

    trim, flipped = orient_ccw(trim)
    if flipped:
        print("    reversed the loop to counter-clockwise")

    transform = uniform_fit(vertices + [p for c in trim for p in c["points"]])
    for curve in trim:
        curve["points"] = [(*transform(p), 0.0) for p in curve["points"]]
    layout = {
        "vertices": [list(transform(v)) for v in vertices],
        "quads": [list(f) for f in faces],
    }

    surface = surface_json(
        coons_control_net(read_boundary_surface(assets / "3" / "surface_1x1_2.txt"))
    )
    write_case(out, "thesis_doublevb", surface, trim, layout)

    print("\nSee the module docstring for the four reconstruction decisions this makes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
