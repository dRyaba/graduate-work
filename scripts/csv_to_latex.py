r"""Generate LaTeX tables for the thesis chapter from a v7-format CSV.

For each (graph, s, t) group emits a table block in the format
the user is already using in the thesis draft:

  \begin{table}[h]
  \centering
  \caption{<human-readable graph name>, пара $(s,t) = (s+1, t+1)$, $d_{s,t} = <dist>$}
  \label{tab:<slug>}
  \begin{tabular}{|c|c|c|c|c|c|c|c|}
  \hline
  $D$ & $\mathrm{gap}$ & m0 & m1 & m2 & m3 & m4 & m5 \\
  \hline
  <d> & <gap> & <t_m0> & ... & <t_m5> \\
  ...
  \end{tabular}
  \end{table}

Vertex indices in the caption are converted to **1-based** to match
the existing thesis convention (CSV stores 0-based). Cell values are
the median wall time (seconds) for OK rows, `TO` for TIMEOUT
(annotated with the budget), `ERR` for ERROR.

Usage:
    python scripts/csv_to_latex.py cross_check_v7.csv [--outdir DIR]

stdlib only; UTF-8 output (LaTeX captions are in Russian).
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


METHOD_IDS = [0, 1, 2, 3, 4, 5]


GRAPH_LABELS = {
    "K4_kao.txt":                                      ("K4",                              "k4"),
    "2_blocks_sausage_s3_kao.txt":                     ("Цепь из 2 решёток $s_3$",          "s3-2"),
    "3_blocks_sausage_3x3_kao.txt":                    ("Цепь из 3 решёток $3 \\times 3$",  "grid3"),
    "4_blocks_sausage_3x3_kao.txt":                    ("Цепь из 4 решёток $3 \\times 3$",  "grid4"),
    "4_blocks_sausage_s3_kao.txt":                     ("Цепь из 4 решёток $s_3$",          "s3-4"),
    "5_blocks_sausage_3x3_kao.txt":                    ("Цепь из 5 решёток $3 \\times 3$",  "grid5"),
    "6_blocks_sausage_3x3_kao.txt":                    ("Цепь из 6 решёток $3 \\times 3$",  "grid6"),
    "6_blocks_sausage_s3_kao.txt":                     ("Цепь из 6 решёток $s_3$",          "s3-6"),
    "Geant2004_kao.txt":                               ("Geant2004",                       "geant2004"),
    "Geant2009_kao.txt":                               ("Geant2009",                       "geant2009"),
    "IEEE-118-node_kao.txt":                           ("IEEE-118",                        "ieee118"),
    "UPS_of_Russia_composed_with_colored_cut_vertices_kao.txt":
                                                       ("UPS of Russia",                   "ups"),
}


def fmt_time(t: float) -> str:
    """Format wall time with a reasonable number of significant digits."""
    if t is None:
        return ""
    if t == 0.0:
        return "0"
    if t < 1e-3:
        return f"{t * 1e6:.1f}\\,$\\mu$s"
    if t < 1.0:
        return f"{t * 1e3:.2f}\\,мс"
    if t < 60.0:
        return f"{t:.2f}\\,с"
    minutes = t / 60.0
    if minutes < 60.0:
        return f"{minutes:.1f}\\,мин"
    return f"{minutes / 60.0:.1f}\\,ч"


def cell_value(row) -> str:
    status = row.get("Status", "")
    if status == "OK":
        try:
            return fmt_time(float(row.get("TimeSec") or 0))
        except ValueError:
            return ""
    if status == "TIMEOUT":
        try:
            budget = int(row.get("TimeoutBudgetSec") or 0)
            return f"\\textit{{TO}}\\,$\\!>\\!${budget // 60}\\,мин"
        except ValueError:
            return "\\textit{TO}"
    if status == "ERROR":
        return "\\textit{ERR}"
    return ""


def load(path: str):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def group_rows(rows):
    """Group by (graph, s, t). Return list of (label, slug, s, t, dist, [(d, {m: row}) sorted by d])."""
    groups = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        try:
            s = int(r["S"]); t = int(r["T"])
            d = int(r["D"]); m = int(r["MethodId"])
        except (KeyError, ValueError):
            continue
        groups[(r["Graph"], s, t)][d][m] = r

    out = []
    for (graph, s, t), per_d in groups.items():
        ds = sorted(per_d.keys())
        if not ds:
            continue
        dist = min(ds)  # d_count auto-ranges from dist upward
        label, slug = GRAPH_LABELS.get(graph, (graph, graph.replace(".", "-")))
        rows_by_d = [(d, per_d[d]) for d in ds]
        out.append({
            "graph":  graph,
            "label":  label,
            "slug":   slug,
            "s":      s,
            "t":      t,
            "dist":   dist,
            "rows":   rows_by_d,
        })
    # Stable order matching the user's existing thesis chapter:
    # sausage chains first (sorted by file name), then real-world.
    def sort_key(g):
        name = g["graph"]
        is_sausage = "sausage" in name.lower()
        return (0 if is_sausage else 1, name)
    out.sort(key=sort_key)
    return out


def render_table(g) -> str:
    """LaTeX table for one (graph, s, t) group."""
    # 1-based vertices for the caption (thesis convention).
    s1 = g["s"] + 1
    t1 = g["t"] + 1
    cols = "|c|c|c|c|c|c|c|c|"
    head = (
        f"\\begin{{table}}[h]\n"
        f"\\centering\n"
        f"\\caption{{{g['label']}, пара $(s,t) = ({s1}, {t1})$, "
        f"$d_{{s,t}} = {g['dist']}$}}\n"
        f"\\label{{tab:{g['slug']}}}\n"
        f"\\begin{{tabular}}{{{cols}}}\n"
        f"\\hline\n"
        f"$D$ & $\\mathrm{{gap}}$ & m0 & m1 & m2 & m3 & m4 & m5 \\\\\n"
        f"\\hline\n"
    )
    body_lines = []
    for d, methods in g["rows"]:
        gap = d - g["dist"]
        cells = []
        for mid in METHOD_IDS:
            r = methods.get(mid)
            cells.append(cell_value(r) if r else "")
        body_lines.append(
            f"{d} & {gap} & " + " & ".join(cells) + " \\\\\n\\hline\n"
        )
    tail = "\\end{tabular}\n\\end{table}\n"
    return head + "".join(body_lines) + tail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--outdir", default=None,
                    help="Where to write tables.tex (default = CSV's dir)")
    args = ap.parse_args()

    if not os.path.exists(args.csv):
        print(f"File not found: {args.csv}", file=sys.stderr)
        return 1

    outdir = args.outdir or os.path.dirname(os.path.abspath(args.csv)) or "."
    os.makedirs(outdir, exist_ok=True)

    rows = load(args.csv)
    groups = group_rows(rows)
    if not groups:
        print("No groups found", file=sys.stderr)
        return 1

    out_path = os.path.join(outdir, "cross_check_tables.tex")
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("% Generated by scripts/csv_to_latex.py\n")
        f.write(f"% Source: {os.path.basename(args.csv)}\n")
        f.write(f"% Groups: {len(groups)}\n\n")
        for g in groups:
            f.write(render_table(g))
            f.write("\n")

    print(f"Wrote {out_path} ({len(groups)} tables)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
