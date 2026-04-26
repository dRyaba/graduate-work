"""Post-process cross_check_results.csv into thesis-ready artefacts.

Reads the raw cross-check CSV (columns:
  Graph, S, T, D, MethodId, Method, Status, Reliability, TimeSec,
  Recursions, Error)
and produces:

  cross_check_agreement_matrix.csv
      6x6 matrix: agreement[i][j] = #cells where both m_i and m_j
      completed OK and |R_i - R_j| <= tolerance. By default, trivial
      cells (all OK methods return R == 0) are excluded — use
      --keep-trivial to include them.

  cross_check_disagreements.csv
      Rows for every (Graph, D, i, j) pair where both completed OK but
      |R_i - R_j| > tolerance.

  cross_check_trivial_cells.csv
      Cells where every OK method returned R == 0 (i.e. d < dist(s,t)
      or no s-t path survives). Useful for auditing the experimental
      range, not for agreement claims.

  cross_check_perf_summary.csv
      Per-method performance: count of OK/TIMEOUT/ERROR, median TimeSec
      and 95th-percentile TimeSec (OK cells only), aggregated globally
      and per graph-family. Excludes trivial cells by default.

Usage:
  python scripts/cross_check_analyze.py cross_check_results.csv
  python scripts/cross_check_analyze.py --tolerance 1e-10 cross_check_results.csv
  python scripts/cross_check_analyze.py --keep-trivial cross_check_results.csv

No external dependencies: uses only the Python 3 standard library.
"""

import argparse
import csv
import os
import statistics
import sys
from collections import defaultdict


METHOD_IDS = [0, 1, 2, 3, 4, 5]


def family_of(graph: str) -> str:
    g = graph.lower()
    if g.startswith("k4"):
        return "K4"
    if "sausage" in g:
        return "sausage"
    if "geant" in g:
        return "Geant2004"
    if "ieee" in g:
        return "IEEE-118"
    return "other"


def load(path: str):
    """Return a list of row dicts with typed fields."""
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            def as_float(x):
                try:
                    return float(x) if x != "" else None
                except ValueError:
                    return None

            def as_int(x):
                try:
                    return int(x) if x != "" else None
                except ValueError:
                    return None

            rows.append({
                "Graph":       r.get("Graph", ""),
                "S":           as_int(r.get("S", "")),
                "T":           as_int(r.get("T", "")),
                "D":           as_int(r.get("D", "")),
                "MethodId":    as_int(r.get("MethodId", "")),
                "Method":      r.get("Method", ""),
                "Status":      r.get("Status", ""),
                "Reliability": as_float(r.get("Reliability", "")),
                "TimeSec":     as_float(r.get("TimeSec", "")),
                "Recursions":  as_int(r.get("Recursions", "")),
                "Error":       r.get("Error", ""),
            })
    return rows


def group_by_cell(rows):
    """Return {(graph, d): {method_id: reliability}} for OK rows only."""
    by_cell = defaultdict(dict)
    for r in rows:
        if r["Status"] != "OK":
            continue
        if r["MethodId"] is None or r["Reliability"] is None:
            continue
        by_cell[(r["Graph"], r["D"])][r["MethodId"]] = r["Reliability"]
    return by_cell


def is_trivial_cell(per_method, tol: float) -> bool:
    """True iff every OK method in this cell returned R == 0 (within tol)."""
    if not per_method:
        return False
    return all(abs(v) <= tol for v in per_method.values())


def split_trivial(by_cell, tol: float):
    """Split {(graph, d): {mid: R}} into (non_trivial, trivial) dicts."""
    non_trivial, trivial = {}, {}
    for key, per_method in by_cell.items():
        (trivial if is_trivial_cell(per_method, tol) else non_trivial)[key] = per_method
    return non_trivial, trivial


def agreement_matrix(by_cell, tol: float):
    counts = [[0] * 6 for _ in range(6)]
    for per_method in by_cell.values():
        for i in METHOD_IDS:
            for j in METHOD_IDS:
                if i in per_method and j in per_method:
                    if abs(per_method[i] - per_method[j]) <= tol:
                        counts[i][j] += 1
    return counts


def pair_disagreements(by_cell, tol: float):
    out = []
    for (graph, d), per_method in by_cell.items():
        methods = sorted(per_method.keys())
        for a_i, i in enumerate(methods):
            for j in methods[a_i + 1:]:
                diff = abs(per_method[i] - per_method[j])
                if diff > tol:
                    out.append({
                        "Graph": graph,
                        "D": d,
                        "MethodA": i,
                        "MethodB": j,
                        "R_A": per_method[i],
                        "R_B": per_method[j],
                        "AbsDiff": diff,
                    })
    out.sort(key=lambda x: (-x["AbsDiff"], x["Graph"], x["D"]))
    return out


def trivial_cells_rows(trivial_by_cell):
    """Flatten trivial cells into CSV rows for audit."""
    out = []
    for (graph, d), per_method in sorted(trivial_by_cell.items()):
        out.append({
            "Graph": graph,
            "D": d,
            "OKMethods": ";".join(f"m{m}" for m in sorted(per_method.keys())),
            "Reliability": 0.0,
        })
    return out


def percentile(xs, q):
    """Linear-interpolation quantile. Works with as few as 1 sample."""
    if not xs:
        return None
    xs = sorted(xs)
    if len(xs) == 1:
        return xs[0]
    idx = q * (len(xs) - 1)
    lo = int(idx)
    hi = min(lo + 1, len(xs) - 1)
    frac = idx - lo
    return xs[lo] * (1 - frac) + xs[hi] * frac


def perf_summary(rows, excluded_keys=None):
    """Per-(scope, method) perf stats. `excluded_keys` filters out rows from
    specified (graph, d) cells (typically trivial cells)."""
    excluded = excluded_keys or set()
    filtered = [r for r in rows if (r["Graph"], r["D"]) not in excluded]
    by_scope_method = defaultdict(lambda: defaultdict(list))

    def add_scope(scope_name, scope_rows):
        for r in scope_rows:
            if r["MethodId"] is None:
                continue
            by_scope_method[scope_name][r["MethodId"]].append(r)

    add_scope("ALL", filtered)
    families = defaultdict(list)
    for r in filtered:
        families[family_of(r["Graph"])].append(r)
    for fam, fam_rows in families.items():
        add_scope(fam, fam_rows)

    out = []
    for scope in ["ALL"] + sorted(k for k in by_scope_method.keys() if k != "ALL"):
        methods = by_scope_method[scope]
        for mid in METHOD_IDS:
            if mid not in methods:
                continue
            sub = methods[mid]
            statuses = [r["Status"] for r in sub]
            ok_times = [r["TimeSec"] for r in sub
                        if r["Status"] == "OK" and r["TimeSec"] is not None]
            median_t = statistics.median(ok_times) if ok_times else None
            p95_t = percentile(ok_times, 0.95) if ok_times else None
            out.append({
                "Scope": scope,
                "MethodId": mid,
                "n_OK":      statuses.count("OK"),
                "n_TIMEOUT": statuses.count("TIMEOUT"),
                "n_ERROR":   statuses.count("ERROR"),
                "median_time_sec": median_t,
                "p95_time_sec":    p95_t,
            })
    return out


def write_matrix_csv(path, mat):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([""] + [f"m{j}" for j in METHOD_IDS])
        for i, row in enumerate(mat):
            w.writerow([f"m{i}"] + row)


def write_rows_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def format_matrix(mat):
    lines = ["    " + "  ".join(f"m{j:>2}" for j in METHOD_IDS)]
    for i, row in enumerate(mat):
        lines.append(f"m{i}  " + "  ".join(f"{v:>3}" for v in row))
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", help="Path to cross_check_results.csv")
    ap.add_argument("--tolerance", type=float, default=1e-10)
    ap.add_argument("--outdir", default=None,
                    help="Where to write *_agreement_matrix.csv etc. "
                         "Defaults to the CSV's directory.")
    ap.add_argument("--keep-trivial", action="store_true",
                    help="Include cells where every OK method returned R=0 "
                         "in the agreement matrix and perf summary. By default "
                         "they are excluded and written to "
                         "cross_check_trivial_cells.csv for audit only.")
    args = ap.parse_args()

    if not os.path.exists(args.csv):
        print(f"File not found: {args.csv}", file=sys.stderr)
        return 1

    outdir = args.outdir or os.path.dirname(os.path.abspath(args.csv)) or "."
    os.makedirs(outdir, exist_ok=True)
    rows = load(args.csv)
    print(f"Loaded {len(rows)} rows from {args.csv}")

    by_cell_all = group_by_cell(rows)
    non_trivial, trivial = split_trivial(by_cell_all, args.tolerance)
    print(f"Non-trivial cells: {len(non_trivial)}, trivial (all-zero): {len(trivial)}")

    by_cell_for_matrix = by_cell_all if args.keep_trivial else non_trivial
    excluded_for_perf = set() if args.keep_trivial else set(trivial.keys())

    mat = agreement_matrix(by_cell_for_matrix, args.tolerance)
    mat_path = os.path.join(outdir, "cross_check_agreement_matrix.csv")
    write_matrix_csv(mat_path, mat)
    print(f"\nWrote {mat_path}"
          f" ({'incl.' if args.keep_trivial else 'excl.'} trivial cells)")
    print(format_matrix(mat))

    dis = pair_disagreements(by_cell_for_matrix, args.tolerance)
    dis_path = os.path.join(outdir, "cross_check_disagreements.csv")
    write_rows_csv(dis_path, dis,
                   ["Graph", "D", "MethodA", "MethodB",
                    "R_A", "R_B", "AbsDiff"])
    print(f"\nWrote {dis_path} ({len(dis)} pair-disagreements > {args.tolerance})")
    for r in dis[:20]:
        print(f"  {r['Graph']} d={r['D']} m{r['MethodA']} vs m{r['MethodB']}"
              f"  R_A={r['R_A']:.15f} R_B={r['R_B']:.15f}"
              f"  |diff|={r['AbsDiff']:.3e}")

    triv_rows = trivial_cells_rows(trivial)
    triv_path = os.path.join(outdir, "cross_check_trivial_cells.csv")
    write_rows_csv(triv_path, triv_rows,
                   ["Graph", "D", "OKMethods", "Reliability"])
    print(f"\nWrote {triv_path} ({len(triv_rows)} all-zero cells)")

    perf = perf_summary(rows, excluded_keys=excluded_for_perf)
    perf_path = os.path.join(outdir, "cross_check_perf_summary.csv")
    write_rows_csv(perf_path, perf,
                   ["Scope", "MethodId", "n_OK", "n_TIMEOUT", "n_ERROR",
                    "median_time_sec", "p95_time_sec"])
    print(f"\nWrote {perf_path}")
    header = f"{'Scope':<12} {'M':>2} {'OK':>4} {'TO':>4} {'ERR':>4} {'median':>10} {'p95':>10}"
    print(header)
    for r in perf:
        med = f"{r['median_time_sec']:.4f}" if r['median_time_sec'] is not None else "—"
        p95 = f"{r['p95_time_sec']:.4f}"    if r['p95_time_sec']    is not None else "—"
        print(f"{r['Scope']:<12} {r['MethodId']:>2} {r['n_OK']:>4} "
              f"{r['n_TIMEOUT']:>4} {r['n_ERROR']:>4} {med:>10} {p95:>10}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
