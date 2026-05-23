"""Cross-method speedup / coverage / win-loss tables from a v7-format CSV.

Reads `cross_check_v*.csv` (rows: Graph, S, T, D, MethodId, Method,
Status, Reliability, TimeSec, Recursions, CompletedRuns, MinTimeSec,
MaxTimeSec, TimeoutBudgetSec, Error) and emits four artefacts next
to the input file (or under --outdir):

1. cross_check_speedup_vs_m4.csv
   For every cell where both `m_i` and the reference `m4` finished OK,
   speedup = m4.TimeSec / m_i.TimeSec. Aggregated per (family, method)
   into median, geomean, min, max plus n cells contributing. Helps
   the thesis answer "how much faster is m5 than m4 on sausage chains
   where both finish?".

2. cross_check_coverage.csv
   % OK cells per (family, method). Rough "how often does the method
   work at all".

3. cross_check_winloss.csv
   For each ordered pair (m_i, m_j), the count of cells where m_i
   finishes strictly faster than m_j (both OK) plus ties (|t_i-t_j|/
   max < 1e-3). Symmetric mostly-OK pair-wise compete.

4. cross_check_recursions_per_sec.csv
   Recursions / TimeSec per method × family. Probes "work intensity"
   — m0 does ~1e5 rec/s, m3 ~3e4 rec/s but each rec is cheaper. A
   complementary lens to wall time.

Usage:
    python scripts/analyze_speedups.py cross_check_v7.csv [--baseline 4] [--outdir DIR]

stdlib only.
"""

import argparse
import csv
import math
import os
import statistics
import sys
from collections import defaultdict


METHOD_IDS = [0, 1, 2, 3, 4, 5]


def family_of(graph: str) -> str:
    g = graph.lower()
    if g.startswith("k4"):
        return "K4"
    if "sausage_s3" in g:
        return "sausage-s3"
    if "sausage" in g:
        return "sausage-3x3"
    if "geant2004" in g:
        return "Geant2004"
    if "geant2009" in g:
        return "Geant2009"
    if "ieee" in g:
        return "IEEE-118"
    if "ups" in g and "russia" in g:
        return "UPS-Russia"
    return "other"


def load(path: str):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            def fnum(k, conv):
                v = r.get(k, "")
                if v is None or v == "":
                    return None
                try:
                    return conv(v)
                except ValueError:
                    return None
            rows.append({
                "Graph":       r.get("Graph", ""),
                "Family":      family_of(r.get("Graph", "")),
                "D":           fnum("D", int),
                "MethodId":    fnum("MethodId", int),
                "Status":      r.get("Status", ""),
                "TimeSec":     fnum("TimeSec", float),
                "Recursions":  fnum("Recursions", int),
                "MinTimeSec":  fnum("MinTimeSec", float),
                "MaxTimeSec":  fnum("MaxTimeSec", float),
            })
    return rows


def cell_key(r):
    return (r["Graph"], r["D"])


def speedup_vs_baseline(rows, baseline_id, outdir):
    """For every cell with both m_i and m_baseline OK, compute t_baseline/t_i."""
    # Build per-cell time per method
    by_cell = defaultdict(dict)  # (graph, d) -> {mid: time}
    for r in rows:
        if r["Status"] != "OK" or r["TimeSec"] is None:
            continue
        by_cell[cell_key(r)][r["MethodId"]] = (r["TimeSec"], r["Family"])

    # Aggregate speedup per (family, method)
    bucket = defaultdict(list)  # (family, mid) -> list of speedups
    for ck, methods in by_cell.items():
        if baseline_id not in methods:
            continue
        t_base, _ = methods[baseline_id]
        if t_base <= 0:
            continue
        for mid, (t_m, fam) in methods.items():
            if t_m <= 0:
                continue
            # Guard against 0-duration rows (m4/m5 sub-microsecond)
            speedup = t_base / max(t_m, 1e-9)
            bucket[(fam, mid)].append(speedup)

    # ALL bucket (across all families)
    all_bucket = defaultdict(list)
    for (fam, mid), vals in bucket.items():
        all_bucket[mid].extend(vals)
    for mid, vals in all_bucket.items():
        bucket[("ALL", mid)] = vals

    out = []
    families = ["ALL"] + sorted(set(fam for (fam, _) in bucket.keys() if fam != "ALL"))
    for fam in families:
        for mid in METHOD_IDS:
            vals = bucket.get((fam, mid))
            if not vals:
                continue
            geomean = math.exp(statistics.fmean(math.log(v) for v in vals))
            out.append({
                "Family":   fam,
                "MethodId": mid,
                "n_cells":  len(vals),
                "median":   round(statistics.median(vals), 6),
                "geomean":  round(geomean, 6),
                "min":      round(min(vals), 6),
                "max":      round(max(vals), 6),
            })

    path = os.path.join(outdir, "cross_check_speedup_vs_m{}.csv".format(baseline_id))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["Family", "MethodId", "n_cells",
                                          "median", "geomean", "min", "max"])
        w.writeheader()
        for row in out:
            w.writerow(row)
    return path, out


def coverage(rows, outdir):
    """% OK per (family, method)."""
    fam_method_total = defaultdict(int)
    fam_method_ok = defaultdict(int)
    for r in rows:
        if r["MethodId"] is None:
            continue
        key = (r["Family"], r["MethodId"])
        fam_method_total[key] += 1
        if r["Status"] == "OK":
            fam_method_ok[key] += 1

    # Also produce ALL row
    method_total = defaultdict(int)
    method_ok = defaultdict(int)
    for (fam, mid), n in fam_method_total.items():
        method_total[mid] += n
    for (fam, mid), n in fam_method_ok.items():
        method_ok[mid] += n
    for mid in METHOD_IDS:
        fam_method_total[("ALL", mid)] = method_total[mid]
        fam_method_ok[("ALL", mid)] = method_ok[mid]

    out = []
    families = ["ALL"] + sorted(set(fam for (fam, _) in fam_method_total.keys()
                                    if fam != "ALL"))
    for fam in families:
        for mid in METHOD_IDS:
            total = fam_method_total.get((fam, mid), 0)
            if total == 0:
                continue
            ok = fam_method_ok.get((fam, mid), 0)
            out.append({
                "Family":     fam,
                "MethodId":   mid,
                "n_total":    total,
                "n_ok":       ok,
                "pct_ok":     round(100.0 * ok / total, 1),
            })

    path = os.path.join(outdir, "cross_check_coverage.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["Family", "MethodId", "n_total",
                                          "n_ok", "pct_ok"])
        w.writeheader()
        for row in out:
            w.writerow(row)
    return path, out


def winloss(rows, outdir, tie_rel=1e-3):
    """Pairwise win/loss/tie count on cells where both methods OK."""
    by_cell = defaultdict(dict)
    for r in rows:
        if r["Status"] != "OK" or r["TimeSec"] is None:
            continue
        by_cell[cell_key(r)][r["MethodId"]] = r["TimeSec"]

    # (mi, mj) -> [wins, losses, ties]
    counts = defaultdict(lambda: [0, 0, 0])
    for methods in by_cell.values():
        present = sorted(methods.keys())
        for i in present:
            for j in present:
                if i == j:
                    continue
                ti, tj = methods[i], methods[j]
                hi = max(ti, tj, 1e-9)
                if abs(ti - tj) / hi < tie_rel:
                    counts[(i, j)][2] += 1
                elif ti < tj:
                    counts[(i, j)][0] += 1
                else:
                    counts[(i, j)][1] += 1

    path = os.path.join(outdir, "cross_check_winloss.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["MethodA", "MethodB", "A_wins", "A_losses", "ties",
                    "comparable_cells"])
        for i in METHOD_IDS:
            for j in METHOD_IDS:
                if i == j:
                    continue
                wins, losses, ties = counts[(i, j)]
                tot = wins + losses + ties
                if tot == 0:
                    continue
                w.writerow([i, j, wins, losses, ties, tot])
    return path


def recursions_per_sec(rows, outdir):
    """Recursions per OK second, aggregated."""
    bucket = defaultdict(list)
    for r in rows:
        if (r["Status"] != "OK" or r["TimeSec"] is None
                or r["Recursions"] is None or r["TimeSec"] <= 0):
            continue
        rps = r["Recursions"] / r["TimeSec"]
        bucket[(r["Family"], r["MethodId"])].append(rps)
        bucket[("ALL", r["MethodId"])].append(rps)

    out = []
    families = ["ALL"] + sorted(set(fam for (fam, _) in bucket.keys()
                                    if fam != "ALL"))
    for fam in families:
        for mid in METHOD_IDS:
            vals = bucket.get((fam, mid))
            if not vals:
                continue
            out.append({
                "Family":   fam,
                "MethodId": mid,
                "n_cells":  len(vals),
                "median_rps": round(statistics.median(vals), 2),
                "geomean_rps": round(math.exp(statistics.fmean(
                    math.log(max(v, 1.0)) for v in vals)), 2),
            })

    path = os.path.join(outdir, "cross_check_recursions_per_sec.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["Family", "MethodId", "n_cells",
                                          "median_rps", "geomean_rps"])
        w.writeheader()
        for row in out:
            w.writerow(row)
    return path


def print_speedup_pretty(stats, baseline_id):
    print(f"\n=== Speedup vs m{baseline_id} (t_m{baseline_id} / t_method, > 1 = faster than m{baseline_id}) ===")
    print(f"{'Family':<14} {'M':>2} {'n':>4} {'median':>10} {'geomean':>10} {'min':>10} {'max':>12}")
    for r in stats:
        print(f"{r['Family']:<14} {r['MethodId']:>2} {r['n_cells']:>4} "
              f"{r['median']:>10.3f} {r['geomean']:>10.3f} "
              f"{r['min']:>10.3f} {r['max']:>12.3f}")


def print_coverage_pretty(stats):
    print(f"\n=== Coverage (% OK per family x method) ===")
    print(f"{'Family':<14} {'M':>2} {'n_total':>8} {'n_ok':>5} {'pct_ok':>7}")
    for r in stats:
        print(f"{r['Family']:<14} {r['MethodId']:>2} {r['n_total']:>8} "
              f"{r['n_ok']:>5} {r['pct_ok']:>6}%")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--baseline", type=int, default=4,
                    help="Method id for the speedup denominator (default 4)")
    ap.add_argument("--outdir", default=None,
                    help="Where to write outputs (default = CSV's dir)")
    args = ap.parse_args()

    if not os.path.exists(args.csv):
        print(f"File not found: {args.csv}", file=sys.stderr)
        return 1

    outdir = args.outdir or os.path.dirname(os.path.abspath(args.csv)) or "."
    os.makedirs(outdir, exist_ok=True)

    rows = load(args.csv)
    print(f"Loaded {len(rows)} rows from {args.csv}")

    sp_path, sp_stats = speedup_vs_baseline(rows, args.baseline, outdir)
    print(f"Wrote {sp_path} ({len(sp_stats)} rows)")
    print_speedup_pretty(sp_stats, args.baseline)

    cov_path, cov_stats = coverage(rows, outdir)
    print(f"\nWrote {cov_path} ({len(cov_stats)} rows)")
    print_coverage_pretty(cov_stats)

    wl_path = winloss(rows, outdir)
    print(f"\nWrote {wl_path}")

    rps_path = recursions_per_sec(rows, outdir)
    print(f"Wrote {rps_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
