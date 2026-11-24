import argparse
import math
import os
import sys

COLS = ["tool", "label", "truth_min_umis", "call_min_umis",
        "n_shared_inf", "n_truthpos_inf", "TP", "FN", "sensitivity", "sens_lo", "sens_hi",
        "n_shared_uninf", "n_truthpos_uninf", "n_truthneg_uninf", "FP", "TN", "FPR", "fpr_lo", "fpr_hi",
        "FP_strict", "FPR_strict", "fprs_lo", "fprs_hi"]


def load(path):
    d = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if header[:2] != ["barcode", "hiv_umis"]:
            sys.exit(f"ERROR: {path} does not look like a hiv_counts.py table (header: {header})")
        for line in f:
            b, v = line.rstrip("\n").split("\t")[:2]
            d[b.split("-")[0]] = int(v)
    return d


def wilson(k, n, z=1.96):
    if n == 0:
        return None, None, None
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, max(0.0, centre - half), min(1.0, centre + half)


def fmt(x):
    return "NA" if x is None else (f"{x:.4f}" if isinstance(x, float) else str(x))

def check_universe(sample, truth, tool, mode, min_cov):
    """Compare the barcode set of one truth table with one tool table."""
    if mode == "off":
        return
    shared = len(set(truth) & set(tool))
    cov = shared / len(truth) if truth else 0.0
    n_zero = sum(1 for v in tool.values() if v == 0)
    print(f"barcode universe ({sample}): truth {len(truth)}, tool {len(tool)}, shared {shared} "
          f"({cov:.1%} of truth), tool-only {len(tool) - shared}, tool rows with 0 UMIs {n_zero}",
          file=sys.stderr)
    problems = []
    if cov < min_cov:
        problems.append(f"tool table covers only {cov:.1%} of the {len(truth)} truth barcodes "
                        f"(minimum {min_cov:.0%}); the {len(truth) - shared} missing truth cells are dropped from scoring")
    if tool and n_zero == 0:
        problems.append("tool table has no zero-count rows, so it probably lists only positive cells "
                        "and its negatives are missing")
    for p in problems:
        print(f"{'ERROR' if mode == 'error' else 'WARNING'} ({sample}): {p}", file=sys.stderr)
    if problems and mode == "error":
        return problems

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tool", required=True)
    ap.add_argument("--label", default="", help="free text, e.g. 25M_rep1")
    ap.add_argument("--truth-inf", required=True)
    ap.add_argument("--truth-uninf", required=True)
    ap.add_argument("--calls-inf", required=True)
    ap.add_argument("--calls-uninf", required=True)
    ap.add_argument("--truth-thresholds", type=int, nargs="+", default=[1, 2, 5])
    ap.add_argument("--call-min", type=int, default=1)
    ap.add_argument("--out", help="append rows to this file")
    ap.add_argument("--universe-check", choices=["error", "warn", "off"], default="error",
                help="what to do if a tool table's barcodes don't cover the truth barcodes")
    ap.add_argument("--min-truth-coverage", type=float, default=0.95,
                help="minimum fraction of truth barcodes the tool table must contain")
    a = ap.parse_args()

    ti, tu = load(a.truth_inf), load(a.truth_uninf)
    ci, cu = load(a.calls_inf), load(a.calls_uninf)
    bad = [check_universe("infected", ti, ci, a.universe_check, a.min_truth_coverage),check_universe("uninfected", tu, cu, a.universe_check, a.min_truth_coverage)]
    if a.universe_check == "error" and any(bad):
        sys.exit("Fix the tool table (see hiv_counts.py) or rerun with --universe-check warn/off to override.")
    shared_i = sorted(set(ti) & set(ci))
    shared_u = sorted(set(tu) & set(cu))
    print(f"shared barcodes: infected {len(shared_i)} (truth {len(ti)}, tool {len(ci)}); "
          f"uninfected {len(shared_u)} (truth {len(tu)}, tool {len(cu)})", file=sys.stderr)
    if not shared_i or not shared_u:
        sys.exit("ERROR: no shared barcodes: check that both tables use the same barcode format")

    fp_strict = sum(1 for b in shared_u if cu[b] >= a.call_min)
    pS, loS, hiS = wilson(fp_strict, len(shared_u))

    rows = []
    for t in a.truth_thresholds:
        tp_cells = [b for b in shared_i if ti[b] >= t]
        tp = sum(1 for b in tp_cells if ci[b] >= a.call_min)
        fn = len(tp_cells) - tp
        neg_u = [b for b in shared_u if tu[b] < t]
        fp = sum(1 for b in neg_u if cu[b] >= a.call_min)
        tn = len(neg_u) - fp
        s, slo, shi = wilson(tp, len(tp_cells))
        f, flo, fhi = wilson(fp, len(neg_u))
        rows.append([a.tool, a.label, t, a.call_min,
                     len(shared_i), len(tp_cells), tp, fn, s, slo, shi,
                     len(shared_u), len(shared_u) - len(neg_u), len(neg_u), fp, tn, f, flo, fhi,
                     fp_strict, pS, loS, hiS])

    lines = ["\t".join(fmt(x) for x in r) for r in rows]
    print("\t".join(COLS))
    print("\n".join(lines))
    if a.out:
        new = not os.path.exists(a.out) or os.path.getsize(a.out) == 0
        with open(a.out, "a") as f:
            if new:
                f.write("\t".join(COLS) + "\n")
            f.write("\n".join(lines) + "\n")
        print(f"appended {len(lines)} rows to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()