import argparse
import itertools
import sys
from barcode_tables import load, check_universes

COLS = ["sample", "tool_a", "tool_b", "truth_min_umis", "call_min_umis", "n_cells",
        "both_pos", "only_a", "only_b", "both_neg", "agreement", "jaccard",
        "err_a", "err_b", "err_both"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth-inf", required=True)
    ap.add_argument("--truth-uninf", required=True)
    ap.add_argument("--tool", nargs=3, action="append", metavar=("NAME", "INF", "UNINF"), required=True)
    ap.add_argument("--truth-threshold", type=int, default=1)
    ap.add_argument("--call-min", type=int, default=1)
    ap.add_argument("--universe-check", choices=["error", "warn", "off"], default="error",
                help="what to do if a tool table's barcodes don't cover the truth barcodes")
    ap.add_argument("--min-truth-coverage", type=float, default=0.95,
                help="minimum fraction of truth barcodes each tool table must contain")
    ap.add_argument("--out")
    a = ap.parse_args()
    if len(a.tool) < 2:
        sys.exit("ERROR: give at least two --tool entries")

    truths = {"infected": load(a.truth_inf), "uninfected": load(a.truth_uninf)}
    calls = {name: {"infected": load(i), "uninfected": load(u)} for name, i, u in a.tool}
    check_universes([(f"{name} {sample}", truths[sample], calls[name][sample])
                 for name in calls for sample in truths],
                a.universe_check, a.min_truth_coverage)
    
    rows = []
    for sample, truth in truths.items():
        for na, nb in itertools.combinations(calls, 2):
            ca, cb = calls[na][sample], calls[nb][sample]
            cells = sorted(set(truth) & set(ca) & set(cb))
            pa = {b for b in cells if ca[b] >= a.call_min}
            pb = {b for b in cells if cb[b] >= a.call_min}
            both, only_a, only_b = pa & pb, pa - pb, pb - pa
            neither = len(cells) - len(pa | pb)
            if sample == "infected":      # error = a truth-positive cell that the tool misses
                universe = {b for b in cells if truth[b] >= a.truth_threshold}
                ea, eb = universe - pa, universe - pb
            else:                         # error = a truth-negative cell that the tool calls
                universe = {b for b in cells if truth[b] < a.truth_threshold}
                ea, eb = universe & pa, universe & pb
            union = len(pa | pb)
            rows.append([sample, na, nb, a.truth_threshold, a.call_min, len(cells),
                         len(both), len(only_a), len(only_b), neither,
                         f"{(len(both) + neither) / len(cells):.4f}" if cells else "NA",
                         f"{len(both) / union:.4f}" if union else "NA",
                         len(ea), len(eb), len(ea & eb)])

    lines = ["\t".join(str(x) for x in r) for r in rows]
    print("\t".join(COLS))
    print("\n".join(lines))
    if a.out:
        with open(a.out, "w") as f:
            f.write("\t".join(COLS) + "\n" + "\n".join(lines) + "\n")
        print(f"wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
