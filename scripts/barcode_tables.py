"""Shared helpers for the scoring scripts: reading hiv_counts.py tables and checking barcode universes."""
import sys


def load(path):
    """Read a hiv_counts.py table into {barcode: hiv_umis}; the -1 suffix is stripped."""
    d = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if header[:2] != ["barcode", "hiv_umis"]:
            sys.exit(f"ERROR: {path} does not look like a hiv_counts.py table (header: {header})")
        for line in f:
            b, v = line.rstrip("\n").split("\t")[:2]
            d[b.split("-")[0]] = int(v)
    return d


def check_universes(items, mode, min_cov):
    """items: [(label, truth, tool), ...] with {barcode: umis} dicts. mode: error, warn or off.
    Tool-only barcodes are reported but never fail; they drop out of scoring by design."""
    if mode == "off":
        return
    failed = False
    for label, truth, tool in items:
        shared = len(set(truth) & set(tool))
        cov = shared / len(truth) if truth else 0.0
        n_zero = sum(1 for v in tool.values() if v == 0)
        print(f"barcode universe ({label}): truth {len(truth)}, tool {len(tool)}, shared {shared} "
              f"({cov:.1%} of truth), tool-only {len(tool) - shared}, tool rows with 0 UMIs {n_zero}",
              file=sys.stderr)
        problems = []
        if cov < min_cov:
            problems.append(f"tool table covers only {cov:.1%} of the {len(truth)} truth barcodes "
                            f"(minimum {min_cov:.0%}); the {len(truth) - shared} missing truth cells "
                            f"are dropped from scoring")
        if tool and n_zero == 0:
            problems.append("tool table has no zero-count rows, so it probably lists only positive "
                            "cells and its negatives are missing")
        for p in problems:
            print(f"{'ERROR' if mode == 'error' else 'WARNING'} ({label}): {p}", file=sys.stderr)
        failed = failed or bool(problems)
    if failed and mode == "error":
        sys.exit("Fix the tool table(s) (see hiv_counts.py) or rerun with --universe-check warn/off to override.")