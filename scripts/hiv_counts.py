import argparse
import gzip
import os
import re
import sys


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def find_file(folder, stem):
    """folder/stem or folder/stem.gz, whichever exists."""
    for name in (stem, stem + ".gz"):
        p = os.path.join(folder, name)
        if os.path.exists(p):
            return p
    sys.exit(f"ERROR: neither {stem} nor {stem}.gz found in {folder}")


def norm(bc):
    return bc.strip().split("-")[0]


def write_table(path, barcodes, counts):
    with open(path, "w") as out:
        out.write("barcode\thiv_umis\n")
        for b in barcodes:
            out.write(f"{b}\t{counts.get(b, 0)}\n")


def cmd_cellranger(a):
    # 1. which genes are HIV: every gene_id on the HIV contig in the GTF the reference was built from
    hiv_ids = set()
    with opener(a.gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[0] != a.contig:
                continue
            m = re.search(r'gene_id "([^"]+)"', c[8])
            if m:
                hiv_ids.add(m.group(1))
    if not hiv_ids:
        sys.exit(f"ERROR: no genes on contig '{a.contig}' in {a.gtf} (check --contig and --gtf)")

    # 2. which matrix rows those are (features.tsv column 1 = gene_id, column 2 = gene name)
    feats = [l.rstrip("\n").split("\t") for l in opener(find_file(a.matrix_dir, "features.tsv"))]
    rows = {i + 1 for i, ft in enumerate(feats) if ft[0] in hiv_ids}
    names = [ft[1] if len(ft) > 1 else ft[0] for i, ft in enumerate(feats) if (i + 1) in rows]
    if not rows:
        sys.exit("ERROR: none of the HIV genes in the GTF are in features.tsv (wrong reference for this matrix?)")
    print(f"HIV genes found in matrix: {len(rows)} of {len(hiv_ids)} in the GTF: {', '.join(names)}", file=sys.stderr)

    # 3. barcodes (column order of the matrix)
    barcodes = [norm(l) for l in opener(find_file(a.matrix_dir, "barcodes.tsv"))]

    # 4. stream matrix.mtx (MatrixMarket coordinate: row col value, 1-based), keep HIV rows only
    counts = {}
    with opener(find_file(a.matrix_dir, "matrix.mtx")) as f:
        dims_seen = False
        for line in f:
            if line.startswith("%"):
                continue
            if not dims_seen:           # first non-comment line = "n_rows n_cols n_entries"
                nr, nc, _ = line.split()
                if int(nc) != len(barcodes) or int(nr) != len(feats):
                    sys.exit(f"ERROR: matrix is {nr}x{nc} but features/barcodes files give {len(feats)}x{len(barcodes)}")
                dims_seen = True
                continue
            r, c, v = line.split()
            if int(r) in rows:
                b = barcodes[int(c) - 1]
                counts[b] = counts.get(b, 0) + int(float(v))
    write_table(a.out, barcodes, counts)
    pos = sum(1 for b in barcodes if counts.get(b, 0) > 0)
    print(f"{a.out}: {len(barcodes)} barcodes, {pos} with >=1 HIV UMI, {sum(counts.values())} HIV UMIs", file=sys.stderr)


def cmd_pathogentrack(a):
    rx = re.compile(a.taxon_regex, re.I)
    with open(a.microbes) as f:
        header = f.readline().rstrip("\n").split("\t")
        barcodes = [norm(b) for b in header[1:]]
        counts, matched = {}, []
        for line in f:
            c = line.rstrip("\n").split("\t")
            if not rx.search(c[0]):
                continue
            matched.append(c[0])
            for b, v in zip(barcodes, c[1:]):
                counts[b] = counts.get(b, 0) + int(v)
    if not matched:
        print("No taxon matching the HIV pattern in microbes.tsv: every cell is a negative.", file=sys.stderr)
    else:
        print(f"Taxa summed: {', '.join(matched)}", file=sys.stderr)
    write_table(a.out, barcodes, counts)
    pos = sum(1 for b in barcodes if counts.get(b, 0) > 0)
    print(f"{a.out}: {len(barcodes)} barcodes, {pos} with >=1 HIV UMI, {sum(counts.values())} HIV UMIs", file=sys.stderr)


def cmd_dvsc(a):
    # DVsc's countbyumi file lists ONLY cells that have virus reads: virus <TAB> barcode <TAB> count.
    # The cell universe is therefore the UMI-tools whitelist DVsc was run with; every whitelist cell
    # that is not listed gets 0. (Otherwise the scorer would silently drop DVsc's negative cells.)
    universe = []
    with open(a.whitelist) as f:
        for line in f:
            if line.strip():
                universe.append(norm(line.split("\t")[0]))
    universe_set = set(universe)
    rx = re.compile(a.virus_regex, re.I)
    counts, other, outside, n_lines = {}, set(), 0, 0
    with open(a.countbyumi) as f:
        for line in f:
            c = line.split()
            if len(c) < 3:
                continue
            n_lines += 1
            if not rx.search(c[0]):
                other.add(c[0])
                continue
            b = norm(c[1])
            if b not in universe_set:
                outside += 1
            counts[b] = counts.get(b, 0) + int(float(c[2]))
    if other:
        print(f"WARNING: rows for other viruses ignored (do not match --virus-regex): {', '.join(sorted(other))}", file=sys.stderr)
    if outside:
        print(f"WARNING: {outside} rows have a barcode that is not in the whitelist {a.whitelist}", file=sys.stderr)
    barcodes = universe + sorted(b for b in counts if b not in universe_set)
    write_table(a.out, barcodes, counts)
    pos = sum(1 for b in barcodes if counts.get(b, 0) > 0)
    print(f"{a.out}: {len(barcodes)} barcodes (whitelist {len(universe)}), {pos} with >=1 HIV UMI, "
          f"{sum(counts.values())} HIV UMIs (sum of column 3: check it is UMIs, not reads)", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("dvsc")
    d.add_argument("--countbyumi", required=True, help="DVsc countbyumi file (virus, barcode, count)")
    d.add_argument("--whitelist", required=True, help="the UMI-tools whitelist DVsc was run with (defines the cells)")
    d.add_argument("--virus-regex", default="immunodeficiency|HIV|pNL43")
    d.add_argument("--out", required=True)
    d.set_defaults(fn=cmd_dvsc)
    c = sub.add_parser("cellranger")
    c.add_argument("--matrix-dir", required=True, help="Cell Ranger outs/filtered_feature_bc_matrix")
    c.add_argument("--gtf", required=True, help="genes.gtf(.gz) of the reference the matrix was made with")
    c.add_argument("--contig", default="HIV-1_pNL4-3")
    c.add_argument("--out", required=True)
    c.set_defaults(fn=cmd_cellranger)
    t = sub.add_parser("pathogentrack")
    t.add_argument("--microbes", required=True, help="PathogenTrack microbes.tsv")
    t.add_argument("--taxon-regex", default="immunodeficiency|HIV")
    t.add_argument("--out", required=True)
    t.set_defaults(fn=cmd_pathogentrack)
    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()