- Official GitHub repo: [github.com/ncrna/PathogenTrack](https://github.com/ncrna/PathogenTrack)
- Quick start guide: [PathogenTrack/blob/master/doc/QUICK_START.md](https://github.com/ncrna/PathogenTrack/blob/master/doc/QUICK_START.md)
- Also installable via Bioconda: [bioconda.github.io/recipes/pathogentrack](https://bioconda.github.io/recipes/pathogentrack/README.html)
- Zhang W, Xu X, Fu Z, Chen J, Chen S, Tan Y. PathogenTrack and Yeskit: tools for identifying intracellular pathogens from single-cell RNA-sequencing datasets as illustrated by application to COVID-19. Frontiers of Medicine. 2022;16(2):251-262. doi:10.1007/s11684-021-0915-9

# What PathogenTrack does

For each cell barcode on a list you supply, PathogenTrack:

1. **umi_tools extract**: reads the 16 bp cell barcode and 10 bp UMI from read 1, keeps only barcodes on the list, and moves them into the read 2 name.
2. **fastp**: removes low-quality and low-complexity read 2 reads.
3. **STAR**: aligns read 2 to the human index. Reads that map to human are discarded.
4. **Kraken2**: classifies the unmapped reads against a database (here, HIV only).
5. Removes duplicates by UMI and writes a taxon-by-cell UMI count matrix (`microbes.tsv`).

It needs raw FASTQs, a host STAR index, a Kraken2 database and a list of valid cell barcodes. We use Cell Ranger's cell list
# Installation

1. Clone this repo:

```bash
git clone https://github.com/MeganMarieMartinez/host_pathogen_scrna_benchmarking.git
cd host_pathogen_scrna_benchmarking
```

2. Create the PathogenTrack and UMI-tools environments (exported with `scripts/export_envs.sh`):

```bash
conda env create -f envs/pathogentrack.yml   # PathogenTrack 0.2.3, STAR 2.7.10b, Kraken2 2.1.3, fastp 0.22.0
conda env create -f envs/umitools.yml        # UMI-tools 1.1.0
```

3. PathogenTrack calls `umi_tools`, which is in its own environment. `scripts/run_pathogentrack.sh` adds it to the `PATH` itself. When running PathogenTrack by hand, do:

```bash
export PATH=$HOME/miniconda3/envs/umitools/bin:$PATH
```
# Running PathogenTrack on the HIV data

The workflow has four steps. Steps 1 through 5 are scripts in `scripts/` and write timing rows to `results/timing.tsv`.

**Inputs needed before starting:** the 25M FASTQ slices (`scripts/download.sh`), the human reference `refdata-gex-GRCh38-2024-A`, and the HIV contig from the human + HIV reference (see the main README for the sources).

### 1. Cell barcode lists and renamed FASTQs (Cell Ranger)

```bash
bash scripts/cellranger_barcodes.sh
```

This builds the Cell Ranger reference (human + HIV + EGFP), renames the FASTQs to Cell Ranger's naming (`<SRR>_S1_L001_R1_001.fastq.gz`) in `~/cellranger_work/fastq_25M/`, runs `cellranger count`, and writes one barcode per line, with the `-1` suffix removed, to `~/cellranger_work/barcodes_<SRR>-25M.tsv`. These barcode files are PathogenTrack's `--barcode` input.

Cell Ranger counts for the same libraries are also the truth (see the main README).

### 2. References (run once)

```bash
conda activate pathogentrack
bash scripts/pathogentrack_setup.sh
```

- **Human STAR index:** built from `refdata-gex-GRCh38-2024-A` (FASTA and GTF) with `--sjdbOverhang 100`. It needs about 60 GB of RAM (see the resource table).
- **HIV-only Kraken2 database:** the `HIV-1_pNL4-3` contig, with its header rewritten to `>HIV-1_pNL4-3|kraken:taxid|11676`, plus the NCBI taxonomy files `names.dmp` and `nodes.dmp`. Without internet, download `taxdump.tar.gz` elsewhere and pass `TAXDUMP=/path/to/taxdump.tar.gz`.

Parameters can be overridden with environment variables (`HOST_REF`, `HIV_SRC`, `HIV_CONTIG`, `HIV_TAXID`, `WORK`, `THREADS`).

### 3. (Optional) Smoke test (1M read pairs)

```bash
SMOKE=1 bash scripts/run_pathogentrack.sh
```

This runs the first 1,000,000 read pairs of each sample to check that every step works before the full run. If the Cell Ranger barcode list is not there yet, it uses the UMI-tools whitelist instead. Smoke-test timings are labelled `pathogentrack_smoke` and are not results.

### 4. Full run

Ran inside `tmux`:
```bash
bash scripts/run_pathogentrack.sh              # both samples, REP=1
REP=2 bash scripts/run_pathogentrack.sh        # a second replicate for timing
ONLY=SRR6825024 bash scripts/run_pathogentrack.sh   # one sample only
```

The script runs this command for each sample (inside `time_wrapper.sh`):

```bash
PathogenTrack count --project_id HIV_inf1_25M_rep1 \
  --pattern CCCCCCCCCCCCCCCCNNNNNNNNNN \
  --thread 8 --min_reads 10 --confidence 0.11 \
  --star_index ~/pathogentrack_work/refs/star_hg38 \
  --kraken_db  ~/pathogentrack_work/refs/kraken_hiv \
  --barcode ~/cellranger_work/barcodes_SRR6825024-25M.tsv \
  --read1 ~/cellranger_work/fastq_25M/SRR6825024_S1_L001_R1_001.fastq.gz \
  --read2 ~/cellranger_work/fastq_25M/SRR6825024_S1_L001_R2_001.fastq.gz
```

| Option | Value | Why |
| --- | --- | --- |
| `--pattern` | `CCCCCCCCCCCCCCCCNNNNNNNNNN` | 10x v2: 16 bp barcode + 10 bp UMI. A v3 dataset (12 bp UMI) needs a different pattern |
| `--thread` | 8 | Passed to fastp, STAR and Kraken2. UMI extraction and counting are single-threaded. Same thread count as Cell Ranger |
| `--min_reads` | 10 | <!-- TODO: say what this means, from the PathogenTrack docs, and whether it is the default --> |
| `--confidence` | 0.11 | Kraken2 confidence threshold <!-- TODO: state whether this is the default / the paper's value --> |

# Output

Each run writes a folder `~/pathogentrack_work/runs_25M/<dataset>_rep<REP>/` containing:

| File                                                 | Content                                                                     |
| ---------------------------------------------------- | --------------------------------------------------------------------------- |
| `microbes.tsv`                                       | **Main result.** One row per taxon, one column per cell barcode, UMI counts |
| `classified_report.txt`, `classified_table.txt`      | Kraken2 report and per-read classifications                                 |
| `reads_barcoded.fq.gz`, `reads_barcoded_clean.fq.gz` | Reads after barcode extraction, and after fastp                             |
| `unmapped_reads.fq`                                  | Reads that did not align to human (the Kraken2 input)                       |

Example (first lines of `~/HIV_inf1_25M_rep1/microbes.tsv`): 

```
gene    AAACCTGAGAATAGGG    AAACCTGAGACACGAC    
AAACCTGAGAGTCTGG    AAACCTGAGCCTTGAT    
AAACCTGAGTACATGA    AAACCTGAGTGTTTGC
```

The script prints the number of cells with at least one UMI per taxon at the end of each run. To turn `microbes.tsv` into the per-cell table used for scoring:

```bash
python scripts/hiv_counts.py pathogentrack \
  --microbes ~/pathogentrack_work/runs_25M/HIV_inf1_25M_rep1/microbes.tsv \
  --out results/scoring/calls_PathogenTrack_SRR6825024_25M.tsv
```

# Results

Run time and memory (`results/timing.tsv`, 8 threads, `rep` 1):

| Step                                 | Dataset                | Elapsed  | Peak memory (kB) |
| ------------------------------------ | ---------------------- | -------- | ---------------- |
| STAR index build (setup, once)       | `star_hg38`            | 1:05:35  | 60,506,884       |
| Kraken2 database build (setup, once) | `kraken_hiv`           | 0:08.99  | 968,960          |
| Smoke test (1M pairs, not a result)  | `smoke_HIV_inf1_25M`   | 2:01.55  | 29,850,012       |
| Smoke test (1M pairs, not a result)  | `smoke_HIV_uninf1_25M` | 1:56.15  | 29,849,248       |
| Full run, 25M                        | `HIV_inf1_25M`         | 34:35.24 | 29,906,164       |
| Full run, 25M                        | `HIV_uninf1_25M`       | 32:33.47 | 29,916,064       |
| Full run, 5M                         | `HIV_inf1_5M`          | 9:07.16  | 29,864,328       |
| Full run, 5M                         | `HIV_uninf1_5M`        | 7:06.77  | 29,862,740       |

The peak memory is dominated by loading the STAR human index, which is why it is almost the same for 1M and 25M read pairs. Cell Ranger's own run time (needed for the barcode list) is not included; see `results/timing.tsv` (`cellranger_count`).

Accuracy against the Cell Ranger truth is in `results/scoring/scores.tsv`:<!-- TODO: add the PathogenTrack row summary (sensitivity, FPR) after scoring -->

# Caveats

- **HIV-only Kraken2 database.** The database contains only HIV, so every non-human read is classified against it alone. Results can differ from a run with a full microbial database. <!-- TODO: finalise with Megan which database goes in the paper. -->
- **Minimum UMI rule.** <!-- TODO: describe exactly how PathogenTrack decides a taxon is present in a cell (the minimum UMI rule you found in PathogenTrack.py), since it sets the lowest count a cell can have and so affects sensitivity at low UMI counts. -->
- **Barcodes come from Cell Ranger.** PathogenTrack only sees the cells on the list, so its cell universe is Cell Ranger's.
- **Timings overlap with other jobs** on the same server for `rep` 1 (see `load1_start`). <!-- TODO: update after the REP=2 run. -->