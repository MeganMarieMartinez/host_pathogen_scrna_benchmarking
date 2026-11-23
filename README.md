# Benchmarking Pathogen Detection Tools for Single-Cell RNA-Seq Data

Skye La Rosa & Megan Martinez 
## Overview

This project aims to benchmark four tools that detect and quantify pathogen reads in 10x single-cell RNA-seq data. We compare them on accuracy (sensitivity and false-positive rate), run time, and peak memory.

The RNA-seq data comes from Bradley et al. 2018, using HIV-infected and uninfected CD4 T cells (infected donor 1 = SRR6825024, uninfected donor 1 = SRR6825026), subsampled to 25 million read pairs each. 

| Tool | Approach in one line | Run by | Docs |
| --- | --- | --- | --- |
| DVsc | fastp trim, HISAT2 against host then virus index, Snakemake pipeline | Skye | `docs/DVsc_README.md` |
| PathogenTrack | STAR host removal, then Kraken2 k-mer classification of unmapped reads | Skye | `docs/PathogenTrack_README.md` |
| scPathoQuant | Cell Ranger unmapped reads, bbmap/bowtie2 to the pathogen, UMI counts | Megan | `docs/scpathoquant_README.md` |
| Viral-Track | STAR host+virus alignment, StringTie, UMI-tools | Megan | `docs/viral-track_README.md` |

# Repository Layout

```
README.md                  this file
data/datasets.tsv          candidate datasets and the metadata used to choose them
data/raw/                  two 200,000-read-pair example FASTQs (install/usage test input only)
data/viraltrack_files/     Viral-Track parameter and target files for the 25M runs
docs/                      one README per tool, plus cellranger_preprocessing.md
envs/                      exported conda environments (one per tool)
scripts/                   every script used for the results (see "Reproducing the Analysis")
results/timing.tsv         run time and memory of every timed step
results/OLD_timing.tsv     first 5M pilot timings (not used in the paper; see Reproducibility Notes)
results/scoring/           truth tables, per-tool call tables, scores, call-threshold sweep, tool agreement
dvsc_changes.patch         edits to upstream DVsc's analysis.snake (see docs/DVsc_README.md)
```

# Data Sources

| Data                                                                   | Source                                                                                                                                                                                                                                                                                                                          |
| ---------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| HIV scRNA-seq reads                                                    | Bradley et al. 2018, [SRP134859](https://www.ncbi.nlm.nih.gov/sra?term=SRP134859) / BioProject [PRJNA437838](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA437838). Infected donor 1 = [SRR6825024](https://www.ncbi.nlm.nih.gov/sra/SRR6825024), uninfected donor 1 = [SRR6825026](https://www.ncbi.nlm.nih.gov/sra/SRR6825026) |
| Human reference (PathogenTrack STAR index, host Cell Ranger baseline)  | 10x `refdata-gex-GRCh38-2024-A`, [10x downloads](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads)                                                                                                                                                                                        |
| Human + HIV-1 pNL4-3 reference (truth, DVsc indexes, Kraken2 database) | `refdata-cellranger-hg38-pEGFP-N1-HIVpNL4-3_HIVonly-4.0.0` <!-- TODO: say where this came from (scPathoQuant authors? a download link?) -->                                                                                                                                                                                     |
| EGFP FASTA/GTF added to that reference                                 | `egfp_ref/EGFP.fa`, `EGFP.gtf` <!-- TODO: say how these were made -->                                                                                                                                                                                                                                                           |
| NCBI taxonomy for Kraken2                                              | [taxdump.tar.gz](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)                                                                                                                                                                                                                                                      |

Chemistry: 10x v2. Read 1 is 26 bp (16 bp cell barcode + 10 bp UMI); read 2 is 98 bp cDNA.
## Downloading and Subsampling

The script `scripts/download.sh` documents how slices were made from each sample. The script makes one 25M slice per sample by combining a series of "chunks". The script works as follows:

1. `prefetch` downloads the `.sra` file.
2. The run is converted in chunks of 20,000,000 spots with `fastq-dump -N <start> -X <end> --split-files`.
3. Each chunk is subsampled with `seqtk sample -s <100 + chunk index>`, using the same seed for both mates so pairs stay in sync. The number kept per chunk is proportional to its size, so the chunks add up to exactly 25,000,000 pairs.
4. The chunks are concatenated into `<SRR>_25M_1.fastq` and `<SRR>_25M_2.fastq`, and the script checks the line count (4 × 25,000,000), that R1 and R2 read names match, and the read lengths (26 and 98).

The infected and uninfected slices are made the same way. Parameters can be changed through environment variables (`N`, `TAG`, `CHUNK`, `JOBS`, `SEED`, `SHARED`) which are explained in the script.

`scripts/make_slices.sh` makes 10M and 5M read-pair slices by sampling **from the 25M slice**
(10M from 25M, then 5M from 10M, `seqtk sample -s 100` on both mates). It also gzips them into
Cell Ranger-style names under `~/cellranger_work/fastq_<TAG>/` and links
`barcodes_<SRR>-<TAG>.tsv` to the 25M barcode list. The 5M and 10M runs therefore use
**Cell Ranger's 25M cell list**, and the 5M/10M slices are used for run time and memory only;
accuracy is scored at 25M.

# Software

| Software      | Version                                                                  | Source                                                             |
| ------------- | ------------------------------------------------------------------------ | ------------------------------------------------------------------ |
| DVsc          | commit `16e60b5eac7a469f80bed91308b99440613ae380` + `dvsc_changes.patch` | https://github.com/madisonglad/DVsc                                |
| PathogenTrack | 0.2.3 (PyPI)                                                             | https://github.com/ncrna/PathogenTrack                             |
| scPathoQuant  | <!-- TODO -->                                                            | https://github.com/galelab/scPathoQuant                            |
| Viral-Track   | <!-- TODO -->                                                            | https://github.com/PierreBSC/Viral-Track                           |
| Cell Ranger   | <!-- TODO: output of `cellranger --version` -->                          | https://www.10xgenomics.com/support/software/cell-ranger/downloads |
| STAR          | 2.7.10b                                                                  | https://github.com/alexdobin/STAR                                  |
| Kraken2       | 2.1.3                                                                    | https://github.com/DerrickWood/kraken2                             |
| fastp         | 0.22.0                                                                   | https://github.com/OpenGene/fastp                                  |
| HISAT2        | 2.2.3                                                                    | https://daehwankimlab.github.io/hisat2/                            |
| samtools      | 1.21                                                                     | https://www.htslib.org/                                            |
| Snakemake     | 9.27.0                                                                   | https://snakemake.github.io/                                       |
| UMI-tools     | 1.1.0                                                                    | https://github.com/CGATOxford/UMI-tools                            |
| seqtk         | 1.4-r122                                                                 | https://github.com/lh3/seqtk                                       |
| SRA Toolkit   | 3.0.3                                                                    | https://github.com/ncbi/sra-tools                                  |
# Installation Guides

1. **Clone the repo and create conda environments.** Separate conda environments were used for different tools. The environments used were exported using `scripts/export_envs.sh` and can easily be recreated:

```bash
git clone https://github.com/MeganMarieMartinez/host_pathogen_scrna_benchmarking.git
cd host_pathogen_scrna_benchmarking
conda env create -f envs/dvsc.yml
conda env create -f envs/umitools.yml
conda env create -f envs/pathogentrack.yml
conda env create -f envs/scpathoquant.yml
conda env create -f envs/viral-track.yml
```
2. **Install Cell Ranger.** See `docs/cellranger_preprocessing.md` for a more in-depth guide.
3. **Install the data tools.** `scripts/download.sh` needs `prefetch`, `srapath`, `sra-stat`, `fastq-dump` and `seqtk` on the `PATH`. <!-- TODO: say which env they live in, or add an envs/download.yml -->

4. **Per-tool setup:**
	- **DVsc**: `docs/DVsc_README.md`
	- **PathogenTrack**: `docs/PathogenTrack_README.md`
	- **scPathoQuant**: `docs/scpathoquant_README.md`
	- **Viral-Track**: `docs/viral-track_README.md`
# Truth Tables and Scoring

All tables have two columns, `barcode` and `hiv_umis`. Barcodes have the `-1` suffix removed.

`scripts/hiv_counts.py` turns each tool's native output into one common two-column table (`barcode`, `hiv_umis`), with one subcommand per tool:
`cellranger` (truth), `dvsc`, `pathogentrack`, `scpathoquant`, `viraltrack`. It needs `pandas` (installed in the `dvsc`, `umitools`, `scpathoquant` and `viral-track` envs, **not** in `pathogentrack`).
`scripts/score_calls.py` compares one tool's tables with the truth tables. `scripts/tool_agreement.py` compares two or more tools with each other.

### `scripts/hiv_counts.py`

```bash
# truth: HIV UMIs per cell from Cell Ranger
python scripts/hiv_counts.py cellranger \
  --matrix-dir ~/cellranger_work/runs_truth/SRR6825024-25M/outs/filtered_feature_bc_matrix \
  --gtf ~/cellranger_work/ref_build/hg38_HIV_EGFP/genes/genes.gtf.gz \
  --out results/scoring/truth_SRR6825024_25M.tsv
# (repeat for SRR6825026)

# DVsc calls (cells = the UMI-tools whitelist DVsc was run with)
python scripts/hiv_counts.py dvsc \
  --countbyumi ~/dvsc_work/runs_25M_rep1/HIV_inf1_25M/out/HIV_inf1_25M.countbyumi.filter.txt \
  --whitelist  ~/dvsc_work/preproc_25M_rep1/HIV_inf1_25M_whitelist.txt \
  --out results/scoring/calls_DVsc_SRR6825024_25M.tsv
# (repeat for HIV_uninf1_25M / SRR6825026)

# PathogenTrack calls
python scripts/hiv_counts.py pathogentrack \
  --microbes ~/pathogentrack_work/runs_25M/HIV_inf1_25M_rep1/microbes.tsv \
  --out results/scoring/calls_PathogenTrack_SRR6825024_25M.tsv
```
### `scripts/score_calls.py`

```bash
# score every tool against the truth (add one --tool block per tool)
python scripts/score_calls.py --tool PathogenTrack --label 25M_rep1 \
  --truth-inf   results/scoring/truth_SRR6825024_25M.tsv \
  --truth-uninf results/scoring/truth_SRR6825026_25M.tsv \
  --calls-inf   results/scoring/calls_PathogenTrack_SRR6825024_25M.tsv \
  --calls-uninf results/scoring/calls_PathogenTrack_SRR6825026_25M.tsv \
  --truth-thresholds 1 2 5 --call-min 1 --out results/scoring/scores.tsv

# call-threshold sweep -> scores_callsweep.tsv  (TEMPLATE: I could not find the command you used; replace with it)
for tool in DVsc PathogenTrack; do for c in 1 2 3 5; do
  python scripts/score_calls.py --tool $tool --label 25M --call-min $c \
    --truth-inf   results/scoring/truth_SRR6825024_25M.tsv \
    --truth-uninf results/scoring/truth_SRR6825026_25M.tsv \
    --calls-inf   results/scoring/calls_${tool}_SRR6825024_25M.tsv \
    --calls-uninf results/scoring/calls_${tool}_SRR6825026_25M.tsv \
    --out results/scoring/scores_callsweep.tsv
done; done
```
### `scripts/tool_agreement.py`

``` bash
# how much two tools agree -> agreement_25M.tsv
python scripts/tool_agreement.py \
  --truth-inf results/scoring/truth_SRR6825024_25M.tsv --truth-uninf results/scoring/truth_SRR6825026_25M.tsv \
  --tool DVsc results/scoring/calls_DVsc_SRR6825024_25M.tsv results/scoring/calls_DVsc_SRR6825026_25M.tsv \
  --tool PathogenTrack results/scoring/calls_PathogenTrack_SRR6825024_25M.tsv results/scoring/calls_PathogenTrack_SRR6825026_25M.tsv \
  --out results/scoring/agreement_25M.tsv
```
## How Truth is Defined + Metrics

The infected and uninfected libraries were both run through Cell Ranger with a human + HIV reference. Each cell's HIV UMI count from Cell Ranger is treated as the truth, and each tool's per-cell HIV UMI count is compared against it.

**Sensitivity:** Of the truth-positive cells in the infected sample, what fraction the tool also calls positive.
**False-Positive Rate (FPR):** Of the truth-negative cells in the uninfected sample, what the tools calls positive.
**Run Time and Peak Memory:** From `usr/bin/time -v` via `scripts/time_wrapper.sh`

- The truth for a cell is the number of HIV UMIs Cell Ranger counted for it: all genes on the `HIV-1_pNL4-3` contig of the combined reference (`hiv_counts.py cellranger`).
- A cell is truth-positive if it has at least `--truth-thresholds` HIV UMIs. The primary threshold is 1, and 2 and 5 are reported as a check.
- Sensitivity uses truth-positive cells in the infected library. FPR uses the uninfected library's cells whose truth count is below the threshold (so uninfected cells that Cell Ranger itself shows HIV in are not counted as false positives).
- `FP_strict` and `FPR_strict` count every called cell in the uninfected library, regardless of Cell Ranger's count.
- Scoring is on cells present in both the truth table and the tool's table (Cell Ranger's cell calls). Confidence intervals are 95% Wilson intervals.
- A tool call is positive when the tool reports at least `--call-min` HIV UMIs (1).

> [!NOTE] 
> Because truth is another tool's output, the scores measure agreement with Cell Ranger, not with a known ground truth. We therefore report sensitivity and FPR, not precision.

# Output Files

### `results/timing.tsv`

One row per timed step, appended by `scripts/time_wrapper.sh`.

| Column                     | Description                                                                                    |
| -------------------------- | ---------------------------------------------------------------------------------------------- |
| `tool`, `dataset`, `stage` | What was run (for example `DVsc`, `HIV_inf1_25M`, `dvsc_snakemake`)                            |
| `elapsed`, `elapsed_s`     | Wall-clock time (as printed, and in seconds)                                                   |
| `user_s`, `sys_s`          | CPU time                                                                                       |
| `max_rss_kb`               | Peak resident memory in kB                                                                     |
| `threads`, `rep`           | Threads given to the step and the replicate number                                             |
| `exit_status`              | 0 = success. Rows with a non-zero status are failed attempts and are excluded from the results |
| `load1_start`              | Machine load when the step started. A high value means a noisy measurement                     |
| `start_time`               | Start time of the step                                                                         |
Rows with `stage` = `cellranger_mkref`, `cellranger_count`, `star_index_build`, `kraken2_build` are **setup costs** shared by the whole workflow. They are reported separately from per-sample run times. `pathogentrack_smoke` rows use 1M read pairs and are not results.

> **Timing caveat:** the 25M `rep 1` runs of DVsc and PathogenTrack overlapped with each other and with Cell Ranger on the same server (compare `start_time` and `elapsed_s`; `load1_start` was 7 to 11 for several). The 5M runs (Nov 23) ran one after another and have a lower load. <!-- TODO: update after the REP=2 25M runs. -->

### `results/scoring/`

| File                         | Description                                                                                                                      |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `truth_<SRR>_25M.tsv`        | Cell Ranger HIV UMIs per cell                                                                                                    |
| `calls_<tool>_<SRR>_25M.tsv` | Tool HIV UMIs per cell                                                                                                           |
| `scores.tsv`                 | One row per tool and truth threshold: TP, FN, sensitivity, FP, TN, FPR, `FP_strict`, `FPR_strict`, with 95% CIs (`*_lo`, `*_hi`) |
| `scores_callsweep.tsv` | Same columns as `scores.tsv`, with the tool's minimum UMI count for a positive call (`call_min_umis`) set to 1, 2, 3 and 5 |
| `agreement_25M.tsv` | Pairwise tool agreement per sample: cells called by both (`both_pos`), only one tool, or neither; `agreement`, `jaccard`; and `err_a`, `err_b`, `err_both` (errors against the truth, and how many are shared) |

## Results map

<!-- TODO: Will eventually map to actual figures in the paper. -->

| Paper section | Figure/table | Script | Input |
| --- | --- | --- | --- |
| Accuracy | Sensitivity and FPR table | `scripts/score_calls.py` | `results/scoring/*.tsv`,`scripts/tool_agreement.py`,`results/scoring/scores_callsweep.tsv` |
| Accuracy | Sensitivity vs FPR across truth thresholds | <!-- TODO: scripts/make_figures.py --> | `results/scoring/scores.tsv` |
| Run time and memory | Time and peak memory table (setup separated) | <!-- TODO: scripts/make_figures.py --> | `results/timing.tsv` |
| Scaling | 5M vs 25M run time and memory | <!-- TODO --> | `results/timing.tsv` |

# Reproducibility Notes

- **Random Seeds**: Read slicing is random; `seqtk sample` with seed `100 + chunk index`. The tools themselves were run without random options.
- **Parameters:** All non-default options are recorded in each tool README (PathogenTrack: `--min_reads 10 --confidence 0.11 --thread 8`; Cell Ranger: `--chemistry SC3Pv2 --localcores 8 --localmem 64 --create-bam true`; DVsc: see its config).

## Paths and Manual Steps

Large files (FASTQ slices, references, etc.) exist in a shared directory `SHARED=/home/skyemegan`. Scripts read it from `SHARED` or otherwise hardcode this path.

`run_dvsc.sh` reads the shared reference from `REF=/home/skyemegan/ref_dvsc`.

 <!-- TODO: Absolute paths for Viral-Track/SPQ? -->

**The following were done by hand and not using a script:**
- Downloading Cell Ranger and the 10x human reference
- Obtaining the human + HIV reference and the EGFP FASTA/GTF <!-- TODO: describe how -->


# References

1. Whitmore LS, Tisoncik-Go J, Gale M Jr. scPathoQuant: a tool for efficient alignment and quantification of pathogen sequence reads from 10× single cell sequencing datasets. Bioinformatics. 2024 Apr 1;40(4):btae145. doi:10.1093/bioinformatics/btae145 
2. Bradley T, Ferrari G, Haynes BF, Margolis DM, Browne EP. Single-Cell Analysis of Quiescent HIV Infection Reveals Host Transcriptional Profiles that Regulate Proviral Latency. Cell Rep. 2018 Oct 2;25(1):107-117.e3. doi: 10.1016/j.celrep.2018.09.020. PMID: 30282021; PMCID: PMC6258175. 
3. Leng F, Mei S, Zhou X, et al. DVsc: An Automated Framework for Efficiently Detecting Viral Infection from Single-cell Transcriptomics Data. Genomics, Proteomics & Bioinformatics. 2024;22(2). doi:10.1093/gpbjnl/qzad007 
4. Zhang W, Xu X, Fu Z, Chen J, Chen S, Tan Y. PathogenTrack and Yeskit: tools for identifying intracellular pathogens from single-cell RNA-sequencing datasets as illustrated by application to COVID-19. Frontiers of Medicine. 2022;16(2):251-262. doi:10.1007/s11684-021-0915-9 
5. Bost P, Giladi A, Liu Y, et al. Host-Viral Infection Maps Reveal Signatures of Severe COVID-19 Patients. Cell. 2020;181(7):1475-1488.e12. doi:10.1016/j.cell.2020.05.006 
6. Ravindra NG, Alfajaro MM, Gasque V, Huston NC, Wan H, Szigeti-Buck K, Yasumoto Y, Greaney AM, Habet V, Chow RD, Chen JS, Wei J, Filler RB, Wang B, Wang G, Niklason LE, Montgomery RR, Eisenbarth SC, Chen S, Williams A, Iwasaki A, Horvath TL, Foxman EF, Pierce RW, Pyle AM, van Dijk D, Wilen CB. Single-cell longitudinal analysis of SARS-CoV-2 infection in human airway epithelium identifies target cells, alterations in gene expression, and cell state changes. PLoS Biol. 2021 Mar 17;19(3):e3001143. doi: 10.1371/journal.pbio.3001143. PMID: 33730024; PMCID: PMC8007021. 
7. Rosin NL, Jaffer A, Sinha S, Mulloy RP, Robinson C, Labit E, et al. SARS-CoV-2 Infection of Circulating Immune Cells is Not Responsible for Virus Dissemination in Severe COVID-19 Patients [SSRN Scholarly Paper] [Internet]. Rochester, NY: Social Science Research Network; 2021 [cited 2026 Sep 22]. Available from: https://papers.ssrn.com/abstract=3771665 doi:10.2139/ssrn.3771665