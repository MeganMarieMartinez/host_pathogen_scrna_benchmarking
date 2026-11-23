# Benchmarking Pathogen Detection Tools for Single-Cell RNA-Seq Data

Skye La Rosa & Megan Martinez 
## Overview

This project aims to benchmark four tools that detect and quantify pathogen reads in 10x single-cell RNA-seq data. We compare them on accuracy (sensitivity and false-positive rate), run time, and peak memory.

The RNA-seq data comes from Bradley et al. 2018, using HIV-infected and uninfected CD4 T cells (infected donor 1 = SRR6825024, uninfected donor 1 = SRR6825026), subsampled to 25 million read pairs each. 

# Repository Layout

```
README.md                  this file
data/datasets.tsv          candidate datasets and the metadata used to choose them
data/raw/                  two 200,000-read-pair example FASTQs (pilot / smoke-test input only)
data/processed/            small scPathoQuant example outputs
docs/                      one README per tool, plus cellranger_preprocessing.md
envs/                      exported conda environments (one per tool)
scripts/                   every script used for the results (see "Reproducing the analysis")
results/timing.tsv         run time and memory of every timed step (current, 25M)
results/OLD_timing.tsv     5M-read pilot timings (kept for reference, not used in the paper)
results/scoring/           truth tables, per-tool call tables and scores.tsv
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

1. **Clone the repo and create condo environments.** Separate conda environments were used for different tools. The environments used were exported using `scripts/export_envs.sh` and can easily be recreated:

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

<!-- TODO: Explain how the hiv_counts.py and score_calls.py scripts work. -->

```bash
# truth: HIV UMIs per cell from Cell Ranger
python scripts/hiv_counts.py cellranger \
  --matrix-dir ~/cellranger_work/runs_truth/SRR6825024-25M/outs/filtered_feature_bc_matrix \
  --gtf ~/cellranger_work/ref_build/hg38_HIV_EGFP/genes/genes.gtf.gz \
  --out results/scoring/truth_SRR6825024_25M.tsv
# (repeat for SRR6825026)

# DVsc calls (cells = the UMI-tools whitelist DVsc was run with)
python scripts/hiv_counts.py dvsc \
  --countbyumi ~/dvsc_work/HIV_inf1_25M/<sample>.countbyumi.filter.txt \
  --whitelist  ~/dvsc_work/preproc/HIV_inf1_25M_whitelist.txt \
  --out results/scoring/calls_DVsc_SRR6825024_25M.tsv
# (repeat for HIV_uninf1_25M / SRR6825026)

# PathogenTrack calls
python scripts/hiv_counts.py pathogentrack \
  --microbes ~/pathogentrack_work/runs_25M/HIV_inf1_25M_rep1/microbes.tsv \
  --out results/scoring/calls_PathogenTrack_SRR6825024_25M.tsv
```

<!-- TODO: check the DVsc --countbyumi path and file name against your ~/dvsc_work output. -->

```bash
python scripts/score_calls.py --tool DVsc --label 25M_rep1 \
  --truth-inf   results/scoring/truth_SRR6825024_25M.tsv \
  --truth-uninf results/scoring/truth_SRR6825026_25M.tsv \
  --calls-inf   results/scoring/calls_DVsc_SRR6825024_25M.tsv \
  --calls-uninf results/scoring/calls_DVsc_SRR6825026_25M.tsv \
  --truth-thresholds 1 2 5 --call-min 1 \
  --out results/scoring/scores.tsv
```
## How Truth is Defined + Metrics

The infected and uninfected libraries were both run through Cell Ranger with a human + HIV reference. Each cell's HIV UMI count from Cell Ranger is treated as the truth, and each tool's per-cell HIV UMI count is compared against it.

**Sensitivity:** Of the truth-positive cells in the infected sample, what fraction the tool also calls positive.
**False-Positive Rate (FPR):** Of the truth-negative cells in the uninfected sample, what fraction the tools calls positive.
**Run Time and Peak Memory:** From `usr/bin/time -v` via `scripts/time_wrapper.sh`

- The truth for a cell is the number of HIV UMIs Cell Ranger counted for it: all genes on the `HIV-1_pNL4-3` contig of the combined reference (`hiv_counts.py cellranger`). <!-- TODO: say whether the EGFP row is included. The code counts only the HIV-1_pNL4-3 contig; the project plan talks about EGFP+ cells. Make the text match the code. -->
- A cell is truth-positive if it has at least `--truth-thresholds` HIV UMIs. The primary threshold is 1, and 2 and 5 are reported as a check.
- Sensitivity uses truth-positive cells in the infected library. FPR uses the uninfected library's cells whose truth count is below the threshold (so uninfected cells that Cell Ranger itself shows HIV in are not counted as false positives).
- `FP_strict` and `FPR_strict` count every called cell in the uninfected library, regardless of Cell Ranger's count.
- Scoring is on cells present in both the truth table and the tool's table (Cell Ranger's cell calls). Confidence intervals are 95% Wilson intervals.
- A tool call is positive when the tool reports at least `--call-min` HIV UMIs (1).

> [!NOTE] 
> Because truth is another tool's output, the scores measure agreement with Cell Ranger, not with a known ground truth. We therefore report sensitivity and FPR, not precision.

# Output Files

### `results/timing.tsv`

Includes one row per timed step, appended by `scripts/time_wrapper.sh`.

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

> **Timing caveat:** the `rep = 1` timings were collected while other jobs shared the server (see `load1_start`, and note the overlapping `start_time` values).

### `results/scoring/`

| File                         | Description                                                                                                                      |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `truth_<SRR>_25M.tsv`        | Cell Ranger HIV UMIs per cell                                                                                                    |
| `calls_<tool>_<SRR>_25M.tsv` | Tool HIV UMIs per cell                                                                                                           |
| `scores.tsv`                 | One row per tool and truth threshold: TP, FN, sensitivity, FP, TN, FPR, `FP_strict`, `FPR_strict`, with 95% CIs (`*_lo`, `*_hi`) |
## Results map

<!-- TODO: Will eventually map to actual figures in the paper. -->

| Paper section | Figure/table | Script | Input |
| --- | --- | --- | --- |
| Accuracy | Sensitivity and FPR table | `scripts/score_calls.py` | `results/scoring/*.tsv` |
| Accuracy | Sensitivity vs FPR across truth thresholds | <!-- TODO: scripts/make_figures.py --> | `results/scoring/scores.tsv` |
| Run time and memory | Time and peak memory table (setup separated) | <!-- TODO: scripts/make_figures.py --> | `results/timing.tsv` |
| Scaling | 5M vs 25M run time and memory | <!-- TODO --> | `results/OLD_timing.tsv`, `results/timing.tsv` |

# Reproducibility Notes

- **Random Seeds**: Read slicing is random; `seqtk sample` with seed `100 + chunk index`. The tools themselves were run without random options.
- **Parameters:** All non-default options are recorded in each tool README (PathogenTrack: `--min_reads 10 --confidence 0.11 --thread 8`; Cell Ranger: `--chemistry SC3Pv2 --localcores 8 --localmem 64 --create-bam true`; DVsc: see its config).

## Paths and Manual Steps

Large files (FASTQ slices, references, etc.) exist in a shared directory `SHARED=/home/skyemegan`. Scripts read it from `SHARED` or otherwise hardcode this path.

`scripts/copy_hiv_config.sh` writes absolute `/home/slarosa/...` paths into the DVsc config.

**The following were done by hand and not using a script:**
- Downloading Cell Ranger and the 10x human reference
- Obtaining the human + HIV reference and the EGFP FASTA/GTF <!-- TODO: describe how -->
- Creating `hiv_uninf1.config` from `hiv_inf1.config` with `sed` (shown in `docs/DVsc_README.md`).

# References

1. Whitmore LS, Tisoncik-Go J, Gale M Jr. scPathoQuant: a tool for efficient alignment and quantification of pathogen sequence reads from 10× single cell sequencing datasets. Bioinformatics. 2024 Apr 1;40(4):btae145. doi:10.1093/bioinformatics/btae145 
2. Bradley T, Ferrari G, Haynes BF, Margolis DM, Browne EP. Single-Cell Analysis of Quiescent HIV Infection Reveals Host Transcriptional Profiles that Regulate Proviral Latency. Cell Rep. 2018 Oct 2;25(1):107-117.e3. doi: 10.1016/j.celrep.2018.09.020. PMID: 30282021; PMCID: PMC6258175. 
3. Leng F, Mei S, Zhou X, et al. DVsc: An Automated Framework for Efficiently Detecting Viral Infection from Single-cell Transcriptomics Data. Genomics, Proteomics & Bioinformatics. 2024;22(2). doi:10.1093/gpbjnl/qzad007 
4. Zhang W, Xu X, Fu Z, Chen J, Chen S, Tan Y. PathogenTrack and Yeskit: tools for identifying intracellular pathogens from single-cell RNA-sequencing datasets as illustrated by application to COVID-19. Frontiers of Medicine. 2022;16(2):251-262. doi:10.1007/s11684-021-0915-9 
5. Bost P, Giladi A, Liu Y, et al. Host-Viral Infection Maps Reveal Signatures of Severe COVID-19 Patients. Cell. 2020;181(7):1475-1488.e12. doi:10.1016/j.cell.2020.05.006 
6. Ravindra NG, Alfajaro MM, Gasque V, Huston NC, Wan H, Szigeti-Buck K, Yasumoto Y, Greaney AM, Habet V, Chow RD, Chen JS, Wei J, Filler RB, Wang B, Wang G, Niklason LE, Montgomery RR, Eisenbarth SC, Chen S, Williams A, Iwasaki A, Horvath TL, Foxman EF, Pierce RW, Pyle AM, van Dijk D, Wilen CB. Single-cell longitudinal analysis of SARS-CoV-2 infection in human airway epithelium identifies target cells, alterations in gene expression, and cell state changes. PLoS Biol. 2021 Mar 17;19(3):e3001143. doi: 10.1371/journal.pbio.3001143. PMID: 33730024; PMCID: PMC8007021. 
7. Rosin NL, Jaffer A, Sinha S, Mulloy RP, Robinson C, Labit E, et al. SARS-CoV-2 Infection of Circulating Immune Cells is Not Responsible for Virus Dissemination in Severe COVID-19 Patients [SSRN Scholarly Paper] [Internet]. Rochester, NY: Social Science Research Network; 2021 [cited 2026 Sep 22]. Available from: https://papers.ssrn.com/abstract=3771665 doi:10.2139/ssrn.3771665