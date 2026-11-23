- Official code repository (NGDC BioCode): [ngdc.cncb.ac.cn/biocode/tools/BT007367](https://ngdc.cncb.ac.cn/biocode/tools/BT007367)
- GitHub mirror: [github.com/madisonglad/DVsc](https://github.com/madisonglad/DVsc)
- Leng F, Mei S, Zhou X, et al. DVsc: An Automated Framework for Efficiently Detecting Viral Infection from Single-cell Transcriptomics Data. Genomics, Proteomics & Bioinformatics. 2024;22(2). doi:10.1093/gpbjnl/qzad007
# What DVsc does

DVsc is a Snakemake pipeline. From a UMI-extracted read 2 FASTQ and a cell whitelist it:

1. trims reads with **fastp**;
2. aligns them to a **human index with HISAT2**, and sets the aligned reads aside;
3. aligns the remaining reads to a **virus index with HISAT2** (here, the HIV contig `NC_pNL43`);
4. sorts and filters the alignments with **samtools** and the DVsc scripts in `src/`, builds transcripts, and counts viral reads and UMIs per cell.

Unlike PathogenTrack, DVsc does not take Cell Ranger's cell list. Its cells come from `umi_tools whitelist`, and the whitelist is intersected with Cell Ranger's cells at scoring time.
# Installation

1. Clone DVsc

```bash
git clone https://github.com/madisonglad/DVsc ~/DVsc
```
2. Create DVsc and UMI-tools environments (exported from `scripts/export_envs.sh):

```bash
conda env create -f envs/dvsc.yml
conda env create -f envs/umitools.yml
```
3. Apply patch to DVsc's `analysis.snake` (see next section for explanation):

```bash
cd DVsc
git apply ../host_pathogen_scrna_benchmarking/dvsc_changes.patch
```


4. Create work directory

```bash
bash scripts/setup_dvsc.sh
```

# Changes made to upstream DVsc

The following edits were made to DVsc's `analysis.snake` file and saved as `dvsc_changes.patch`:

1. `qc` **rule: declare `count_ref.txt` as an input.** `src/qc.py` reads `<sample>/count_ref.txt`, which is produced by `bamdeal`, but the rule did not list it. Snakemake therefore started `qc` at the same time as the alignment, and it crashed with `FileNotFoundError: mini/count_ref.txt`. Adding the input makes `qc` wait for `bamdeal`.
2. **`alignment_hisat2` rule: `-P {params.threads}` → `-p {params.threads}`.**
3. **`count_by_umi` rule: pass the cell-barcode whitelist.** `src/count_by_umi.py` needs the whitelist (the file made by `umi_tools whitelist`) in scRNA-seq mode, but the rule passed none (`whitelistfile: None` gave a `TypeError`). The whitelist path is in the config (`"whitelist"`, written by `copy_hiv_config.sh`) and is passed to the script in the rule (`w {config[whitelist]}`).
# Workflow and scripts:

All listed scripts used are in `scripts/`. They require the environment shown to be activated (`conda activate {env_name}`) before running.

| Script                    | Env      | Description                                                                                                                           |
| ------------------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `download.sh`             | n/a      | Records how the 5M-read slices were made                                                                                              |
| `setup_dvsc.sh`           | any      | Creates `~/dvsc_work`, which includes `preproc/`, `4/`, and  the `src/` symlink (to `src/` in original DVsc repo)                     |
| `preprocess_umi_tools.sh` | umitools | UMI-tools whitelist + extract for infected and uninfected samples                                                                     |
| `build_dvsc_index.sh`     | dvsc     | Splits the genome, renames HIV contig to `NC_pNL43`, builds the host index and the host+HIV index, tests HIV alignment. Ran in `tmux` |
| `update_dvsc_anno.sh`     | any      | Makes `anno_hiv.txt` from DVsc's `anno.txt`. Appends a tab-separated row for the HIV contig `NC_pNL43` (length taken from the `.fai`) |
| `copy_hiv_config.sh`      | any      | Writes `hiv_inf1.config`.                                                                                                             |
| `time_wrapper.sh`         | dvsc     | Wraps a command with `/usr/bin/time -v` and appends elapsed time and peak memory to `results/timing.tsv`                              |
| `export_envs.sh`          | base     | Exports the conda environments used, including `dvsc` and`umitools`, to `envs/*.yml` (--no-builds)`                                   |

# Running DVsc

Run inside `tmux` (jobs survive an SSH disconnect), with the `dvsc` environment active, from `~/dvsc_work`.

**1. Preprocess (umitools env):**

```bash
bash ~/host_pathogen_scrna_benchmarking/scripts/preprocess_umi_tools.sh   # SIZE=25M by default
```

**2. Make the configs.** `copy_hiv_config.sh` writes the infected config. The uninfected config is the same with the sample name changed:

```bash
bash ~/host_pathogen_scrna_benchmarking/scripts/copy_hiv_config.sh
cd ~/dvsc_work
sed 's#HIV_inf1#HIV_uninf1#g' hiv_inf1.config > hiv_uninf1.config
grep -E "fastq|whitelist|sample|outpath" hiv_uninf1.config    # should all say uninf1_25M
```

Key config values: `readlength: 98`, `method: scRNAseq`, `alignment: hisat2`, tool threads 4. The config also points to the host index, the host+HIV index and `anno_hiv.txt`.

**3. Run Snakemake with timing:**

```bash
conda activate dvsc
cd ~/dvsc_work
T=~/host_pathogen_scrna_benchmarking/scripts/time_wrapper.sh

STAGE=dvsc_snakemake THREADS=8 REP=1 bash $T DVsc HIV_inf1_25M \
  snakemake -s ~/DVsc/analysis.snake --configfile hiv_inf1.config --cores 8
STAGE=dvsc_snakemake THREADS=8 REP=1 bash $T DVsc HIV_uninf1_25M \
  snakemake -s ~/DVsc/analysis.snake --configfile hiv_uninf1.config --cores 8
```

<!-- TODO: check these match the commands you actually ran (cores, REP). -->

Delete `~/dvsc_work/HIV_inf1_25M/` before a timed re-run: if the outputs exist, Snakemake finishes in about a second.

## Output

Results for a sample are in `~/dvsc_work/<sample>/`:

| File | Content |
| --- | --- |
| `*.virus_count.filter.txt` | Reads per virus (filtered) |
| `*.countbyumi.filter.txt` | Per virus and cell barcode: virus, barcode, UMI count. Only cells with virus reads are listed |
| `*.gtf` | Assembled viral transcripts |
| `count_ref.txt` | Reads per reference from the alignment |

Example (first lines of `*.countbyumi.filter.txt`): <!-- TODO: paste 2 to 3 real lines from HIV_inf1_25M -->

To turn this into the per-cell table used for scoring, see "Step 8" in the main README (`hiv_counts.py dvsc`).

# Results

25M read pairs per sample (random slices). Time and memory are for the `dvsc_snakemake` step with 8 threads; preprocessing is listed separately. From `results/timing.tsv`.

| Dataset | Elapsed | Peak memory (kB) | Cells with ≥1 HIV UMI | Total HIV UMIs |
| --- | --- | --- | --- | --- |
| HIV_inf1_25M | 10:05.85 | 4,539,200 | 11,110 of 12,611 whitelist cells (88.1%) | 170,699 |
| HIV_uninf1_25M | 9:33.55 | 4,533,440 | 394 of 13,740 whitelist cells (2.9%) | 2,817 |

UMI-tools preprocessing (single thread):

| Dataset | Step | Elapsed | Peak memory (kB) |
| --- | --- | --- | --- |
| HIV_inf1_25M | whitelist | 16:44.76 | 461,752 |
| HIV_inf1_25M | extract | 37:02.20 | 89,280 |
| HIV_uninf1_25M | whitelist | 17:53.17 | 425,580 |
| HIV_uninf1_25M | extract | 29:08.72 | 89,376 |

<!-- TODO: add "HIV reads mapped" from count_ref.txt (the HIV row) for the two 25M samples. -->

**Accuracy** against the Cell Ranger truth (`results/scoring/scores.tsv`, truth threshold 1, tool call threshold 1):

| Sensitivity (infected) | FPR (uninfected, truth-negative cells) | FPR_strict (all uninfected cells) |
| --- | --- | --- |
| 0.9924 (9,809 of 9,884; 95% CI 0.9905 to 0.9939) | 0.0057 (70 of 12,307) | 0.0286 (361 of 12,602) |

At truth thresholds of 2 and 5, sensitivity is 1.000 and FPR rises to 0.0132 and 0.0211. The uninfected sample is not free of HIV reads: Cell Ranger also finds HIV UMIs in 295 of its 12,602 cells (2.3%), so some DVsc "false positives" are background that Cell Ranger sees too. That is why FPR counts only cells that are truth-negative, and `FPR_strict` is also reported.

**5M pilot** (random infected slice, first-5M-reads uninfected slice, `--cores 4`; not comparable to 25M because the uninfected slice was sampled differently):

| Dataset | Elapsed | Peak memory (kB) | HIV reads mapped | Cells with ≥1 HIV UMI | Total HIV UMIs |
| --- | --- | --- | --- | --- | --- |
| HIV_inf1_5M | 5:46 | 4,501,760 | 37,367 | 8,299 of 12,616 (65.8%) | 36,550 |
| HIV_uninf1_5M | 4:40 | 4,503,680 | 668 | 237 of 13,738 (1.7%) | 644 |

# Caveats

- `results/timing.tsv` contains two failed DVsc attempts (`exit_status` 1, about 1 to 3 seconds, at 16:06). Ignore them; the successful runs are the ones above. <!-- TODO: say in a few words what failed (or delete those two rows). -->
- The rep 1 timings were collected while other jobs ran on the same server. <!-- TODO: update after the REP=2 run. -->
- The HIV contig is renamed to `NC_pNL43` so DVsc's annotation lookup works. The virus index also contains the human genome, so a looser alignment score can pull some human reads into the HIV hits; the HIV row of `count_ref.txt` is the real signal.
