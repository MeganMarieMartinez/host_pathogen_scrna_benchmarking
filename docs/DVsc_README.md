- Official code repository (NGDC BioCode): [ngdc.cncb.ac.cn/biocode/tools/BT007367](https://ngdc.cncb.ac.cn/biocode/tools/BT007367)
- GitHub mirror: [github.com/madisonglad/DVsc](https://github.com/madisonglad/DVsc)
- Leng F, Mei S, Zhou X, et al. DVsc: An Automated Framework for Efficiently Detecting Viral Infection from Single-cell Transcriptomics Data. Genomics, Proteomics & Bioinformatics. 2024;22(2). doi:10.1093/gpbjnl/qzad007
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
3. **`count_by_umi` rule: pass the cell-barcode whitelist.** `src/count_by_umi.py` needs the whitelist (the file made by `umi_tools whitelist`) in scRNA-seq mode, but the rule passed none (`whitelistfile: None` → `TypeError`). The whitelist path is in the config (`"whitelist"`, written by `copy_hiv_config.sh`) and is passed to the script in the rule (`w {config[whitelist]}`)
# Scripts and how they were run:

All scripts used are in `scripts`. They require the environment shown to be activated `conda activate {env_name}` before running.

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

Ran inside `tmux` (jobs survive an SSH disconnect), with the `dvsc` env active, from `~/dvsc_work`:

```bash
conda activate dvsc
cd ~/dvsc_work

bash ~/host_pathogen_scrna_benchmarking/scripts/time_wrapper.sh DVsc HIV_inf1_5M \
  snakemake -s ~/DVsc/analysis.snake --configfile hiv_inf1.config --cores 4

```

Outputs are in `~/dvsc_work/HIV_inf1_5M/`

For the uninfected sample, a new config was made first:

```bash
sed 's#HIV_inf1#HIV_uninf1#g' hiv_inf1.config > hiv_uninf1.config
grep -E "fastq|whitelist|sample|outpath" hiv_uninf1.config
```
And then the full 5M run with timing was done again, but with the uninfected sample this time:

```bash
bash ~/host_pathogen_scrna_benchmarking/scripts/time_wrapper.sh DVsc HIV_uninf1_5M \
  snakemake -s ~/DVsc/analysis.snake --configfile hiv_uninf1.config --cores 4
```

# Results

| Dataset       | Elapsed | Peak memory (kB) | HIV reads mapped | Cells with ≥1 HIV UMI   | Total HIV UMIs |
| ------------- | ------- | ---------------- | ---------------- | ----------------------- | -------------- |
| HIV_inf1_5M   | 5:46    | 4,501,760        | 37,367           | 8,299 of 12,616 (65.8%) | 36,550         |
| HIV_uninf1_5M | 4:40    | 4,503,680        | 668              | 237 of 13,738 (1.7%)*   | 644            |
\*The uninfected sample shows the presence of HIV infected cells, more is being done to look into this at the moment
