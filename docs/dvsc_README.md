1. DVsc cloned and installed in own conda environment (dvsc)
2. HIV Infected 1 and Uninfected 2: 5M Read-pair slices
3. UMI tools barcode/UMI extraction, Infected 1 and 2
4. Human and HIV HISAT2 indexes + annotation file
5. 100k-read smoke test
6. 5M-read run on infected 1 with timing
7. Uninfected 1 run

# Changes made to upstream DVsc

# Scripts and the order to run them:
1. download.sh
2. setup_dvsc.sh
3. preprocess_umi_tools.sh (umitools)
4. build_dvsc_index.sh (dvsc)
5. update_dvsc_anno.sh
6. copy_hiv_config.sh
7. time_wrapper.sh (dvsc)
8. export_envs.sh (base)

# Running DVsc

Run inside tmux (jobs survive an SSH disconnect), with the dvsc env active, from ~/dvsc_work:

```bash
conda activate dvsc
cd ~/dvsc_work

# 1. dry run: lists the 10 jobs, nothing executes
snakemake -s ~/DVsc/analysis.snake --configfile hiv_inf1.config -n

# 2. mini test: first 100,000 reads
zcat preproc/HIV_inf1_R2_extracted.fastq.gz | head -n 400000 | gzip > preproc/mini.fastq.gz
sed 's#HIV_inf1_R2_extracted#mini#; s#HIV_inf1_5M#mini#' hiv_inf1.config > mini.config
snakemake -s ~/DVsc/analysis.snake --configfile mini.config --cores 4

# 3. full 5M run with timing
bash ~/host_pathogen_scrna_benchmarking/scripts/time_wrapper.sh DVsc HIV_inf1_5M \
  snakemake -s ~/DVsc/analysis.snake --configfile hiv_inf1.config --cores 4

```

Outputs are in ~/dvsc_work/HIV_inf1_5M/

