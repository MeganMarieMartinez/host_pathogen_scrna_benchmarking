# Benchmarking Pathogen Detection Tools for Single-Cell RNA-Seq Data

Skye La Rosa & Megan Martinez

## Overview


# Data Sources

## HIV

**HIV Dataset** from Bradley et al. 2018 (SRP134859).
 Read 1 is 26 bp (16 bp cell barcode + 10 bp UMI); read 2 is 98 bp cDNA.
 Infected donor 1 = SRR6825024; uninfected donor 1 = SRR6825026

### Extracting and Subsetting Reads

Slices of 5 million read pairs per sample were extracted and stored in a shared directory. The script `scripts/download.sh` documents how these slices were made. This script is for documentation purposes and should not be re-run (currently)

**Note:** The infected and uninfected slices were made differently. The infected slices were made with a random subset of 5 million reads whereas the the uninfected slices were made with the first 5 million reads. (Possible TODO: re-slice the uninfected sample)

Two 200,000 read pair slices of the uninfected sample are provided in `data/raw` as an example.
Reads were extracted via fasterq-dump (3.0.3)

```
fasterq-dump SRRXXXXXX.sra --split-files
```

Samples were subsetted to 200,000 reads via seqtk (1.4-r122)

```
seqtk sample -s100 SRRXXXXXX_1.fastq 200000 > SRRXXXXXX-200k_S1_L001_R1_0001.fastq
seqtk sample -s100 SRRXXXXXX_2.fastq 200000 > SRRXXXXXX-200k_S1_L001_R2_0001.fastq
```

# Installation Guides

Separate conda environments were used for different tools. The environments used were exported using `scripts/export_envs.sh` and can easily be recreated:

```bash
git clone https://github.com/MeganMarieMartinez/host_pathogen_scrna_benchmarking.git
cd host_pathogen_scrna_benchmarking
conda env create -f envs/pathogentrack.yml
conda env create -f envs/dvsc.yml
conda env create -f envs/umitools.yml
conda env list
```

## DVsc

Source: https://github.com/madisonglad/DVsc (commit used: `16e60b5eac7a469f80bed91308b99440613ae380`)

For an in-depth installation and usage guide, see `docs/DVsc_README.md`

Several changes were made to DVsc's `analysis.snake`, and are saved in `dvsc_changes.patch`. These changes need to be applied before running DVsc:

```bash 
conda activate dvsc
git clone https://github.com/madisonglad/DVsc.git
cd DVsc
git apply ../host_pathogen_scrna_benchmarking/dvsc_changes.patch
snakemake -s analysis.snake --configfile test.config -n
```

## PathogenTrack

Source: https://github.com/ncrna/PathogenTrack (installed from PyPI)

```bash
conda activate pathogentrack
```
PathogenTrack calls `umi_tools`, which is in its own environment. Add it to your PATH before running:

```bash
export PATH=$HOME/miniconda3/envs/umitools/bin:$PATH
```
# References

1. Whitmore LS, Tisoncik-Go J, Gale M Jr. scPathoQuant: a tool for efficient alignment and quantification of pathogen sequence reads from 10× single cell sequencing datasets. Bioinformatics. 2024 Apr 1;40(4):btae145. doi:10.1093/bioinformatics/btae145 
2. Bradley T, Ferrari G, Haynes BF, Margolis DM, Browne EP. Single-Cell Analysis of Quiescent HIV Infection Reveals Host Transcriptional Profiles that Regulate Proviral Latency. Cell Rep. 2018 Oct 2;25(1):107-117.e3. doi: 10.1016/j.celrep.2018.09.020. PMID: 30282021; PMCID: PMC6258175. 
3. Leng F, Mei S, Zhou X, et al. DVsc: An Automated Framework for Efficiently Detecting Viral Infection from Single-cell Transcriptomics Data. Genomics, Proteomics & Bioinformatics. 2024;22(2). doi:10.1093/gpbjnl/qzad007 
4. Zhang W, Xu X, Fu Z, Chen J, Chen S, Tan Y. PathogenTrack and Yeskit: tools for identifying intracellular pathogens from single-cell RNA-sequencing datasets as illustrated by application to COVID-19. Frontiers of Medicine. 2022;16(2):251-262. doi:10.1007/s11684-021-0915-9 
5. Bost P, Giladi A, Liu Y, et al. Host-Viral Infection Maps Reveal Signatures of Severe COVID-19 Patients. Cell. 2020;181(7):1475-1488.e12. doi:10.1016/j.cell.2020.05.006 
6. Ravindra NG, Alfajaro MM, Gasque V, Huston NC, Wan H, Szigeti-Buck K, Yasumoto Y, Greaney AM, Habet V, Chow RD, Chen JS, Wei J, Filler RB, Wang B, Wang G, Niklason LE, Montgomery RR, Eisenbarth SC, Chen S, Williams A, Iwasaki A, Horvath TL, Foxman EF, Pierce RW, Pyle AM, van Dijk D, Wilen CB. Single-cell longitudinal analysis of SARS-CoV-2 infection in human airway epithelium identifies target cells, alterations in gene expression, and cell state changes. PLoS Biol. 2021 Mar 17;19(3):e3001143. doi: 10.1371/journal.pbio.3001143. PMID: 33730024; PMCID: PMC8007021. 
7. Rosin NL, Jaffer A, Sinha S, Mulloy RP, Robinson C, Labit E, et al. SARS-CoV-2 Infection of Circulating Immune Cells is Not Responsible for Virus Dissemination in Severe COVID-19 Patients [SSRN Scholarly Paper] [Internet]. Rochester, NY: Social Science Research Network; 2021 [cited 2026 Sep 22]. Available from: https://papers.ssrn.com/abstract=3771665 doi:10.2139/ssrn.3771665
## Software Sources
