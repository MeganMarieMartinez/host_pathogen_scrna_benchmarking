# Cellranger Pre-Processing

## Installing CellRanger
1. Install and unpack cell ranger tar file: https://www.10xgenomics.com/support/software/cell-ranger/downloads#download-links. Unpack tar file in a convenient directory
```
# all files will be unpacked in a directory named cellranger_setup
cd /cellranger_setup
tar -xzvf cellranger-x.y.z.tar.gz
```

2. Install and unpack Human reference (GRCh38) - 2024-A: https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads
```
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz
```

3. Export to PATH
```
export PATH=/cellranger_setup/cellranger-x.y.z:$PATH
```

## Running CellRanger Count
To create cellranger output, run the cellranger_alignment.py file in the ./scripts directory

```
# assuming we are running within cellranger_setup directory
python cellranger_alignment.py -i path/to/fastq/files -g ./refdata-gex-GRCh38-2020-A
```
# Truth Dataset - Cell Ranger Count
Path to reference genome (human + HIV with integrated EGFP)
```
/home/skyemegan/ref_genomes_w_pathogen/refdata-cellranger-hg38-pEGFP-N1-HIVpNL4-3_HIVonly-4.0.0
```
