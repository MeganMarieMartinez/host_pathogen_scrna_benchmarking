# Benchmarking Pathogen Detection Tools for Single-Cell RNA-Seq Data

Skye la Rosa & Megan Martinez

## Overview


## Data Sources


### Extracting and Subsetting Reads

Reads were extracted via fasterq-dump (3.0.3)

```
fasterq-dump SRRXXXXXX.sra --split-files
```

Samples were subsetted to 200,000 reads via seqtk (1.4-r122)

```
seqtk sample -s100 SRRXXXXXX_1.fastq 200000 > SRRXXXXXX-200k_S1_L001_R1_0001.fastq
seqtk sample -s100 SRRXXXXXX_2.fastq 200000 > SRRXXXXXX-200k_S1_L001_R2_0001.fastq
```


## Software Sources
