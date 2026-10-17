# Benchmarking Pathogen Detection Tools for Single-Cell RNA-Seq Data

Skye la Rosa & Megan Martinez

## Overview


## Data Sources


### Extracting and Subsetting Reads

Reads were extracted via fasterq-dump (3.0.3)

```
fasterq-dump SRRXXXXXX.sra --split-files
```

Samples were subsetted to 5 million reads using seqtk (1.4-r122)

```
seqtk sample -s100 SRRXXXXXX_1.fastq 5000000 > SRRXXXXXX_5M_1.fastq
seqtk sample -s100 SRRXXXXXX_2.fastq 5000000 > SRRXXXXXX_5M_2.fastq
```


## Software Sources
