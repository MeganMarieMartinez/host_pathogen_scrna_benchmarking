# Viral-Track Installation and Usage

Build Viral-Track environment
```
conda env create -f envs/viral-track.yml
```

## Creation of Index and Annotation File
### STAR Index
Create STAR index for both host and virus reference genomes:
```
mkdir /path/to/index
STAR --runThreadN N --runMode genomeGenerate --genomeDir /path/to/index --genomeFastaFiles /path/to/Virusite_file.fa  path/to/Host_genome.fa 
```

Command to create STAR index for our analysis:
```
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles data/hiv_genome/pNL4-3renamed_HIVonly.fa /home/skyemegan/refdata-gex-GRCh38-2024-A/fasta/genome.fa
```
### Annotation File
Viral-track requires an annotation file that lists all the viruses included in the index and their genome length. The annotation file including the EFGP-tagged HIV genome can be found at ```data/viraltrack_files/viral_track_annotation.txt```

## Data Preprocessing and Required Files
UMI Tools is used to begin demultiplexing the raw fastq data. Note: the cell barcode pattern argument below is specific to our 10x style read structure.

Step 1: Identify correct cell barcodes
```
umi_tools whitelist --stdin /path/to/forward_read \
                    --bc-pattern CCCCCCCCCCCCCCCCNNNNNNNNNN \
                    --log2stderr > whitelist.txt
```
Step 2: Extract barcdoes and UMIs and add to read names
```         
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin /path/to/forward_read \
                  --stdout /path/to/forward_read_extracted \
                  --read2-in /path/to/reverse_read \
                  --read2-out /path/to/reverse_read_extracted \
                  --whitelist whitelist.txt
```
A test-run of the preprocessing on the 25M subset of SRR6825024 and SRR6825026 can be found in ```scripts/viral_track_preprocessing.sh``` to serve as an example file. 

### Parameter file
Viral-Track requires a parameter file with the required format below:
```
N_thread=N
Output_directory="/path/to/output/directory"
Index_genome="/path/to/index/genome"
Viral_annotation_file="/path/to/viral/annotationfile"
Name_run="run name"
Load_STAR_module=FALSE
Load_samtools_module=FALSE
Load_stringtie_module=FALSE
Minimal_read_mapped=50
```
Example parameter files for 25M test-run can be found in ```data/viraltrack_files```

### Target file
Viral-Track requires a target file including the path of each reverse read resulting from the UMI tools preprocessing. Example target files for 25M test-run can be found in ```data/viraltrack_files```.

## Edits to Original Viral-Track Files
The original Viral-Track files were edited for usability. 

| Original Script | Edited Script | Edit Reason |
| --- | --- | --- |
| Viral_Track_scanning.R | edited_Viral_Track_scanning.R | Data structure mismatched: removed unecessary conditional check |
|  Viral_Track_cell_demultiplexing.R | edited_VT_demux.R | Improper gtf file creation: fixed overwriting issue |

## Running Viral-Track
Step 1: Scanning
```
 Rscript edited_Viral_Track_scanning.R /path/to/parameters.txt /path/to/target_file.txt
```
Step 2: Demultiplexing
```
 Rscript edited_VT_demux.R /path/to/parameters.txt /path/to/target_file.txt
```

A test-run running both Viral-Track steps on the 25M subset of SRR6825024 and SRR6825026 can be found in ```scripts/run_viraltrack_testdata.sh``` to serve as an example file. 

## Output
For the purposes of this analysis, the only output of interest are:
* Expression_table.tsv: UMI table with associated pathogen counts