#!/usr/bin/env bash

forward_read=$1
reverse_read=$2

mkdir ViralTrackResults
cd ViralTrackResults

# changing the naming conventions to rename output files later

# keeping just the filename, not the rest of the path 
fwd_base=$(basename "$forward_read")

# extracting just the sample name (SRR) and the suffix 
sample=${fwd_base%%_S*} # splits at the _S 
suffix=${fwd_base#"$sample"}

# Build output names
out_r1="${sample}_extracted${suffix}"                  # output names for umi_tools extract, R1
r2_in="$(dirname "$forward_read")/${fwd_base/_R1_/_R2_}"   # R2 input in same dir as R1
out_r2="${out_r1/_R1_/_R2_}"                           # output names for umi_tools extract, R2

## running UMI tools ##

# identify the correct cell barcodes
umi_tools whitelist --stdin "$forward_read" \
                    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                    --log2stderr > "${sample}_whitelist.txt"

# extract barcodes and UMIS -> add to readnames
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin "$forward_read" \
                  --stdout "$out_r1" \
                  --read2-in "$r2_in" \
                  --read2-out="$out_r2" \
                  --whitelist "${sample}_whitelist.txt"

mkdir ../human_ref_genome/viraltrack_star_index

# running STAR
STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir ../human_ref_genome/viraltrack_star_index \
--genomeFastaFiles /home/mmartinez30/ Whitmore_scPathoQuant_testSets/Bradley_et_al_testSet/genomefiles/pNL4-3renamed_HIVonly.fa /home/skyemegan/refdata-gex-GRCh38-2024-A/fasta/genome.fa

wait

# finding the parameter file
parameter_file="${sample}_parameters.txt$" # this is made outside of this script

# makiung the target file
target_file="${sample}_targetfile.txt$"

echo "./$out_r2" >> $target_file # i'm not sure if this is providing the full path this might be problematic

Rscript edited_Viral_Track_scanning.R ./Viral-Trackfiles/$parameter_file ./Viral-Trackfiles/$target_file

wait

Rscript edited_VT_demux.R ./Viral-Trackfiles/$parameter_file ./Viral-Trackfiles/$target_file

wait