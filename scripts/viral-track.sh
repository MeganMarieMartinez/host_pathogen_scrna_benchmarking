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

# STAR was already run, so to save time I will not include it in this script, but write it into another script for future use.


# finding the parameter file
parameter_file="${sample}_parameters.txt" # this is made outside of this script

# makiung the target file
target_file="${sample}_targetfile.txt"

echo "./Viral-Trackfiles/$out_r2" >> $target_file

Rscript edited_Viral_Track_scanning.R ./Viral-Trackfiles/$parameter_file ./Viral-Trackfiles/$target_file

wait

Rscript edited_VT_demux.R ./Viral-Trackfiles/$parameter_file ./Viral-Trackfiles/$target_file

wait

echo "Viral-Track analysis done for sample $sample"