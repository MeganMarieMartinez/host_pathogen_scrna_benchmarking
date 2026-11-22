#!/usr/bin/env bash

mkdir -p ViralTrack_Preprocessing
cd ViralTrack_Preprocessing

##### running analysis on infected hiv sample (SRR6825024-25M) ####

infected_r1="/home/skyemegan/formatted_SRR6825024_25M/SRR6825024-25M_S1_L001_R1_0001.fastq"
infected_r2="/home/skyemegan/formatted_SRR6825024_25M/SRR6825024-25M_S1_L001_R2_0001.fastq"
uninfected_r1="/home/skyemegan/formatted_SRR6825026_25M/SRR6825026-25M_S1_L001_R1_0001.fastq"
uninfected_r2="/home/skyemegan/formatted_SRR6825026_25M/SRR6825026-25M_S1_L001_R2_0001.fastq"


### umi tools preprocessing for infected samples ###
echo "UMI Tools preprocessing for SRR6825024_25M"

# identify the correct cell barcodes
umi_tools whitelist --stdin "$infected_r1" \
                    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                    --log2stderr > SRR6825024-25M_whitelist.txt;

# extract barcodes and UMIS -> add to readnames
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin "$infected_r1" \
                  --stdout=SRR6825024-25M_extracted_S1_L001_R1_0001.fastq \
                  --read2-in "$infected_r2" \
                  --read2-out=SRR6825024-25M_extracted_S1_L001_R2_0001.fastq  \
                  --whitelist=SRR6825024-25M_whitelist.txt;



### umi tools preprocessing for uninfected samples ###
echo "UMI Tools preprocessing for SRR6825026_25M"

umi_tools whitelist --stdin "$uninfected_r1" \
                    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                    --log2stderr > SRR6825024-25M_whitelist.txt;

# extract barcodes and UMIS -> add to readnames
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin "$uninfected_r1" \
                  --stdout=SRR6825026-25M_extracted_S1_L001_R1_0001.fastq \
                  --read2-in "$uninfected_r2" \
                  --read2-out=SRR6825026-25M_extracted_S1_L001_R1_0001.fastq  \
                  --whitelist=SRR6825026-25M_whitelist.txt;