mkdir -p ~/DVsc/preproc && cd ~/DVsc/preproc
SH=/home/skyemegan/HIV_Infected_1_SRR6825024/
R1=$SH/SRR6825024_5M_1.fastq; R2=$SH/SRR6825024_5M_2.fastq

umi_tools whitelist --stdin $R1 --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
  --plot-prefix=HIV_inf1 --log2stderr > HIV_inf1_whitelist.txt

umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
  --stdin $R1 --stdout HIV_inf1_R1_extracted.fastq.gz \
  --read2-in $R2 --read2-out HIV_inf1_R2_extracted.fastq.gz \
  --filter-cell-barcode --whitelist=HIV_inf1_whitelist.txt --log2stderr