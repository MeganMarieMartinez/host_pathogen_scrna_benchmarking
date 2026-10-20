cd /home/skyemegan/ref_dvsc
HG=/home/skyemegan/fasta/genome.fa

#### Split FASTA file into Human only + HIV ####
grep 'HIV' $HG.fai # shows contig name + length
cut -f1 $HG.fai | grep -v -x 'HIV-1_pNL4-3' > human_contigs.txt
wc -l human_contigs.txt # expect 194

samtools faidx $HG -r human_contigs.txt > human_only.fa
samtools faidx $HG 'HIV-1_pNL4-3' > hiv_pNL43.fa
grep -c '>' human_only.fa hiv_pNL43.fa    # expect 194 and 1

#### Build HISAT2 Indexes ####
conda activate DVsc

#host only
hisat2-build -p 8 human_only.fa host_hisat2_grch38 > host_build.log 2>&1

# virus + host (virus FASTA first)
hisat2-build -p 8 hiv_pNL43.fa,human_only.fa virus_hisat2_grch38 > virus_build.log 2>&1

ls host_hisat2_grch38.*.ht2 | wc -l #expect 8


#### Add HIV Row to Annotation File ####


#### Alignment Test With New ID ####
zcat ~/HIV_inf1_R2_extracted.fastq.gz | head -n 4000000 > /tmp/test1M.fastq
hisat2 --dta -q -x /home/skyemegan/ref_dvsc/virus_hisat2_grch38 -U /tmp/test1M.fastq --sensitive -p 4 \
  | samtools view -F 4 - | cut -f3 | grep -c 'HIV-1_pNL4-3'

#expect a nonzero count 