#Pull 5M slice of infected reads
fastq-dump -X 5000000 --split-files SRR6825024
awk 'FNR==2{print FILENAME, length($0)}' SRR6825024_1.fastq SRR6825024_2.fastq
mv SRR6825026_1.fastq SRR6825024_5M_1.fastq; mv SRR6825024_2.fastq SRR6825024_5M_2.fastq

#Run smoke test on infected slice
cd /home/skyemegan/HIV_Infected_1_SRR6825024
awk 'FNR==2{print FILENAME, length($0)}' SRR6825024_5M_1.fastq SRR6825024_5M_2.fastq   # expect 26 and 98
head -1 SRR6825024_5M_1.fastq SRR6825024_5M_2.fastq 

#Pull 5M slice of uninfected reads
fastq-dump -X 5000000 --split-files SRR6825026
awk 'FNR==2{print FILENAME, length($0)}' SRR6825026_1.fastq SRR6825026_2.fastq
mv SRR6825026_1.fastq SRR6825026_5M_1.fastq; mv SRR6825026_2.fastq SRR6825026_5M_2.fastq

#Run smoke test on uninfected slice
cd /home/skyemegan/HIV_uninfected_1_SRR6825026
awk 'FNR==2{print FILENAME, length($0)}' SRR6825026_5M_1.fastq SRR6825026_5M_2.fastq   # expect 26 and 98
head -1 SRR6825026_5M_1.fastq SRR6825026_5M_2.fastq 