SHARED=/home/skyemegan
REPO=$HOME/host_pathogen_scrna_benchmarking

N=25000000; TAG=25M

#Pull 25M read pairs slice of infected reads (SRR6825024)
mkdir -p $SHARED/HIV_Infected_1_SRR6825024
cd $SHARED/HIV_Infected_1_SRR6825024
fasterq-dump SRR6825024 --split-files
seqtk sample -s100 SRR6825024_1.fastq $N > SRR6825024_${TAG}_1.fastq
seqtk sample -s100 SRR6825024_2.fastq $N > SRR6825024_${TAG}_2.fastq

#Pull 25M slice of uninfected reads (SRR6825026)
mkdir -p $SHARED/HIV_Uninfected_1_SRR6825026
cd $SHARED/HIV_Uninfected_1_SRR6825026
fasterq-dump SRR6825026 --split-files
seqtk sample -s100 SRR6825026_1.fastq $N > SRR6825026_${TAG}_1.fastq
seqtk sample -s100 SRR6825026_2.fastq $N > SRR6825026_${TAG}_2.fastq

#Run smoke tests on uninfected and infected slices
#expect: 100000000 lines
cd $SHARED/HIV_Infected_1_SRR6825024
wc -l SRR6825024_5M_1.fastq SRR6825024_5M_2.fastq
awk 'FNR==2{print FILENAME, length($0)}' SRR6825024_5M_1.fastq SRR6825024_5M_2.fastq
head -1 SRR6825024_5M_1.fastq SRR6825024_5M_2.fastq

cd $SHARED/HIV_Uninfected_1_SRR6825026
wc -l SRR6825026_5M_1.fastq SRR6825026_5M_2.fastq
awk 'FNR==2{print FILENAME, length($0)}' SRR6825026_5M_1.fastq SRR6825026_5M_2.fastq
head -1 SRR6825026_5M_1.fastq SRR6825026_5M_2.fastq