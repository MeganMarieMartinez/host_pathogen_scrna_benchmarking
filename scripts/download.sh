SHARED=/home/skyemegan
REPO=$HOME/host_pathogen_scrna_benchmarking

#Pull 5M read pairs slice of infected reads (SRR6825024)
#NOTE: This gets RANDOM 5M reads from entire run
#For the uninfected slice, it gets the FIRST 5M reads
#WILL NEED RESLICING LATER
mkdir -p $SHARED/HIV_Infected_1_SRR6825024
cd $SHARED/HIV_Infected_1_SRR6825024
fasterq-dump SRR6825024 --split-files
seqtk sample -s100 SRR6825024_1.fastq 5000000 > SRR6825024_5M_1.fastq
seqtk sample -s100 SRR6825024_2.fastq 5000000 > SRR6825024_5M_2.fastq

#Pull 5M slice of uninfected reads (SRR6825026)
mkdir -p $SHARED/HIV_uninfected_1_SRR6825026
cd $SHARED/HIV_uninfected_1_SRR6825026
fastq-dump -X 5000000 --split-files SRR6825026
mv SRR6825026_1.fastq SRR6825026_5M_1.fastq
mv SRR6825026_2.fastq SRR6825026_5M_2.fastq

#Run smoke tests on uninfected and infected slices
#expect: 20000000 lines and
# read lengths 26 (_1) and 98 (_2)
cd $SHARED/HIV_Infected_1_SRR6825024
wc -l SRR6825024_5M_1.fastq SRR6825024_5M_2.fastq
awk 'FNR==2{print FILENAME, length($0)}' SRR6825024_5M_1.fastq SRR6825024_5M_2.fastq
head -1 SRR6825024_5M_1.fastq SRR6825024_5M_2.fastq

cd $SHARED/HIV_uninfected_1_SRR6825026
wc -l SRR6825026_5M_1.fastq SRR6825026_5M_2.fastq
awk 'FNR==2{print FILENAME, length($0)}' SRR6825026_5M_1.fastq SRR6825026_5M_2.fastq
head -1 SRR6825026_5M_1.fastq SRR6825026_5M_2.fastq