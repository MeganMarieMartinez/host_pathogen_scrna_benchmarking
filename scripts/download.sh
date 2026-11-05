SHARED=/home/skyemegan
REPO=$HOME/host_pathogen_scrna_benchmarking

N=25000000; TAG=25M

log() { echo "[$(date '+%F %T')] $*"; }
SECONDS=0 # Log seconds since starting script

log "Starting 25M slice generation (N=$N, TAG=$TAG)"


#Pull 25M read pairs slice of infected reads (SRR6825024)
log "=== Infected sample (SRR6825024) ==="

mkdir -p $SHARED/HIV_Infected_1_SRR6825024
cd $SHARED/HIV_Infected_1_SRR6825024

log "Downloading + converting full run with fasterq-dump"
fasterq-dump --progress SRR6825024 --split-files
log "fasterq-dump finished (elapsed: ${SECONDS}s)"

log "Sampling R1 with seqtk"
seqtk sample -s100 SRR6825024_1.fastq $N > SRR6825024_${TAG}_1.fastq
log "R1 done (elapsed: ${SECONDS}s)"

log "Sampling R2 with seqtk"
seqtk sample -s100 SRR6825024_2.fastq $N > SRR6825024_${TAG}_2.fastq
log "R2 done (elapsed: ${SECONDS}s)"

ls -lh SRR6825024_${TAG}_1.fastq SRR6825024_${TAG}_2.fastq

#Pull 25M slice of uninfected reads (SRR6825026)
log "=== Uninfected sample (SRR6825026) ==="

mkdir -p $SHARED/HIV_Uninfected_1_SRR6825026
cd $SHARED/HIV_Uninfected_1_SRR6825026

log "Downloading + converting full run with fasterq-dump"
fasterq-dump --progress SRR6825026 --split-files
log "fasterq-dump finished (elapsed: ${SECONDS}s)"

log "Sampling R1 with seqtk"
seqtk sample -s100 SRR6825026_1.fastq $N > SRR6825026_${TAG}_1.fastq
log "R1 done (elapsed: ${SECONDS}s)"

log "Sampling R2 with seqtk"
seqtk sample -s100 SRR6825026_2.fastq $N > SRR6825026_${TAG}_2.fastq
log "R2 done (elapsed: ${SECONDS}s)"

ls -lh SRR6825026_${TAG}_1.fastq SRR6825026_${TAG}_2.fastq

#Run smoke tests on uninfected and infected slices
#expect: 100000000 lines, read lengths 26 (_1) and 98 (_2)
log "=== Smoke tests (expect 100000000 lines per file) ==="

cd $SHARED/HIV_Infected_1_SRR6825024

log "Infected:"
wc -l SRR6825024_25M_1.fastq SRR6825024_25M_2.fastq
awk 'FNR==2{print FILENAME, length($0)}' SRR6825024_25M_1.fastq SRR6825024_25M_2.fastq
head -1 SRR6825024_25M_1.fastq SRR6825024_25M_2.fastq

cd $SHARED/HIV_Uninfected_1_SRR6825026

log "Uninfected:"
wc -l SRR6825026_25M_1.fastq SRR6825026_25M_2.fastq
awk 'FNR==2{print FILENAME, length($0)}' SRR6825026_25M_1.fastq SRR6825026_25M_2.fastq
head -1 SRR6825026_25M_1.fastq SRR6825026_25M_2.fastq

log "All done (total elapsed: ${SECONDS}s)"