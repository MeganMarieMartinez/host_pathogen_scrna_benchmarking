set -eu
 
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate umitools

SHARED=/home/skyemegan
OUT=$HOME/dvsc_work/preproc
PATTERN=CCCCCCCCCCCCCCCCNNNNNNNNNN
mkdir -p "$OUT" && cd "$OUT"

SIZE=25M

log() { echo "[$(date '+%F %T')] $*"; }
SECONDS=0

run_sample() {   # run_sample <label> <folder in $SHARED> <SRR id>
  label=$1; dir=$SHARED/$2; srr=$3
  R1=$dir/${srr}_${SIZE}_1.fastq
  R2=$dir/${srr}_${SIZE}_2.fastq
 
  if [ -s "${label}_R2_extracted.fastq.gz" ]; then
    echo "Skipping $label: ${label}_R2_extracted.fastq.gz already exists"
    return
  fi
  
  log "[$label] umi_tools whitelist starting"
  umi_tools whitelist --stdin "$R1" --bc-pattern=$PATTERN \
    --plot-prefix="$label" --log2stderr > "${label}_whitelist.txt"
  log "[$label] whitelist done: $(wc -l < "${label}_whitelist.txt") barcodes (elapsed ${SECONDS}s)"
 
  log "[$label] umi_tools extract starting"
  umi_tools extract --bc-pattern=$PATTERN \
    --stdin "$R1" --stdout "${label}_R1_extracted.fastq.gz" \
    --read2-in "$R2" --read2-out "${label}_R2_extracted.fastq.gz" \
    --filter-cell-barcode --whitelist="${label}_whitelist.txt" --log2stderr
  log "[$label] extract done (elapsed ${SECONDS}s)"
}
 
run_sample HIV_inf1_${SIZE}   HIV_Infected_1_SRR6825024   SRR6825024
run_sample HIV_uninf1_${SIZE} HIV_Uninfected_1_SRR6825026 SRR6825026