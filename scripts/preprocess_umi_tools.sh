set -eu
 
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate umitools

SHARED=/home/skyemegan
OUT=$HOME/dvsc_work/preproc
PATTERN=CCCCCCCCCCCCCCCCNNNNNNNNNN
mkdir -p "$OUT" && cd "$OUT"
 
run_sample() {   # run_sample <label> <folder in $SHARED> <SRR id>
  label=$1; dir=$SHARED/$2; srr=$3
  R1=$dir/${srr}_5M_1.fastq
  R2=$dir/${srr}_5M_2.fastq
 
  if [ -s "${label}_R2_extracted.fastq.gz" ]; then
    echo "Skipping $label: ${label}_R2_extracted.fastq.gz already exists"
    return
  fi
 
  umi_tools whitelist --stdin "$R1" --bc-pattern=$PATTERN \
    --plot-prefix="$label" --log2stderr > "${label}_whitelist.txt"
 
  umi_tools extract --bc-pattern=$PATTERN \
    --stdin "$R1" --stdout "${label}_R1_extracted.fastq.gz" \
    --read2-in "$R2" --read2-out "${label}_R2_extracted.fastq.gz" \
    --filter-cell-barcode --whitelist="${label}_whitelist.txt" --log2stderr
}
 
run_sample HIV_inf1   HIV_Infected_1_SRR6825024   SRR6825024
run_sample HIV_uninf1 HIV_uninfected_1_SRR6825026 SRR6825026