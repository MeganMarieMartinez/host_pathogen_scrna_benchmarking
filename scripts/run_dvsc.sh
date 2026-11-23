set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dvsc

TAG=${TAG:-25M}
REP=${REP:-1}
ONLY=${ONLY:-}
THREADS=${THREADS:-8}
WORKROOT=${WORKROOT:-$HOME/dvsc_work}
PRE=${PRE:-$WORKROOT/preproc_${TAG}_rep${REP}}
RUN_DIR=${RUN_DIR:-$WORKROOT/runs_${TAG}_rep${REP}}
SCRIPTS=${SCRIPTS:-$HOME/host_pathogen_scrna_benchmarking/scripts}
WRAP=${WRAP:-$SCRIPTS/time_wrapper.sh}
DVSC_DIR=${DVSC_DIR:-$HOME/DVsc}                      
REF=${REF:-/home/skyemegan/ref_dvsc}                  
ANNO=${ANNO:-$WORKROOT/anno_hiv.txt}                  
STAR_TINY=${STAR_TINY:-$HOME/tiny_star_index}          
SRC_LINK=${SRC_LINK:-$WORKROOT/src}                    

SAMPLES=("SRR6825024:HIV_inf1" "SRR6825026:HIV_uninf1")

log() { echo "[$(date '+%F %T')] $*"; }
die() { log "ERROR: $*"; exit 1; }
SECONDS=0

#### Step 0: checks ####
for t in snakemake hisat2 fastp samtools; do
  command -v "$t" >/dev/null || die "$t is not on PATH (is the dvsc env active?)"
done
[ "$(ls "$REF"/host_hisat2_grch38.*.ht2 2>/dev/null | wc -l)" -eq 8 ]  || die "host HISAT2 index incomplete in $REF"
[ "$(ls "$REF"/virus_hisat2_grch38.*.ht2 2>/dev/null | wc -l)" -eq 8 ] || die "virus HISAT2 index incomplete in $REF"
[ -s "$ANNO" ] || die "annotation file $ANNO not found (update_dvsc_anno.sh makes it)"
[ -s "$DVSC_DIR/analysis.snake" ] || die "$DVSC_DIR/analysis.snake not found"
[ -d "$SRC_LINK" ] || die "$SRC_LINK (DVsc's src/ folder) not found (setup_dvsc.sh makes it)"
[ -s "$SCRIPTS/preprocess_umi_tools.sh" ] || die "$SCRIPTS/preprocess_umi_tools.sh not found"
log "TAG=$TAG REP=$REP THREADS=$THREADS ONLY=${ONLY:-both} PRE=$PRE RUN_DIR=$RUN_DIR"
log "snakemake $(snakemake --version) | hisat2 $(hisat2 --version | head -n 1 | awk '{print $NF}') | fastp $(fastp --version 2>&1)"

#### Step 1: UMI-tools whitelist + extract (both timed inside the script) ####
log "=== Step 1: UMI-tools preprocessing ==="
OUT="$PRE" SIZE="$TAG" REP="$REP" ONLY="$ONLY" WRAP="$WRAP" bash "$SCRIPTS/preprocess_umi_tools.sh"

#### Step 2 and 3: config + snakemake for each sample ####
mkdir -p "$RUN_DIR"
ln -sfn "$(readlink -f "$SRC_LINK")" "$RUN_DIR/src"
mkdir -p "$RUN_DIR/$THREADS"

for entry in "${SAMPLES[@]}"; do
  IFS=: read -r srr label <<< "$entry"
  if [ -n "$ONLY" ] && [ "$ONLY" != "$srr" ]; then continue; fi
  dataset=${label}_${TAG}
  marker=$RUN_DIR/${dataset}.done
  if [ -e "$marker" ]; then log "[$dataset] already finished, skipping"; continue; fi

  fq=$PRE/${dataset}_R2_extracted.fastq.gz
  wl=$PRE/${dataset}_whitelist.txt
  [ -s "$fq" ] && [ -s "$wl" ] || die "[$dataset] preprocessing output missing: $fq / $wl"

  cfg=$RUN_DIR/${dataset}.config
  cat > "$cfg" <<EOF
{
  "resources":{
    "humdb":"$REF/host_hisat2_grch38",
    "virusdb":"$REF/virus_hisat2_grch38",
    "virusdb_star":"$STAR_TINY",
    "annotatedfile":"$ANNO"
  },
  "tools":{
    "fastp":{"call":"fastp","threads":$THREADS},
    "hisat2":{"call":"hisat2","threads":$THREADS},
    "star":{"call":"STAR","threads":$THREADS},
    "samtools":{"call":"samtools","threads":$THREADS}
  },
  "fastq":"$fq",
  "whitelist":"$wl",
  "sample":"$dataset",
  "method":"scRNAseq",
  "readlength":98,
  "outpath":"$RUN_DIR/${dataset}_output",
  "alignment":"hisat2"
}
EOF
  log "[$dataset] config written: $cfg"

  cd "$RUN_DIR"
  if [ -e "$dataset" ]; then
    aside=${dataset}.incomplete.$(date +%s)
    log "[$dataset] moving unfinished earlier output to $aside"
    mv "$dataset" "$aside"
  fi

  log "[$dataset] snakemake starting"
  STAGE=dvsc_snakemake THREADS="$THREADS" REP="$REP" bash "$WRAP" DVsc "$dataset" \
    snakemake -s "$DVSC_DIR/analysis.snake" --configfile "$cfg" --cores "$THREADS"

  final=$RUN_DIR/$dataset/out/${dataset}.countbyumi.filter.txt
  [ -e "$final" ] || die "[$dataset] snakemake finished but $final was not written"
  log "[$dataset] result: $(wc -l < "$final") cells with >=1 HIV UMI, $(awk '{s+=$3} END{print s+0}' "$final") UMIs (column 3), file: $final"
  touch "$marker"
  log "[$dataset] done (elapsed ${SECONDS}s)"
done

log "All done (total elapsed: ${SECONDS}s)."