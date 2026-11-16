#   1. umi_tools extract: reads the 16 bp cell barcode + 10 bp UMI from R1, keeps only
#      barcodes on the supplied list, and moves barcode+UMI into the R2 read name.
#   2. fastp: removes low-quality / low-complexity R2 reads.
#   3. STAR: aligns R2 to the human index. Reads that map to human are discarded.
#   4. Kraken2: classifies the unmapped reads against the HIV-only database.
#   5. Deduplicates by UMI and writes a taxon-by-cell UMI count matrix.

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate pathogentrack
export PATH="$(conda info --base)/envs/umitools/bin:$PATH"

TAG=${TAG:-25M}                
REP=${REP:-1}                   
SMOKE=${SMOKE:-0}
SMOKE_READS=${SMOKE_READS:-1000000}
ONLY=${ONLY:-}                  
WORK=${WORK:-$HOME/pathogentrack_work}
CR_WORK=${CR_WORK:-$HOME/cellranger_work}
FQ_DIR=${FQ_DIR:-$CR_WORK/fastq_${TAG}}      
BC_DIR=${BC_DIR:-$CR_WORK}                   
PREPROC_DIR=${PREPROC_DIR:-$HOME/dvsc_work/preproc}  
STAR_INDEX=${STAR_INDEX:-$WORK/refs/star_hg38}
KRAKEN_DB=${KRAKEN_DB:-$WORK/refs/kraken_hiv}

MIN_READS=${MIN_READS:-10}
CONFIDENCE=${CONFIDENCE:-0.11}
PATTERN=${PATTERN:-CCCCCCCCCCCCCCCCNNNNNNNNNN}   
WRAP=${WRAP:-$HOME/host_pathogen_scrna_benchmarking/scripts/time_wrapper.sh}


SAMPLES=("SRR6825024:HIV_inf1" "SRR6825026:HIV_uninf1")

log() { echo "[$(date '+%F %T')] $*"; }
die() { log "ERROR: $*"; exit 1; }
SECONDS=0

#### Step 0: checks ####
for t in PathogenTrack STAR kraken2 fastp umi_tools; do
  command -v "$t" >/dev/null || die "$t is not on PATH"
done
log "STAR $(STAR --version 2>&1 | head -n 1) | $(kraken2 --version 2>&1 | head -n 1) | umi_tools $(umi_tools --version 2>&1 | head -n 1) | fastp $(fastp --version 2>&1 | head -n 1)"
[ -s "$STAR_INDEX/SAindex" ]    || die "STAR index not found at $STAR_INDEX (run setup_pathogentrack_refs.sh)"
[ -s "$KRAKEN_DB/hash.k2d" ]    || die "Kraken2 database not found at $KRAKEN_DB (run setup_pathogentrack_refs.sh)"
if [ ! -d "$FQ_DIR" ]; then
  die "FASTQ folder $FQ_DIR not found. cellranger_barcodes.sh creates it in its Step 2 (older versions used $CR_WORK/fastq: set FQ_DIR if so)"
fi
if [ "$SMOKE" = 1 ]; then
  RUN_DIR=$WORK/runs_smoke; STAGE_NAME=pathogentrack_smoke; DS_PREFIX=smoke_
else
  RUN_DIR=$WORK/runs_${TAG}; STAGE_NAME=pathogentrack_run; DS_PREFIX=
fi
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
log "TAG=$TAG REP=$REP SMOKE=$SMOKE MIN_READS=$MIN_READS CONFIDENCE=$CONFIDENCE RUN_DIR=$RUN_DIR"

make_subset() {  
  local in=$1 out=$2 n=$3 lines
  if [ -s "$out" ]; then return 0; fi

  ( set +o pipefail; zcat "$in" | head -n $(( 4 * n )) | gzip > "$out.tmp" )
  lines=$(zcat "$out.tmp" | wc -l)
  [ "$lines" -gt 0 ] && [ $(( lines % 4 )) -eq 0 ] || die "subset of $in has $lines lines (not a multiple of 4)"
  mv "$out.tmp" "$out"
}

#### Step 1: run each sample ####
for entry in "${SAMPLES[@]}"; do
  IFS=: read -r srr label <<< "$entry"
  if [ -n "$ONLY" ] && [ "$ONLY" != "$srr" ]; then continue; fi
  dataset=${DS_PREFIX}${label}_${TAG}
  id=${dataset}_rep${REP}                
  marker=$RUN_DIR/${id}.done

  if [ -e "$marker" ]; then
    log "[$id] already finished, skipping"
    continue
  fi

  r1=$FQ_DIR/${srr}_S1_L001_R1_001.fastq.gz
  r2=$FQ_DIR/${srr}_S1_L001_R2_001.fastq.gz
  [ -s "$r1" ] && [ -s "$r2" ] || die "[$id] missing FASTQs: $r1 / $r2"

  # barcode list
  bc=$BC_DIR/barcodes_${srr}-${TAG}.tsv
  if [ ! -s "$bc" ]; then
    if [ "$SMOKE" = 1 ]; then
      wl=$PREPROC_DIR/${label}_${TAG}_whitelist.txt
      [ -s "$wl" ] || die "[$id] no Cell Ranger barcodes ($bc) and no UMI-tools whitelist ($wl) for the smoke test"
      bc=$RUN_DIR/smoke_barcodes_${srr}.tsv
      cut -f1 "$wl" > "$bc"
      log "[$id] smoke test: using UMI-tools whitelist barcodes (not for results)"
    else
      die "[$id] barcode file $bc not found. Wait for cellranger_barcodes.sh to finish (Step 4 writes it)"
    fi
  fi

  if [ "$SMOKE" = 1 ]; then
    make_subset "$r1" "$RUN_DIR/smoke_${srr}_R1.fastq.gz" "$SMOKE_READS"
    make_subset "$r2" "$RUN_DIR/smoke_${srr}_R2.fastq.gz" "$SMOKE_READS"
    r1=$RUN_DIR/smoke_${srr}_R1.fastq.gz
    r2=$RUN_DIR/smoke_${srr}_R2.fastq.gz
    log "[$id] smoke test on the first $SMOKE_READS read pairs"
  fi
  log "[$id] $(wc -l < "$bc") barcodes in $bc"


  if [ -e "$id" ]; then
    aside=${id}.incomplete.$(date +%s)
    log "[$id] moving unfinished earlier output to $aside"
    mv "$id" "$aside"
  fi


  runlog=$RUN_DIR/${id}.log
  log "[$id] PathogenTrack count starting, THREADS=$THREADS (README: 4-6 hours per sample at full size)"
  STAGE=$STAGE_NAME THREADS="$THREADS" REP="$REP" bash "$WRAP" PathogenTrack "$dataset" \
    PathogenTrack count --project_id "$id" --pattern "$PATTERN" \
    --thread "$THREADS" --min_reads "$MIN_READS" --confidence "$CONFIDENCE" \
    --star_index "$STAR_INDEX" --kraken_db "$KRAKEN_DB" \
    --barcode "$bc" --read1 "$r1" --read2 "$r2" 2>&1 | tee "$runlog"

  
  for f in reads_barcoded.fq.gz reads_barcoded_clean.fq.gz unmapped_reads.fq \
           classified_report.txt classified_table.txt microbes.tsv; do
    [ -e "$id/$f" ] || die "[$id] expected output $id/$f was not written; see $runlog"
  done
  [ -s "$id/microbes.tsv" ] || die "[$id] $id/microbes.tsv is empty; see $runlog"
  if grep -qE '^status is [1-9]' "$runlog"; then
    log "[$id] WARNING: these PathogenTrack steps reported a non-zero status:"
    grep -E '^status is [1-9]' "$runlog" | cut -c1-300 | sed 's/^/    /'
  fi
  log "[$id] microbes.tsv: $(( $(head -n 1 "$id/microbes.tsv" | awk -F'\t' '{print NF-1}') )) barcodes, $(( $(wc -l < "$id/microbes.tsv") - 1 )) taxa kept"
  awk -F'\t' 'NR>1 { n=0; u=0; for(i=2;i<=NF;i++) if($i>0){n++; u+=$i} printf "    %s: %d cells with >=1 UMI, %d UMIs\n", $1, n, u }' "$id/microbes.tsv"
  touch "$marker"
  log "[$id] done (elapsed ${SECONDS}s)"
done

log "All done (total elapsed: ${SECONDS}s). Outputs are under $RUN_DIR"