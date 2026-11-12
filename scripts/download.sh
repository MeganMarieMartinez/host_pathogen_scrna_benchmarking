set -euo pipefail
 
SHARED=${SHARED:-/home/skyemegan}
N=${N:-25000000}           # read pairs wanted per sample
TAG=${TAG:-25M}            # used in output file names
CHUNK=${CHUNK:-20000000}   # spots converted per chunk (about 7 GB of FASTQ)
JOBS=${JOBS:-4}            # chunks processed in parallel
SEED=${SEED:-100}          # chunk i uses seed SEED+i for both mates
 
log() { echo "[$(date '+%F %T')] $*"; }
die() { log "ERROR: $*"; exit 1; }
SECONDS=0
 
# Check the tools BEFORE anything slow, so a missing one is found now, not hours later.
for tool in prefetch srapath sra-stat fastq-dump seqtk; do
  command -v "$tool" >/dev/null || die "$tool not found on PATH (is the right conda env active?)"
done
log "Starting $TAG slice generation (N=$N, CHUNK=$CHUNK, JOBS=$JOBS)"
 
# Total number of spots (read pairs) in the run. Taken from SPOTS_<srr> if set,
# otherwise read from sra-stat.
total_spots() {   # total_spots <srr> <local .sra path>; prints the count, or returns 1
  local var="SPOTS_$1" n
  if [ -n "${!var:-}" ]; then echo "${!var}"; return 0; fi
  n=$(sra-stat --xml --quick "$2" 2>/dev/null | grep -oE 'spot_count="[0-9]+"' | head -n 1 | grep -oE '[0-9]+' || true)
  if [ -z "$n" ]; then
    # message goes to stderr: stdout of this function is captured by $(...) in the caller
    log "could not read the spot count from sra-stat; set $var=<spots> and re-run" >&2
    return 1
  fi
  echo "$n"
}
 
# Convert one chunk of spots and sample from it. Safe to re-run: skipped if finished.
# Uses (exported): SRA WORK CHUNK TOTAL N SEED
do_chunk() {   # do_chunk <chunk index>
  local i=$1 start end prev_end n cdir r f lines
  start=$(( i * CHUNK + 1 ))
  end=$(( start + CHUNK - 1 ))
  if [ "$end" -gt "$TOTAL" ]; then end=$TOTAL; fi
  prev_end=$(( i * CHUNK ))
  # Pairs to keep from this chunk. Using floor(N*end/TOTAL) - floor(N*prev_end/TOTAL)
  # makes the per-chunk counts add up to exactly N over all chunks.
  n=$(( N * end / TOTAL - N * prev_end / TOTAL ))
  if [ -e "$WORK/chunk_$i.done" ]; then return 0; fi
 
  cdir=$WORK/tmp_$i
  rm -rf "$cdir"; mkdir -p "$cdir"
  log "chunk $i: spots $start-$end, keeping $n pairs"
  fastq-dump -N "$start" -X "$end" --split-files -O "$cdir" "$SRA" >/dev/null \
    || { log "chunk $i: fastq-dump failed"; return 1; }
 
  for r in 1 2; do
    local files=("$cdir"/*_"$r".fastq)       # fastq-dump names these <accession>_<read>.fastq
    f=${files[0]}
    [ -s "$f" ] || { log "chunk $i: no R$r FASTQ produced"; return 1; }
    seqtk sample -s "$(( SEED + i ))" "$f" "$n" > "$WORK/chunk_${i}_R$r.fastq.tmp" \
      || { log "chunk $i: seqtk failed on R$r"; return 1; }
    lines=$(wc -l < "$WORK/chunk_${i}_R$r.fastq.tmp")
    [ "$lines" -eq $(( 4 * n )) ] \
      || { log "chunk $i R$r: $lines lines, expected $(( 4 * n ))"; return 1; }
    mv "$WORK/chunk_${i}_R$r.fastq.tmp" "$WORK/chunk_${i}_R$r.fastq"
  done
 
  # R1 and R2 must contain the same reads in the same order
  cmp -s <(awk 'NR%4==1{print $1}' "$WORK/chunk_${i}_R1.fastq") \
         <(awk 'NR%4==1{print $1}' "$WORK/chunk_${i}_R2.fastq") \
    || { log "chunk $i: R1/R2 read names differ (pairs out of sync)"; return 1; }
 
  rm -rf "$cdir"                  # free the big temporary FASTQs
  touch "$WORK/chunk_$i.done"     # marker: this chunk never needs redoing
  log "chunk $i done"
}
export -f do_chunk log die
 
make_slice() {   # make_slice <SRR id> <folder name under $SHARED>
  local srr=$1 dir=$SHARED/$2 nchunks r i out lines
  mkdir -p "$dir"
  cd "$dir"
  log "=== $srr ($dir) ==="
  log "Free space here: $(df -h . | tail -n 1 | awk '{print $4}')"
 
  # ---- Step 1: get the .sra file (resumes / skips if already complete) ----
  log "prefetch $srr into $dir (re-run this script to resume if it is interrupted)"
  prefetch --max-size 500G -O "$dir" "$srr" || die "prefetch $srr failed"
  # prefetch -O <dir> normally writes <dir>/<srr>/<srr>.sra. Look there first, then
  # anywhere under <dir>, then ask srapath, so a different layout is found right away
  # instead of after hours of work.
  SRA=$dir/$srr/$srr.sra
  if [ ! -s "$SRA" ]; then
    SRA=$(find "$dir" -maxdepth 2 -name "$srr.sra" -size +0 -print -quit)
    [ -n "$SRA" ] || SRA=$(srapath "$srr")
  fi
  case $SRA in
    /*) ;;
    *) die "could not find the downloaded $srr.sra (srapath returned '$SRA'); prefetch may not have completed" ;;
  esac
  [ -s "$SRA" ] || die "local file $SRA is missing or empty"
  log "Local SRA file: $SRA ($(du -h "$SRA" | cut -f1))"
 
  TOTAL=$(total_spots "$srr" "$SRA") || die "cannot determine the number of spots for $srr"
  [ "$N" -lt "$TOTAL" ] || die "$srr has only $TOTAL spots, which is fewer than the N=$N requested"
  nchunks=$(( (TOTAL + CHUNK - 1) / CHUNK ))
  WORK=$dir/work_${TAG}
  mkdir -p "$WORK"
  export SRA TOTAL WORK CHUNK N SEED
  log "$srr: $TOTAL spots -> $nchunks chunks of up to $CHUNK, $JOBS at a time"
 
  # ---- Step 2: convert + sample each chunk (finished chunks are skipped) ----
  seq 0 $(( nchunks - 1 )) \
    | xargs -P "$JOBS" -I{} bash -c 'set -euo pipefail; do_chunk "$1"' _ {} \
    || die "some chunks failed (see messages above); fix the cause and re-run: finished chunks are skipped"
  log "all $nchunks chunks done (elapsed: ${SECONDS}s)"
 
  # ---- Step 3: combine, verify, publish ----
  for r in 1 2; do
    out=${srr}_${TAG}_${r}.fastq
    if [ -s "$out" ]; then
      log "$out already exists; skipping"
      continue
    fi
    : > "$out.tmp"
    for i in $(seq 0 $(( nchunks - 1 ))); do
      cat "$WORK/chunk_${i}_R$r.fastq" >> "$out.tmp"
    done
    lines=$(wc -l < "$out.tmp")
    [ "$lines" -eq $(( 4 * N )) ] || die "$out.tmp has $lines lines, expected $(( 4 * N ))"
    mv "$out.tmp" "$out"       # only a verified slice gets the final name
    log "$out done (elapsed: ${SECONDS}s)"
  done
 
  # ---- Smoke tests ----
  cmp -s <(awk 'NR%4==1{print $1} NR>=400000{exit}' "${srr}_${TAG}_1.fastq") \
         <(awk 'NR%4==1{print $1} NR>=400000{exit}' "${srr}_${TAG}_2.fastq") \
    || die "read names in ${srr}_${TAG}_1/2.fastq differ: pairs out of sync"
  log "$srr: line counts (expect $(( 4 * N )) each), read lengths (expect 26 and 98), first headers"
  wc -l "${srr}_${TAG}_1.fastq" "${srr}_${TAG}_2.fastq"
  awk 'FNR==2{print FILENAME, length($0)}' "${srr}_${TAG}_1.fastq" "${srr}_${TAG}_2.fastq"
  head -n 1 "${srr}_${TAG}_1.fastq" "${srr}_${TAG}_2.fastq"
  log "$srr: slices verified. Chunk files in $WORK and the .sra can be deleted to free space."
}
 
make_slice SRR6825024 HIV_Infected_1_SRR6825024
make_slice SRR6825026 HIV_Uninfected_1_SRR6825026
 
log "All done (total elapsed: ${SECONDS}s)"