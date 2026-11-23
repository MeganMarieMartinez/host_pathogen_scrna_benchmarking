set -euo pipefail

SHARED=${SHARED:-/home/skyemegan}
CR_WORK=${CR_WORK:-$HOME/cellranger_work}
SEED=${SEED:-100}
SAMPLES=("SRR6825024:${DIR_SRR6825024:-HIV_Infected_1_SRR6825024}"
         "SRR6825026:${DIR_SRR6825026:-HIV_Uninfected_1_SRR6825026}")

read -r -a STEPS <<< "${STEPS_SPEC:-10M:10000000:25M 5M:5000000:10M}"   

log() { echo "[$(date '+%F %T')] $*"; }
die() { log "ERROR: $*"; exit 1; }
SECONDS=0

command -v seqtk >/dev/null || die "seqtk is not on PATH (conda install -c bioconda seqtk)"
log "$(seqtk 2>&1 | grep -i version) | SEED=$SEED | SHARED=$SHARED | CR_WORK=$CR_WORK"

for entry in "${SAMPLES[@]}"; do
  IFS=: read -r srr folder <<< "$entry"
  dir=$SHARED/$folder
  [ -d "$dir" ] || die "folder $dir not found (set DIR_$srr if it has another name)"
  [ -s "$CR_WORK/barcodes_${srr}-25M.tsv" ] || die "$CR_WORK/barcodes_${srr}-25M.tsv not found (run the Cell Ranger script first)"

  for step in "${STEPS[@]}"; do
    IFS=: read -r tag n parent <<< "$step"
    in1=$dir/${srr}_${parent}_1.fastq;  in2=$dir/${srr}_${parent}_2.fastq
    out1=$dir/${srr}_${tag}_1.fastq;    out2=$dir/${srr}_${tag}_2.fastq
    marker=$dir/.${srr}_${tag}.verified

    if [ -e "$marker" ]; then
      log "[$srr $tag] already made and verified, skipping"
    else
      [ -s "$in1" ] && [ -s "$in2" ] || die "[$srr $tag] parent slice missing: $in1 / $in2"
    
      if [ -e "$out1" ] || [ -e "$out2" ]; then
        die "[$srr $tag] $out1 or $out2 already exists and was not made by this script. Not overwriting. Rename the old files first, e.g.: mv $out1 $dir/${srr}_${tag}_old_1.fastq; mv $out2 $dir/${srr}_${tag}_old_2.fastq"
      fi
      log "[$srr $tag] sampling $n pairs from the $parent slice (seed $SEED)"
     
      seqtk sample -s "$SEED" "$in1" "$n" > "$out1.tmp" || die "[$srr $tag] seqtk failed on R1"
      seqtk sample -s "$SEED" "$in2" "$n" > "$out2.tmp" || die "[$srr $tag] seqtk failed on R2"
      for f in "$out1.tmp" "$out2.tmp"; do
        lines=$(wc -l < "$f")
        [ "$lines" -eq $(( 4 * n )) ] || die "[$srr $tag] $f has $lines lines, expected $(( 4 * n ))"
      done
      cmp -s <(awk 'NR%4==1{print $1}' "$out1.tmp") <(awk 'NR%4==1{print $1}' "$out2.tmp") \
        || die "[$srr $tag] R1 and R2 read names differ: pairs are out of sync, do not use these files"
      mv "$out1.tmp" "$out1"; mv "$out2.tmp" "$out2"
      touch "$marker"
      log "[$srr $tag] done and verified: $n pairs, R1/R2 names identical"
    fi

  
    fqdir=$CR_WORK/fastq_${tag}
    mkdir -p "$fqdir"
    for r in 1 2; do
      dst=$fqdir/${srr}_S1_L001_R${r}_001.fastq.gz
      if [ ! -s "$dst" ]; then
        gzip -c "$dir/${srr}_${tag}_${r}.fastq" > "$dst.tmp" && mv "$dst.tmp" "$dst"
        log "[$srr $tag] wrote $dst"
      fi
    done
 
    ln -sfn "barcodes_${srr}-25M.tsv" "$CR_WORK/barcodes_${srr}-${tag}.tsv"
  done
done

log "All done (elapsed ${SECONDS}s)."
log "Check: ls -lh $CR_WORK/fastq_10M $CR_WORK/fastq_5M ; ls -l $CR_WORK/barcodes_*"
log "These slices can be deleted after the results are final; the 25M slices are the source."