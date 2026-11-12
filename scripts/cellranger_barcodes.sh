set -euo pipefail
 
SHARED=/home/skyemegan
SRC=$SHARED/ref_genomes_w_pathogen/refdata-cellranger-hg38-pEGFP-N1-HIVpNL4-3_HIVonly-4.0.0
EGFP_DIR=$SHARED/egfp_ref # only used if SRC has no EGFP
WORK=${WORK:-$HOME/cellranger_work}
THREADS=${THREADS:-8}
MEM=${MEM:-64}
CHEMISTRY=${CHEMISTRY:-SC3Pv2}
WRAP=${WRAP:-$HOME/host_pathogen_scrna_benchmarking/scripts/time_wrapper.sh} # timing wrapper
TAG=${TAG:-25M}
 
SAMPLES=(SRR6825024 SRR6825026)
declare -A DIRS=(
  [SRR6825024]=${DIR_SRR6825024:-HIV_Infected_1_SRR6825024}
  [SRR6825026]=${DIR_SRR6825026:-HIV_Uninfected_1_SRR6825026}
)
 
log() { echo "[$(date '+%F %T')] $*"; }
die() { log "ERROR: $*"; exit 1; }
SECONDS=0
 
#### Step 0: Preliminary Checks  ####
command -v cellranger >/dev/null || die "cellranger is not on PATH"
log "Using $(cellranger --version 2>&1 | head -n 1)"
log "WORK=$WORK THREADS=$THREADS MEM=${MEM}GB CHEMISTRY=$CHEMISTRY"
[ -s "$SRC/fasta/genome.fa" ]    || die "Missing $SRC/fasta/genome.fa"
[ -s "$SRC/genes/genes.gtf.gz" ] || die "Missing $SRC/genes/genes.gtf.gz"
 
#### Step 1: Truth Reference ####
log "=== Step 1: Truth Reference ==="

# Check if existing reference already annotates EGFP
# Count GTF lines that contain it
N_EGFP=$(zcat "$SRC/genes/genes.gtf.gz" | grep -ci egfp || true)
 
if [ "$N_EGFP" -gt 0 ]; then
  REF=$SRC
  log "EGFP is annotated in $SRC ($N_EGFP GTF lines); using it as-is"

else
  REF=$WORK/ref_build/hg38_HIV_EGFP

  # When mkref finishes, reference.json is created
  if [ -s "$REF/reference.json" ]; then
    log "Combined reference already built: $REF"

  else
    [ -s "$EGFP_DIR/EGFP.fa" ] && [ -s "$EGFP_DIR/EGFP.gtf" ] \
      || die "No EGFP in $SRC and $EGFP_DIR/EGFP.fa / EGFP.gtf not found"

    # the GTF seqname must match a FASTA header or mkref will fail or drop EGFP
    for n in $(grep -v '^#' "$EGFP_DIR/EGFP.gtf" | cut -f1 | sort -u); do
      grep -qE "^>$n([[:space:]]|$)" "$EGFP_DIR/EGFP.fa" \
        || die "GTF seqname '$n' not found in $EGFP_DIR/EGFP.fa headers"
    done

    mkdir -p "$WORK/ref_build"
    cd "$WORK/ref_build"

    # An already existing folder but with no reference.json means there was an incomplete build
    if [ -d "$REF" ]; then
      log "Removing incomplete earlier build at $REF"
      rm -rf "$REF"
    fi

    log "Concatenating human+HIV and EGFP FASTA/GTF"
    # Append EGFP sequence and annotation to human + HIV ones
    sed -e '$a\' "$SRC/fasta/genome.fa" > genome_all.fa     
    sed -e '$a\' "$EGFP_DIR/EGFP.fa" >> genome_all.fa
    zcat "$SRC/genes/genes.gtf.gz" | sed -e '$a\' > genes_all.gtf
    sed -e '$a\' "$EGFP_DIR/EGFP.gtf" >> genes_all.gtf
    # 194 human contigs + 1 HIV + 1 EGFP
    N_CONTIGS=$(grep -c '>' genome_all.fa)
    log "Combined FASTA has $N_CONTIGS contigs (expect 196: 194 human + HIV + EGFP)"
    [ "$N_CONTIGS" -eq 196 ] || log "WARNING: unexpected contig count, check EGFP.fa"

    # mkref builds STAR genome index and packages reference
    log "cellranger mkref starting"
    STAGE=cellranger_mkref THREADS="$THREADS" bash "$WRAP" CellRanger hg38_HIV_EGFP \
      cellranger mkref --genome=hg38_HIV_EGFP --fasta=genome_all.fa --genes=genes_all.gtf \
      --nthreads="$THREADS" --memgb="$MEM"
    log "mkref done (elapsed ${SECONDS}s)"
    log "Last contigs in new reference:" # Should end with HIV-1_pNL4-3 and EGFP
    cut -f1,2 "$REF/fasta/genome.fa.fai" | tail -3
  fi
fi
log "Reference for counting: $REF"
 
#### Step 2: FASTQ Prep ####
#CellRanger expects format: <sample>_S<sample#>_L<lane>_R<read>_001.fastq.gz

log "=== Step 2: FASTQ prep ==="

FQ=$WORK/fastq_${TAG} # one folder per slice size
mkdir -p "$FQ"
GZ=$(command -v pigz || command -v gzip)

for s in "${SAMPLES[@]}"; do
  dir=$SHARED/${DIRS[$s]}

  for r in 1 2; do
    in=$dir/${s}_${TAG}_${r}.fastq
    out=$FQ/${s}_S1_L001_R${r}_001.fastq.gz
    [ -s "$in" ] || die "Missing input $in (run download.sh first)"
    # Skip finished output if its exists already
    if [ -s "$out" ]; then
      log "Already exists: $out"
      continue

    fi
    log "Compressing $in -> $out"
    # Write to temporary file, only rename when gzip finishes successfully
    "$GZ" -c "$in" > "$out.tmp"
    mv "$out.tmp" "$out"
  done
done

log "FASTQ prep done (elapsed ${SECONDS}s)"
 
#### Step 3: CellRanger Count ####
# turns raw reads into a cell x gene count matrix
log "=== Step 3: cellranger count ==="

mkdir -p "$WORK/runs_truth"
cd "$WORK/runs_truth"

for s in "${SAMPLES[@]}"; do
  id=$s-$TAG

  # A completed run creates barcodes.tsv.gz
  if [ -s "$id/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" ]; then
    log "[$id] already finished, skipping"
    continue

  fi
  log "[$id] cellranger count starting"
  # --id          name of the output folder
  # --sample      which FASTQ files to use
  # --transcriptome  the full reference from Step 1
  # --fastqs      folder with the FASTQs for BOTH samples
  # --chemistry   barcode/UMI layout
  # --localcores/--localmem  resource usage
  # --create-bam  also write the alignments
  STAGE=cellranger_count THREADS="$THREADS" bash "$WRAP" CellRanger "$id" \
    cellranger count --id="$id" --sample="$s" \
    --transcriptome="$REF" --fastqs="$FQ" \
    --chemistry="$CHEMISTRY" \
    --localcores="$THREADS" --localmem="$MEM" \
    --create-bam=true

  log "[$id] cellranger count done (elapsed ${SECONDS}s)"
done
 
#### Step 4: Barcode Lists ####
# Used for PathogenTrack
log "=== Step 4: Barcode Lists ==="

for s in "${SAMPLES[@]}"; do
  bc=$WORK/barcodes_${s}-${TAG}.tsv
  zcat "$WORK/runs_truth/$s-$TAG/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" \
    | sed 's/-1$//' > "$bc"
  log "$bc: $(wc -l < "$bc") barcodes" # Number of cells CellRanger called

  # Check reference worked: EGFP should appear as a row in features.tsv.gz 
  n_egfp=$(zcat "$WORK/runs_truth/$s-$TAG/outs/filtered_feature_bc_matrix/features.tsv.gz" | grep -ci egfp || true)
  if [ "$n_egfp" -gt 0 ]; then
    log "[$s] EGFP feature present in the count matrix"

  else
    log "[$s] WARNING: no EGFP feature in the count matrix; check the reference"
  fi
done
 
log "All done (total elapsed: ${SECONDS}s)"
log "Cell Ranger outputs: $WORK/runs_truth/<sample>-$TAG/outs/"
log "PathogenTrack --barcode files: $WORK/barcodes_<sample>-$TAG.tsv"