set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate pathogentrack

SHARED=/home/skyemegan
HOST_REF=${HOST_REF:-$SHARED/refdata-gex-GRCh38-2024-A}   
HIV_SRC=${HIV_SRC:-$SHARED/ref_genomes_w_pathogen/refdata-cellranger-hg38-pEGFP-N1-HIVpNL4-3_HIVonly-4.0.0}
HIV_CONTIG=${HIV_CONTIG:-HIV-1_pNL4-3} 
HIV_TAXID=${HIV_TAXID:-11676} 
WORK=${WORK:-$HOME/pathogentrack_work}
STAR_INDEX=${STAR_INDEX:-$WORK/refs/star_hg38}
KRAKEN_DB=${KRAKEN_DB:-$WORK/refs/kraken_hiv}
THREADS=${THREADS:-8}
TAXDUMP=${TAXDUMP:-}
WRAP=${WRAP:-$HOME/host_pathogen_scrna_benchmarking/scripts/time_wrapper.sh}

log() { echo "[$(date '+%F %T')] $*"; }
die() { log "ERROR: $*"; exit 1; }
SECONDS=0
trap 'rm -f "$WORK/refs/host_genes.gtf.tmp"' EXIT

#### Step 0: checks ####
command -v STAR >/dev/null          || die "STAR is not on PATH (is the pathogentrack env active?)"
command -v kraken2-build >/dev/null || die "kraken2-build is not on PATH (is the pathogentrack env active?)"
[ -s "$HOST_REF/fasta/genome.fa" ] || die "Missing $HOST_REF/fasta/genome.fa"
HOST_GTF=""
for g in "$HOST_REF/genes/genes.gtf" "$HOST_REF/genes/genes.gtf.gz"; do
  if [ -s "$g" ]; then HOST_GTF=$g; break; fi
done
[ -n "$HOST_GTF" ] || die "No genes.gtf or genes.gtf.gz under $HOST_REF/genes"
[ -s "$HIV_SRC/fasta/genome.fa" ] || die "Missing $HIV_SRC/fasta/genome.fa"
log "STAR $(STAR --version) | $(kraken2-build --version 2>&1 | head -n 1)"
log "WORK=$WORK THREADS=$THREADS"
mkdir -p "$WORK/refs"

#### Step 1: human STAR index ####

log "=== Step 1: STAR index ==="
if [ -s "$STAR_INDEX/SAindex" ] && [ -s "$STAR_INDEX/genomeParameters.txt" ]; then
  log "STAR index already built: $STAR_INDEX"
else
  rm -rf "$STAR_INDEX" # a folder without SAindex is an interrupted build
  mkdir -p "$STAR_INDEX"
  GTF=$HOST_GTF
  case $HOST_GTF in
    *.gz) GTF=$WORK/refs/host_genes.gtf.tmp; zcat "$HOST_GTF" > "$GTF" ;;
  esac
  log "STAR genomeGenerate starting (slow; needs about 30 GB RAM)"
  STAGE=star_index_build THREADS="$THREADS" bash "$WRAP" PathogenTrack star_hg38 \
    STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$STAR_INDEX" \
    --genomeFastaFiles "$HOST_REF/fasta/genome.fa" --sjdbGTFfile "$GTF" --sjdbOverhang 100 \
    --limitGenomeGenerateRAM 64000000000 --outFileNamePrefix "$STAR_INDEX/"
  [ -s "$STAR_INDEX/SAindex" ] || die "STAR finished but $STAR_INDEX/SAindex is missing"
  log "STAR index done (elapsed ${SECONDS}s): $(du -sh "$STAR_INDEX" | cut -f1)"
fi

#### Step 2: HIV-only Kraken2 database ####

log "=== Step 2: Kraken2 HIV database ==="
if [ -s "$KRAKEN_DB/hash.k2d" ] && [ -s "$KRAKEN_DB/opts.k2d" ] && [ -s "$KRAKEN_DB/taxo.k2d" ]; then
  log "Kraken2 database already built: $KRAKEN_DB"
else
  rm -rf "$KRAKEN_DB"
  mkdir -p "$KRAKEN_DB/taxonomy"

  
  awk -v c=">$HIV_CONTIG" '/^>/{p=($1==c)} p' "$HIV_SRC/fasta/genome.fa" > "$WORK/refs/hiv_raw.fa"
  [ "$(grep -c '>' "$WORK/refs/hiv_raw.fa")" -eq 1 ] \
    || die "expected exactly one contig named $HIV_CONTIG in $HIV_SRC/fasta/genome.fa"
  { echo ">${HIV_CONTIG}|kraken:taxid|${HIV_TAXID}"; grep -v '>' "$WORK/refs/hiv_raw.fa"; } > "$WORK/refs/hiv_kraken.fa"
  log "HIV sequence: $(grep -v '>' "$WORK/refs/hiv_kraken.fa" | tr -d '\n' | wc -c) bp, taxid $HIV_TAXID"

  # taxonomy files
  if [ -n "$TAXDUMP" ]; then
    [ -s "$TAXDUMP" ] || die "TAXDUMP=$TAXDUMP does not exist"
    cp "$TAXDUMP" "$WORK/refs/taxdump.tar.gz"
  else
    log "Downloading NCBI taxonomy"
    curl -fsSL --retry 3 -o "$WORK/refs/taxdump.tar.gz" https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz \
      || die "could not download taxdump.tar.gz (no route to NCBI?). Download it elsewhere, copy it to the server, and rerun with TAXDUMP=/path/to/taxdump.tar.gz"
  fi
  tar -xzf "$WORK/refs/taxdump.tar.gz" -C "$KRAKEN_DB/taxonomy" names.dmp nodes.dmp \
    || die "could not extract names.dmp and nodes.dmp from the taxdump archive"
  grep -qE "^${HIV_TAXID}[[:space:]]" "$KRAKEN_DB/taxonomy/names.dmp" \
    || die "taxid $HIV_TAXID not found in names.dmp"

  kraken2-build --add-to-library "$WORK/refs/hiv_kraken.fa" --db "$KRAKEN_DB"
  log "kraken2-build --build starting"
  STAGE=kraken2_build THREADS="$THREADS" bash "$WRAP" PathogenTrack kraken_hiv \
    kraken2-build --build --db "$KRAKEN_DB" --threads "$THREADS"
  for f in hash.k2d opts.k2d taxo.k2d; do
    [ -s "$KRAKEN_DB/$f" ] || die "kraken2-build finished but $KRAKEN_DB/$f is missing"
  done
  log "Kraken2 database done (elapsed ${SECONDS}s): $(du -sh "$KRAKEN_DB" | cut -f1)"
fi

log "All done (total elapsed: ${SECONDS}s)"
log "STAR index:  $STAR_INDEX"
log "Kraken2 DB:  $KRAKEN_DB"