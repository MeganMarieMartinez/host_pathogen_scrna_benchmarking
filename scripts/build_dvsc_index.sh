set -eu
 
REFDIR=/home/skyemegan/ref_dvsc
HG=/home/skyemegan/fasta/genome.fa
HIV_CONTIG='HIV-1_pNL4-3'
THREADS=8
FQ="$HOME/HIV_inf1_R2_extracted.fastq.gz"   # output from preprocess_umi_tools.sh
 
mkdir -p "$REFDIR"
cd "$REFDIR"
 
#load conda shell hook
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate DVsc

#### Split FASTA file into Human only + HIV ####
grep 'HIV' $HG.fai # shows contig name + length
cut -f1 $HG.fai | grep -v -x 'HIV-1_pNL4-3' > human_contigs.txt
wc -l human_contigs.txt # expect 194

samtools faidx $HG -r human_contigs.txt > human_only.fa
samtools faidx $HG 'HIV-1_pNL4-3' > hiv_pNL43.fa
grep -c '>' human_only.fa hiv_pNL43.fa    # expect 194 and 1

#### Build HISAT2 Indexes ####

# host only
hisat2-build -p "$THREADS" human_only.fa host_hisat2_grch38 > host_build.log 2>&1 &
pid_host=$!

# virus + host (virus FASTA first according to DVsc README)
hisat2-build -p "$THREADS" hiv_pNL43.fa,human_only.fa virus_hisat2_grch38 > virus_build.log 2>&1 &
pid_virus=$!

wait "$pid_host"
wait "$pid_virus"
 
for idx in host_hisat2_grch38 virus_hisat2_grch38; do
  echo "$idx: $(ls ${idx}.*.ht2 | wc -l) index files (expect 8)"
  du -ch ${idx}.*.ht2 | tail -1       
done

#### Alignment Test With New HIV ID ####

TMP=$(mktemp -p "$HOME" test1M.XXXXXX.fastq)
zcat "$FQ" | head -n 4000000 > "$TMP"
N_HIV=$(hisat2 --dta -q -x "$REFDIR/virus_hisat2_grch38" -U "$TMP" --sensitive -p 4 \
  | samtools view -F 4 - | cut -f3 | grep -c -x "$HIV_CONTIG" || true)
rm -f "$TMP"
echo "Reads aligned to $HIV_CONTIG in first 1M reads: $N_HIV (expect > 0)"