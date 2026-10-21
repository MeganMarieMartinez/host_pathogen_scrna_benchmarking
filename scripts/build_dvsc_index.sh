set -eu
 
REFDIR=/home/skyemegan/ref_dvsc
HG=/home/skyemegan/ref_genomes_w_pathogen/refdata-cellranger-hg38-pEGFP-N1-HIVpNL4-3_HIVonly-4.0.0/fasta/genome.fa
HIV_CONTIG='HIV-1_pNL4-3'
DVSC_ID='NC_pNL43'
THREADS=8
FQ="$HOME/dvsc_work/preproc/HIV_inf1_R2_extracted.fastq.gz"   # output from preprocess_umi_tools.sh
 
mkdir -p "$REFDIR"
cd "$REFDIR"
 
#load conda shell hook
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dvsc

#### Split FASTA file into Human only + HIV ####
grep 'HIV' $HG.fai # shows contig name + length
cut -f1 $HG.fai | grep -v -x 'HIV-1_pNL4-3' > human_contigs.txt
wc -l human_contigs.txt # expect 194

samtools faidx $HG -r human_contigs.txt > human_only.fa
samtools faidx $HG 'HIV-1_pNL4-3' > hiv_pNL43.fa
sed -i "1s/.*/>${DVSC_ID}/" hiv_pNL43.fa  
grep '>' hiv_pNL43.fa # expect >NC_pNL43
grep -c '>' human_only.fa hiv_pNL43.fa # expect 194 and 1

#### Build HISAT2 Indexes ####

# rebuild the virus index 
rm -f virus_hisat2_grch38.*.ht2

# host only
pid_host=""
if [ "$(ls host_hisat2_grch38.*.ht2 2>/dev/null | wc -l)" -eq 8 ]; then
  echo "Host index already complete; skipping host build"
else
  hisat2-build -p "$THREADS" human_only.fa host_hisat2_grch38 > host_build.log 2>&1 &
  pid_host=$!
fi

# virus + host (virus FASTA first according to DVsc README)
hisat2-build -p "$THREADS" hiv_pNL43.fa,human_only.fa virus_hisat2_grch38 > virus_build.log 2>&1 &
pid_virus=$!
 
[ -z "$pid_host" ] || wait "$pid_host"
wait "$pid_virus"
 
for idx in host_hisat2_grch38 virus_hisat2_grch38; do
  echo "$idx: $(ls ${idx}.*.ht2 | wc -l) index files (expect 8)"
  du -ch ${idx}.*.ht2 | tail -1
done
 

#### Alignment Test With New HIV ID ####

TMP=$(mktemp -p "$HOME" test1M.XXXXXX.fastq)
zcat "$FQ" | head -n 4000000 > "$TMP"
N_HIV=$(hisat2 --dta -q -x "$REFDIR/virus_hisat2_grch38" -U "$TMP" --sensitive -p 4 \
  | samtools view -F 4 - | cut -f3 | grep -c -x "$DVSC_ID" || true)
rm -f "$TMP"
echo "Reads aligned to $DVSC_ID in first 1M reads: $N_HIV (expect > 0)"