- Official GitHub repo: [github.com/ncrna/PathogenTrack](https://github.com/ncrna/PathogenTrack)
- Quick start guide: [PathogenTrack/blob/master/doc/QUICK_START.md](https://github.com/ncrna/PathogenTrack/blob/master/doc/QUICK_START.md)
- Also installable via Bioconda: [bioconda.github.io/recipes/pathogentrack](https://bioconda.github.io/recipes/pathogentrack/README.html)
- Zhang W, Xu X, Fu Z, Chen J, Chen S, Tan Y. PathogenTrack and Yeskit: tools for identifying intracellular pathogens from single-cell RNA-sequencing datasets as illustrated by application to COVID-19. Frontiers of Medicine. 2022;16(2):251-262. doi:10.1007/s11684-021-0915-9

# Installation

1. Clone this repo:

```bash
git clone https://github.com/MeganMarieMartinez/host_pathogen_scrna_benchmarking.git
cd host_pathogen_scrna_benchmarking
```
2. Create PathogenTrack and UMI-tools environments (exported from `scripts/export_envs.sh):

```bash
conda env create -f envs/pathogentrack.yml
conda env create -f envs/umitools.yml
```
3. Add `umi_tools` to your PATH whenever you run PathogenTrack:

```bash
export PATH=$HOME/miniconda3/envs/umitools/bin:$PATH
```
# Example Run

1. Human STAR index (Here using the human mitochondrial genome)

```bash
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_012920.1&rettype=fasta&retmode=text" -o chrM.fa
mkdir -p tiny_star_index
STAR --runMode genomeGenerate --genomeDir tiny_star_index/ \ --genomeFastaFiles chrM.fa --genomeSAindexNbases 6
```

2. Kraken2 database using SARS-CoV-2 as pathogen

```bash
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_045512.2&rettype=fasta&retmode=text" -o pathogen_genome.fa
sed -i '1s/.*/>NC_045512.2|kraken:taxid|2697049 SARS-CoV-2/' pathogen_genome.fa

mkdir -p tiny_kraken_db/taxonomy
cd tiny_kraken_db/taxonomy
curl -O https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz tar -xzf taxdump.tar.gz
cd ../..
kraken2-build --add-to-library pathogen_genome.fa --db tiny_kraken_db
kraken2-build --build --db tiny_kraken_db
```

3. Simulate 1,000 read pairs

```bash
dwgsim -N 1000 chrM.fa test_reads
zcat test_reads.bwa.read1.fastq.gz | awk 'NR%4==1{print "@read"++i} NR%4!=1{print}' | gzip > R1_renamed.fastq.gz
zcat test_reads.bwa.read2.fastq.gz | awk 'NR%4==1{print "@read"++i} NR%4!=1{print}' | gzip > R2_renamed.fastq.gz
```
4. Build barcode list from the reads

```bash
zcat R1_renamed.fastq.gz | awk 'NR%4==2 {print substr($0,1,16)}' | sort -u > barcodes.tsv
```

5. Run
```bash
export PATH=$HOME/miniconda3/envs/umitools/bin:$PATH
PathogenTrack count --project_id test_out4 \
  --pattern CCCCCCCCCCCCCCCCNNNNNNNNNN --min_reads 1 --confidence 0.11 \
  --star_index tiny_star_index/ --kraken_db tiny_kraken_db/ \
  --barcode barcodes.tsv \
  --read1 R1_renamed.fastq.gz --read2 R2_renamed.fastq.gz
```
