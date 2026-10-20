set -eu
 
cd ~/DVsc
HG=/home/skyemegan/ref_genomes_w_pathogen/refdata-cellranger-hg38-pEGFP-N1-HIVpNL4-3_HIVonly-4.0.0/fasta/genome.fa
HIV_CONTIG='HIV-1_pNL4-3'
 
cp anno.txt anno_hiv.txt  # re-copies each run, so the HIV row is only added once
grep -n -i -E "immunodeficiency|NC_001802" anno_hiv.txt || echo "No existing HIV rows in anno.txt"
 
# make sure the file ends in a newline so the new row isn't added onto the last one
if [ -n "$(tail -c1 anno_hiv.txt)" ]; then
  echo >> anno_hiv.txt
fi
 
LEN=$(awk -v c="$HIV_CONTIG" '$1==c{print $2}' "$HG.fai")
if [ -z "$LEN" ]; then
  echo "ERROR: $HIV_CONTIG not found in $HG.fai" >&2
  exit 1
fi
echo "length = $LEN"
 
printf '%s\t%s\tHuman immunodeficiency virus 1\tcomplete genome\trefseq\n' "$HIV_CONTIG" "$LEN" >> anno_hiv.txt
tail -2 anno_hiv.txt | cat -A | cut -c1-120 
 