cd ~/DVsc
ls -lh $HOME/dvsc_work/preproc/HIV_inf1_R2_extracted.fastq.gz

cat > hiv_inf1.config << 'EOF'
{
  "resources":{
    "humdb":"/home/skyemegan/ref_dvsc/host_hisat2_grch38",
    "virusdb":"/home/skyemegan/ref_dvsc/virus_hisat2_grch38",
    "virusdb_star":"/home/slarosa/tiny_star_index",
    "annotatedfile":"/home/slarosa/DVsc/anno_hiv.txt"
  },
  "tools":{
    "fastp":{"call":"fastp","threads":4},
    "hisat2":{"call":"hisat2","threads":4},
    "star":{"call":"STAR","threads":4},
    "samtools":{"call":"samtools","threads":4}
  },
  "fastq":"/home/slarosa/HIV_inf1_R2_extracted.fastq.gz",
  "sample":"HIV_inf1_5M",
  "method":"scRNAseq",
  "readlength":98,
  "outpath":"/home/slarosa/DVsc/HIV_inf1_5M_output",
  "alignment":"hisat2"
}
EOF