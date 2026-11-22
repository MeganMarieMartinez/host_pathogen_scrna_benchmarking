#!/usr/bin/env bash
conda activate viral-track

### analysis on infected sample ###
echo 'Running Viral Track for Infected Sample (SRR6825024-25M)'

echo 'Viral Track Scanning'

Rscript ./edited_Viral_Track_scanning.R ../data/viraltrack_files/SRR6825024-25M_parameters.txt ../data/viraltrack_files/SRR6825024-25M_targetfile.txt 

wait

echo 'Viral Track Demultiplexing'

Rscript ./edited_VT_demux.R ../data/viraltrack_files/SRR6825024-25M_parameters.txt ../data/viraltrack_files/SRR6825024-25M_targetfile.txt 

wait

echo 'Viral Track Analysis Done for (SRR6825024-25M) '

### analysis on uninfected sample ###
echo 'Running Viral Track for Infected Sample (SRR6825026-25M)'

echo 'Viral Track Scanning'

Rscript ./edited_Viral_Track_scanning.R ../data/viraltrack_files/SRR6825026-25M_parameters.txt ../data/viraltrack_files/SRR6825026-25M_targetfile.txt 

wait

echo 'Viral Track Demultiplexing'

Rscript ./edited_VT_demux.R ../data/viraltrack_files/SRR6825026-25M_parameters.txt ../data/viraltrack_files/SRR6825026-25M_targetfile.txt 

wait