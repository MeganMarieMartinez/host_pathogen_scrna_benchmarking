REPO=$HOME/host_pathogen_scrna_benchmarking

mkdir -p $REPO/envs

conda env export -n pathogentrack --no-builds > $REPO/envs/pathogentrack.yml
conda env export -n dvsc --no-builds > $REPO/envs/dvsc.yml
conda env export -n umitools --no-builds > $REPO/envs/umitools.yml