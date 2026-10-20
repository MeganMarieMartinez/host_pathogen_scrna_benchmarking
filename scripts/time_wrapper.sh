#!/usr/bin/env bash
# usage: time_wrapper.sh <tool> <dataset> <command...>
#Logs elapsed time and peak memory to a TSV

set -euo pipefail
tool=$1; dataset=$2; shift 2
log=$(mktemp)
/usr/bin/time -v -o "$log" "$@"
elapsed=$(grep 'Elapsed (wall clock)' "$log" | awk '{print $NF}')
maxrss=$(grep 'Maximum resident set size' "$log" | awk '{print $NF}')

RESULTS_DIR=${RESULTS_DIR:-$HOME/host_pathogen_scrna_benchmarking/results}
mkdir -p "$RESULTS_DIR"
[ -f "$RESULTS_DIR/timing.tsv" ] || printf 'tool\tdataset\telapsed\tmax_rss_kb\n' > "$RESULTS_DIR/timing.tsv"
printf '%s\t%s\t%s\t%s\n' "$tool" "$dataset" "$elapsed" "$maxrss" >> "$RESULTS_DIR/timing.tsv"