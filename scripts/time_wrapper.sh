#!/usr/bin/env bash
# usage: time_wrapper.sh <tool> <dataset> <command...>
set -euo pipefail
tool=$1; dataset=$2; shift 2
log=$(mktemp)
/usr/bin/time -v -o "$log" "$@"
elapsed=$(grep 'Elapsed (wall clock)' "$log" | awk '{print $NF}')
maxrss=$(grep 'Maximum resident set size' "$log" | awk '{print $NF}')
mkdir -p results
[ -f results/timing.tsv ] || printf 'tool\tdataset\telapsed\tmax_rss_kb\n' > results/timing.tsv
printf '%s\t%s\t%s\t%s\n' "$tool" "$dataset" "$elapsed" "$maxrss" >> results/timing.tsv