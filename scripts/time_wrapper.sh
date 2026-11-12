tool=$1; dataset=$2; shift 2
stage=${STAGE:-all}; rep=${REP:-1}; threads=${THREADS:-NA}
load1=$(cut -d' ' -f1 /proc/loadavg)    # machine load when the call started: a high value means a noisy measurement
start=$(date '+%F %T')
log=$(mktemp)
 
set +e
/usr/bin/time -v -o "$log" "$@"
status=$?
set -eu -o pipefail
 
elapsed=$(grep 'Elapsed (wall clock)' "$log" | awk '{print $NF}' || true)
maxrss=$(grep 'Maximum resident set size' "$log" | awk '{print $NF}' || true)
user=$(grep 'User time (seconds)' "$log" | awk '{print $NF}' || true)
sys=$(grep 'System time (seconds)' "$log" | awk '{print $NF}' || true)
# elapsed is m:ss.xx or h:mm:ss; convert to plain seconds
elapsed_s=$(awk -v t="$elapsed" 'BEGIN{n=split(t,a,":"); s=0; for(i=1;i<=n;i++) s=s*60+a[i]; print s}')
rm -f "$log"
 
RESULTS_DIR=${RESULTS_DIR:-$HOME/host_pathogen_scrna_benchmarking/results}
mkdir -p "$RESULTS_DIR"
f="$RESULTS_DIR/timing.tsv"
header=$'tool\tdataset\telapsed\tmax_rss_kb\tstage\trep\tthreads\telapsed_s\tuser_s\tsys_s\texit_status\tload1_start\tstart_time'
if [ ! -f "$f" ]; then
  printf '%s\n' "$header" > "$f"
elif [ "$(head -n 1 "$f" | awk -F'\t' '{print NF}')" -eq 4 ]; then
  # older 4-column file: extend the header and pad the old rows with NA, keeping their values
  { printf '%s\n' "$header"
    tail -n +2 "$f" | awk -F'\t' -v OFS='\t' '{print $0,"NA","NA","NA","NA","NA","NA","NA","NA","NA"}'
  } > "$f.new" && mv "$f.new" "$f"
fi
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$tool" "$dataset" "${elapsed:-NA}" "${maxrss:-NA}" "$stage" "$rep" "$threads" "${elapsed_s:-NA}" "${user:-NA}" "${sys:-NA}" "$status" "$load1" "$start" >> "$f"
echo "[time_wrapper] $tool $dataset $stage rep$rep: ${elapsed_s}s wall, ${user}s user, ${maxrss} kB peak, exit $status" >&2
exit "$status"