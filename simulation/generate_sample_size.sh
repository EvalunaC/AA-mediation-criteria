#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")"
RSCRIPT=${RSCRIPT:-Rscript}
N_CORES=${1:-16}; CHUNK=${2:-5}
N_REPS=500

mkdir -p logs ../data/sample_size
STAMP=$(date +%Y%m%d_%H%M%S)
TASKS=logs/tasks_ss_${STAMP}.txt; : > "$TASKS"
for n in 100 200 300; do
  for sc in sc_null sc_snull sc_modl sc_med sc_direct sc_mixed; do
    for ((from=1; from<=N_REPS; from+=CHUNK)); do
      to=$((from+CHUNK-1)); (( to > N_REPS )) && to=$N_REPS
      echo "$n $sc $from $to" >> "$TASKS"
    done
  done
done
echo "queued $(wc -l < "$TASKS") tasks on $N_CORES workers"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
xargs -a "$TASKS" -P "$N_CORES" -L 1 \
  sh -c "$RSCRIPT worker_sample_size.R \$0 \$1 \$2 \$3 >> logs/${STAMP}_n\$0_\$1.log 2>&1 || echo \"FAILED: \$0 \$1 \$2 \$3\" >> logs/${STAMP}_failures.log"

echo "on disk: $(ls ../data/sample_size/*.rds 2>/dev/null | wc -l) / $((N_REPS * 18)) files"
[ -s "logs/${STAMP}_failures.log" ] && { echo "FAILURES:"; cat "logs/${STAMP}_failures.log"; }
