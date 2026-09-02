#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")"
RSCRIPT=${RSCRIPT:-Rscript}
N_CORES=${1:-4}

mkdir -p logs ../data/case_interim
STAMP=$(date +%Y%m%d_%H%M%S)
TASKS=logs/tasks_case_${STAMP}.txt; : > "$TASKS"
for cohort in "309 all" "309 KRAS wt"; do
  for cut in 6 9 12 18 24 final; do
    echo "${cohort}|${cut}" >> "$TASKS"
  done
done
echo "queued $(wc -l < "$TASKS") cells on $N_CORES workers"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export RSCRIPT STAMP
xargs -a "$TASKS" -P "$N_CORES" -d '\n' -I@LINE@ bash -c '
  line="@LINE@"; cohort="${line%%|*}"; cut="${line##*|}"
  "$RSCRIPT" worker_case.R "$cohort" "$cut" >> "logs/${STAMP}_case.log" 2>&1 ||
    echo "FAILED: $line" >> "logs/${STAMP}_failures.log"'

echo "on disk: $(ls ../data/case_interim/*.rds 2>/dev/null | wc -l) / 12 cells"
[ -s "logs/${STAMP}_failures.log" ] && { echo "FAILURES:"; cat "logs/${STAMP}_failures.log"; }
