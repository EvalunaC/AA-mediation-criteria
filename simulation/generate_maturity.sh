#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")"
RSCRIPT=${RSCRIPT:-Rscript}
N_SIMS=${1:-500}; N_CORES=${2:-16}; CHUNK=${3:-5}
shift 3 2>/dev/null || true

ALL=(sc_null sc_snull sc_modl sc_med sc_direct sc_mixed)
if [ "$#" -gt 0 ]; then SCEN=("$@"); else SCEN=("${ALL[@]}"); fi
for sc in "${SCEN[@]}"; do
  case " ${ALL[*]} " in *" $sc "*) ;; *) echo "unknown scenario '$sc'"; exit 2;; esac
done

mkdir -p logs ../data/maturity
STAMP=$(date +%Y%m%d_%H%M%S)
TASKS=logs/tasks_${STAMP}.txt; : > "$TASKS"
for ((from=1; from<=N_SIMS; from+=CHUNK)); do
  to=$((from+CHUNK-1)); (( to > N_SIMS )) && to=$N_SIMS
  for sc in "${SCEN[@]}"; do echo "$sc $from $to" >> "$TASKS"; done
done
echo "queued $(wc -l < "$TASKS") tasks on $N_CORES workers"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
xargs -a "$TASKS" -P "$N_CORES" -L 1 \
  sh -c "$RSCRIPT worker_maturity.R \$0 \$1 \$2 >> logs/${STAMP}_\$0.log 2>&1 || echo \"FAILED: \$0 \$1 \$2\" >> logs/${STAMP}_failures.log"

echo "on disk: $(ls ../data/maturity/*.rds 2>/dev/null | wc -l) / $((N_SIMS * 6)) files"
[ -s "logs/${STAMP}_failures.log" ] && { echo "FAILURES:"; cat "logs/${STAMP}_failures.log"; }
