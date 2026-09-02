#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")"
RSCRIPT=${RSCRIPT:-Rscript}
N_SIMS=${1:-500}; N_CORES=${2:-16}; CHUNK=${3:-5}
shift 3 2>/dev/null || true

ALL=(sc_snull sc_med sc_null sc_modl sc_direct sc_mixed)
if [ "$#" -gt 0 ]; then SCEN=("$@"); else SCEN=("${ALL[@]}"); fi
for s in "${SCEN[@]}"; do
  case " ${ALL[*]} " in *" $s "*) ;; *) echo "unknown scenario '$s'"; exit 2;; esac
done
NU=4

mkdir -p logs ../data/confounding
STAMP=$(date +%Y%m%d_%H%M%S)
TASKS=logs/tasks_conf_${STAMP}.txt; : > "$TASKS"
for ((from=1; from<=N_SIMS; from+=CHUNK)); do
  to=$((from+CHUNK-1)); (( to > N_SIMS )) && to=$N_SIMS
  for ui in $(seq 1 $NU); do for s in "${SCEN[@]}"; do echo "$s $ui $from $to" >> "$TASKS"; done; done
done
echo "queued $(wc -l < "$TASKS") tasks on $N_CORES workers"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
xargs -a "$TASKS" -P "$N_CORES" -L 1 \
  sh -c "$RSCRIPT worker_confound.R \$0 \$1 \$2 \$3 >> logs/${STAMP}_\$0_u\$1.log 2>&1 || echo \"FAILED: \$0 \$1 \$2 \$3\" >> logs/${STAMP}_failures.log"

echo "on disk: $(ls ../data/confounding/*.rds 2>/dev/null | wc -l) / $((N_SIMS * NU * 6)) files"
[ -s "logs/${STAMP}_failures.log" ] && { echo "FAILURES:"; cat "logs/${STAMP}_failures.log"; }
