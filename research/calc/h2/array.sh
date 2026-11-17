#!/bin/bash
# Slurm のジョブ配列の 1 本(submit.sh から sbatch で投げる)。 array.sh <task> <n_tasks> <out_dir>
# 配列の番号 SLURM_ARRAY_TASK_ID(1 から)で受け持ちを決め、run.py に渡す。出力は out_dir/<task>_<番号>.json。
set -u
task=$1; ntask=$2; out=$3
source "$HOME/Superradiance_LVKpop/research/calc/h2/env.sh"
cd "$HOME/Superradiance_LVKpop/research/calc"
echo "=== $task $SLURM_ARRAY_TASK_ID/$ntask job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} on $(hostname) ($(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs)) start $(date '+%F %T')"
python h2/run.py "$task" "$SLURM_ARRAY_TASK_ID" "$ntask" "$out"
rc=$?
echo "=== end $(date '+%F %T') exit $rc"
exit $rc
