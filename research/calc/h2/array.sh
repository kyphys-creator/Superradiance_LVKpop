#!/bin/bash
# Slurm のジョブ配列の 1 本(submit.sh から sbatch で投げる)。 array.sh <task> <n_tasks> <out_dir>
# 配列の番号 SLURM_ARRAY_TASK_ID(1 から)で受け持ちを決め、run.py に渡す。出力は out_dir/<task>_<番号>.json。
# 終了の合図(SIGTERM。時間の上限や scancel で Slurm が送る)を受けたら python に渡し、python が保存し終わるのを待ってから
# 終わる。bash が合図で先に終わると、Slurm は残った python を保存の前に強制終了する(2026-11-17、止めた本で見つけた)。
set -u
task=$1; ntask=$2; out=$3
source "$HOME/Superradiance_LVKpop/research/calc/h2/env.sh"
cd "$HOME/Superradiance_LVKpop/research/calc"
echo "=== $task $SLURM_ARRAY_TASK_ID/$ntask job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} on $(hostname) ($(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs)) start $(date '+%F %T')"
python h2/run.py "$task" "$SLURM_ARRAY_TASK_ID" "$ntask" "$out" &
pid=$!
trap 'echo "=== SIGTERM $(date "+%F %T"): passed to python, waiting for it to save"; kill -TERM "$pid" 2>/dev/null' TERM
wait "$pid"
rc=$?
while kill -0 "$pid" 2>/dev/null; do          # wait が合図で中断されたら、python が終わるまで待ち直す
    wait "$pid"
    rc=$?
done
echo "=== end $(date '+%F %T') exit $rc"
exit $rc
