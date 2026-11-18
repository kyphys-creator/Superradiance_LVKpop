#!/bin/bash
# ジョブ配列を投げる(Hoffman2 の Slurm のログインノード hoffman2.oarc.ucla.edu で)。
#   [SRLVK_CAP=<呼び出しの上限>] h2/submit.sh <task> <n_tasks> [time=4:00:00] [mem=2G] [cpus=1]
# 1 本は cpus コア・メモリ mem(1 本の運動論は 0.1-0.2 GB。cpus 並列なら mem はその倍に余裕を足す。--mem を書かないと
# ノードのメモリを全部取るので必ず書く)。環境変数(SRLVK_CAP など)はジョブに引き継がれる。
# 公式の作法: キャンパスのジョブは 24 時間まで、配列の番号は 1000 まで、1 分以下で終わる本を短時間にたくさん流さない
# (n_tasks を減らして 1 本に複数の例をまとめる)。出力とログは $SCRATCH/srlvk/<task>/(14 日さわらないと消える。
# 消さないために日付だけ更新するのは禁止。要る結果は手元か $HOME に移す)。
set -euo pipefail
task=$1; ntask=$2; time=${3:-4:00:00}; mem=${4:-2G}; cpus=${5:-1}
out=$SCRATCH/srlvk/$task
mkdir -p "$out/log"
here=$(cd "$(dirname "$0")" && pwd)
sbatch --job-name="srlvk_$task" --array="1-$ntask" --time="$time" --mem="$mem" --cpus-per-task="$cpus" --ntasks=1 \
       --output="$out/log/%A_%a.out" "$here/array.sh" "$task" "$ntask" "$out"
