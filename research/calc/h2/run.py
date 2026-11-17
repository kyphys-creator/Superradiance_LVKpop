"""Hoffman2 のジョブ配列の中身。 python h2/run.py <task> <番号(1 から)> <本数> <out_dir>
task:
  corrgrid  kinn.corr_grid の 55 例(kinn._corr_job、補正のあり・なしを evolveN_checked で)。数秒で終わる例が多いので、
            本数を 11 ほどにして 1 本に 5 例ずつまとめる(1 分以下の本をたくさん流さない、Hoffman2 の作法)。
            Mac の kinn.py corrcheck の cases と同じ形の結果を出すので、Mac と照らせる。
  merge     python h2/run.py merge <task> <out_dir>: <task>_*.json を一つにまとめて <out_dir>/<task>_all.json に。
出力: <out_dir>/<task>_<番号>.json(すでにあれば飛ばす。途中で切れても投げ直せば続きから)。"""
import glob
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def cases(task):
    import kinn
    if task == "corrgrid":
        return kinn.corr_grid(), kinn._corr_job
    raise ValueError(f"unknown task {task}")


def main():
    if sys.argv[1] == "merge":
        task, out = sys.argv[2], sys.argv[3]
        rows = []
        for p in sorted(glob.glob(os.path.join(out, f"{task}_[0-9]*.json"))):
            rows.extend(json.load(open(p)))
        json.dump(rows, open(os.path.join(out, f"{task}_all.json"), "w"))
        print(f"{len(rows)} rows -> {task}_all.json")
        return
    task, tid, ntask, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    path = os.path.join(out, f"{task}_{tid:05d}.json")
    if os.path.exists(path):
        print("already done:", path)
        return
    allc, fn = cases(task)
    mine = allc[tid - 1::ntask]                     # 番号 tid の受け持ち(本数で割った余り)
    rows = []
    for c in mine:
        t0 = time.time()
        r = fn(c)
        r["host_secs"] = time.time() - t0
        rows.append(r)
    tmp = path + ".tmp"
    json.dump(rows, open(tmp, "w"), default=float)
    os.replace(tmp, path)
    print(f"{len(rows)} cases -> {path}")


if __name__ == "__main__":
    main()
