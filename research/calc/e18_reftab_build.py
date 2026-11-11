"""E18: chi_SR の値の表を、基準の運動論(fp3.evolve3m)から直接作る(表の形は calc/reftab.py の levels3)。

節(細かい格子。偶数番目だけを取ると粗い格子になる):
  alpha: 0.008--0.5 の対数で 47 点
  s(初期スピンの臨界に対する位置): 322 から育つ側 SA、211 から育つ側 SB(s - 1)。粗い 13 点と、その対数の中点で 25 点
  f: 1e8--1e20.5 GeV を 0.25 dex で 51 点
  M: 1、300 Msun(L = ln(N0 M^2) を通してだけ効く)
211 の臨界が 0.995 を越える alpha では 211 側の節は作らない。途中で止まった節は、止まった tau までの準位を持つ。
粗い格子の節を先に解き、10000 節ごとに途中の結果を保存する(--resume で続きから)。
出力: out/reftab.npz(途中は out/reftab_partial.npz)
実行: cd research && ~/miniforge3/bin/python3 calc/e18_reftab_build.py [--resume] [--coarse-only]
"""
import os
import sys
import time
import warnings
from multiprocessing import Pool

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import reftab as R  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
PARTIAL = os.path.join(OUT, "reftab_partial.npz")


def _fine(coarse):
    c = np.asarray(coarse, float)
    return np.sort(np.concatenate([c, np.sqrt(c[1:] * c[:-1])]))


A_NODES = np.geomspace(0.008, 0.5, 47)
SA = _fine([1e-3, 3e-3, 0.01, 0.03, 0.07, 0.15, 0.25, 0.4, 0.55, 0.7, 0.85, 0.95, 0.999])     # s(322 側)
SB = _fine([1e-3, 3e-3, 0.01, 0.03, 0.07, 0.15, 0.25, 0.4, 0.55, 0.7, 0.85, 0.95, 1.0])       # s - 1(211 側)
LF_NODES = np.round(np.arange(8.0, 20.51, 0.25), 6)
M_NODES = np.array([1.0, 300.0])
_TAB = None


def _node(args):
    global _TAB
    if _TAB is None:
        _TAB = chisr.build(use_reftab=False)
    piece, ia, js, kf, il = args
    a = A_NODES[ia]
    s = SA[js] if piece == 0 else 1 + SB[js]
    c1 = float(chisr.chic_vec(a, 1, tab=_TAB))
    if piece == 1 and c1 >= chisr.CHI_TOP:
        return args, None, None, 2, 0.0, 0.0, 0.0
    t0 = time.time()
    tau, chi, ok, tlast, c1, c2, chi0 = R.node_curve(a, s, 10 ** LF_NODES[kf], M_NODES[il], _TAB, max_nfev=600000)
    lt, rz, pm = R.levels3(tau, chi, s, c1, c2)
    return args, lt.astype(np.float32), rz.astype(np.float32), (0 if ok else 1), tlast, pm, time.time() - t0


def _save(path, arrays, done):
    np.savez(path, A_NODES=A_NODES, SA=SA, SB=SB, LF_NODES=LF_NODES, M_NODES=M_NODES, P_LEV=R.P_LEV2, G_LT=R.G_LT,
             TAU_MAX=R.TAU_MAX, done=done, **arrays)


def main():
    shp = (2, len(A_NODES), len(SA), len(LF_NODES), len(M_NODES))
    arrays = {"LT": np.full(shp + (len(R.P_LEV2),), np.nan, np.float32),
              "RG": np.full(shp + (len(R.G_LT),), np.nan, np.float32),
              "status": np.full(shp, -1, np.int8), "tau_valid": np.zeros(shp), "pmax": np.zeros(shp),
              "secs": np.zeros(shp)}
    done = np.zeros(shp, bool)
    if "--resume" in sys.argv and os.path.exists(PARTIAL):
        z = np.load(PARTIAL)
        for k in arrays:
            arrays[k] = z[k]
        done = z["done"]
        print(f"resume: {done.sum()} nodes already done", flush=True)
    idx = [t for t in np.ndindex(shp) if not done[t]]
    coarse = [t for t in idx if t[1] % 2 == 0 and t[2] % 2 == 0 and t[3] % 2 == 0]
    rest = [t for t in idx if not (t[1] % 2 == 0 and t[2] % 2 == 0 and t[3] % 2 == 0)]
    tasks = coarse + ([] if "--coarse-only" in sys.argv else rest)
    print(f"{len(tasks)} nodes to run ({len(coarse)} coarse first)", flush=True)
    t0 = time.time()
    n_coarse_left = len(coarse)
    with Pool(10) as p:
        for n, (args, lt, rz, st, tl, pm, sec) in enumerate(p.imap(_node, tasks, chunksize=4)):
            if lt is not None:
                arrays["LT"][args], arrays["RG"][args] = lt, rz
            arrays["status"][args], arrays["tau_valid"][args] = st, tl
            arrays["pmax"][args], arrays["secs"][args] = pm, sec
            done[args] = True
            if n + 1 == n_coarse_left or (n + 1) % 10000 == 0:
                _save(PARTIAL, arrays, done)
                st_all = arrays["status"][done]
                print(f"{n + 1}/{len(tasks)} {time.time() - t0:.0f}s: ok {np.sum(st_all == 0)}, stopped {np.sum(st_all == 1)}, "
                      f"skipped {np.sum(st_all == 2)}" + ("  [coarse subset complete]" if n + 1 == n_coarse_left else ""),
                      flush=True)
    _save(PARTIAL, arrays, done)
    if done.all():
        _save(os.path.join(OUT, "reftab.npz"), arrays, done)
        st = arrays["status"]
        print(f"done in {time.time() - t0:.0f}s: ok {np.sum(st == 0)}, stopped {np.sum(st == 1)}, skipped {np.sum(st == 2)}; "
              f"max p {arrays['pmax'].max():.3f}; slowest {arrays['secs'].max():.0f}s; saved out/reftab.npz", flush=True)


if __name__ == "__main__":
    main()
