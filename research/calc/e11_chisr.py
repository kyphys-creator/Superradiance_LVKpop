"""E11: 集団解析に渡す chi_SR(準平衡と自由場をつないだもの)の照合と計算時間。

1. E5 の ODE(GW241011 の主星、mu = 0.71e-12 eV、chi0 = 0.78、f = 1e11--1e15 GeV、T = 1e6, 1e7, 1e8 yr)と比べる。
2. 自由場の領域(f = 1e19 GeV)で、E2 の半解析(sr.chi_final_seq、相対論的な率、211 と 322)と比べる。
3. f を 1e10--1e19 GeV で走査し、領域の境目でスピンが飛ばないかを見る(図のデータ)。
4. 1e6 個の標本の計算時間。
出力: out/e11.json、../figs/data/e11_scan.dat
実行: cd research && ~/miniforge3/bin/python3 calc/e11_chisr.py
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")


def main():
    res = {}
    tab = chisr.build()
    mu = 0.71e-12

    # ---------- 1. E5 の ODE と ----------
    d = np.genfromtxt(os.path.join(FIG, "e05_chi_vs_f.dat"), names=True)
    sel = d["f"] <= 1e15
    diffs = {}
    for T, col in [(1e6, "chi_t6"), (1e7, "chi_t7"), (1e8, "chi_t8")]:
        c, r = chisr.chi_sr(tab, 19.6, 0.78, T, mu, d["f"][sel])
        diffs[f"T={T:.0e}"] = {"max_abs_diff": float(np.max(np.abs(c - d[col][sel]))),
                               "regimes": sorted(set(int(x) for x in r))}
    res["vs_E5_ODE"] = diffs

    # ---------- 2. 自由場で E2 の半解析と ----------
    rng = np.random.default_rng(11)
    n = 400
    Ms = rng.uniform(5, 60, n)
    mus = 10 ** rng.uniform(np.log10(0.2e-12), np.log10(2e-12), n)
    c0s = rng.uniform(0.3, 0.99, n)
    Ts = 10 ** rng.uniform(3, 10, n)
    cu, ru = chisr.chi_sr(tab, Ms, c0s, Ts, mus, 1e19)
    cs = np.array([sr.chi_final_seq(M, c0, m_, [T], [(2, 1, 1), (3, 2, 2)], rel=True)[0]
                   for M, c0, m_, T in zip(Ms, c0s, mus, Ts)])
    dd = np.abs(cu - cs)
    res["vs_E2_seq_free"] = {"n": n, "median_abs_diff": float(np.median(dd)),
                             "p95_abs_diff": float(np.percentile(dd, 95)), "max_abs_diff": float(dd.max()),
                             "frac_agree_1e-2": float(np.mean(dd < 1e-2)),
                             "note": "残りの差は、臨界を水素原子の主要項の E で解く(chi_sr)か alpha^5 までの E で解く(seq)かの違い"}

    # ---------- 3. f の走査 ----------
    fs = np.logspace(10, 19, 181)
    rows = [fs]
    hdr = ["f"]
    for (M, c0, lab) in [(19.6, 0.78, "GW241011"), (10.0, 0.7, "M10"), (15.0, 0.7, "M15")]:
        for T in [1e7, 1e9]:
            c, r = chisr.chi_sr(tab, M, c0, T, mu, fs)
            rows += [c, r]
            hdr += [f"chi_{lab}_T{int(np.log10(T))}", f"reg_{lab}_T{int(np.log10(T))}"]
            jumps = np.abs(np.diff(c))
            res.setdefault("scan_max_step", {})[f"{lab}_T{T:.0e}"] = float(jumps.max())
    np.savetxt(os.path.join(FIG, "e11_scan.dat"), np.column_stack(rows), header=" ".join(hdr),
               comments="", fmt="%.6e")

    # ---------- 4. 計算時間 ----------
    nev = 1_000_000
    M = rng.uniform(5, 80, nev)
    c0 = rng.uniform(0.0, 0.99, nev)
    T = 10 ** rng.uniform(6, 10, nev)
    mus = rng.uniform(0.5e-12, 1.0e-12, nev)
    fs = 10 ** rng.uniform(11, 18, nev)
    t0 = time.time()
    c, r = chisr.chi_sr(tab, M, c0, T, mus, fs)
    res["timing_1e6_s"] = time.time() - t0
    res["regime_fractions_random"] = {str(k): float(np.mean(r == k)) for k in range(4)}

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "e11.json"), "w") as fo:
        json.dump(res, fo, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
