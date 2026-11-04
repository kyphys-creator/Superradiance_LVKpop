"""E8: 窓を (mu, f) 面に描く。X 線連星の制限をこのノートの模型で計算し、GW241011 と 1 世代目の条件と重ねる。

模型: E3 の準平衡(si.py)、E4 の相対論的な超放射の率、gamma_BH は主要項(補正 C_BH ありも出す、E7)。
すべて「準平衡で chi_i から chi_f まで回転を失う時間」と年齢を比べる(si.t_spindown、f^-2 に比例)。
  X 線連星(2412.03655 Table III、Refs/notes/window_consistency.md):
    GRS 1915+105: 12.4 Msun、a~ > 0.95、tau = tau_Sal ~ 4.5e7 yr  -> 0.995 から 0.95 まで落ちない
    Cyg X-1(控えめなスピン): 21.2 Msun、a~ = 0.92 (+0.05 -0.07)、4.8e6 yr -> 0.995 から 0.85 まで落ちない
    Cyg X-1(fiducial スピン): 21.2 Msun、a~ > 0.9985、4.8e6 yr -> 0.9995 から 0.9985 まで落ちない
      (chi > 0.995 は連分数法の表の外。表の 3 次の外挿は、chi = 0.9995 の直接の解と 2e-3 で合う(E9)ので、
       相対論的な率を使う。E8 の最初の版は解析式の率で、上限を 29% 小さく出していた)
  GW241011 の主星: 19.6 Msun、0.78 から 0.69 まで、年齢 T2 のあいだ落ちない(f の上限)
  1 世代目: M1、chi0 = 0.7 から (211 の臨界 + delta) まで、年齢 T1 のうちに落ちる(f の下限)
出力: out/e08.json、../figs/data/e08_*.dat
実行: cd research && ~/miniforge3/bin/python3 calc/e08_window.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import si  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")
FREF = 1e13
MUS = np.geomspace(0.4e-12, 1.2e-12, 17)


def f_limit(M, chi_i, chi_f, T, mu, rel=True):
    """準平衡で chi_i -> chi_f にかかる時間が T になる f [GeV](時間は f^-2 に比例)。"""
    si.REL = rel
    t = si.t_spindown(M, chi_i, chi_f, mu, FREF)
    si.REL = False
    if not np.isfinite(t):
        return np.nan
    return FREF * np.sqrt(t / T)


def main():
    os.makedirs(OUT, exist_ok=True)
    res = {}
    curves = {}
    for cb in [False, True]:
        si.CBH_ON = cb
        key = "cbh" if cb else "lo"
        def row_at(mu):
            row = [mu,
                   f_limit(12.4, 0.995, 0.95, 4.5e7, mu),            # GRS 1915+105
                   f_limit(21.2, 0.995, 0.85, 4.8e6, mu),            # Cyg X-1 控えめ
                   f_limit(21.2, 0.9995, 0.9985, 4.8e6, mu),             # Cyg X-1 fiducial
                   f_limit(19.6, 0.78, 0.69, 1e7, mu),               # GW241011, T2 = 1e7
                   f_limit(19.6, 0.78, 0.69, 3e7, mu),               # GW241011, T2 = 3e7
                   f_limit(19.6, 0.78, 0.69, 1e8, mu)]               # GW241011, T2 = 1e8
            for M1 in [10.0, 15.0, 20.0]:
                thr = si.chi_threshold_211(M1, mu)
                for dl in [0.05, 0.18]:
                    row.append(f_limit(M1, 0.7, thr + dl, 6e9, mu) if 0.7 > thr + dl else np.nan)
            return row
        rows = np.array([row_at(mu) for mu in MUS])
        curves[key] = rows
        np.savetxt(os.path.join(FIG, f"e08_window_{key}.dat"), rows, comments="",
                   header="mu GRS CygCons CygFid GW241011_T1e7 GW241011_T3e7 GW241011_T1e8 "
                          "M10_d05 M10_d18 M15_d05 M15_d18 M20_d05 M20_d18", fmt="%.5e")
        names = ["GRS", "CygCons", "CygFid", "GW241011_T1e7", "GW241011_T3e7", "GW241011_T1e8",
                 "1G_M10_d05_T6e9", "1G_M10_d18_T6e9", "1G_M15_d05_T6e9", "1G_M15_d18_T6e9",
                 "1G_M20_d05_T6e9", "1G_M20_d18_T6e9"]
        # 要約は格子点ではなく、ちょうどその mu で計算する(E8 の最初の版は最も近い格子点 0.69e-12 を使った)
        for tag, mu in [("0.60", 0.6e-12), ("0.69", 0.69e-12), ("0.71", 0.71e-12), ("0.89", 0.89e-12)]:
            vals = row_at(mu)
            res[f"mu={tag}e-12_{key}"] = {n: float(vals[k + 1]) for k, n in enumerate(names)}
    si.CBH_ON = False
    with open(os.path.join(OUT, "e08.json"), "w") as f:
        json.dump(res, f, indent=1)
    for k, v in res.items():
        print(k)
        for kk, vv in v.items():
            print(f"   {kk:22s} {vv:.3e}")


if __name__ == "__main__":
    main()
