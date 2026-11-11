"""E17: m = 3(433)の超放射の印と、433 の相対論的な表をやめた理由。

1. 433 の連分数法: 倍精度と多倍長(mpmath)の比較。Im omega / Re omega が 1e-16 ほどで、倍精度では alpha = 0.30 付近の
   高いスピンで 4 倍ずれる(表を補間すると率の符号まで誤る)。
2. 433 の臨界スピンの表が alpha について単調か。
3. mu = 0.71e-12 eV で、433 が 1 個の量子から(自由場で)飽和するまでの時間。
4. 印の割合: 1 世代目の自然スピン Beta(0.69, 2.28)、質量 8--50 と 8--80 Msun、mu = 0.71e-12 と 0.89e-12 eV、6 Gyr。
出力: out/e17.json
実行: cd research && ~/miniforge3/bin/python3 calc/e17_m3_flag.py
"""
import json
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import leaver  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


def main():
    res = {}
    # 1. 433 の精度
    prec = []
    for a, c in [(0.30, 0.995), (0.32, 0.95), (0.35, 0.95), (0.40, 0.95)]:
        wd = leaver.bound_state(a, c, 4, 3, 3)
        wm = leaver.bound_state_mp(a, c, 4, 3, 3, dps=40)
        prec.append({"alpha": a, "chi": c, "im_double": float(wd.imag), "im_mp": float(wm.imag),
                     "im_analytic": float(sr.gamma_over_mu(a, c, 4, 3, 3) * a),
                     "im_over_re": float(abs(wm.imag / wm.real))})
    res["leaver_433_precision"] = prec

    tab = chisr.build()
    # 2. 臨界の表の単調性
    res["thr_nonmonotone_alpha"] = {k: [float(x) for x in tab["thr_alpha"][np.flatnonzero(np.diff(tab[k]) < -1e-9)]]
                                    for k in ("thr1", "thr2", "thr3")}

    # 3. 433 の飽和の時間
    mu = 0.71e-12
    mus = mu / sr.HBAR_EVS
    rows = []
    for M in [30, 35, 40, 45, 50, 55, 60]:
        a = sr.alpha_of(M, mu)
        c2 = float(chisr.chic_vec(a, 2, tab=tab))
        c3 = float(chisr.chic_vec(a, 3, tab=tab))
        x = 0.5 * (c2 + c3)
        g = float(chisr.g_vec(a, x, 4, 3, 3)) * mus
        t = np.log((x - c3) / 3 * sr.N0 * M ** 2) / g / sr.YR
        cc, rr = chisr.chi_sr(tab, M, x, 6e9, mu, 1e12)
        rows.append({"M": M, "alpha": round(a, 4), "chi_c433": round(c3, 4), "chi_c322": round(c2, 4), "chi0": round(x, 4),
                     "t_sat_433_yr": t, "regime_6Gyr_f1e12": int(rr)})
    res["t_sat_433"] = rows

    # 4. 印の割合
    rng = np.random.default_rng(12)
    n = 200000
    out = {}
    for Mmax in (50, 80):
        Mp = rng.uniform(8, Mmax, n)
        cp = np.minimum(rng.beta(0.69, 2.28, n), 0.999)
        for muv in (0.71e-12, 0.89e-12):
            ap = sr.alpha_of(Mp, muv)
            below322 = (cp <= chisr.chic_vec(ap, 2, tab=tab)) & (ap <= chisr.ALPHA_MAX)
            for f in (1e12, 1e16, 1e19):
                _, r = chisr.chi_sr(tab, Mp, cp, 6e9, muv, f)
                d = {str(k): round(float(np.mean(r == k)), 4) for k in range(5)}
                d["3_from_m3"] = round(float(np.mean((r == 3) & below322)), 4)   # 211・322 の臨界より下で 433 の印がついたもの
                out[f"M8-{Mmax}_mu{muv:.2e}_f{f:.0e}"] = d
    res["regime_fractions_T6Gyr"] = out

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "e17.json"), "w") as fo:
        json.dump(res, fo, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
