"""E9: 独立な検証(理論・コード・観測)で見つかった点を反映し、窓を楽観的な仮定と現実的な仮定で描き直す。

仮定の組(すべて準平衡、相対論的な超放射の率、gamma_BH は主要項。C_BH を入れた値も出す):
  楽観的(E8 と同じ): GRS 1915+105 は 0.995 -> 0.95 を tau_Sal = 4.5e7 yr で保つ。
                      GW241011 の主星は 0.78 で生まれ、T2 = 3e7 yr のあいだ 0.69 を保つ。
  現実的:             GRS は 0.995 -> 0.87 を BH の年齢 3e9 yr で保つ(Reid+14 の系統誤差込み 3 sigma の下限 0.92 から、
                      伴星が軽いので降着で戻せるのは高々 ~0.05 として差し引く)。
                      GW241011 の主星は 0.72(2 世代目の予言の上端)で生まれ、T2 = 5e7 yr(重い 2 世代目の遅れの最小値の目安)。
  1 世代目: chi0 = 0.7 から (211 の臨界 + delta) まで、T1 = 6e9 yr で落ちる(f の下限)。
ほかに: 質量の変化の補正(si.MASS_CORR)、Cyg X-1 fiducial の別の判定、E6 の訂正の chi_mu - chi_c、ずれ機構の暗黒物質。
出力: out/e09.json、../figs/data/e09_window_{opt,real}.dat
実行: cd research && ~/miniforge3/bin/python3 calc/e09_review.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import si  # noqa: E402
import sr  # noqa: E402
from e08_window import f_limit, MUS  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")

SCEN = {
    "opt": {"grs": (0.995, 0.95, 4.5e7), "gw": (0.78, 3e7)},
    "real": {"grs": (0.995, 0.87, 3e9), "gw": (0.72, 5e7)},
}


def bounds(mu, scen):
    g0, g1, gT = SCEN[scen]["grs"]
    c0, T2 = SCEN[scen]["gw"]
    row = [mu, f_limit(12.4, g0, g1, gT, mu), f_limit(19.6, c0, 0.69, T2, mu)]
    for M1 in [10.0, 15.0, 20.0, 25.0]:
        thr = si.chi_threshold_211(M1, mu)
        for dl in [0.18, 0.05]:
            row.append(f_limit(M1, 0.7, thr + dl, 6e9, mu) if 0.7 > thr + dl else np.nan)
    return row


def chi_at_omegaH(target):
    lo, hi = 0.0, 0.999999
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if sr.omegaH(mid) > target else (mid, hi)
    return 0.5 * (lo + hi)


def main():
    res = {}
    hdr = "mu GRS GW241011 M10_d18 M10_d05 M15_d18 M15_d05 M20_d18 M20_d05 M25_d18 M25_d05"
    names = hdr.split()[1:]
    for cb in [False, True]:
        si.CBH_ON = cb
        key = "cbh" if cb else "lo"
        for scen in SCEN:
            rows = np.array([bounds(mu, scen) for mu in MUS])
            if not cb:
                np.savetxt(os.path.join(FIG, f"e09_window_{scen}.dat"), rows, header=hdr, comments="", fmt="%.5e")
            v = bounds(0.71e-12, scen)
            res[f"{scen}_{key}_mu0.71"] = {n: float(v[k + 1]) for k, n in enumerate(names)}
    si.CBH_ON = False

    # 仮定を一つずつ変えたときの上限(mu = 0.71e-12、主要項)
    mu = 0.71e-12
    res["variants_mu0.71"] = {
        "GRS_0.95_tauSal": f_limit(12.4, 0.995, 0.95, 4.5e7, mu),
        "GRS_0.92_tauSal": f_limit(12.4, 0.995, 0.92, 4.5e7, mu),
        "GRS_0.95_3Gyr": f_limit(12.4, 0.995, 0.95, 3e9, mu),
        "GRS_0.87_3Gyr": f_limit(12.4, 0.995, 0.87, 3e9, mu),
        "GRS_0.87_5Gyr": f_limit(12.4, 0.995, 0.87, 5e9, mu),
        **{f"GW241011_born{c0}_T2={T:.0e}": f_limit(19.6, c0, 0.69, T, mu)
           for c0 in [0.78, 0.75, 0.72, 0.70] for T in [1e7, 3e7, 5e7, 3e8]},
        "CygFid_0.9995to0.9985_rel": f_limit(21.2, 0.9995, 0.9985, 4.8e6, mu),
        "CygFid_0.9995to0.9985_analytic": f_limit(21.2, 0.9995, 0.9985, 4.8e6, mu, rel=False),
        "CygFid_0.998to0.997_rel": f_limit(21.2, 0.998, 0.997, 4.8e6, mu),
    }
    # 質量の変化の補正: 回転を失う時間の比(補正あり / なし)
    si.REL = True
    mc = {}
    for name, (M, a, b) in {"GW241011": (19.6, 0.78, 0.69), "GRS": (12.4, 0.995, 0.95),
                            "1G_M10": (10.0, 0.7, si.chi_threshold_211(10.0, mu) + 0.18),
                            "1G_M15": (15.0, 0.7, si.chi_threshold_211(15.0, mu) + 0.18)}.items():
        si.MASS_CORR = False
        t0 = si.t_spindown(M, a, b, mu, 1e13)
        si.MASS_CORR = True
        t1 = si.t_spindown(M, a, b, mu, 1e13)
        mc[name] = t1 / t0
    si.MASS_CORR = False
    si.REL = False
    res["mass_correction_time_ratio_mu0.71"] = mc
    # E6 の訂正: m=1 の吸収で上がりうるスピンの上限 chi_mu(Omega_H = mu)と 211 の臨界の差
    res["chi_mu_minus_chi_c211_mu0.71"] = {}
    for name, M in {"GW241011": 19.6, "GRS": 12.4, "CygX-1": 21.2, "1G_10Msun": 10.0}.items():
        a = sr.alpha_of(M, mu)
        c211 = chi_at_omegaH(a * sr.omega_over_mu(a, 2, 1, 1, 0.5))
        cmu = chi_at_omegaH(a)
        res["chi_mu_minus_chi_c211_mu0.71"][name] = {"alpha": a, "chi_c211": c211, "chi_mu": cmu, "diff": cmu - c211}
    # ずれ機構: 放射優勢で H = mu(または 3H = mu)のとき振動が始まり、その後 n/s が保存するとして
    # Omega h^2 = 0.12 となる f theta(g* = 10.75、M_pl = 1.22e19 GeV、rho_DM/s0 = 4.37e-10 GeV)
    mis = {}
    for lab, k in [("H=mu", 1.0), ("3H=mu", 3.0)]:
        m = 0.71e-12 * 1e-9                          # GeV
        gs = 10.75
        Tosc = np.sqrt(m / k * si.MPL_GEV / (1.66 * np.sqrt(gs)))
        s = 2 * np.pi ** 2 / 45 * gs * Tosc ** 3
        mis[lab] = {"T_osc_MeV": Tosc * 1e3, "f_theta_GeV": float(np.sqrt(4.37e-10 * s / (0.5 * m ** 2)))}
    res["misalignment_mu0.71"] = mis

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "e09.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
