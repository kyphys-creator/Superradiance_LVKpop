"""E3: 自己相互作用の準平衡で、GW241011 の主星を残しつつ 1 世代目の回転を奪う f はあるか(見積もり)。

準平衡では回転を失う速さが f^2 に比例するので、二つの BH が同じ f で回転を失う時間の比は f によらない。
  t_2G: GW241011 の主星(19.6 Msun)が chi = 0.78 から 90% 区間の下端 0.69 まで落ちる時間
  t_1G: 質量 M1 の 1 世代目が chi0 から (211 の臨界スピン + delta) まで落ちる時間
窓があるには、1 世代目の年齢 T1 と 2 世代目の年齢 T2 が T1/T2 >= R = t_1G/t_2G を満たす必要がある。
出力: out/e03.json、../figs/data/e03_ratio_*.dat
実行: cd research && ~/miniforge3/bin/python3 calc/e03_quasieq.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import si  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")

M2G, CHI2G_OBS, CHI2G_LOW = 19.6, 0.78, 0.69        # GW241011 主星(LVK 文献レビュー T04-28)
MUS = {"0.60": 0.60e-12, "0.71": 0.71e-12, "0.89": 0.89e-12}   # 2609.02678 の 90% 区間
FREF = 1e13                                          # 時間は f^-2 に比例するので基準の f で計算する
res = {}

# ---------- 1. 照合: BGLS の準平衡(主要項の率) ----------
mu_s = 1e-12 / sr.HBAR_EVS
e211, e322 = si.eps_eq(0.01, 0.9, mu_s, 1e15, leading=True)
res["check_BGLS_alpha0.01_a0.9_f1e15"] = {
    "eps211": e211, "eps211_BGLS_55a": 0.25,
    "eps322": e322, "eps322_BGLS_55b": 6.9e-6,
    "tau_sd_yr_mu1e-12": 1.0 / abs(si.chidot(0.01, 0.9, mu_s, 1e15, leading=True)) / sr.YR,
    "tau_sd_BGLS_60_printed_yr": 1e7,
}

# ---------- 2. 年齢の比 R(M1) ----------
M1s = np.linspace(6.0, 40.0, 69)
cases = [(0.5, 0.05), (0.5, 0.18), (0.7, 0.05), (0.7, 0.18)]
res["ratio"] = {}
for tag, mu in MUS.items():
    t2 = si.t_spindown(M2G, CHI2G_OBS, CHI2G_LOW, mu, FREF)
    d = {"t2G_yr_at_f1e13": t2, "alpha_2G": sr.alpha_of(M2G, mu)}
    rows = []
    for M1 in M1s:
        thr = si.chi_threshold_211(M1, mu)
        row = [M1, sr.alpha_of(M1, mu), thr]
        for chi0, dl in cases:
            row.append(si.t_spindown(M1, chi0, thr + dl, mu, FREF) / t2 if chi0 > thr + dl else np.nan)
        rows.append(row)
    rows = np.array(rows)
    np.savetxt(os.path.join(FIG, f"e03_ratio_mu{tag}.dat"), rows, comments="",
               header="M1 alpha thr R_c05_d05 R_c05_d18 R_c07_d05 R_c07_d18", fmt="%.6g")
    for M1 in [8.0, 10.0, 15.0, 20.0, 30.0]:
        i = np.argmin(abs(M1s - M1))
        d[f"M1={M1:g}"] = {"alpha": rows[i, 1], "thr": rows[i, 2],
                           "R(chi0=0.5,d=0.05)": rows[i, 3], "R(0.5,0.18)": rows[i, 4],
                           "R(0.7,0.05)": rows[i, 5], "R(0.7,0.18)": rows[i, 6]}
    res["ratio"][f"mu={tag}e-12"] = d

# ---------- 3. f の目安: GW241011 が T2 のあいだ 0.69 を保つ上限 ----------
mu = MUS["0.71"]
t2 = si.t_spindown(M2G, CHI2G_OBS, CHI2G_LOW, mu, FREF)
res["f_max_GW241011_mu0.71e-12"] = {f"T2={T:.0e}yr": FREF * np.sqrt(t2 / T) for T in [1e6, 1e7, 1e8, 1e9]}

# ---------- 4. 準平衡の近似が成り立つか(GW241011、f = 上限 at T2 = 1e7 yr) ----------
a = sr.alpha_of(M2G, mu)
mu_s = mu / sr.HBAR_EVS
f = res["f_max_GW241011_mu0.71e-12"]["T2=1e+07yr"]
thr = si.chi_threshold_211(M2G, mu)
emax = CHI2G_OBS - thr                                # 自由場で飽和するときの eps211(質量の変化を無視)
e211, e322 = si.eps_eq(a, CHI2G_OBS, mu_s, f)
gs = si.g_sr(a, CHI2G_OBS, mu_s)
eps_crit = 1024 * np.pi * f ** 2 / (9 * a ** 3 * si.MPL_GEV ** 2)      # BGLS 式 88(bosenova)
lnGM2 = np.log(M2G ** 2 * sr.N0)
res["validity_GW241011_f_at_T2=1e7"] = {
    "f_GeV": f, "alpha": a,
    "eps211_eq_over_eps_max": e211 / emax,          # << 1 なら自己相互作用が成長を止める
    "eps211_eq_over_eps_bosenova": e211 / eps_crit,
    "t_to_equilibrium_yr": np.log(e211 * M2G ** 2 * sr.N0) / gs / sr.YR,
    "lnGM2": lnGM2,
}

os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "e03.json"), "w") as fo:
    json.dump(res, fo, indent=1)
print(json.dumps(res, indent=1))
