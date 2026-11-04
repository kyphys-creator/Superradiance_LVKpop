"""E2: 自由場の超放射で、2 世代目の BH(生まれたときのスピン 0.69)がどこまで回転を失うか。

出力:
  out/e02.json                 照合の数値
  ../figs/data/e02_chif_t{6,8,10}.dat   chi_f(M, mu) の格子(pgfplots の surf/contour 用)
  ../figs/data/e02_regge_*.dat  mu を固定したときの chi_f(M)
実行: cd research && ~/miniforge3/bin/python3 calc/e02_freefield.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

NORM = sr.NORM_N                                 # 占有数の成長率 = 2 x 1908.10370 の Gamma
CHI0 = 0.69                                      # 2 世代目の BH の生まれたときのスピン
LEVELS = [(2, 1, 1), (3, 2, 2), (4, 3, 3)]
T_YR = [1e6, 1e8, 1e10]

res = {"norm": NORM, "chi0": CHI0, "levels": LEVELS}

# ---------- 1. 照合: 半解析と ODE(単一の準位)、見積もり ----------
a = sr.alpha_of(15.0, 7e-13)
g = sr.gamma_over_mu(a, CHI0, 2, 1, 1, NORM)      # 占有数の率 / mu
mu_s = 7e-13 / sr.HBAR_EVS
chic = sr.chi_crit(a, 1)
res["check_15Msun_7e-13eV"] = {
    "alpha": a,
    "Gamma_over_mu_chi0.69": g,
    "efold_yr": 1.0 / (g * mu_s) / sr.YR,
    "n_efolds_lnNmax": float(np.log((CHI0 - chic) * 15.0 ** 2 * sr.N0)),
    "chi_crit_211_omega=mu": chic,
}
# 半解析(飽和で区切る)と ODE を、211 だけで比べる。
# 比べる量: 飽和したあとのスピンと、スピンが中点まで落ちる時刻。
tt = np.logspace(-2, 10, 4001)
comp = []
for M, mu in [(15.0, 7e-13), (5.0, 1e-12), (30.0, 3e-13), (60.0, 1e-13), (10.0, 2e-12)]:
    s1 = sr.evolve_free(M, CHI0, mu, LEVELS[:1], 1e10, NORM)
    c_ode = sr.chi_at(s1, tt, 1)
    c_seq = sr.chi_final_seq(M, CHI0, mu, tt, LEVELS[:1], NORM)
    mid = 0.5 * (CHI0 + c_ode[-1])
    t_ode = tt[np.argmax(c_ode < mid)] if c_ode[-1] < CHI0 - 1e-3 else None
    t_seq = tt[np.argmax(c_seq < mid)] if c_seq[-1] < CHI0 - 1e-3 else None
    comp.append({"M": M, "mu": mu, "alpha": sr.alpha_of(M, mu),
                 "chif_ode": float(c_ode[-1]), "chif_seq": float(c_seq[-1]),
                 "t_mid_ode_yr": t_ode, "t_mid_seq_yr": t_seq})
res["seq_vs_ode_211"] = comp

# GW241011 の主星(m1 = 19.6 Msun、chi1 = 0.78、LVK 文献レビュー T04-28)を mu = 7.1e-13 eV の雲に置いたとき、
# スピンが 90% 区間の下端 0.69 を下回るまでの時間。
tt2 = np.logspace(-2, 10, 6001)
c241 = sr.chi_final_seq(19.6, 0.78, 7.1e-13, tt2, LEVELS)
res["GW241011_primary_mu7.1e-13"] = {
    "alpha": sr.alpha_of(19.6, 7.1e-13),
    "t_below_0.69_yr": float(tt2[np.argmax(c241 < 0.69)]) if c241[-1] < 0.69 else None,
    "chif_1e10yr": float(c241[-1]),
}

# ---------- 2. chi_f(M, mu) の格子 ----------
Ms = np.logspace(np.log10(3.0), np.log10(150.0), 121)
mus = np.logspace(-14, -11, 121)
chif = np.zeros((len(T_YR), len(mus), len(Ms)))
for j, mu in enumerate(mus):
    for i, M in enumerate(Ms):
        chif[:, j, i] = sr.chi_final_seq(M, CHI0, mu, T_YR, LEVELS, NORM)
for it, t in enumerate(T_YR):
    tag = int(round(np.log10(t)))
    with open(os.path.join(FIG, f"e02_chif_t{tag}.dat"), "w") as f:
        f.write("lgM lgmu chif\n")
        for j, mu in enumerate(mus):
            for i, M in enumerate(Ms):
                f.write(f"{np.log10(M):.5f} {np.log10(mu):.5f} {chif[it, j, i]:.5f}\n")
            f.write("\n")

# 2 世代目が 0.5 を下回る mu の範囲(質量帯ごと)
bands = {"10-20": (10.0, 20.0), "45-100": (45.0, 100.0)}
res["mu_range_chif_below_0.5"] = {}
for it, t in enumerate(T_YR):
    d = {}
    for name, (lo, hi) in bands.items():
        sel = (Ms >= lo) & (Ms <= hi)
        # 帯の中のすべての質量で 0.5 を下回る mu
        allbelow = np.all(chif[it][:, sel] < 0.5, axis=1)
        d[name] = [float(mus[allbelow].min()), float(mus[allbelow].max())] if allbelow.any() else None
    res["mu_range_chif_below_0.5"][f"t={t:.0e}yr"] = d

# chi_f = 0.5 の等高線(pgfplots 用、線ごとに空行で区切る)
import contourpy  # noqa: E402
lgM, lgmu = np.meshgrid(np.log10(Ms), np.log10(mus))
for it, t in enumerate(T_YR):
    tag = int(round(np.log10(t)))
    lines = contourpy.contour_generator(lgM, lgmu, chif[it]).lines(0.5)
    with open(os.path.join(FIG, f"e02_c05_t{tag}.dat"), "w") as f:
        f.write("lgM lgmu\n")
        for ln in lines:
            for x, y in ln:
                f.write(f"{x:.5f} {y:.5f}\n")
            f.write("\n")

# ---------- 3. Regge 面の軌跡: mu を固定した chi_f(M) ----------
Mf = np.logspace(np.log10(3.0), np.log10(150.0), 400)
for mu in [1e-13, 7.1e-13, 3e-12]:
    rows = []
    for M in Mf:
        rows.append([M] + list(sr.chi_final_seq(M, CHI0, mu, T_YR, LEVELS, NORM)))
    np.savetxt(os.path.join(FIG, f"e02_regge_mu{mu:.1e}.dat"), np.array(rows),
               header="M chif_t6 chif_t8 chif_t10", comments="", fmt="%.5f")

with open(os.path.join(OUT, "e02.json"), "w") as f:
    json.dump(res, f, indent=1)
print(json.dumps(res, indent=1))
