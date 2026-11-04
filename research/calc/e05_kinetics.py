"""E5: 自己相互作用のある二準位(211・322)の運動論を ODE で解き、E3 の準平衡と照合する。

式: BGLS(arXiv:2011.11646)の式 34--36 の主な項(Refs/notes/theory_rates.md §4.3)。BH の質量は固定。
  de211/dt = g_SR211 e211 - 2 g_BH e211^2 e322 + g_inf e322^2 e211 - 2 g_GW e211^2 + g_tr e211 e322
  de322/dt = g_SR322 e322 + g_BH e211^2 e322 - 2 g_inf e322^2 e211 - g_tr e211 e322
  dchi/dt  = -g_SR211 e211 - 2 g_SR322 e322
  g_SR は E4 の相対論的な率(占有数の率)。g_GW = 1e-2 alpha^14 mu(211 の消滅、Table I)、
  g_tr = 5e-6 alpha^10 mu(322 -> 211 の重力波の遷移、式 101・Table IV)。
変数は ln e211、ln e322、chi。陰的な解法(Radau)。初期値は各準位 1 個の量子(e = 1/(G M^2))。
出力: out/e05.json、../figs/data/e05_*.dat
実行: cd research && ~/miniforge3/bin/python3 calc/e05_kinetics.py
"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(__file__))
import si  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")
YR = sr.YR
M2G, CHI0, MU = 19.6, 0.78, 0.71e-12


def evolve(M, chi0, mu_eV, f_gev, t_end_yr, t_eval_yr=None):
    mu_s = mu_eV / sr.HBAR_EVS
    a = sr.alpha_of(M, mu_eV)
    gGW = 1e-2 * a ** 14 * mu_s
    gtr = 5e-6 * a ** 10 * mu_s
    ginf = si.g_inf(a, mu_s, f_gev)
    e0 = 1.0 / (M ** 2 * sr.N0)

    def rhs(t, y):
        e211, e322, chi = np.exp(y[0]), np.exp(y[1]), y[2]
        g1 = sr.gamma_rel_over_mu(a, chi, 2, 1, 1, sr.NORM_N) * mu_s
        g2 = sr.gamma_rel_over_mu(a, chi, 3, 2, 2, sr.NORM_N) * mu_s
        gbh = si.g_bh(a, chi, mu_s, f_gev)
        d1 = g1 - 2 * gbh * e211 * e322 + ginf * e322 ** 2 - 2 * gGW * e211 + gtr * e322
        d2 = g2 + gbh * e211 ** 2 - 2 * ginf * e322 * e211 - gtr * e211
        dc = -g1 * e211 - 2 * g2 * e322
        return [d1, d2, dc]

    y0 = [np.log(e0), np.log(e0), chi0]
    te = None if t_eval_yr is None else np.asarray(t_eval_yr) * YR
    return solve_ivp(rhs, (0.0, t_end_yr * YR), y0, method="Radau", rtol=1e-8, atol=1e-10,
                     t_eval=te, dense_output=t_eval_yr is None)


def chi_vs_f(f):
    T = [1e6, 1e7, 1e8]
    s = evolve(M2G, CHI0, MU, f, 1e8, t_eval_yr=T)
    if not s.success or s.y.shape[1] < len(T):
        return [np.nan] * len(T)
    return list(s.y[2])


def main():
    os.makedirs(OUT, exist_ok=True)
    res = {}
    si.REL = True
    # ---------- 1. 時間発展(f = 4e12 GeV)と準平衡の比較 ----------
    f = 4e12
    tt = np.logspace(-1, 8, 181)
    s = evolve(M2G, CHI0, MU, f, 1e8, t_eval_yr=tt)
    mu_s = MU / sr.HBAR_EVS
    a = sr.alpha_of(M2G, MU)
    e211_eq = np.array([si.eps_eq(a, c, mu_s, f)[0] for c in s.y[2]])
    e322_eq = np.array([si.eps_eq(a, c, mu_s, f)[1] for c in s.y[2]])
    # 準平衡だけでのスピン(E3 の式を時間積分)
    chi_qe = [CHI0]
    for k in range(1, len(tt)):
        c = chi_qe[-1]
        dt = (tt[k] - tt[k - 1]) * YR
        chi_qe.append(c + si.chidot(a, c, mu_s, f) * dt)
    np.savetxt(os.path.join(FIG, "e05_evol_f4e12.dat"),
               np.column_stack([tt, np.exp(s.y[0]), np.exp(s.y[1]), s.y[2], e211_eq, e322_eq, chi_qe]),
               header="t_yr e211 e322 chi e211_eq e322_eq chi_qe", comments="", fmt="%.6e")
    late = tt > 1e3
    res["evol_f4e12"] = {
        "max_rel_dev_e211_from_eq_t>1e3yr": float(np.max(abs(np.exp(s.y[0][late]) / e211_eq[late] - 1))),
        "max_rel_dev_e322_from_eq_t>1e3yr": float(np.max(abs(np.exp(s.y[1][late]) / e322_eq[late] - 1))),
        "chi_1e7yr_ode": float(np.interp(1e7, tt, s.y[2])),
        "chi_1e7yr_qe": float(np.interp(1e7, tt, chi_qe)),
        "chi_1e8yr_ode": float(s.y[2][-1]),
        "chi_1e8yr_qe": float(chi_qe[-1]),
    }
    # ---------- 2. 年齢 1e6, 1e7, 1e8 yr でのスピン vs f ----------
    # f >~ 1e17 GeV では 211 が約 34 年で自由に飽和し(chi <= 0.40)、そのあと 322 が 1 個の量子から
    # 急に増える段階で陰的解法が止まる。0.69 を保てるかの判定には効かないので、走査は 1e16 GeV までにする。
    fs = np.logspace(11, 16, 41)
    with Pool(10) as p:
        out = np.array(p.map(chi_vs_f, fs))
    free = sr.chi_final_seq(M2G, CHI0, MU, [1e6, 1e7, 1e8], [(2, 1, 1), (3, 2, 2)], rel=True)
    qe = []
    for fv in fs:
        t2 = si.t_spindown(M2G, CHI0, 0.69, MU, fv)
        qe.append(t2)
    np.savetxt(os.path.join(FIG, "e05_chi_vs_f.dat"),
               np.column_stack([fs, out, np.array(qe)]), comments="",
               header="f chi_t6 chi_t7 chi_t8 tqe_to_069_yr", fmt="%.6e")
    # 0.69 を保つ f の上限(ODE、t = 1e7 yr)
    ok = out[:, 1] >= 0.69
    res["chi_vs_f"] = {
        "free_field_chi_t1e6_1e7_1e8": list(map(float, free)),
        "f_max_keep_0.69_t1e7_ode_GeV": float(fs[ok].max()) if ok.any() else None,
        "f_max_keep_0.69_t1e7_qe_GeV": float(1e13 * np.sqrt(si.t_spindown(M2G, CHI0, 0.69, MU, 1e13) / 1e7)),
        "n_failed": int(np.isnan(out[:, 0]).sum()),
    }
    si.REL = False
    with open(os.path.join(OUT, "e05.json"), "w") as fo:
        json.dump(res, fo, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
