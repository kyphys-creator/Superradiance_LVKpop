"""E15: 三準位の外の準位(411、433、211 の臨界より上の 544)が育つかの判定と、印 4 の付け直し。

1. 固定点での 411・433・544 の増加率の符号(主要項の率と、BH が吸収する過程すべてに E7 の C_BH を掛けた率)。
2. 411 の相対論的な率: Leaver の連分数法で 411 の表を作る(out/rates_411.npz)。相対論的/解析式の比も 211 と比べる。
3. 五準位(211・322・411・433・544)の ODE: alpha = 0.22(BGLS の Fig. 411 の例、10 Msun、chi0 = 0.9)で、
   主要項と C_BH 込みの両方を三準位と比べる。
4. 新しい印の割合(1 世代目の自然スピン Beta(0.69, 2.28))。
5. 図のデータ: 育つ場所の (alpha, chi) の地図と、五準位の ODE の時間変化。
出力: out/e15.json、../figs/data/e15_*.dat
実行: cd research && ~/miniforge3/bin/python3 calc/e15_extra_levels.py
"""
import json
import os
import subprocess
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import fp3  # noqa: E402
import leaver  # noqa: E402
import si  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")
MU = 0.71e-12


def _solve_411(args):
    """411 の準束縛状態を 1 点解く。隣の 311・511 に飛んでいないかを、実部が n = 4 の水素原子の値にいちばん近いかで確かめる。"""
    alpha, chi = args
    from e04_rates import ncf
    try:
        w = leaver.bound_state(alpha, chi, 4, 1, 1, N=ncf(alpha))
    except Exception:
        return np.nan, np.nan
    e = {n: 1 - alpha ** 2 / (2 * n ** 2) for n in (3, 4, 5)}
    x = w.real / alpha
    if min(e, key=lambda n: abs(x - e[n])) != 4:
        return np.nan, np.nan
    return w.real, w.imag


def make_rates_411():
    """411 の相対論的な率の表(E4 と同じ形式、out/rates_411.npz)。alpha 0.05--0.45、chi は E4 と同じ格子。"""
    from multiprocessing import Pool
    from e04_rates import CHIS
    alphas = np.geomspace(0.05, 0.45, 64)
    with Pool(10) as p:
        out = p.map(_solve_411, [(al, chi) for al in alphas for chi in CHIS], chunksize=40)
    out = np.array(out).reshape(len(alphas), len(CHIS), 2)
    np.savez(os.path.join(OUT, "rates_411.npz"), alpha=alphas, chi=CHIS, re=out[..., 0], im=out[..., 1], n=4, l=1, m=1)
    return {"n_alpha": len(alphas), "n_chi": len(CHIS), "n_nan": int(np.isnan(out[..., 1]).sum())}


def growth_map(use_c):
    rows = []
    for a in [0.12, 0.15, 0.17, 0.20, 0.22, 0.25, 0.28, 0.30]:
        C = float(np.interp(a, si.CBH_ALPHA, si.CBH_VALUE)) if use_c else 1.0
        c1, c2 = chisr._zero(a, 2, 1, 1), chisr._zero(a, 3, 2, 2)
        pts = {}
        for lab, chi in [("c211+0.2", min(c1 + 0.2, 0.99)), ("c211+0.05", c1 + 0.05), ("c211+0.005", c1 + 0.005),
                         ("c211-0.01", c1 - 0.01), ("mid", 0.5 * (c1 + c2)), ("c322+0.01", c2 + 0.01)]:
            g, k = fp3.coeffs(a, chi, C, cbh_all=True)
            p3, p2 = fp3.fp_P3(g, k), fp3.fp_P2(g, k)
            if g[0] >= 0 or fp3.percap(p3, g, k)[2] < 0:
                x, tag = p3, "P3"
            else:
                x, tag = p2, "P2"
            l411, l433 = fp3.extra_growth(a, chi, x, C)
            l544 = fp3.percap(x, g, k)[2] if tag == "P3" else 0.0
            pts[lab] = {"fp": tag, "grows": [n for n, v in (("411", l411), ("433", l433), ("544", l544)) if v > 0]}
        rows.append({"alpha": a, "C_BH": C, "c211": round(c1, 4), "c322": round(c2, 4), "points": pts})
    return rows


def leaver_411():
    out = []
    for a, chi in [(0.1, 0.9), (0.2, 0.9), (0.25, 0.95), (0.3, 0.95)]:
        r = {}
        for n in (2, 4):
            w = leaver.bound_state(a, chi, n, 1, 1)
            r[f"{n}11"] = float(w.imag / a / sr.gamma_over_mu(a, chi, n, 1, 1))
        out.append({"alpha": a, "chi": chi, "rel_over_analytic_211": r["211"], "rel_over_analytic_411": r["411"]})
    return out


def ode5(mode, M, c0, mu, f, Ts):
    """mode: 5c / 3c / 5n(五準位か三準位か、C_BH を掛けるか)。別のプロセスで呼ぶ。"""
    a = sr.alpha_of(M, mu)
    C = float(np.interp(a, si.CBH_ALPHA, si.CBH_VALUE)) if mode.endswith("c") else 1.0
    if mode.startswith("3"):
        fp3.PROC5 = [p for p in fp3.PROC5 if 2 not in p[:3] and 3 not in p[:3]]
        orig = fp3.sr_rates5
        fp3.sr_rates5 = lambda al, ch: orig(al, ch) * np.array([1, 1, 0, 0, 1])
    s = fp3.evolve5(M, c0, mu, f, Ts[-1], Ts, cbh=C)
    print(json.dumps({"status": int(s.status), "chi": [float(v) for v in s.y[5]],
                      "log10eps_peak": [float(v) / np.log(10) for v in s.y[:5].max(axis=1)],
                      "log10eps": [[float(v) / np.log(10) for v in row] for row in s.y[:5]]}))


def run_ode5(mode, M, c0, mu, f, Ts, timeout=600):
    try:
        r = subprocess.run([sys.executable, os.path.abspath(__file__), "--ode5", mode, str(M), str(c0), str(mu), str(f),
                            ",".join(str(t) for t in Ts)], capture_output=True, text=True, timeout=timeout)
        return json.loads(r.stdout.strip().splitlines()[-1])
    except (subprocess.TimeoutExpired, IndexError, json.JSONDecodeError):
        return {"status": -9, "chi": []}


T_MAP_YR = 1e10   # 地図の基準: 固定点にとどまったとき、1 個の量子から 10 Gyr のうちに育つか(mu = 0.71e-12 eV)


def first_growth(a, chis, use_c, which):
    """chis(昇順)のうち、固定点で which の準位が 10 Gyr のうちに 1 個の量子から育つ最小の chi(なければ nan)。
    増加率の符号だけでは、211 の臨界のすぐ上で固定点の占有数が 0 に近づき、433 自身のごく遅い超放射が残る細い帯が出るので、
    lambda mu T_MAP > ln(GM^2) を条件にする(印 4 の判定と同じ考え方)。"""
    C = float(np.interp(a, si.CBH_ALPHA, si.CBH_VALUE)) if use_c else 1.0
    M = a / sr.alpha_of(1.0, MU)
    lam_min = np.log(sr.N0 * M ** 2) / (MU / sr.HBAR_EVS * T_MAP_YR * sr.YR)
    for chi in chis:
        g, k = fp3.coeffs(a, chi, C, cbh_all=True)
        p3, p2 = fp3.fp_P3(g, k), fp3.fp_P2(g, k)
        if g[0] >= 0 or fp3.percap(p3, g, k)[2] < 0:
            x, tag = p3, "P3"
        else:
            x, tag = p2, "P2"
        l411, l433 = fp3.extra_growth(a, chi, x, C)
        l544 = fp3.percap(x, g, k)[2] if tag == "P3" else -1.0
        lam = {"411": l411, "433": l433, "544": l544}
        if max(lam[w] for w in which) > lam_min:
            return float(chi)
    return np.nan


def polygon(alphas, lo, top=0.995):
    """lo(alpha) から top までの領域の多角形(lo が nan の alpha は除く)。"""
    m = np.isfinite(lo)
    if not m.any():
        return np.zeros((0, 2))
    xa, ya = alphas[m], lo[m]
    return np.vstack([np.column_stack([xa, ya]), [[xa[-1], top], [xa[0], top]]])


def figure_data(res, M=10.0, c0=0.9, f=1e14):
    """図のデータ: (a) 育つ場所の地図、(b)(c) 五準位の ODE の時間変化。"""
    alphas = np.linspace(0.10, 0.30, 81)
    c1s = np.array([chisr._zero(a, 2, 1, 1) for a in alphas])
    c2s = np.array([chisr._zero(a, 3, 2, 2) for a in alphas])
    lo_any, lo411, lo_lead, below_any = [], [], [], []
    for a, c1, c2 in zip(alphas, c1s, c2s):
        above = c1 + np.concatenate([np.geomspace(1e-5, 1e-2, 12, endpoint=False), np.linspace(1e-2, 0.995 - c1, 150)])
        lo_any.append(first_growth(a, above, True, ("433", "544")))
        lo411.append(first_growth(a, above, True, ("411",)))
        lo_lead.append(first_growth(a, above, False, ("411", "433", "544")))
        below = np.linspace(c2 + 1e-3, c1 - 2e-3, 40)
        below_any.append(np.isfinite(first_growth(a, below, True, ("411", "433", "544"))))
    lo_any, lo411, lo_lead = map(np.array, (lo_any, lo411, lo_lead))
    np.savetxt(os.path.join(FIG, "e15_map.dat"), np.column_stack([alphas, c1s, c2s, lo_any, lo411, lo_lead]),
               header="alpha c211 c322 lo_cbh lo411_cbh lo_lead", comments="", fmt="%.6e")
    for name, lo in (("e15_reg_cbh.dat", lo_any), ("e15_reg411.dat", lo411), ("e15_reg_lead.dat", lo_lead)):
        np.savetxt(os.path.join(FIG, name), polygon(alphas, lo), header="x y", comments="", fmt="%.6e")
    res["map_summary"] = {"first_alpha_433_544_cbh": float(alphas[np.isfinite(lo_any)].min()),
                          "first_alpha_411_cbh": float(alphas[np.isfinite(lo411)].min()),
                          "first_alpha_leading": float(alphas[np.isfinite(lo_lead)].min()),
                          "any_growth_below_c211_cbh": bool(np.any(below_any))}
    # (b)(c) 時間変化
    mu = 0.22 / sr.alpha_of(M, 1.0)
    Ts = list(np.logspace(-2, 3, 61))
    o5c = run_ode5("5c", M, c0, mu, f, Ts)
    o3c = run_ode5("3c", M, c0, mu, f, Ts)
    o5n = run_ode5("5n", M, c0, mu, f, Ts)
    tab = chisr.build()
    csr, _ = chisr.chi_sr(tab, M, c0, np.array(Ts), mu, f)
    np.savetxt(os.path.join(FIG, "e15_ode5_chi.dat"),
               np.column_stack([Ts, o5c["chi"], o3c["chi"], o5n["chi"], csr]),
               header="t chi5c chi3c chi5n chisr", comments="", fmt="%.6e")
    np.savetxt(os.path.join(FIG, "e15_ode5_eps.dat"), np.column_stack([Ts] + [o5c["log10eps"][i] for i in range(5)]),
               header="t e211 e322 e411 e433 e544", comments="", fmt="%.6e")
    res["ode5_c211"] = chisr._zero(0.22, 2, 1, 1)


def main():
    res = {"rates_411_table": make_rates_411()}
    res.update({"growth_map_leading": growth_map(False), "growth_map_cbh": growth_map(True), "leaver_411": leaver_411()})
    tab = chisr.build()

    # ---------- 3. 五準位の ODE ----------
    M, c0, f = 10.0, 0.9, 1e14
    mu = 0.22 / sr.alpha_of(M, 1.0)
    Ts = [1, 10, 100, 300, 1000]
    c, r = chisr.chi_sr(tab, M, c0, np.array(Ts, float), mu, f)
    res["ode5_alpha0.22"] = {"M": M, "chi0": c0, "f": f, "mu_eV": mu, "T_yr": Ts,
                             "five_levels_cbh": run_ode5("5c", M, c0, mu, f, Ts),
                             "three_levels_cbh": run_ode5("3c", M, c0, mu, f, Ts),
                             "five_levels_leading": run_ode5("5n", M, c0, mu, f, Ts),
                             "chi_sr": [float(v) for v in c], "chi_sr_regime": [int(v) for v in r]}

    # ---------- 4. 印の割合 ----------
    rng = np.random.default_rng(12)
    n = 200000
    Mpop = rng.uniform(8, 50, n)
    c0pop = np.minimum(rng.beta(0.69, 2.28, n), 0.999)
    res["natal_beta_regimes_T6Gyr"] = {}
    for fv in [1e12, 1e14, 1e16]:
        _, rr = chisr.chi_sr(tab, Mpop, c0pop, 6e9, MU, fv)
        res["natal_beta_regimes_T6Gyr"][f"{fv:.0e}"] = {str(k): round(float(np.mean(rr == k)), 3) for k in range(5)}
    res["examples_T6Gyr_f1e12"] = {}
    for Mx, cx in [(19.6, 0.78), (25, 0.8), (30, 0.8), (35, 0.8), (40, 0.8), (45, 0.7), (50, 0.6)]:
        cc, rg = chisr.chi_sr(tab, Mx, cx, 6e9, MU, 1e12)
        res["examples_T6Gyr_f1e12"][f"M{Mx}_chi{cx}"] = {"alpha": round(sr.alpha_of(Mx, MU), 3), "chi": float(cc), "regime": int(rg)}

    figure_data(res)

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "e15.json"), "w") as fo:
        json.dump(res, fo, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--ode5":
        ode5(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]),
             [float(x) for x in sys.argv[7].split(",")])
    else:
        main()
