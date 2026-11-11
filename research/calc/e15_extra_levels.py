"""E15: 三準位の外の準位(411、433、211 の臨界より上の 544)が育つかの判定と、印 4 の付け直し。

1. 固定点での 411・433・544 の増加率の符号(主要項の率と、BH が吸収する過程すべてに E7 の C_BH を掛けた率)。
2. 411 の相対論的な率: Leaver の連分数法で 411 を解き、相対論的/解析式の比が 211 と同じ程度であることを確かめる。
3. 五準位(211・322・411・433・544)の ODE: alpha = 0.22(BGLS の Fig. 411 の例、10 Msun、chi0 = 0.9)で、
   主要項と C_BH 込みの両方を三準位と比べる。
4. 新しい印の割合(1 世代目の自然スピン Beta(0.69, 2.28))。
出力: out/e15.json
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
MU = 0.71e-12


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
                      "log10eps_peak": [float(v) / np.log(10) for v in s.y[:5].max(axis=1)]}))


def run_ode5(mode, M, c0, mu, f, Ts, timeout=600):
    try:
        r = subprocess.run([sys.executable, os.path.abspath(__file__), "--ode5", mode, str(M), str(c0), str(mu), str(f),
                            ",".join(str(t) for t in Ts)], capture_output=True, text=True, timeout=timeout)
        return json.loads(r.stdout.strip().splitlines()[-1])
    except (subprocess.TimeoutExpired, IndexError, json.JSONDecodeError):
        return {"status": -9, "chi": []}


def main():
    res = {"growth_map_leading": growth_map(False), "growth_map_cbh": growth_map(True), "leaver_411": leaver_411()}
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
