"""E10: 方針の見直し。GW のスピンだけで閉じた窓、2 世代目の生まれたときのスピンの分布、集団解析用の chi_SR の表。

1. GW241011 の主星だけで f の上限を決め(X 線連星を外す)、回転を失える 1 世代目の最小質量 M_min(mu) を出す。
   仮定の組は E9 と同じ(楽観的: 0.78 で生まれ T2 = 3e7 yr、現実的: 0.72 で生まれ T2 = 5e7 yr)。
   1 世代目: chi0 = 0.7 から (211 の臨界 + delta)、T1 = 6e9 yr。比べるため、現実的な GRS を入れた場合も出す。
2. 2 世代目の生まれたときのスピン: hbr16.py(1605.01938 の一般スピン版)で、親の分布を変えて計算する。
3. chisr.py の表を、直接の積分(臨界からの距離の対数で変数変換)と E5 の ODE と照合し、計算時間を測る。
出力: out/e10.json、out/chisr_table.npz、../figs/data/e10_*.dat
実行: cd research && ~/miniforge3/bin/python3 calc/e10_gwonly.py
"""
import json
import os
import sys
import time

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import hbr16  # noqa: E402
import si  # noqa: E402
import sr  # noqa: E402
from e08_window import f_limit  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")
SCEN = {"opt": (0.78, 3e7), "real": (0.72, 5e7)}


def f_1g(M1, mu, dl):
    thr = si.chi_threshold_211(M1, mu)
    return f_limit(M1, 0.7, thr + dl, 6e9, mu) if 0.7 > thr + dl else np.inf


def m_cap(mu, dl):
    """211 の臨界 + delta が初期スピン 0.7 に等しくなる質量(これより重い 1 世代目は判定の対象外)。"""
    return brentq(lambda M: si.chi_threshold_211(M, mu) + dl - 0.7, 1.0, 200.0, xtol=1e-3)


def m_min(f_up, mu, dl, lo=5.0):
    """f_1G(M) <= f_up となる最小の M。f_1G は M について減少し、m_cap で 0 に近づく。"""
    hi = m_cap(mu, dl) - 0.01
    g = lambda M: np.log(f_1g(M, mu, dl)) - np.log(f_up)   # noqa: E731
    if g(hi) > 0:
        return np.nan
    if g(lo) < 0:
        return lo
    return brentq(g, lo, hi, xtol=0.05)


def main():
    res = {}
    # ---------- 1. GW だけの窓と最小質量 ----------
    mus = np.geomspace(0.45e-12, 1.1e-12, 14)
    rows = []
    for mu in mus:
        row = [mu]
        for scen, (c0, T2) in SCEN.items():
            fgw = f_limit(19.6, c0, 0.69, T2, mu)
            for dl in [0.18, 0.05]:
                row.append(m_min(fgw, mu, dl))
        fup = min(f_limit(19.6, 0.72, 0.69, 5e7, mu), f_limit(12.4, 0.995, 0.87, 3e9, mu))   # 現実的 + GRS
        row.append(m_min(fup, mu, 0.18))
        rows.append(row)
    rows = np.array(rows)
    np.savetxt(os.path.join(FIG, "e10_mmin.dat"), rows, comments="",
               header="mu opt_d18 opt_d05 real_d18 real_d05 realGRS_d18", fmt="%.5e")
    res["m_min"] = {}
    for mu in [0.60e-12, 0.71e-12, 0.89e-12]:
        d = {}
        for scen, (c0, T2) in SCEN.items():
            fgw = f_limit(19.6, c0, 0.69, T2, mu)
            d[f"{scen}_f_up_GeV"] = fgw
            for dl in [0.18, 0.05]:
                d[f"{scen}_d{dl}"] = m_min(fgw, mu, dl)
        fup = min(f_limit(19.6, 0.72, 0.69, 5e7, mu), f_limit(12.4, 0.995, 0.87, 3e9, mu))
        d["real_with_GRS_d0.18"] = m_min(fup, mu, 0.18)
        d["m_cap_d0.18"] = m_cap(mu, 0.18)
        d["m_cap_d0.05"] = m_cap(mu, 0.05)
        res["m_min"][f"mu={mu:.2e}"] = d

    # ---------- 2. 2 世代目の生まれたときのスピン ----------
    rng = np.random.default_rng(20261107)
    n = 400000
    pops = {
        "nonspin_q0.5-1": dict(a="zero", qlo=0.5),
        "U(0,0.5)_iso_q0.5-1": dict(a="uni", qlo=0.5),
        "Beta(2,5)_iso_q0.5-1": dict(a="beta", qlo=0.5),
        "Beta(2,5)_iso_q0.25-1": dict(a="beta", qlo=0.25),
    }
    res["natal_2G"] = {}
    hist_edges = np.linspace(0.3, 0.95, 66)
    hists = [0.5 * (hist_edges[1:] + hist_edges[:-1])]
    for name, p in pops.items():
        q = rng.uniform(p["qlo"], 1.0, n)
        if p["a"] == "zero":
            a1 = a2 = np.zeros(n)
        elif p["a"] == "uni":
            a1, a2 = rng.uniform(0, 0.5, n), rng.uniform(0, 0.5, n)
        else:
            a1, a2 = rng.beta(2, 5, n), rng.beta(2, 5, n)
        cb, cg = rng.uniform(-1, 1, n), rng.uniform(-1, 1, n)
        dphi = rng.uniform(0, 2 * np.pi, n)
        cal = cb * cg + np.sqrt(1 - cb ** 2) * np.sqrt(1 - cg ** 2) * np.cos(dphi)
        af = hbr16.a_final(q, a1, a2, cb, cg, cal)
        res["natal_2G"][name] = {
            "median": float(np.median(af)), "p05": float(np.percentile(af, 5)), "p95": float(np.percentile(af, 95)),
            "P(>=0.72)": float(np.mean(af >= 0.72)), "P(>=0.78)": float(np.mean(af >= 0.78)),
        }
        h, _ = np.histogram(af, bins=hist_edges, density=True)
        hists.append(h)
    np.savetxt(os.path.join(FIG, "e10_natal.dat"), np.column_stack(hists), comments="",
               header="chi " + " ".join(f"p{i}" for i in range(len(pops))), fmt="%.5e")
    res["natal_2G_checks"] = {"q1_nonspin": float(hbr16.a_final(1, 0, 0, 1, 1, 1)),
                              "q1/3_nonspin": float(hbr16.a_final(1 / 3, 0, 0, 1, 1, 1)),
                              "q1_aligned_0.5": float(hbr16.a_final(1, 0.5, 0.5, 1, 1, 1))}

    # ---------- 3. chi_SR の表の照合と計算時間 ----------
    t0 = time.time()
    tab = chisr.build()
    res["chisr_build_s"] = time.time() - t0
    np.savez(os.path.join(OUT, "chisr_table.npz"), alphas=chisr.ALPHAS, U=chisr.U, **tab)
    si.REL = True
    mu = 0.71e-12
    mu_s = mu / sr.HBAR_EVS

    def chic_rel(a):
        return brentq(lambda c: sr.gamma_rel_over_mu(a, c, 2, 1, 1), 0.01, 0.99, xtol=1e-13)

    def t_exact(a, cr, c0, c1, f):
        s0, s1 = np.log(c1 - cr), np.log(c0 - cr)
        return quad(lambda s: np.exp(s) / abs(si.chidot(a, cr + np.exp(s), mu_s, f)), s0, s1,
                    limit=400, epsrel=1e-10)[0] / sr.YR

    rng = np.random.default_rng(1)
    d = []
    for _ in range(300):
        M = rng.uniform(8, 30)
        f = 10 ** rng.uniform(11.5, 13)
        T = 10 ** rng.uniform(6, 9.8)
        c0 = rng.uniform(0.3, 0.95)
        a = sr.alpha_of(M, mu)
        cr = chic_rel(a)
        if c0 <= cr + 1e-3:
            continue
        cT = float(chisr.chi_after(tab, a, c0, chisr.tau_of(mu, f, T)))
        lo = cr + 1e-9
        g = lambda x: t_exact(a, cr, c0, x, f) - T   # noqa: E731
        ce = lo if g(lo) < 0 else brentq(g, lo, c0 - 1e-12, xtol=1e-12)
        d.append(abs(cT - ce))
    d = np.array(d)
    a241 = sr.alpha_of(19.6, mu)
    nev = 1_000_000
    A = rng.uniform(0.03, 0.25, nev)
    C = rng.uniform(0.2, 0.95, nev)
    TAU = 10 ** rng.uniform(-6, 0, nev) * chisr.tau_of(mu, 1e12, 1e9)
    t0 = time.time()
    chisr.chi_after(tab, A, C, TAU)
    res["chisr_check"] = {
        "n_random": int(len(d)), "abs_err_median": float(np.median(d)),
        "abs_err_p95": float(np.percentile(d, 95)), "abs_err_max": float(d.max()),
        "GW241011_f4e12_T1e7": float(chisr.chi_after(tab, a241, 0.78, chisr.tau_of(mu, 4e12, 1e7))),
        "GW241011_f4e12_T1e8": float(chisr.chi_after(tab, a241, 0.78, chisr.tau_of(mu, 4e12, 1e8))),
        "E5_ODE_T1e7": 0.6951128, "E5_ODE_T1e8": 0.4488703,
        "time_1e6_evals_s": time.time() - t0,
    }
    si.REL = False
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "e10.json"), "w") as fo:
        json.dump(res, fo, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
