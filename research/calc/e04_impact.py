"""E4: 相対論的な超放射の率に替えたとき、E2・E3 の主な数値がどれだけ動くか。

出力: out/e04_impact.json、../figs/data/e04_R_rel.dat
実行: cd research && ~/miniforge3/bin/python3 calc/e04_rates.py && ~/miniforge3/bin/python3 calc/e04_impact.py
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
LEVELS = [(2, 1, 1), (3, 2, 2), (4, 3, 3)]
res = {}

# ---------- E2 の数値(自由場) ----------
tt = np.logspace(-2, 10, 6001)
e2 = {}
for rel in [False, True]:
    tag = "rel" if rel else "analytic"
    c = sr.chi_final_seq(15.0, 0.69, 7e-13, tt, LEVELS[:1], rel=rel)
    mid = 0.5 * (0.69 + c[-1])
    c241 = sr.chi_final_seq(19.6, 0.78, 7.1e-13, tt, LEVELS, rel=rel)
    # 帯の中のすべての質量で chi_f < 0.5 になる mu(t = 1e6, 1e8 yr)
    Ms = np.logspace(np.log10(3.0), np.log10(150.0), 121)
    mus = np.logspace(-14, -11, 121)
    T = [1e6, 1e8]
    chif = np.array([[sr.chi_final_seq(M, 0.69, mu, T, LEVELS, rel=rel) for M in Ms] for mu in mus])
    bands = {}
    for it, t in enumerate(T):
        for name, (lo, hi) in {"10-20": (10, 20), "45-100": (45, 100)}.items():
            sel = (Ms >= lo) & (Ms <= hi)
            ok = np.all(chif[:, sel, it] < 0.5, axis=1)
            bands[f"{name}_t{t:.0e}"] = [float(mus[ok].min()), float(mus[ok].max())] if ok.any() else None
    e2[tag] = {"t_mid_15Msun_7e-13_yr": float(tt[np.argmax(c < mid)]),
               "GW241011_t_below_0.69_yr": float(tt[np.argmax(c241 < 0.69)]) if c241[-1] < 0.69 else None,
               "bands_mu_range": bands}
res["E2"] = e2

# ---------- E3 の数値(準平衡、gamma_SR だけ相対論的に替える) ----------
mu = 0.71e-12
e3 = {}
rows = {}
M1s = np.linspace(6.0, 28.0, 45)
for rel in [False, True]:
    si.REL = rel
    tag = "rel" if rel else "analytic"
    t2 = si.t_spindown(19.6, 0.78, 0.69, mu, 1e13)
    d = {"f_max_T2=1e7yr_GeV": 1e13 * np.sqrt(t2 / 1e7)}
    r = []
    for M1 in M1s:
        thr = si.chi_threshold_211(M1, mu)
        r.append([si.t_spindown(M1, c0, thr + dl, mu, 1e13) / t2 if c0 > thr + dl else np.nan
                  for c0, dl in [(0.7, 0.05), (0.7, 0.18)]])
    rows[tag] = np.array(r)
    for M1 in [8.0, 10.0, 15.0, 20.0]:
        i = np.argmin(abs(M1s - M1))
        d[f"R(chi0=0.7,d=0.05)_M1={M1:g}"] = rows[tag][i, 0]
        d[f"R(chi0=0.7,d=0.18)_M1={M1:g}"] = rows[tag][i, 1]
    e3[tag] = d
si.REL = False
res["E3"] = e3
np.savetxt(os.path.join(FIG, "e04_R_rel.dat"),
           np.column_stack([M1s, rows["analytic"], rows["rel"]]), comments="",
           header="M1 Ran_d05 Ran_d18 Rrel_d05 Rrel_d18", fmt="%.6g")

with open(os.path.join(OUT, "e04_impact.json"), "w") as f:
    json.dump(res, f, indent=1)
print(json.dumps(res, indent=1))
