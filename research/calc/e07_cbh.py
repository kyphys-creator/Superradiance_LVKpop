"""E7: gamma_BH(211x211 -> 322xBH)の相対論的な補正が、窓の数値をどれだけ動かすか(感度の検定)。

補正は 2412.03655 の Fig. 3 右から目で読んだ値(si.CBH_*、a~ = 0.95)。超放射の率は E4 の相対論的な率。
出力: out/e07.json、../figs/data/e07_R_mu*.dat
実行: cd research && ~/miniforge3/bin/python3 calc/e07_cbh.py
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
M2G, C2G, C2G_LOW = 19.6, 0.78, 0.69
M1s = np.linspace(6.0, 28.0, 45)
res = {}
si.REL = True
for tag, mu in {"0.60": 0.60e-12, "0.71": 0.71e-12, "0.89": 0.89e-12}.items():
    d = {}
    cols = []
    for on in [False, True]:
        si.CBH_ON = on
        key = "cbh" if on else "nocbh"
        t2 = si.t_spindown(M2G, C2G, C2G_LOW, mu, 1e13)
        d[f"f_max_T2=1e7_{key}_GeV"] = 1e13 * np.sqrt(t2 / 1e7)
        d[f"f_max_T2=1e8_{key}_GeV"] = 1e13 * np.sqrt(t2 / 1e8)
        R = []
        for M1 in M1s:
            thr = si.chi_threshold_211(M1, mu)
            R.append([si.t_spindown(M1, c0, thr + dl, mu, 1e13) / t2 if c0 > thr + dl else np.nan
                      for c0, dl in [(0.7, 0.05), (0.7, 0.18)]])
        R = np.array(R)
        cols.append(R)
        for M1 in [8.0, 10.0, 15.0, 20.0]:
            i = np.argmin(abs(M1s - M1))
            d[f"R(0.7,0.05)_M1={M1:g}_{key}"] = R[i, 0]
            d[f"R(0.7,0.18)_M1={M1:g}_{key}"] = R[i, 1]
    np.savetxt(os.path.join(FIG, f"e07_R_mu{tag}.dat"), np.column_stack([M1s, cols[0], cols[1]]),
               comments="", header="M1 R_d05 R_d18 Rc_d05 Rc_d18", fmt="%.6g")
    res[f"mu={tag}e-12"] = d
si.CBH_ON = False
si.REL = False
with open(os.path.join(OUT, "e07.json"), "w") as f:
    json.dump(res, f, indent=1)
print(json.dumps(res, indent=1))
