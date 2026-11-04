"""論文の図: (mu, f) 面の窓(楽観的な仮定と現実的な仮定)。研究ノート E9 の図を細かい格子で描き直すためのデータ。

仮定と模型は research/calc/e09_review.py と同じ(準平衡、相対論的な超放射の率、gamma_BH は主要項)。
出力: paper/figs/data/fig_window_{opt,real}.dat
実行: cd paper/figs && ~/miniforge3/bin/python3 make_fig_window.py && lualatex fig_window.tex
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "research", "calc"))
from e09_review import bounds  # noqa: E402

MUS = np.geomspace(0.4e-12, 1.2e-12, 61)
HDR = "mu GRS GW241011 upper M10_d18 M10_d05 M15_d18 M15_d05 M20_d18 M20_d05 M25_d18 M25_d05"

os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
for scen in ["opt", "real"]:
    rows = []
    for mu in MUS:
        b = bounds(mu, scen)
        rows.append([b[0], b[1], b[2], min(b[1], b[2])] + b[3:])
    np.savetxt(os.path.join(HERE, "data", f"fig_window_{scen}.dat"), np.array(rows),
               header=HDR, comments="", fmt="%.5e")
    print(scen, "done")
