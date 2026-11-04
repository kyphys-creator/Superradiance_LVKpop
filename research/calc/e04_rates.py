"""E4: 相対論的な超放射の率を Leaver の連分数法で表にする(211・322・433)。

照合:
  1. Dolan, arXiv:0705.2880 の Table I(211 の最大成長率)
  2. 多倍長(mpmath、40 桁)との比較
  3. 小さい alpha で 1908.10370 の式 2.28(実部)・2.32(虚部)に近づくこと
出力:
  out/rates_{211,322,433}.npz   alpha、chi の格子と M omega(実部・虚部)
  out/e04.json                  照合の数値
  ../figs/data/e04_ratio_*.dat  数値/解析式の比(図用)
実行: cd research && ~/miniforge3/bin/python3 calc/e04_rates.py
"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import leaver  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")
NCF = 4000

LEVELS = {"211": (2, 1, 1, 0.03, 0.55), "322": (3, 2, 2, 0.08, 1.0), "433": (4, 3, 3, 0.15, 1.5)}
CHIS = np.concatenate([np.linspace(0.0, 0.9, 91), np.linspace(0.905, 0.995, 19)])


def solve_point(args):
    alpha, chi, n, l, m = args
    try:
        w = leaver.bound_state(alpha, chi, n, l, m, N=NCF)
    except Exception:
        return np.nan, np.nan
    # 別の節の根に飛んでいないか(実部が水素原子の値から 5% 以上ずれたら捨てる)
    if abs(w.real / alpha - (1 - alpha ** 2 / (2 * n ** 2))) > 0.05:
        return np.nan, np.nan
    return w.real, w.imag


def main():
    os.makedirs(OUT, exist_ok=True)
    res = {}
    # ---------- 1. Dolan Table I ----------
    dolan = [(0.7, 0.187, 3.33e-10), (0.8, 0.231, 2.16e-9), (0.9, 0.293, 1.55e-8),
             (0.95, 0.343, 4.88e-8), (0.98, 0.393, 1.11e-7), (0.99, 0.421, 1.50e-7)]
    res["dolan_tableI"] = [{"a": a, "Mmu": al, "Im_ref": ref,
                            "Im_this": leaver.bound_state(al, a, 2, 1, 1, N=NCF).imag}
                           for a, al, ref in dolan]
    # ---------- 2. 多倍長との比較 ----------
    mpc = []
    for (n, l, m, al, chi) in [(2, 1, 1, 0.05, 0.9), (2, 1, 1, 0.293, 0.9), (3, 2, 2, 0.15, 0.9), (4, 3, 3, 0.4, 0.95)]:
        wd = leaver.bound_state(al, chi, n, l, m, N=NCF)
        wm = leaver.bound_state_mp(al, chi, n, l, m, N=NCF, dps=40)
        mpc.append({"level": f"{n}{l}{m}", "alpha": al, "chi": chi,
                    "Im_double": wd.imag, "Im_mp": float(wm.imag), "rel_diff": abs(wd.imag / float(wm.imag) - 1)})
    res["double_vs_mp"] = mpc

    # ---------- 3. 格子 ----------
    for name, (n, l, m, amin, amax) in LEVELS.items():
        alphas = np.geomspace(amin, amax, 64)
        tasks = [(al, chi, n, l, m) for al in alphas for chi in CHIS]
        with Pool(10) as p:
            out = p.map(solve_point, tasks, chunksize=20)
        out = np.array(out).reshape(len(alphas), len(CHIS), 2)
        np.savez(os.path.join(OUT, f"rates_{name}.npz"), alpha=alphas, chi=CHIS,
                 re=out[..., 0], im=out[..., 1], n=n, l=l, m=m)
        res[f"grid_{name}"] = {"n_alpha": len(alphas), "n_chi": len(CHIS),
                               "n_nan": int(np.isnan(out[..., 1]).sum())}
        # 数値/解析式の比(chi = 0.5, 0.7, 0.9, 0.99)
        rows = []
        for i, al in enumerate(alphas):
            row = [al]
            for chi in [0.5, 0.7, 0.9, 0.99]:
                j = np.argmin(abs(CHIS - chi))
                an = sr.gamma_over_mu(al, CHIS[j], n, l, m) * al
                row.append(out[i, j, 1] / an if an > 0 and out[i, j, 1] > 0 else np.nan)
            rows.append(row)
        np.savetxt(os.path.join(FIG, f"e04_ratio_{name}.dat"), np.array(rows), comments="",
                   header="alpha r05 r07 r09 r099", fmt="%.6g")
    with open(os.path.join(OUT, "e04.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
