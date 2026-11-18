"""1 つの節で崩壊定数 f だけを変え、スピン・BH の質量・雲の大きさの時間変化を比べる図(2026-11-18、説明用)。

節は reftab8 の (区分 1, alpha 番号 45, 位置 22, f 番号, M 番号 0): alpha = 0.457、chi0 = 0.995、M = 1 Msun(M は
ln(N0 M^2) を通してだけ効く)。運動論は本番の表と同じ kinn.evolveN_checked(8 準位、相対論的な補正、多倍長の率、
計算の上限 1e7 回と平均の確かめ)。曲線は有効な年齢(tau_valid)まで。
実行(Mac Studio、f の数だけ並列。f = 1e16・3e14 GeV は 1 本約 30 分):
  cd research && ~/miniforge3/bin/python3 calc/fscan_example.py solve     -> calc/out/work/fscan/kfNN.json
  cd research && ~/miniforge3/bin/python3 calc/fscan_example.py plot      -> calc/out/work/fscan_alpha0457.pdf
"""
import json
import os
import sys
import time
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "out", "work", "fscan")
KFS = [46, 38, 32, 26, 20, 14, 8]                      # f = 10^(8 + 0.25 kf) GeV: 3e19 ... 1e10
KEY = (1, 45, 22, None, 0)


def _thin(lt, ys, dl=0.04, dy=(2e-3, 2e-3, 0.05)):
    """図用に間引く: log tau で dl ごと、またはどれかの量が dy 以上動いたら点を残す。"""
    keep, last = [0], 0
    for i in range(1, len(lt)):
        if lt[i] - lt[last] >= dl or any(abs(y[i] - y[last]) >= d for y, d in zip(ys, dy)):
            keep.append(i)
            last = i
    if keep[-1] != len(lt) - 1:
        keep.append(len(lt) - 1)
    return np.array(keep)


def solve_one(kf):
    warnings.simplefilter("ignore")
    import kinn
    import reftab8 as R
    import sr
    key = (KEY[0], KEY[1], KEY[2], kf, KEY[4])
    a, sig, f, M0, cs, chi0 = R.node_params(*key)
    mu = a / sr.alpha_of(M0, 1.0)
    t0 = time.time()
    sol = kinn.evolveN_checked(R._system(), M0, chi0, mu, f, R.TAU_MAX / mu, R.TAU_GRID / mu, max_nfev=10 ** 7)
    nv = sol.n_valid
    N = sol.y.shape[0] - 2
    lt = np.log10(R.TAU_GRID[:nv])
    M = sol.y[N, :nv]
    chi = kinn.chi_of(sol)[:nv]
    cloud = (sr.alpha_of(1.0, mu) * np.exp(sol.y[:N, :nv]) / M).sum(0)     # 雲の質量 / BH の質量
    lc = np.log10(np.maximum(cloud, 1e-30))
    k = _thin(lt, [chi, M / M0, lc])
    d = dict(kf=kf, f=f, alpha=a, chi0=chi0, M0=M0, mu=mu, secs=time.time() - t0, nv=int(nv), reason=sol.reason,
             crit=[float(c) for c in cs], lt=np.round(lt[k], 4).tolist(), chi=np.round(chi[k], 5).tolist(),
             m=np.round(M[k] / M0, 5).tolist(), lc=np.round(lc[k], 3).tolist())
    json.dump(d, open(os.path.join(OUT, f"kf{kf:02d}.json"), "w"))
    return kf, f, d["secs"], nv, sol.reason


def solve(nproc=len(KFS)):
    from multiprocessing import Pool
    os.makedirs(OUT, exist_ok=True)
    with Pool(nproc) as p:
        for r in p.imap_unordered(solve_one, KFS):
            print("kf=%d f=%.2e %.0f s valid %d reason '%s'" % r, flush=True)


def plot(path=os.path.join(HERE, "out", "work", "fscan_alpha0457.pdf")):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {46: "#2a78d6", 38: "#eb6834", 32: "#1baf7a", 26: "#eda100", 20: "#e87ba4", 14: "#008300", 8: "#6250d6"}
    dash = {46: "-", 38: (0, (6, 3)), 32: (0, (2, 2)), 26: (0, (8, 3, 2, 3)), 20: "-", 14: (0, (6, 3)), 8: (0, (2, 2))}
    ds = [json.load(open(os.path.join(OUT, f"kf{kf:02d}.json"))) for kf in KFS
          if os.path.exists(os.path.join(OUT, f"kf{kf:02d}.json"))]
    d0 = ds[0]
    mu30 = d0["alpha"] / __import__("sr").alpha_of(30.0, 1.0)              # alpha = 0.457 の 30 Msun の BH の mu
    age = np.log10(mu30 * 1e10)
    fig, ax = plt.subplots(3, 1, figsize=(6.8, 8.2), sharex=True, gridspec_kw=dict(height_ratios=[1.4, 1, 1.15]))
    for d in ds:
        e = int(np.floor(np.log10(d["f"]) + 1e-9))
        m = d["f"] / 10 ** e
        lab = (r"$f = 10^{%d}$" % e if abs(m - 1) < 0.05 else r"$f = %.0f\times10^{%d}$" % (m, e)) + " GeV"
        if d["nv"] < 6500:
            lab += r" (to $\tau$ = %.1g)" % 10 ** d["lt"][-1]
        kw = dict(color=col[d["kf"]], ls=dash[d["kf"]], lw=1.6)
        ax[0].plot(d["lt"], d["chi"], label=lab, **kw)
        ax[1].plot(d["lt"], d["m"], **kw)
        ax[2].plot(d["lt"], d["lc"], **kw)
    for c, n in zip(d0["crit"], ["211", "322", "433", "544"]):
        ax[0].axhline(c, color="0.6", lw=0.7, ls=(0, (3, 3)))
        ax[0].text(-15.8, c + (0.006 if n != "211" else -0.028), f"{n} critical", fontsize=7.5, color="0.4")
    for a in ax:
        a.axvline(age, color="0.6", lw=0.7, ls=(0, (3, 3)))
    ax[0].text(age - 0.15, 0.365, r"30 $M_\odot$, 10 Gyr", fontsize=7.5, color="0.4", ha="right")
    ax[0].set_ylim(0.35, 1.0)
    ax[1].set_ylim(0.85, 1.005)
    ax[2].set_ylim(-18, 0)
    ax[2].set_xlim(-16, 0.2)
    ax[0].set_ylabel(r"BH spin $\chi$")
    ax[1].set_ylabel(r"$M/M_0$")
    ax[2].set_ylabel(r"$\log_{10}(M_{\rm cloud}/M)$")
    ax[2].set_xlabel(r"$\log_{10}\,\tau$   ($\tau = \mu\,t$  [eV yr])")
    ax[0].legend(fontsize=7.5, loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=4, frameon=False, columnspacing=1.2,
                 handlelength=2.6)
    fig.suptitle(r"$\alpha = %.3f$, $\chi_0 = %.3f$: 8-level kinetics with relativistic corrections"
                 % (d0["alpha"], d0["chi0"]), fontsize=9)
    for a in ax:
        a.tick_params(labelsize=8)
        a.grid(axis="y", color="0.92", lw=0.6)
    fig.align_ylabels(ax)
    fig.tight_layout(h_pad=0.6)
    fig.savefig(path)
    fig.savefig(path.replace(".pdf", ".png"), dpi=110)
    print(path, [d["kf"] for d in ds])


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "plot"
    if cmd == "solve":
        solve(int(sys.argv[2]) if len(sys.argv) > 2 else len(KFS))
    else:
        plot()
