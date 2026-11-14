"""E19: 独立な査読(2026-11-13)の所見の確かめと、修正の前後の比較(照合の本体は calc/e19_validate.py)。

修正の前の単体版は git の OLD_COMMIT(E18)の share/chi_sr/chi_sr.py。表(chi_sr_tables.npz)は前後で同じ。
  ode      : 査読のときに E18 の式だけから別に書いた三準位の運動論を、fp3.evolve3m と比べる(12 例)。
  rates    : 超放射の率が連分数法の表と解析式で切り替わる alpha(322 は 0.1、211 は 0.03)での率の飛びと、
             322・211 の解析式を多倍長の連分数法と比べた比。図のデータ figs/data/e19_jump.dat も書く。
  stepscan : 322 から育つ BH で mu を変えて alpha を動かし、段差の時刻(表の補間 / 基準の運動論)と、印 0・1・2 の年齢での
             誤差の最大を出す(図のデータ figs/data/e19_step.dat)。
  kink     : 段差の時刻の、alpha の方向の折れ(reftab.compute_ks の d1、d2)を値の表から測る。
  noise322 : 322 の率の表(倍精度の連分数法、E4)を、表の節の上で多倍長と比べる(alpha 0.10--0.165、低いスピン。並列 6)。
  e18      : E18 の照合の入力(seed 1801、1802)を修正の後の版で引き直し、変わった印と、約束の誤差を数える。
  window   : 窓の数値(E12 の表: GW241011 の主星からの f の上限と、回転を失う 1 世代目の最小質量)を表の版で計算し直す。
  coverage : 約束の外(印 3・4・5)の標本の割合を、mu・質量・スピンの分布ごとに数える。
  mass     : 関数の M は生まれたときの質量。観測された(最後の)質量を渡したときのずれ。
  grown322 : 322 から育つ BH(f 1e15--10^20.5 GeV)を alpha のセルごとに 260 点ずつとり、細かい年齢で修正の前後を照合する
             (所見 1 の折れが効くのはこの種類の BH。並列 6)。
  plateau  : 段のあとの平らな所の値が、進みのいちばん遅い角の節で決まる例(seed 1901 の thresh の層の 1 点)。
  m4stage  : 自由場で 544(m = 4)の超放射まで進む点の例(独立な照合 seed 1902 で、折れの印のあとも約束を越えた 1 点)。
             433 の印が返すスピンしか見ないために印 2 に戻っていたことと、修正のあとの印。
  lowalpha : 値の表の下(alpha < 0.008)の点の扱い(雲が育つ時刻の LOW_MARGIN e-fold 前より若ければ chi0、そうでなければ
             印 3)を、雲が年齢のうちに育ちうる軽い BH(1--3.5 Msun)で、基準の運動論と細かい年齢で比べる。
  level422 : (211, 322) の固定点での 422 の増加率(主要項)と、(322, 544) の固定点での 422 の超放射の率。
  usability: pickle、メモリ、速さ、警告、NaN を整数に直す所(x86 の挙動を模して確かめる)。
出力: out/e19.json(部分ごとに書き足す)
実行: cd research && ~/miniforge3/bin/python3 calc/e19_review.py [部分の名前 ...]   (名前なしで全部。noise322・mass・grown322・lowalpha は並列 6)
"""
import importlib.util
import json
import os
import subprocess
import sys
import time
import warnings
from multiprocessing import Pool

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import fp3  # noqa: E402
import leaver  # noqa: E402
import reftab  # noqa: E402
import si  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")
PKG = os.path.join(HERE, "..", "..", "share", "chi_sr")
TABLES = os.path.join(PKG, "chi_sr_tables.npz")
OLD_COMMIT = "7846eae"
MU = 0.71e-12


def old_source():
    return subprocess.run(["git", "show", f"{OLD_COMMIT}:share/chi_sr/chi_sr.py"], cwd=HERE, capture_output=True,
                          text=True, check=True).stdout


def module_from(src, name):
    path = os.path.join(OUT, f"_{name}.py")
    with open(path, "w") as fo:
        fo.write(src)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    os.remove(path)
    return mod


def models():
    """(修正の後の単体版, 修正の前の単体版)。"""
    sys.path.insert(0, PKG)
    import chi_sr
    return chi_sr.ChiSR(TABLES), module_from(old_source(), "chi_sr_old").ChiSR(TABLES)


def tol_err(c, ref, lt, delta=0.02):
    """年齢の格子 lt = ln T の上の曲線どうしの、時刻の delta のずれを許した誤差。"""
    e = np.zeros(len(lt))
    for k in range(len(lt)):
        w = (lt >= lt[k] - delta) & (lt <= lt[k] + delta)
        e[k] = max(0.0, c[k] - ref[w].max(), ref[w].min() - c[k])
    return e


# ---------- ode ----------
def ref_ode(alpha0, chi0, f_gev, L, t_eval, gfun, rtol=1e-9):
    """査読のときに、E18 の式だけから別に書いた三準位の運動論。fp3.evolve3m とは変数と単位が違う:
    y = (ln nu_211, ln nu_322, ln nu_544, m, j)。nu = N / (G M0^2/hbar c)、m = M/M0、j = J/(G M0^2/c)、時間は 1/mu。
    L = ln(G M0^2/(hbar c))。gfun(alpha, chi) は三つの準位の超放射の率 / mu(占有数の率)。"""
    from scipy.integrate import solve_ivp
    KA, KB, KC, KD, KE, KF = 4.3e-7, 1.1e-8, 1.9e-9, 4.4e-11, 6.1e-11, 4.2e-13       # BGLS Table V・VI
    GW1, GW2, GWT = 1e-2, 3e-8, 5e-6                                                 # BGLS Table IV
    lev = ((2, 1, 1), (3, 2, 2), (5, 4, 4))
    nq = np.exp(-L)
    s4 = (si.MPL_GEV / f_gev) ** 4

    def w_h(a, n, l, m, chi):                             # omega/mu(1908.10370、alpha^5 まで)
        return (1 - a ** 2 / (2 * n * n) - a ** 4 / (8 * n ** 4) + (-6.0 / (2 * l + 1) + 2.0 / n) * a ** 4 / n ** 3
                + 16.0 / (2 * l * (2 * l + 1) * (2 * l + 2)) * chi * m * a ** 5 / n ** 3)

    def rhs(t, y):
        nu = np.exp(np.clip(y[:3], -900, 50))
        m = min(max(y[3], 0.3), 1.0 + 1e-9)
        j = max(y[4], 0.0)
        chi = min(max(j / m ** 2, 0.0), 0.9999)
        al = alpha0 * m
        g = gfun(al, chi)
        rp = 1 + np.sqrt(1 - chi * chi)
        e = nu / m ** 2
        eq = nq / m ** 2
        m2 = m * m
        A = KA * al ** 11 * rp * s4 * e[0] ** 2 * (e[1] + eq) * m2
        B = KB * al ** 8 * s4 * e[1] ** 2 * (e[0] + eq) * m2
        C = KC * al ** 11 * rp * s4 * e[1] ** 2 * (e[2] + eq) * m2
        D = KD * al ** 8 * s4 * e[2] ** 2 * (e[1] + eq) * m2
        E = KE * al ** 8 * s4 * e[1] * e[2] * (e[0] + eq) * m2
        F = KF * al ** 8 * s4 * e[2] ** 2 * (e[0] + eq) * m2
        S = np.array([g[i] * (nu[i] + nq) if g[i] > 0 else g[i] * (nu[i] - 1e-2 * nq) for i in range(3)])
        a1 = GW1 * al ** 14 * e[0] ** 2 * m2
        a2 = GW2 * al ** 18 * e[1] ** 2 * m2
        tr = GWT * al ** 10 * e[1] * (e[0] + eq) * m2
        d1 = S[0] - 2 * A + B + E + F - 2 * a1 + tr
        d2 = S[1] + A - 2 * B - 2 * C + D - E - 2 * a2 - tr
        d3 = S[2] + C - 2 * D - E - 2 * F
        w = [w_h(al, n, l, mm, chi) for (n, l, mm) in lev]
        dm = -alpha0 * (w[0] * S[0] + w[1] * S[1] + w[2] * S[2]) + alpha0 * (A * (2 * w[0] - w[1]) + C * (2 * w[1] - w[2]))
        dj = -(S[0] + 2 * S[1] + 4 * S[2])
        out = np.array([d1 / nu[0], d2 / nu[1], d3 / nu[2], dm, dj])
        return np.where(np.isfinite(out), out, 0.0)

    sol = solve_ivp(rhs, (0.0, t_eval[-1]), [np.log(nq)] * 3 + [1.0, chi0], method="LSODA", rtol=rtol, atol=1e-13, t_eval=t_eval)
    return sol.y[4] / sol.y[3] ** 2 if sol.y.shape[1] == len(t_eval) else np.full(len(t_eval), np.nan)


def part_ode():
    new, _ = models()

    def gfun(al, chi):
        c = min(chi, 0.9999)
        return (float(chisr.g_vec(al, c, 2, 1, 1)), float(chisr.g_vec(al, c, 3, 2, 2)), float(sr.gamma_over_mu(al, c, 5, 4, 4, sr.NORM_N)))

    cases = [(19.6, 0.78, MU, 1e12), (19.6, 0.78, MU, 4e12), (19.6, 0.72, MU, 1e12), (19.6, 0.78, MU, 1e14), (19.6, 0.78, MU, 1e16),
             (19.6, 0.78, MU, 1e19), (10., 0.7, MU, 2e12), (15., 0.7, MU, 1e12), (30., 0.5, MU, 1e13), (45., 0.7, MU, 1e15),
             (20., 0.3, MU, 1e14), (8., 0.9, 1e-12, 1e13)]
    T = np.array([1e5, 1e6, 1e7, 5e7, 1e8, 3e8, 1e9, 6e9, 1.4e10])
    rows = []
    for (M, c0, mu, f) in cases:
        al = sr.alpha_of(M, mu)
        mine = ref_ode(al, c0, f, np.log(sr.N0 * M ** 2), T * sr.YR * mu / sr.HBAR_EVS, gfun)
        s3 = fp3.evolve3m(M, c0, mu, f, T[-1], T)
        ev = s3.y[4] / s3.y[3] ** 2 if s3.y.shape[1] == len(T) else np.full(len(T), np.nan)
        c, r = new(M, c0, T, mu, f)
        un = r <= 2
        rows.append({"M": M, "chi0": c0, "mu": mu, "f": f, "alpha": float(al),
                     "max_abs_ode_minus_evolve3m": float(np.nanmax(np.abs(mine - ev))),
                     "max_abs_table_minus_ode_regime012": float(np.nanmax(np.abs(c - mine)[un])) if un.any() else None,
                     "regimes": sorted(set(r.tolist()))})
    return {"cases": rows, "max_abs_ode_minus_evolve3m": max(r["max_abs_ode_minus_evolve3m"] for r in rows),
            "max_abs_table_minus_ode_regime012": max(r["max_abs_table_minus_ode_regime012"] for r in rows
                                                     if r["max_abs_table_minus_ode_regime012"] is not None)}


# ---------- kink ----------
def part_kink():
    """段差の時刻の折れ(E19 の式の d1、d2)を値の表から測る。切り替えを挟むセルと、比べるためのセル。"""
    tab = chisr.build()
    rt = tab["rt"]
    an = rt["A_NODES"]
    x = np.log(an)
    lfn = rt["LF_NODES"]

    def d12(LT, i):
        lt = np.where(np.isfinite(LT[i - 1:i + 3]), LT[i - 1:i + 3], np.nan)
        d1 = np.abs(lt[2] - (lt[1] + (lt[1] - lt[0]) * (x[i + 1] - x[i]) / (x[i] - x[i - 1])))
        d2 = np.abs(lt[1] - (lt[2] + (lt[2] - lt[3]) * (x[i] - x[i + 1]) / (x[i + 1] - x[i + 2])))
        return d1, d2, np.isfinite(lt[1]) & np.isfinite(lt[2])

    res = {"KS_MIN_dex": reftab.KS_MIN, "cells": []}
    sw = [i for i in range(len(an) - 1) if chisr.switch_cell(an, np.array([np.sqrt(an[i] * an[i + 1])]))[0]]
    for i, label in [(sw[1], "switch322"), (sw[0], "switch211"), (sw[1] + 1, "above322"), (sw[1] + 2, "above322_2"),
                     (sw[1] - 1, "below322"), (sw[1] - 3, "analytic322"), (sw[1] + 5, "table322")]:
        for piece, pn in ((0, "322-grown"), (1, "211-grown")):
            d1, d2, both = d12(rt["LT"][piece], i)
            for flo, fhi in [(8, 14), (14, 20.5)]:
                q = both & ((lfn >= flo) & (lfn <= fhi))[None, :, None, None]
                dmax, dmin = np.fmax(d1, d2)[q], np.fmin(d1, d2)[q]
                g = np.isfinite(dmax)
                if g.sum() < 1000:
                    continue
                res["cells"].append({"cell": [float(an[i]), float(an[i + 1])], "label": label, "piece": pn, "log10f": [flo, fhi],
                                     "levels": int(g.sum()), "max_d_median": float(np.median(dmax[g])),
                                     "max_d_p90": float(np.percentile(dmax[g], 90)), "max_d_p99": float(np.percentile(dmax[g], 99)),
                                     "min_d_median": float(np.median(dmin[g])), "frac_kink": float(np.mean(dmin[g] > reftab.KS_MIN))})
    ks = [reftab.compute_ks(rt["LT"][pc], x) for pc in (0, 1)]
    res["frac_blocks_with_kink_by_cell"] = [{"cell": [float(an[i]), float(an[i + 1])],
                                             "322-grown": float(np.mean(ks[0][i] > 0)), "211-grown": float(np.mean(ks[1][i] > 0))}
                                            for i in range(len(an) - 1)]
    return res


# ---------- rates ----------
def part_rates():
    res = {}
    # 1. 切り替えの alpha での率の飛び(連分数法の表 / 解析式。占有数の率)
    rows = []
    for chi in [0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]:
        r322 = sr.gamma_rel_over_mu(0.1 * (1 + 1e-9), chi, 3, 2, 2, 2.0) / sr.gamma_over_mu(0.1, chi, 3, 2, 2, 2.0)
        r211 = sr.gamma_rel_over_mu(0.03 * (1 + 1e-9), chi, 2, 1, 1, 2.0) / sr.gamma_over_mu(0.03, chi, 2, 1, 1, 2.0)
        rows.append({"chi": chi, "322_at_0.1_table_over_analytic": float(r322), "211_at_0.03_table_over_analytic": float(r211)})
    res["jump_at_switch"] = rows
    tab = chisr.build()
    an = tab["rt"]["A_NODES"]
    res["switch_cells"] = [[float(an[i]), float(an[i + 1])] for i in range(len(an) - 1)
                           if chisr.switch_cell(an, np.array([np.sqrt(an[i] * an[i + 1])]))[0]]
    res["switch_cell_mass_at_mu0.71"] = [[float(a / sr.alpha_of(1.0, MU)) for a in c] for c in res["switch_cells"]]

    # 2. 322 の解析式と、多倍長の連分数法(alpha < 0.1 では連分数法の表がなく、解析式を使っている)
    rows = []
    for a in [0.05, 0.07, 0.09, 0.10]:
        c1 = float(chisr.chic_vec(a, 1, tab=tab))
        c2 = float(chisr.chic_vec(a, 2, tab=tab))
        for label, chi in [("c2+0.25(c1-c2)", c2 + 0.25 * (c1 - c2)), ("c2+0.75(c1-c2)", c2 + 0.75 * (c1 - c2)),
                           ("0.5", 0.5), ("0.7", 0.7), ("0.9", 0.9)]:
            w = leaver.bound_state_mp(a, chi, 3, 2, 2, N=int(max(4000, 500 / a)), dps=50)
            im = float(w.imag)
            ana = float(sr.gamma_over_mu(a, chi, 3, 2, 2) * a)
            rows.append({"alpha": a, "chi": float(chi), "where": label, "im_mp": im, "im_analytic": ana,
                         "analytic_over_mp": ana / im})
    res["rate322_analytic_vs_mp"] = rows
    # alpha = 0.1 のすぐ上の表の値は多倍長と合うか(表は正しく、ずれているのは解析式のほう)
    rows = []
    for chi in [0.3, 0.5, 0.9]:
        w = leaver.bound_state_mp(0.1, chi, 3, 2, 2, N=5000, dps=50)
        t = float(sr.gamma_rel_over_mu(0.1 * (1 + 1e-9), chi, 3, 2, 2) * 0.1)
        rows.append({"chi": chi, "im_mp": float(w.imag), "im_table": t, "table_over_mp": t / float(w.imag)})
    res["rate322_table_vs_mp_at_0.1"] = rows
    # 211 の解析式と多倍長(alpha < 0.03 では解析式)
    rows = []
    for a in [0.02, 0.03]:
        for chi in [0.25, 0.6, 0.9, 0.99]:
            w = leaver.bound_state_mp(a, chi, 2, 1, 1, N=int(500 / a), dps=40)
            ana = float(sr.gamma_over_mu(a, chi, 2, 1, 1) * a)
            rows.append({"alpha": a, "chi": chi, "im_mp": float(w.imag), "im_analytic": ana,
                         "analytic_over_mp": ana / float(w.imag)})
    res["rate211_analytic_vs_mp"] = rows

    # 3. 図のデータ (a): 切り替えの alpha での率の飛び
    chis = np.linspace(0.2, 0.99, 80)
    j322 = [sr.gamma_rel_over_mu(0.1 * (1 + 1e-9), c, 3, 2, 2, 2.0) / sr.gamma_over_mu(0.1, c, 3, 2, 2, 2.0) for c in chis]
    j211 = [sr.gamma_rel_over_mu(0.03 * (1 + 1e-9), c, 2, 1, 1, 2.0) / sr.gamma_over_mu(0.03, c, 2, 1, 1, 2.0) for c in chis]
    np.savetxt(os.path.join(FIG, "e19_jump.dat"), np.c_[chis, j322, j211], comments="", header="chi r322 r211", fmt="%.6g")
    return res


# ---------- stepscan ----------
def part_stepscan():
    """322 から育つ BH(60 Msun、chi0 = 0.29、f = 5e18 GeV)で mu を変えて alpha を動かす。スピンが 0.27 を切る時刻の、
    表の補間 / 基準の運動論と、印 0・1・2 の年齢での誤差の最大(修正の前後)。図のデータ figs/data/e19_step.dat。"""
    new, old = models()
    M, c0, f, lev = 60.0, 0.29, 5e18, 0.27
    Tw = np.logspace(7.6, 10.1, 2400)
    ltw = np.log(Tw)

    def tcross(x):
        j = int(np.argmax(x < lev))
        return float(Tw[j]) if x[j] < lev else np.nan

    rows = []
    for a in np.geomspace(0.0915, 0.128, 56):
        mu = a / sr.alpha_of(M, 1.0)
        s = fp3.evolve3m(M, c0, mu, f, Tw[-1], Tw)
        ref = s.y[4] / s.y[3] ** 2
        cn, rn = new(M, c0, Tw, mu, f)
        co, ro = old(M, c0, Tw, mu, f)
        e = tol_err(cn, ref, ltw)
        rows.append([a, tcross(ref), tcross(cn), e[ro <= 2].max() if (ro <= 2).any() else 0.0,
                     e[rn <= 2].max() if (rn <= 2).any() else 0.0, np.mean(rn == 5), np.mean(ro == 5)])
    rows = np.array(rows)
    ratio = rows[:, 2] / rows[:, 1]
    np.savetxt(os.path.join(FIG, "e19_step.dat"), np.c_[rows[:, 0], ratio, rows[:, 3], rows[:, 4]],
               comments="", header="alpha ratio err_old err_new", fmt="%.6g")
    an = chisr.build()["rt"]["A_NODES"]
    sw = chisr.switch_cell(an, rows[:, 0])
    up = (rows[:, 0] > 0.10846) & (rows[:, 0] < 0.12982)
    return {"M": M, "chi0": c0, "f": f, "level": lev, "alpha_range": [float(rows[0, 0]), float(rows[-1, 0])], "n": len(rows),
            "t_ratio_switch_cell": [float(np.nanmin(ratio[sw])), float(np.nanmax(ratio[sw]))],
            "t_ratio_0.108-0.130": [float(np.nanmin(ratio[up])), float(np.nanmax(ratio[up]))],
            "t_ratio_elsewhere": [float(np.nanmin(ratio[~sw & ~up])), float(np.nanmax(ratio[~sw & ~up]))],
            "max_err_unflagged_old": float(rows[:, 3].max()), "max_err_unflagged_new": float(rows[:, 4].max()),
            "max_err_unflagged_old_0.108-0.130": float(rows[up, 3].max()),
            "share_regime5_new_mean": float(rows[:, 5].mean()), "share_regime5_old_mean": float(rows[:, 6].mean())}


# ---------- noise322 ----------
def _mp_node(args):
    a, c = args
    try:
        return float(leaver.bound_state_mp(a, c, 3, 2, 2, N=int(max(4000, 500 / a)), dps=50).imag)
    except Exception:
        return np.nan


def part_noise322(nproc=6):
    """表の節の上の値(out/rates_322.npz の im)と多倍長の連分数法の比。322 の臨界より上で chi <= 0.6 の節を 3 つおきに。"""
    z = np.load(os.path.join(OUT, "rates_322.npz"))
    AL, CH, IM = z["alpha"], z["chi"], z["im"]
    jobs, idx = [], []
    for i, a in enumerate(AL):
        if not (0.099 <= a <= 0.165):
            continue
        c2 = 8 * a / (4 + 4 * a * a)                    # 322 の臨界の目安(節を選ぶだけ)
        for j, c in enumerate(CH):
            if c2 + 0.005 < c <= 0.60 and j % 3 == 0:
                jobs.append((float(a), float(c)))
                idx.append((i, j))
    with Pool(nproc) as p:
        mpv = p.map(_mp_node, jobs, chunksize=2)
    res = {"n_nodes": len(jobs), "by_alpha": [], "nodes_off_by_more_than_1pct": []}
    rat = {}
    for (i, j), m in zip(idx, mpv):
        rat.setdefault(i, []).append((float(CH[j]), float(IM[i, j] / m)))
        if abs(IM[i, j] / m - 1) > 0.01:
            res["nodes_off_by_more_than_1pct"].append({"alpha": float(AL[i]), "chi": float(CH[j]), "im_table": float(IM[i, j]),
                                                       "im_mp": m, "table_over_mp": float(IM[i, j] / m),
                                                       "im_over_re": float(abs(m) / AL[i])})
    for i in sorted(rat):
        r = np.array([x[1] for x in rat[i]])
        res["by_alpha"].append({"alpha": float(AL[i]), "n": len(r), "min": float(np.nanmin(r)), "max": float(np.nanmax(r))})
    # 率の滑らかさ: 値の表のセルの真ん中での log10(1/g322) と、セルの両端を結ぶ直線の差(s = 臨界のあいだの位置)
    tab = chisr.build()
    an = tab["rt"]["A_NODES"]
    rows = []
    for sfrac in (0.25, 0.3162):
        def lr(a):
            c1 = float(chisr.chic_vec(a, 1, tab=tab))
            c2 = float(chisr.chic_vec(a, 2, tab=tab))
            return -np.log10(float(chisr.g_vec(a, c2 + sfrac * (c1 - c2), 3, 2, 2)))
        for i in range(26, 33):
            dev = lr(np.sqrt(an[i] * an[i + 1])) - 0.5 * (lr(an[i]) + lr(an[i + 1]))
            rows.append({"s": sfrac, "cell": [float(an[i]), float(an[i + 1])], "midpoint_minus_chord_dex": float(dev)})
    res["rate_smoothness"] = rows
    return res


# ---------- e18 ----------
def part_e18():
    new, old = models()
    res = {}
    for seed in (1801, 1802):
        v = np.load(os.path.join(OUT, f"e18_validate2_{seed}.npz"))
        M, c0, mu, f, T = v["M"], v["chi0"], v["mu"], v["f"], v["T"]
        args = (M[:, None], c0[:, None], T[None, :], mu[:, None], f[:, None])
        cn, rn = new(*args)
        co, ro = old(*args)
        d = {"cells": int(cn.size), "old_code_reproduces_stored": bool(np.array_equal(co, v["chisr"]) and np.array_equal(ro, v["regime"]))}
        ch = (rn != ro) | (cn != co)
        d["changed_cells"] = int(ch.sum())
        d["changed_by_regime_old_new"] = {f"{i}->{j}": int(np.sum(ch & (ro == i) & (rn == j)))
                                         for i in range(6) for j in range(6) if np.any(ch & (ro == i) & (rn == j))}
        d["value_changed_cells"] = int(np.sum(cn != co))
        # 値が変わらないセルは保存してある約束の誤差、変わったセル(印 3 から値を返すようになった点)は |新しい値 - 基準|
        # (時刻のずれを許さない誤差で、約束の誤差の上限)
        err = np.where(cn == co, v["err_tol"], np.abs(cn - v["ref"]))
        for tag, r, e in (("old", ro, v["err_tol"]), ("new", rn, err)):
            un = r <= 2
            g = un & np.isfinite(e)
            d[tag] = {"regime_share": {str(k): float(np.mean(r == k)) for k in range(6)},
                      "unflagged_cells": int(un.sum()), "unflagged_unchecked": int(np.sum(un & ~np.isfinite(e))),
                      "tol_max": float(e[g].max()), "violations": int(np.sum(e[g] > 0.03))}
        to1 = ch & (ro == 3) & (rn <= 2)
        d["regime3_to_unflagged"] = {"cells": int(to1.sum()),
                                     "max_abs_new_minus_ref": float(np.nanmax(np.abs(cn - v["ref"])[to1])) if to1.any() else None,
                                     "max_abs_new_minus_chi0": float(np.max(np.abs(cn - c0[:, None])[to1])) if to1.any() else None}
        res[str(seed)] = d
    return res


# ---------- window ----------
def part_window():
    from scipy.optimize import brentq
    new, _ = models()
    e12 = {(5e7, 0.78): 1.85e12, (5e7, 0.72): 1.12e12, (1e8, 0.78): 1.31e12, (1e8, 0.72): 7.93e11,
           (3e8, 0.78): 7.55e11, (3e8, 0.72): 4.58e11}             # E12 の表(半解析の版)
    rows = []
    for (T2, c0), ref in e12.items():
        lf = brentq(lambda x: float(new(19.6, c0, T2, MU, 10 ** x)[0]) - 0.69, 10, 14, xtol=1e-6)
        row = {"T2_yr": T2, "chi0": c0, "f_max_GeV": float(10 ** lf), "f_max_E12": ref, "ratio": float(10 ** lf / ref),
               "regime_at_f_max": int(new(19.6, c0, T2, MU, 10 ** lf)[1])}
        for delta in (0.02, 0.05):
            def h(M):
                c1 = float(new.chic(new.alpha(M, MU), 1))
                return float(new(M, 0.7, 6e9, MU, 10 ** lf)[0]) - (c1 + delta)
            row[f"M_min_1G_delta{delta}"] = float(brentq(h, 6, 30, xtol=1e-4))
        rows.append(row)
    res = {"gw241011_primary": rows}
    sc = {}
    for c0 in (0.78, 0.72):
        vals = [10 ** brentq(lambda x: float(new(19.6, c0, T2, MU, 10 ** x)[0]) - 0.69, 10, 14, xtol=1e-6) * np.sqrt(T2)
                for T2 in (3e7, 1e8, 3e8, 1e9)]
        sc[str(c0)] = [float(x / vals[0]) for x in vals]
    res["f_max_sqrtT_relative"] = sc
    return res


# ---------- coverage ----------
def part_coverage():
    new, old = models()
    rng = np.random.default_rng(11)
    n = 40000
    bins = [(5, 10), (10, 20), (20, 40), (40, 80), (80, 150)]
    res = {"note": "share of samples with regime 3 / 4 / 5 [%]; T log-uniform 1e7-1e10 yr; n = 40000 per cell"}
    for sname, draw in [("gen1_beta", lambda k: rng.beta(0.69, 2.28, k)), ("fast_0.6-0.9", lambda k: rng.uniform(0.6, 0.9, k))]:
        for fname, (flo, fhi) in [("f1e11-1e13", (11, 13)), ("f1e16-1e19", (16, 19))]:
            block = {}
            for mu12 in [0.1, 0.3, 0.71, 1.0, 1.5, 2.0, 3.0]:
                row = {}
                for (lo, hi) in bins:
                    M = np.exp(rng.uniform(np.log(lo), np.log(hi), n))
                    c0 = np.minimum(draw(n), 0.9949)
                    T = 10 ** rng.uniform(7, 10, n)
                    f = 10 ** rng.uniform(flo, fhi, n)
                    _, rn = new(M, c0, T, mu12 * 1e-12, f)
                    _, ro = old(M, c0, T, mu12 * 1e-12, f)
                    row[f"M{lo}-{hi}"] = {"new": [round(100 * float(np.mean(rn == k)), 2) for k in (3, 4, 5)],
                                          "old": [round(100 * float(np.mean(ro == k)), 2) for k in (3, 4, 5)]}
                block[f"mu{mu12}e-12"] = row
            res[f"{sname}_{fname}"] = block
    return res


# ---------- mass ----------
T_MASS = np.array([1e7, 1e8, 1e9, 6e9])


def _mass_work(args):
    M, c0, mu, f = args
    try:
        s = fp3.evolve3m(M, c0, mu, f, T_MASS[-1], T_MASS, max_nfev=600000)
        if s.y.shape[1] != len(T_MASS):
            return None
        return s.y[3], s.y[4] / s.y[3] ** 2
    except Exception:
        return None


def part_mass(nproc=6):
    new, _ = models()
    rng = np.random.default_rng(77)
    n = 360
    res = {"note": "true map (M_birth, chi0, T) -> (M_final, chi_final) from fp3.evolve3m; T = 1e7, 1e8, 1e9, 6e9 yr; "
                   "cells where both calls have regime 0-2"}
    for name, flo, fhi in [("f1e11-1e13", 11, 13), ("f1e14-1e16", 14, 16), ("f1e17-1e20", 17, 20)]:
        M = rng.uniform(8, 60, n)
        mu = rng.uniform(0.3e-12, 1.2e-12, n)
        c0 = rng.uniform(0.3, 0.95, n)
        f = 10 ** rng.uniform(flo, fhi, n)
        keep = sr.alpha_of(M, mu) <= 0.3
        jobs = [x for x, k in zip(zip(M, c0, mu, f), keep) if k]
        with Pool(nproc) as p:
            out = p.map(_mass_work, jobs, chunksize=2)
        dm, de, dt, da = [], [], [], []
        for (Mi, ci, mui, fi), r in zip(jobs, out):
            if r is None:
                continue
            Mf, chif = r
            c_i, r_i = new(Mi, ci, T_MASS, mui, fi)       # 正しい使い方: 生まれたときの質量
            c_f, r_f = new(Mf, ci, T_MASS, mui, fi)       # 観測された(最後の)質量を渡した場合
            ok = (r_i <= 2) & (r_f <= 2)
            dm += list((1 - Mf / Mi)[ok])
            de += list(np.abs(c_f - chif)[ok])
            dt += list(np.abs(c_i - chif)[ok])
            da += [float(sr.alpha_of(Mi, mui))] * int(ok.sum())
        dm, de, dt, da = map(np.array, (dm, de, dt, da))
        big = de > 0.03
        res[name] = {"cells": int(len(de)), "mass_loss_median": float(np.median(dm)), "mass_loss_p95": float(np.percentile(dm, 95)),
                     "mass_loss_max": float(dm.max()),
                     "err_birth_mass_p95": float(np.percentile(dt, 95)), "err_birth_mass_max": float(dt.max()),
                     "err_observed_mass_p95": float(np.percentile(de, 95)), "err_observed_mass_p99": float(np.percentile(de, 99)),
                     "err_observed_mass_max": float(de.max()), "err_observed_mass_frac_gt_0.03": float(np.mean(big)),
                     "alpha_range_gt_0.03": [float(da[big].min()), float(da[big].max())] if big.any() else None,
                     "mass_loss_median_gt_0.03": float(np.median(dm[big])) if big.any() else None}
    return res


# ---------- lowalpha ----------
def part_lowalpha(nproc=6):
    """値の表の下(alpha 0.005--0.008)で、雲が年齢のうちに育ちうる軽い BH。基準の曲線は E18 と同じ年齢の格子(2400 点)で
    解き、格子の年齢すべてで、時刻の 2% のずれを許した誤差を見る(e19_validate.py の dense と同じ)。"""
    from numpy.lib.stride_tricks import sliding_window_view
    import e18_validate2 as v18
    new, old = models()
    rng = np.random.default_rng(81)
    n = 500
    a = np.exp(rng.uniform(np.log(0.005), np.log(0.008), n))
    M = np.exp(rng.uniform(np.log(1.0), np.log(3.5), n))
    mu = a / sr.alpha_of(M, 1.0)
    c0 = rng.uniform(0.5, 0.995, n)
    f = 10 ** rng.uniform(6, 22, n)
    keep = mu * 1.4e10 <= 0.125                         # 値の表の時間の範囲(mu T <= 0.13 eV yr)の中
    M, c0, mu, f = M[keep], c0[keep], mu[keep], f[keep]
    with Pool(nproc) as p:
        refs = np.array(p.map(v18.ref_point, list(zip(M, c0, mu, f)), chunksize=2), float)
    half = 4
    Td = v18.TG[half:v18.NG - half]
    w = sliding_window_view(refs, 2 * half + 1, axis=1)
    lo, hi, okw = np.nanmin(w, 2), np.nanmax(w, 2), np.all(np.isfinite(w), 2)
    moved = np.abs(refs[:, half:v18.NG - half] - c0[:, None])
    res = {"points": int(len(M)), "ages": int(len(Td)), "ode_incomplete_points": int(np.sum(~np.isfinite(refs[:, -1]))),
           "note": "alpha 0.005-0.008, M 1-3.5 Msun, chi0 0.5-0.995, f 1e6-1e22 GeV",
           "points_ref_moved_gt_0.01": int(np.sum(np.nanmax(moved, 1) > 0.01)),
           "points_ref_moved_gt_0.01_M_ge_2": int(np.sum((np.nanmax(moved, 1) > 0.01) & (M >= 2)))}
    for tag, m in (("new", new), ("old", old)):
        c, r = m(M[:, None], c0[:, None], Td[None, :], mu[:, None], f[:, None])
        tol = np.where(okw, np.maximum(0, np.maximum(c - hi, lo - c)), np.nan)
        d = {"regime_share": {str(k): float(np.mean(r == k)) for k in range(6)}}
        for k in (1, 2, 3):
            q = (r == k) & np.isfinite(tol)
            if q.any():
                bad = q & (tol > 0.03)
                d[f"regime{k}"] = {"cells": int(q.sum()), "tol_max": float(tol[q].max()), "violations": int(bad.sum()),
                                   "points_with_violation": int(np.any(bad, 1).sum()),
                                   "points_with_violation_M_ge_2": int((np.any(bad, 1) & (M >= 2)).sum()),
                                   "ref_moved_max": float(moved[q].max()), "cells_ref_moved_gt_0.01": int(np.sum(moved[q] > 0.01))}
        res[tag] = d
    return res


# ---------- grown322 ----------
def part_grown322(nproc=6):
    """322 から育つ BH(chi_c322 < chi0 < chi_c211)の照合を、値の表の alpha のセルごとに。年齢は e19_validate.py の dense と同じ。"""
    from numpy.lib.stride_tricks import sliding_window_view
    import e18_validate2 as v18
    new, old = models()
    an = chisr.build()["rt"]["A_NODES"]
    rng = np.random.default_rng(5201)
    n = 260
    isw = int(np.searchsorted(an, 0.1) - 1)
    cells = [isw - 1, isw, isw + 1, isw + 2, isw + 3, isw + 5, isw + 8]
    pts, lab = [], []
    for i in cells:
        a = np.exp(rng.uniform(np.log(an[i]), np.log(an[i + 1]), n))
        M = np.exp(rng.uniform(np.log(2), np.log(200), n))
        mu = a / sr.alpha_of(M, 1.0)
        ok = (mu >= 1e-13) & (mu <= 1e-11)
        c1, c2 = new.chic(a, 1), new.chic(a, 2)
        c0 = c2 + rng.uniform(0.03, 0.98, n) * (c1 - c2)
        f = 10 ** rng.uniform(15, 20.5, n)
        pts += list(zip(M[ok], c0[ok], mu[ok], f[ok]))
        lab += [i] * int(ok.sum())
    pts, lab = np.array(pts), np.array(lab)
    with Pool(nproc) as p:
        refs = np.array(p.map(v18.ref_point, [tuple(x) for x in pts], chunksize=2), float)
    half = 4
    Td = v18.TG[half:v18.NG - half]
    w = sliding_window_view(refs, 2 * half + 1, axis=1)
    lo, hi, okw = np.nanmin(w, 2), np.nanmax(w, 2), np.all(np.isfinite(w), 2)
    M, c0, mu, f = pts.T
    res = {"points": int(len(M)), "ages": int(len(Td)), "ode_incomplete_points": int(np.sum(~np.isfinite(refs[:, -1]))), "cells": []}
    out = {}
    for tag, m in (("new", new), ("old", old)):
        tol = np.full((len(M), len(Td)), np.nan)
        reg = np.zeros((len(M), len(Td)), int)
        for i0 in range(0, len(M), 100):
            sl = slice(i0, i0 + 100)
            c, r = m(M[sl, None], c0[sl, None], Td[None, :], mu[sl, None], f[sl, None])
            tol[sl] = np.where(okw[sl], np.maximum(0, np.maximum(c - hi[sl], lo[sl] - c)), np.nan)
            reg[sl] = r
        out[tag] = (tol, reg)
    for i in cells:
        q = lab == i
        row = {"cell": [float(an[i]), float(an[i + 1])], "points": int(q.sum())}
        for tag in ("old", "new"):
            tol, reg = out[tag]
            un = (reg[q] <= 2) & np.isfinite(tol[q])
            bad = un & (tol[q] > 0.03)
            row[tag] = {"unflagged_cells": int(un.sum()), "tol_max": float(np.nanmax(np.where(un, tol[q], 0.0))),
                        "violations": int(bad.sum()), "points_with_violation": int(np.any(bad, 1).sum()),
                        "share_regime5": float(np.mean(reg[q] == 5))}
        res["cells"].append(row)
    for tag in ("old", "new"):
        res[f"{tag}_points_with_violation"] = sum(r[tag]["points_with_violation"] for r in res["cells"])
        res[f"{tag}_tol_max"] = max(r[tag]["tol_max"] for r in res["cells"])
    return res


# ---------- plateau ----------
def part_plateau():
    """seed 1901 の thresh の層の 1 点(alpha = 0.261、f = 10^17.498 GeV)。211 と 322 が飽和したあとの平らな所で、表の値は
    f の節 10^17.25 GeV の側の角(重み 1% 未満)の進み具合で止まる。角の節ごとの、最後に届く進み具合を見る。"""
    new, old = models()
    tab = chisr.build()
    rt = tab["rt"]
    an, lfn = rt["A_NODES"], rt["LF_NODES"]
    M, c0, mu, f = 24.782731016985288, 0.8617252152392237, 1.4095933531901494e-12, 3.151058871853566e+17
    a = float(sr.alpha_of(M, mu))
    c1, c2 = float(chisr.chic_vec(a, 1, tab=tab)), float(chisr.chic_vec(a, 2, tab=tab))
    s = float(chisr.w_of_chi(c0, c1, c2))
    T = np.array([1e6, 1e8, 2e9])
    sol = fp3.evolve3m(M, c0, mu, f, T[-1], T)
    ref = sol.y[4] / sol.y[3] ** 2
    cn, rn = new(M, c0, T, mu, f)
    co, ro = old(M, c0, T, mu, f)
    b = lambda chi: float(reftab.b_of_w(chisr.w_of_chi(chi, c1, c2), s))   # noqa: E731
    ia = int(np.searchsorted(an, a) - 1)
    js = int(np.searchsorted(rt["SB"], s - 1) - 1)
    lf = np.log10(f)
    kf = int(np.floor((lf - lfn[0]) / (lfn[1] - lfn[0])))
    tf = (lf - lfn[kf]) / (lfn[kf + 1] - lfn[kf])
    LT = rt["LT"][1]
    fin = {}
    for df in (0, 1):
        v = [reftab.P_LEV2[int(np.sum(~np.isnan(LT[ia + da, js + ds, kf + df, il]))) - 1]
             for da in (0, 1) for ds in (0, 1) for il in (0, 1)]
        fin[df] = [float(min(v)), float(max(v))]
    return {"M": M, "chi0": c0, "mu": mu, "f": f, "alpha": a, "log10f": float(lf), "chi_c211": c1, "chi_c322": c2,
            "ref_plateau": float(ref[1]), "b_ref": b(ref[1]), "table_value": float(cn[1]), "b_table": b(cn[1]),
            "abs_table_minus_ref": float(abs(cn[1] - ref[1])),
            "log10f_nodes": [float(lfn[kf]), float(lfn[kf + 1])], "weight_lower_f_node": float(1 - tf),
            "b_final_corners_lower_f_node": fin[0], "b_final_corners_upper_f_node": fin[1],
            "regime_old": sorted(set(ro.tolist())), "regime_new": sorted(set(rn.tolist()))}


# ---------- m4stage ----------
def part_m4stage():
    """seed 1902 の thresh の層の 1 点(alpha = 0.491、322 の臨界のすぐ上、f = 6e20 GeV)。322 が飽和して 322 の臨界で止まり、
    三準位の運動論では数 Gyr あとに 544 が超放射で飽和して m = 4 の臨界へ下がる。その前に 433 が飽和するまでの時間も出す。"""
    new, old = models()
    M, c0, mu, f = 32.50163894378716, 0.7929942191547834, 2.018668903614814e-12, 5.999180699266122e+20
    a = float(sr.alpha_of(M, mu))
    c2, c3 = float(new.chic(a, 2)), float(new.chic(a, 3))
    T = np.logspace(5, np.log10(1.4e10), 600)
    s = fp3.evolve3m(M, c0, mu, f, T[-1], T)
    ref = s.y[4] / s.y[3] ** 2
    cn, rn = new(M, c0, T, mu, f)
    co, ro = old(M, c0, T, mu, f)
    plateau = float(ref[np.searchsorted(T, 1e8)])
    j = int(np.argmax(ref < 0.5 * (plateau + ref[-1])))
    mus = mu / sr.HBAR_EVS
    g433 = float(chisr.g_vec(a, plateau, 4, 3, 3)) * mus
    g544 = float(sr.gamma_over_mu(a, plateau, 5, 4, 4, sr.NORM_N)) * mus
    lnN = np.log(sr.N0 * M ** 2)
    late = T > 1.2e10
    T2 = np.array([1.43e10, 2e10])                       # 照合の年齢の範囲(1.4e10 yr まで)のすぐ先
    cn2, rn2 = new(M, c0, T2, mu, f)
    co2, ro2 = old(M, c0, T2, mu, f)
    return {"M": M, "chi0": c0, "mu": mu, "f": f, "alpha": a, "chi_c322": c2, "chi_c433": c3,
            "beyond_range_T_yr": T2.tolist(), "regime_old_beyond_range": ro2.tolist(), "regime_new_beyond_range": rn2.tolist(),
            "value_beyond_range": [float(x) for x in cn2],
            "ref_plateau_after_322": plateau, "ref_final": float(ref[-1]), "t_544_step_yr": float(T[j]),
            "t_sat_433_from_plateau_yr": float((np.log(max(plateau - c3, 1e-300) / 3) + lnN) / g433 / sr.YR),
            "g433_over_g544_at_plateau": g433 / g544,
            "table_final": float(cn[-1]), "abs_table_minus_ref_final": float(abs(cn[-1] - ref[-1])),
            "regime_old_after_step": sorted(set(ro[late].tolist())), "regime_new_after_step": sorted(set(rn[late].tolist())),
            "regime_new_before_step": sorted(set(rn[(T > 1e6) & (T < 1e9)].tolist()))}


# ---------- level422 ----------
K_422_BH = 1.5e-7     # 211 x 211 -> 422 x BH(BGLS Table V。Refs/notes/theory_rates.md 4.2 節の抜き書き)。alpha^11 r+
K_422_INF = 1.6e-8    # 322 x 422 -> 211 x inf(Witte & Mummery Table I。同 5 節の抜き書き)。alpha^8


def part_level422():
    """(211, 322) の固定点での 422 の 1 粒子あたりの増加率 / mu(f によらない単位):
         lam422 = g422 + K_422_BH alpha^11 r+ C' e1^2 - K_422_INF alpha^8 e1 e2
       固定点では ka e1^2 = u、e2 = (u + g322)/(2 kb e1) なので、C' = C_a(211 x 211 -> 322 x BH と同じ補正)のとき
         lam422 = g422 - (K_422_INF/(2 K_B) - K_422_BH/K_A) u - (K_422_INF/(2 K_B)) g322。
       422 の超放射の率 g422 は、322 の率との比を多倍長の連分数法で求めて使う(解析式は高いスピンで 20--30% 大きい)。"""
    tab = chisr.build(use_reftab=False)
    cu = K_422_INF / (2 * fp3.K_B) - K_422_BH / fp3.K_A
    cg = K_422_INF / (2 * fp3.K_B)
    res = {"K_422_BH": K_422_BH, "K_422_INF": K_422_INF,
           "leading_order_pump_over_loss": K_422_BH / (K_422_INF * fp3.K_A / (2 * fp3.K_B)),
           "coefficients": {"u": cu, "g322": cg}}
    # 422 と 322 の超放射の率の比(多倍長の連分数法)。211 の臨界より上と、211 と 322 の臨界のあいだ
    rows = []
    for a in [0.05, 0.1, 0.2, 0.3]:
        c1 = float(chisr.chic_vec(a, 1, tab=tab))
        c2 = float(chisr.chic_vec(a, 2, tab=tab))
        for where, chi in [("between", c2 + 0.25 * (c1 - c2)), ("between", c2 + 0.75 * (c1 - c2)), ("above", c1 + 0.02),
                           ("above", 0.5 * (c1 + 0.95)), ("above", 0.95)]:
            N = int(max(4000, 500 / a))
            w4 = leaver.bound_state_mp(a, chi, 4, 2, 2, N=N, dps=50)
            w3 = leaver.bound_state_mp(a, chi, 3, 2, 2, N=N, dps=50)
            ana = float(sr.gamma_over_mu(a, chi, 4, 2, 2) * a)
            rows.append({"alpha": a, "chi": float(chi), "where": where, "im422_mp": float(w4.imag), "im322_mp": float(w3.imag),
                         "im422_over_im322_mp": float(w4.imag) / float(w3.imag), "im422_analytic_over_mp": ana / float(w4.imag)})
    res["rate422_mp"] = rows
    ratio = max(r["im422_over_im322_mp"] for r in rows if r["where"] == "above")
    res["g422_over_g322_above_c211_max"] = ratio
    res["g422_over_g322_between_range"] = [min(r["im422_over_im322_mp"] for r in rows if r["where"] == "between"),
                                           max(r["im422_over_im322_mp"] for r in rows if r["where"] == "between")]
    # (211, 322) の固定点: g422 = ratio x g322(控えめに、比の最大)
    rows = []
    for a in [0.03, 0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3]:
        c1 = float(chisr.chic_vec(a, 1, tab=tab))
        chis = np.unique(np.concatenate([c1 + np.geomspace(1e-6, 0.99 - c1, 60), [0.99]]))
        r = chisr._fixed(np.full_like(chis, a), chis)
        g422 = ratio * r["g2"]
        rp = sr.rplus(chis)
        pump = K_422_BH * a ** 11 * rp * r["e1"] ** 2
        loss = K_422_INF * a ** 8 * r["e1"] * r["e2"]
        lam = g422 + pump - loss
        ccrit = (loss - g422) / pump                    # BH の過程 211 x 211 -> 422 x BH だけにこの倍率が掛かると育ちはじめる
        k = int(np.argmax(lam / r["g2"]))
        rows.append({"alpha": a, "chi_c211": c1, "lam422_positive_anywhere": bool(np.any(lam > 0)),
                     "max_lam422_over_g322": float((lam / r["g2"]).max()), "at_chi": float(chis[k]),
                     "pump_over_loss_range": [float((pump / loss).min()), float((pump / loss).max())],
                     "min_extra_factor_on_BH_process_to_grow": float(ccrit.min()),
                     "check_formula_max_rel_diff": float(np.max(np.abs(lam - (g422 - cu * fp3.K_A * a ** 11 * rp * r["e1"] ** 2 - cg * r["g2"]))
                                                                / np.abs(lam)))})
    res["P3_211_322"] = rows
    res["P3_lam422_over_g322_max"] = max(r["max_lam422_over_g322"] for r in rows)
    res["P3_min_extra_factor"] = min(r["min_extra_factor_on_BH_process_to_grow"] for r in rows)
    # 211 の臨界に近づくと g211 -> 0、u -> g322/3 で、lam422/g322 は次の値に近づく(格子の最大はこの値より小さい)
    res["P3_lam422_over_g322_limit_at_c211"] = ratio - cu / 3 - cg
    return res


# ---------- usability ----------
_PROBE = r'''
import sys, time, pickle, resource, warnings
import numpy as np
sys.path.insert(0, sys.argv[1])
import chi_sr
t0 = time.time(); m = chi_sr.ChiSR(sys.argv[2]); t_load = time.time() - t0
rss_load = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6          # macOS: bytes -> MB(10^6 バイト)


def nbytes(x):
    """値の表として持っている配列の大きさ。"""
    if isinstance(x, np.ndarray):
        return x.nbytes
    if isinstance(x, (list, tuple)):
        return sum(nbytes(y) for y in x)
    if isinstance(x, dict):
        return sum(nbytes(y) for y in x.values())
    return nbytes(x.q) if hasattr(x, "q") else 0


out = {"load_s": t_load, "peak_rss_after_load_MB": rss_load, "table_arrays_MB": nbytes(m.rt) / 1e6}
try:
    out["pickle_bytes"] = len(pickle.dumps(m)); out["pickle"] = "ok"
except Exception as e:
    out["pickle"] = f"{type(e).__name__}: {str(e)[:60]}"
rng = np.random.default_rng(3); n = 1000000
M = rng.uniform(8, 50, n); mu = rng.uniform(0.6e-12, 1.0e-12, n); f = 10 ** rng.uniform(10, 19, n); T = 10 ** rng.uniform(5, 10, n)
for name, c0 in (("gen1", np.minimum(rng.beta(0.69, 2.28, n), 0.9949)), ("uniform", rng.uniform(0, 0.99, n))):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        t0 = time.time(); m(M, c0, T, mu, f); out[f"time_1e6_{name}_s"] = time.time() - t0
    out[f"warnings_{name}"] = sorted({str(x.message)[:50] for x in w})
out["peak_rss_after_calls_MB"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
import json; print(json.dumps(out))
'''

_X86 = '''

def _x86_cast(x):
    """x86-64 で NaN や範囲の外の実数を int64 に直した結果(INT64_MIN)を模す。Apple Silicon では 0 になる。"""
    bad = ~np.isfinite(x)
    out = np.where(bad, 0.0, x).astype(np.int64)
    out[bad] = np.iinfo(np.int64).min
    _x86_cast.nbad += int(bad.sum())
    return out


_x86_cast.nbad = 0
'''


def part_usability():
    res = {}
    work = os.path.join(OUT, "work", "e19_old_pkg")
    os.makedirs(work, exist_ok=True)
    with open(os.path.join(work, "chi_sr.py"), "w") as fo:
        fo.write(old_source())
    probe = os.path.join(OUT, "work", "e19_probe.py")
    with open(probe, "w") as fo:
        fo.write(_PROBE)
    for tag, d in (("old", work), ("new", PKG)):
        r = subprocess.run([sys.executable, probe, d, TABLES], capture_output=True, text=True, check=True)
        res[tag] = json.loads(r.stdout.strip().splitlines()[-1])
    # NaN を整数に直す所: x86 の挙動を模した版で、落ちるか
    rng = np.random.default_rng(1)
    n = 200000
    sim = {}
    with open(os.path.join(PKG, "chi_sr.py")) as fo:
        new_src = fo.read()
    for tag, src in (("old", old_source()), ("new", new_src)):
        assert src.count("np.floor(j).astype(int)") == 2
        src = src.replace("np.floor(j).astype(int)", "_x86_cast(np.floor(j))").replace("import numpy as np\n", "import numpy as np\n" + _X86, 1)
        mod = module_from(src, f"chi_sr_x86_{tag}")
        m = mod.ChiSR(TABLES)
        d = {}
        for name, (flo, fhi) in {"f1e10-1e19": (10, 19), "f1e11-1e13": (11, 13), "f1e17-1e20": (17, 20)}.items():
            M = rng.uniform(5, 60, n)
            c0 = rng.uniform(0, 0.99, n)
            T = 10 ** rng.uniform(5, 10, n)
            mu = rng.uniform(0.6e-12, 1.0e-12, n)
            f = 10 ** rng.uniform(flo, fhi, n)
            mod._x86_cast.nbad = 0
            try:
                m(M, c0, T, mu, f)
                d[name] = {"result": "ok", "nan_casts": mod._x86_cast.nbad}
            except Exception as e:
                d[name] = {"result": f"{type(e).__name__}: {str(e)[:70]}", "nan_casts": mod._x86_cast.nbad}
        sim[tag] = d
    res["x86_nan_cast_simulation"] = sim
    return res


PARTS = {"ode": part_ode, "rates": part_rates, "stepscan": part_stepscan, "noise322": part_noise322, "kink": part_kink, "e18": part_e18, "window": part_window, "coverage": part_coverage, "mass": part_mass,
         "grown322": part_grown322, "plateau": part_plateau, "m4stage": part_m4stage, "lowalpha": part_lowalpha, "level422": part_level422, "usability": part_usability}


def main():
    names = sys.argv[1:] or list(PARTS)
    path = os.path.join(OUT, "e19.json")
    res = json.load(open(path)) if os.path.exists(path) else {}
    for nm in names:
        t0 = time.time()
        res[nm] = PARTS[nm]()
        print(f"--- {nm} ({time.time() - t0:.0f} s)")
        print(json.dumps(res[nm], indent=1))
        with open(path, "w") as fo:
            json.dump({k: res[k] for k in PARTS if k in res}, fo, indent=1)


if __name__ == "__main__":
    main()
