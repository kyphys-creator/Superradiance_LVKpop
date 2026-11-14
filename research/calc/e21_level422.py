"""E21: 422 は (322, 544) の固定点で育つか。育ったあと何が起きるか。印 4 の判定に入れる。

  rates    : 水素型の率の計算(calc/hyd.py)を、公表された表の値(BGLS 2011.11646 の Table V・VI、Witte & Mummery
             2412.03655 の Table I・II)と比べる。n00 の準位の減衰の率、角度の積分(数値積分)、Coulomb 波の動径の積分
             (動径の方程式を直接解く)でも確かめる。四準位(211・322・422・544)のあいだの過程を総当たりして、
             fp3.PROC4 がそれで全部であることと、fp3 の 422 の定数が hyd の値を 2 桁に丸めたものであることを確かめる。
             n00 の準位に共鳴する過程(alpha^7)の式も、WM の表の値と比べる。
  sr422    : 422 と 322 の超放射の率の比(多倍長の連分数法)。322 と 211 の臨界のあいだで fp3.RHO_422 以下であること。
  growth   : 固定点での 422 の増加率。(211, 322) では負、(322, 544) では alpha がある値より上で正。育ったあとの四準位の
             固定点 (322, 422, 544) と、降下の速さの比(率を動かした場合も)。図のデータ figs/data/e21_*.dat。
  ode4     : 四準位の運動論(fp3.evolve4)を、422 を外した三準位と比べる。422 のある値が、三準位の年齢 T の値と
             年齢 (1 + S_422) T の値のあいだにあること(印 4 の規則の前提)を確かめる。
  flags    : 422 の規則で変わる印の数(E19 の coverage と同じ標本と、f = 1e13--1e16 GeV)、速さ。
  validate : E19 の照合の入力(seed 1901・1902・1903)で、規則を外した印が保存した E19 の版の出力と一致すること、
             規則を入れた版でも印 0・1・2 が約束を満たすこと(保存した誤差と、細かい年齢での照合)。並列 3。
  survey   : ほかの n <= 5 の準位(311・411・433・511・522・533)が固定点で育つか。占有数の一次の項を作る過程を総当たりで
             拾い、増加率を alpha 7 点 x スピン 8 点で評価する(超放射の率はその点での多倍長の連分数法)。並列 6。
  x86      : 単体版(E21 の版)を x86-64 の命令で走らせ、この計算機(arm64)の出力と比べる。方法と子プロセスは E20
             (calc/e20_x86.py)と同じ。x86-64 版の numpy・scipy を展開したフォルダを環境変数 E21_X86_SITE で渡す(なければ
             飛ばす)。作業ファイル(約 0.2 GB)は環境変数 E21_WORK のフォルダ(なければ out/work/e21)に置く。
  versions : ほかの版の numpy・scipy の入った Python(環境変数 E21_PYTHONS にコロンで区切って渡す)で share/chi_sr/selftest.py を
             走らせる(E20 の versions と同じ)。
出力: out/e21.json(部分ごとに書き足す)、../figs/data/e21_*.dat
実行: cd research && ~/miniforge3/envs/python313/bin/python calc/e21_level422.py [部分の名前 ...]   (名前なしで全部)
      (Mac Studio では ~/miniforge3/bin/python3。mpmath と sympy が要る)
"""
import json
import os
import sys
import time
import warnings
from multiprocessing import Pool

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import fp3  # noqa: E402
import si  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")
PKG = os.path.join(HERE, "..", "..", "share", "chi_sr")
MU = 0.71e-12
SEEDS = (1901, 1902, 1903)


def L(s):
    return (int(s[0]), int(s[1]), int(s[2]))


def name(lv):
    return "".join(str(x) for x in lv)


# ---------- rates ----------
# 公表された値。BGLS = 2011.11646 の Table V(BH)・VI(放出)、WM = 2412.03655 の Table I(n <= 4)・II(n = 5)。
# BGLS の値はプロジェクトの fp3.py と Refs/notes/theory_rates.md にあるもの、WM の値は arXiv の HTML 版の表から 2026-11-14 に読んだもの。
PUB_BH = [("211", "211", "322", 4.3e-7, "BGLS"), ("211", "411", "322", 2.5e-8, "BGLS, WM I"), ("411", "411", "322", 1.6e-11, "WM I"),
          ("211", "322", "433", 9.1e-8, "BGLS, WM I"), ("211", "211", "422", 1.5e-7, "BGLS, WM I"), ("211", "311", "422", 2.7e-7, "WM I"),
          ("311", "311", "422", 1.7e-11, "WM I"), ("322", "322", "544", 1.9e-9, "BGLS, WM II"), ("322", "422", "544", 1.3e-11, "WM II"),
          ("422", "422", "544", 2.3e-9, "WM II"), ("211", "211", "522", 7.5e-8, "WM II"), ("322", "411", "533", 2e-8, "WM II"),
          ("411", "411", "522", 9.0e-11, "WM II"), ("211", "311", "522", 1.0e-7, "WM II"), ("211", "433", "544", 1.1e-9, "BGLS"),
          ("211", "411", "522", 9.9e-8, "WM II")]
PUB_INF = [("322", "322", "211", 1.1e-8, "BGLS"), ("311", "311", "211", 5.1e-8, "WM I"), ("311", "322", "211", 1.2e-8, "WM I"),
           ("322", "411", "211", 3.8e-9, "BGLS, WM I"), ("411", "411", "211", 1.7e-9, "WM I"), ("322", "433", "211", 2.6e-9, "BGLS, WM I"),
           ("433", "433", "211", 9.2e-11, "BGLS, WM I"), ("322", "422", "211", 1.6e-8, "WM I"), ("422", "422", "211", 1.6e-9, "WM I"),
           ("422", "433", "211", 6.1e-10, "WM I"), ("411", "422", "211", 2.2e-9, "WM I"), ("311", "422", "211", 7.0e-9, "WM I"),
           ("311", "411", "211", 1.9e-8, "WM I"), ("311", "433", "211", 2.2e-10, "WM I"), ("411", "433", "211", 1.1e-10, "WM I"),
           ("544", "544", "322", 4.4e-11, "BGLS"), ("544", "544", "322", 4.3e-11, "WM II"), ("422", "544", "322", 2.5e-11, "WM II"),
           ("433", "544", "322", 7.8e-10, "BGLS, WM II"), ("322", "544", "211", 6.1e-11, "BGLS"), ("544", "544", "211", 4.2e-13, "BGLS"),
           ("433", "544", "211", 1.9e-11, "BGLS"), ("422", "533", "322", 1.2e-9, "WM II"), ("422", "522", "322", 1.6e-9, "WM II"),
           ("422", "522", "211", 3.7e-9, "WM II"), ("522", "544", "322", 2.2e-11, "WM II"), ("533", "544", "322", 1.8e-10, "WM II"),
           ("433", "533", "322", 2.8e-9, "WM II"), ("522", "533", "322", 4.4e-10, "WM II"), ("533", "533", "322", 1.3e-9, "WM II"),
           ("433", "522", "322", 6.3e-10, "WM II"), ("522", "522", "322", 1.6e-10, "WM II"), ("322", "511", "211", 1.7e-9, "WM II")]
BH_ALSO = {("411", "411", "322"): ("BGLS", 9.8e-11)}       # 同じ過程の、もう一方の表の値(合わないもの)
# n00 の準位に共鳴する過程(alpha^7 (M_pl/f)^4 r+ mu)。WM の Table I・II
PUB_RES = [("211", "311", "322", 3.1e-10, "WM I"), ("311", "311", "322", 1.6e-10, "WM I"), ("411", "411", "422", 2.3e-7, "WM I"),
           ("211", "411", "422", 3.2e-11, "WM I"), ("311", "411", "422", 3.8e-13, "WM I"), ("211", "422", "433", 7.8e-11, "WM I"),
           ("411", "422", "433", 2.3e-11, "WM I"), ("311", "422", "433", 7.7e-12, "WM I"), ("322", "522", "544", 3.4e-12, "WM II"),
           ("422", "522", "544", 3.7e-14, "WM II"), ("522", "522", "544", 2.7e-13, "WM II"), ("211", "511", "522", 6.6e-12, "WM II"),
           ("211", "522", "533", 2.6e-11, "WM II"), ("211", "533", "544", 4.6e-13, "WM II")]


def two_digits(x):
    """2 桁に丸めた値。"""
    return float(f"{float(x):.1e}")


def ang_numeric(ys, conj):
    """int dOmega prod Y(conj[i] が True なら複素共役)。cos theta は Gauss-Legendre、phi は等間隔の和(どちらも厳密)。"""
    from scipy.special import sph_harm_y
    x, w = np.polynomial.legendre.leggauss(80)
    ph = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    TH, PH = np.meshgrid(np.arccos(x), ph, indexing="ij")
    f = np.ones_like(TH, dtype=complex)
    for (l, m), c in zip(ys, conj):
        y = sph_harm_y(l, m, TH, PH)
        f *= np.conj(y) if c else y
    return complex(np.sum(w[:, None] * f) * (2 * np.pi / len(ph)))


def radial_ode(i, j, k, l, xmax=12000.0):
    """J_l = int x F_l(-1/q, q x) R_i R_j R_k dx を、動径の方程式 u'' = [l(l+1)/x^2 - 2/x - q^2] u を直接解いて求める
    (mpmath の Coulomb 波を使わない確かめ)。u は原点で x^(l+1)。規格化は遠方の WKB の振幅 sqrt(q/k(x))。"""
    import hyd
    from scipy.integrate import solve_ivp
    q = float(hyd.kappa_inf(i, j, k)["q"])
    P, s = hyd.triple(i, j, k)
    Pf = {p: float(c) for p, c in P.items()}
    sf = float(s)

    def rhs(x, y):
        prod = sum(c * x ** p for p, c in Pf.items()) * np.exp(-sf * x)
        return [y[1], (l * (l + 1) / x ** 2 - 2 / x - q * q) * y[0], x * y[0] * prod]

    x0 = 1e-3
    y0 = [x0 ** (l + 1) * (1 - x0 / (l + 1)), (l + 1) * x0 ** l - (l + 2) * x0 ** (l + 1) / (l + 1), 0.0]
    sol = solve_ivp(rhs, (x0, xmax), y0, method="DOP853", rtol=1e-12, atol=1e-300, first_step=1e-4)
    u, du, I = sol.y[:, -1]
    kl = np.sqrt(q * q + 2 / xmax - l * (l + 1) / xmax ** 2)
    return float(I / np.sqrt(kl / q * (u * u + du * du / kl ** 2)))


def part_rates():
    import mpmath as mp
    import hyd
    res = {"decay_n00_times_n3": [float(hyd.decay_n00(n) * n ** 3) for n in (1, 2, 3, 4)]}
    assert all(abs(v - 8) < 1e-20 for v in res["decay_n00_times_n3"])
    rows = []
    for a, b, c, pub, src in PUB_BH:
        r = hyd.kappa_bh(L(a), L(b), L(c))
        row = {"process": f"{a}x{b}->{c}xBH", "mine": float(r["kappa"]), "published": pub, "source": src,
               "ratio": float(r["kappa"]) / pub, "nu": float(r["nu"])}
        if (a, b, c) in BH_ALSO:
            row["other_table"] = {"source": BH_ALSO[(a, b, c)][0], "published": BH_ALSO[(a, b, c)][1],
                                  "ratio": float(r["kappa"]) / BH_ALSO[(a, b, c)][1]}
        rows.append(row)
    res["bh"] = rows
    rows = []
    for a, b, c, pub, src in PUB_INF:
        r = hyd.kappa_inf(L(a), L(b), L(c))
        ps = sorted(r["partial"].items())
        rows.append({"process": f"{a}x{b}->{c}xinf", "mine": float(r["kappa"]), "published": pub, "source": src,
                     "ratio": float(r["kappa"]) / pub, "ratio_lowest_l_only": float(ps[0][1]) / pub, "q": float(r["q"]),
                     "partial": {str(l): float(v) for l, v in ps}})
    res["inf"] = rows
    rows = []
    for a, b, c, pub, src in PUB_RES:
        r = hyd.kappa_res(L(a), L(b), L(c))
        rows.append({"process": f"{a}x{b}->{c}xBH", "mine": float(r["kappa"]), "published": pub, "source": src,
                     "ratio": float(r["kappa"]) / pub, "n00": r["n"], "delta": float(r["delta"])})
    res["res"] = rows
    rr = np.array([x["ratio"] for x in rows])
    good = np.abs(rr - 1) <= 0.05
    res["res_summary"] = {"n": len(rr), "n_within_5pct": int(good.sum()), "ratio_min_of_those": float(rr[good].min()),
                          "ratio_max_of_those": float(rr[good].max()),
                          "outliers": [{"process": x["process"], "mine": x["mine"], "published": x["published"]}
                                       for x, g_ in zip(rows, good) if not g_]}
    # まとめ: 公表値は 2 桁なので、2 桁に丸めて一致するか(放出は、最低の部分波だけの値が一致する場合も数える)
    rb = np.array([x["ratio"] for x in res["bh"]])
    res["bh_summary"] = {"n": len(rb), "ratio_min": float(rb.min()), "ratio_max": float(rb.max()),
                         "n_same_two_digits": int(sum(two_digits(x["mine"]) == x["published"] for x in res["bh"])),
                         "n_within_5pct": int(np.sum(np.abs(rb - 1) <= 0.05))}
    ri = np.array([x["ratio"] for x in res["inf"]])
    rl = np.array([x["ratio_lowest_l_only"] for x in res["inf"]])
    best = np.where(np.abs(rl - 1) < np.abs(ri - 1), rl, ri)
    res["inf_summary"] = {"n": len(ri), "n_total_within_5pct": int(np.sum(np.abs(ri - 1) <= 0.05)),
                          "n_total_or_lowest_l_within_5pct": int(np.sum(np.abs(best - 1) <= 0.05)),
                          "lowest_l_matches_but_total_does_not": [x["process"] for x, a_, b_ in zip(res["inf"], ri, rl)
                                                                  if abs(b_ - 1) <= 0.05 < abs(a_ - 1)],
                          "outliers": [{"process": x["process"], "ratio": x["ratio"], "source": x["source"]}
                                       for x, v in zip(res["inf"], best) if abs(v - 1) > 0.05]}
    # 独立な数値での確かめ: 角度の積分と、Coulomb 波の動径の積分
    chk = {}
    i, j, k = L("422"), L("544"), L("322")
    M = i[2] + j[2] - k[2]
    chk["ang4_422x544_322"] = {str(l): {"sympy": float(hyd.ang4(i, j, k, l)),
                                        "numeric": ang_numeric([i[1:], j[1:], k[1:], (l, M)], [False, False, True, True]).real}
                               for l in (4, 6, 8)}
    chk["ang3"] = {f"{a}x{b}->{c}": {"sympy": float(hyd.ang3(L(a), L(b), L(c))),
                                     "numeric": ang_numeric([L(a)[1:], L(b)[1:], L(c)[1:]], [False, False, True]).real}
                   for a, b, c in [("322", "422", "544"), ("211", "211", "322"), ("422", "422", "544")]}
    amax = max(abs(v["sympy"] - v["numeric"]) for d in (chk["ang4_422x544_322"], chk["ang3"]) for v in d.values())
    chk["ang_max_abs_diff"] = amax
    assert amax < 1e-10, chk
    P, s = hyd.triple(i, j, k)
    q = hyd.kappa_inf(i, j, k)["q"]
    rad = {}
    for l in (4, 6):
        Jm = mp.quad(lambda x: x * mp.coulombf(l, -1 / q, q * x) * hyd._poly(P, x) * mp.exp(-s * x), hyd._breaks(s))
        Jo = radial_ode(i, j, k, l)
        rad[str(l)] = {"mpmath": float(Jm), "ode": Jo, "rel_diff": abs(Jo) / abs(float(Jm)) - 1}
    chk["radial_422x544_322"] = rad
    assert max(abs(v["rel_diff"]) for v in rad.values()) < 1e-4, chk
    res["cross_checks"] = chk
    # 422 を含む過程の値(このエントリで使うもの)
    key = {}
    for nm, (a, b, c) in {"K_422A": ("211", "211", "422"), "K_422B": ("322", "422", "211"), "K_422C": ("322", "422", "544"),
                          "K_422D": ("422", "544", "322"), "K_422E": ("422", "544", "211"), "K_422Q": ("422", "422", "544"),
                          "K_422R": ("422", "422", "211")}.items():
        kap, kind = hyd.kappa(L(a), L(b), L(c))
        key[nm] = {"process": f"{a}x{b}->{c}x{kind}", "hyd": float(kap), "fp3": getattr(fp3, nm)}
        assert two_digits(kap) == getattr(fp3, nm), (nm, float(kap))
    res["constants_422"] = key
    res["K_422D_partial"] = {str(l): float(v) for l, v in sorted(hyd.kappa_inf(L("422"), L("544"), L("322"))["partial"].items())}
    # 四準位のあいだの過程の総当たり: l = 0 の強制振動を BH が吸収する過程と、放出のすべて
    found = {}
    for a in range(4):
        for b in range(a, 4):
            for c in range(4):
                if c in (a, b):
                    continue
                kap, kind = hyd.kappa(fp3.LEVELS4[a], fp3.LEVELS4[b], fp3.LEVELS4[c])
                if kind is not None and kap != 0:
                    found[(a, b, c)] = (float(kap), kind == "BH")
    proc = {(i_, j_, k_): (kp, bh) for i_, j_, k_, kp, bh in fp3.PROC4}
    assert set(found) == set(proc), (sorted(set(found) ^ set(proc)))
    res["proc4"] = [{"process": f"{name(fp3.LEVELS4[a])}x{name(fp3.LEVELS4[b])}->{name(fp3.LEVELS4[c])}x{'BH' if found[(a, b, c)][1] else 'inf'}",
                     "hyd": found[(a, b, c)][0], "fp3": proc[(a, b, c)][0], "ratio": found[(a, b, c)][0] / proc[(a, b, c)][0]}
                    for (a, b, c) in sorted(found)]
    assert all(found[k_][1] == proc[k_][1] for k_ in found)
    assert all(abs(x["ratio"] - 1) < 0.04 for x in res["proc4"]), res["proc4"]
    return res


# ---------- sr422 ----------
def _leaver_pair(args):
    import leaver
    a, chi = args
    N = int(max(4000, 500 / a))
    w4 = leaver.bound_state_mp(a, chi, 4, 2, 2, N=N, dps=50)
    w3 = leaver.bound_state_mp(a, chi, 3, 2, 2, N=N, dps=50)
    return float(w4.imag), float(w3.imag)


def part_sr422(nproc=6):
    tab = chisr.build(use_reftab=False)
    alphas = [0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.30]
    ws = [0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.98]
    jobs, meta = [], []
    for a in alphas:
        c1 = float(chisr.chic_vec(a, 1, tab=tab))
        c2 = float(chisr.chic_vec(a, 2, tab=tab))
        for w in ws:
            jobs.append((a, c2 + w * (c1 - c2)))
            meta.append((a, w))
    with Pool(nproc) as p:
        out = p.map(_leaver_pair, jobs, chunksize=1)
    rows = []
    for (a, w), (a_, chi), (i4, i3) in zip(meta, jobs, out):
        g3p = float(chisr.g_vec(a, chi, 3, 2, 2))                        # プロジェクトの 322 の率(占有数の率 / mu)
        rows.append({"alpha": a, "w": w, "chi": chi, "im422": i4, "im322": i3, "ratio": i4 / i3,
                     "project_g322_over_mp": g3p / (2 * i3 / a)})
    r = np.array([x["ratio"] for x in rows])
    res = {"rows": rows, "ratio_min": float(r.min()), "ratio_max": float(r.max()), "RHO_422": fp3.RHO_422,
           "ratio_by_alpha": {str(a): [min(x["ratio"] for x in rows if x["alpha"] == a), max(x["ratio"] for x in rows if x["alpha"] == a)]
                              for a in alphas},
           "project_g322_over_mp_range": [min(x["project_g322_over_mp"] for x in rows), max(x["project_g322_over_mp"] for x in rows)]}
    assert res["ratio_max"] <= fp3.RHO_422, res["ratio_max"]
    return res


# ---------- growth ----------
def cbh_of(a):
    return float(np.interp(a, si.CBH_ALPHA, si.CBH_VALUE))


def p2_state(a, chi, C=1.0):
    """(322, 544) の固定点の率と占有数(f によらない et)。"""
    r = chisr._fixed(np.array([a]), np.array([chi]), c_a=C, c_c=C)
    return {k: float(r[k][0]) for k in ("g1", "g2", "g3", "e1", "e2", "f2", "f3")} | {"use2": bool(r["use2"][0]), "use3": bool(r["use3"][0])}


def lam422_p2(a, chi, C=1.0):
    """(322, 544) の固定点での 422 の増加率 / 322 の超放射の率。fp3.rate_422 と、閉じた式の両方。"""
    s = p2_state(a, chi, C)
    rp = sr.rplus(chi)
    lam = float(fp3.rate_422(a, rp * C, 0.0, s["f2"], s["f3"], s["g2"])) / s["g2"]
    closed = fp3.RHO_422 - (2 / 3) * (fp3.K_422C / fp3.K_C + fp3.K_422D / (fp3.K_C * a ** 3 * rp * C))
    return lam, closed


def fp4_state(a, chi, C=1.0, sQ=1.0, sloss=1.0):
    """四準位の固定点 (322, 422, 544)。返り値: dict(e422/e322、e322・e544 の P2 に対する比、降下の速さの比)か None。
    sQ は 422 x 422 -> 544 x BH、sloss は 422 を減らす二つの過程(322 x 422 -> 544 x BH、422 x 544 -> 322 x inf)の率に掛ける倍率。"""
    s = p2_state(a, chi, C)
    rp = sr.rplus(chi)
    g2, g4 = s["g2"], s["g3"]
    g3 = fp3.RHO_422 * g2
    kc, kd = fp3.K_C * a ** 11 * rp * C, fp3.K_D * a ** 8
    e = fp3.fixed_point_422(g2, g3, g4, kc, kd, sloss * fp3.K_422C * a ** 11 * rp * C, sloss * fp3.K_422D * a ** 8,
                            sQ * fp3.K_422Q * a ** 11 * rp * C)
    if e is None:
        return None
    H3 = 2 * g2 * s["f2"] + 4 * g4 * s["f3"]
    H4 = 2 * g2 * e[0] + 2 * g3 * e[1] + 4 * g4 * e[2]
    return {"e422_over_e322": float(e[1] / e[0]), "e322_ratio": float(e[0] / s["f2"]), "e544_ratio": float(e[2] / s["f3"]),
            "H_ratio": float(H4 / H3)}


VARIANTS = {"LO": dict(useC=False), "CBH": dict(useC=True), "CBH_Qhalf": dict(useC=True, sQ=0.5), "CBH_Qthird": dict(useC=True, sQ=1 / 3),
            "CBH_noloss": dict(useC=True, sloss=0.0), "CBH_noloss_Qhalf": dict(useC=True, sloss=0.0, sQ=0.5), "LO_Qhalf": dict(useC=False, sQ=0.5)}


def part_growth():
    tab = chisr.build(use_reftab=False)
    ws = np.linspace(0.01, 0.99, 99)
    res = {"RHO_422": fp3.RHO_422, "S_422": chisr.S_422, "D_422": chisr.D_422}

    def crit(a):
        return float(chisr.chic_vec(a, 1, tab=tab)), float(chisr.chic_vec(a, 2, tab=tab))

    # (322, 544) の固定点: 422 が育つ alpha の下限(主要項と、C_BH 込み)
    grid = np.round(np.arange(0.100, 0.3001, 0.001), 3)
    rows = []
    worst = 0.0
    for a in grid:
        c1, c2 = crit(a)
        row = {"alpha": float(a), "c1": c1, "c2": c2}
        for tag, C in (("LO", 1.0), ("CBH", cbh_of(a))):
            v = np.array([lam422_p2(a, c2 + w * (c1 - c2), C) for w in ws])
            worst = max(worst, float(np.max(np.abs(v[:, 0] - v[:, 1]))))
            row[tag] = {"min": float(v[:, 0].min()), "max": float(v[:, 0].max()), "at_w0.5": float(v[49, 0]),
                        "share_positive": float(np.mean(v[:, 0] > 0))}
        rows.append(row)
    res["closed_form_max_abs_diff"] = worst                         # 式 (2/3)[...] と fp3.rate_422 の差(/ g322)
    assert worst < 1e-9, worst
    res["alpha_min_grow"] = {tag: float(min(r["alpha"] for r in rows if r[tag]["max"] > 0)) for tag in ("LO", "CBH")}
    res["alpha_all_w_grow"] = {tag: float(min(r["alpha"] for r in rows if r[tag]["min"] > 0)) for tag in ("LO", "CBH")}
    res["lam_over_g322_at_0.30"] = {tag: [rows[-1][tag]["min"], rows[-1][tag]["max"]] for tag in ("LO", "CBH")}
    res["lam_table"] = [r for r in rows if round(r["alpha"] * 100) in (12, 13, 14, 15, 16, 18, 20, 21, 22, 24, 26, 28, 30)
                        and abs(r["alpha"] * 100 - round(r["alpha"] * 100)) < 1e-9]
    res["J422_alpha_min"] = float(chisr.ALPHAS[np.any(chisr.build(use_reftab=False)["J422"] > 0, axis=1)][0])
    # (211, 322) の固定点: 主要項でも C_BH 込みでも負(E19 の結論を、判定に使う関数で確かめ直す)
    pmax = {"LO": -np.inf, "CBH": -np.inf}
    for a in [0.03, 0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3]:
        c1, _ = crit(a)
        chis = np.unique(np.concatenate([c1 + np.geomspace(1e-6, 0.99 - c1, 80), [0.99]]))
        for tag, C in (("LO", 1.0), ("CBH", cbh_of(a))):
            r = chisr._fixed(np.full_like(chis, a), chis, c_a=C, c_c=C)
            lam = fp3.rate_422(a, sr.rplus(chis) * C, r["e1"], r["e2"], 0.0, r["g2"]) / r["g2"]
            pmax[tag] = max(pmax[tag], float(lam.max()))
    res["P3_max_lam_over_g322"] = pmax
    assert max(pmax.values()) < 0
    # 四準位の固定点
    A = np.round(np.arange(0.13, 0.3001, 0.01), 2)
    W4 = np.linspace(0.02, 0.98, 25)
    fp = {}
    for vn, kw in VARIANTS.items():
        kw = dict(kw)
        useC = kw.pop("useC")
        best = {"H_ratio": 1.0}
        tabv = []
        r2, r4 = [], []
        for a in A:
            c1, c2 = crit(a)
            hs, es = [], []
            for w in W4:
                st = fp4_state(a, c2 + w * (c1 - c2), cbh_of(a) if useC else 1.0, **kw)
                if st is None:
                    continue
                hs.append(st["H_ratio"])
                es.append(st["e422_over_e322"])
                r2.append(st["e322_ratio"])
                r4.append(st["e544_ratio"])
                if st["H_ratio"] > best["H_ratio"]:
                    best = dict(st, alpha=float(a), w=float(w))
            tabv.append({"alpha": float(a), "H_ratio_max": max(hs) if hs else 1.0, "e422_over_e322_max": max(es) if es else 0.0})
        fp[vn] = {"max": best, "by_alpha": tabv, "e422_over_e322_max": max(x["e422_over_e322_max"] for x in tabv),
                  "e322_ratio_range": [min(r2), max(r2)], "e544_ratio_range": [min(r4), max(r4)]}
    res["fp4"] = fp
    assert fp["CBH_Qhalf"]["max"]["H_ratio"] < 1 + chisr.S_422, fp["CBH_Qhalf"]["max"]
    # 図のデータ: alpha ごとの増加率(臨界のあいだの最小と最大)と、降下の速さの比
    with open(os.path.join(FIG, "e21_lam.dat"), "w") as fo:
        fo.write("alpha lo_min lo_max cbh_min cbh_max\n")
        for r in rows:
            fo.write(f"{r['alpha']:.3f} {r['LO']['min']:.5f} {r['LO']['max']:.5f} {r['CBH']['min']:.5f} {r['CBH']['max']:.5f}\n")
    with open(os.path.join(FIG, "e21_fp4.dat"), "w") as fo:
        fo.write("alpha h_lo h_cbh h_cbh_qhalf e_lo e_cbh e_cbh_qhalf\n")
        for n, a in enumerate(A):
            fo.write(f"{a:.2f} " + " ".join(f"{fp[v]['by_alpha'][n]['H_ratio_max']:.4f}" for v in ("LO", "CBH", "CBH_Qhalf")) + " "
                     + " ".join(f"{fp[v]['by_alpha'][n]['e422_over_e322_max']:.4f}" for v in ("LO", "CBH", "CBH_Qhalf")) + "\n")
    return res


# ---------- ode4 ----------
def _g422_true(a):
    """422 の超放射の率 / mu を、多倍長の連分数法で表にしたもの(この alpha でスピンについて 3 次の補間)。"""
    import leaver
    from scipy.interpolate import CubicSpline
    chis = np.linspace(0.05, 0.995, 28)
    g = [2 * float(leaver.bound_state_mp(a, float(c), 4, 2, 2, N=int(max(4000, 500 / a)), dps=50).imag) / a for c in chis]
    return CubicSpline(chis, g)


def _ode_case(args):
    alpha, chi0, f, variant, true_g = args
    kw = dict(VARIANTS[variant])
    useC = kw.pop("useC")
    cbh = cbh_of(alpha) if useC else 1.0
    ks = {}
    if "sQ" in kw:
        ks[11] = kw["sQ"]                                   # PROC4[11] = 422 x 422 -> 544 x BH
    if "sloss" in kw:
        ks[8] = ks[9] = kw["sloss"]                         # 322 x 422 -> 544 x BH、422 x 544 -> 322 x inf
    M = alpha / sr.alpha_of(1.0, MU)
    S = chisr.S_422
    te = np.geomspace(1e2, 1.4e10, 500)
    tt = np.unique(np.concatenate([te, te * (1 + S)]))
    g4 = _g422_true(alpha) if true_g else None

    s4 = fp3.evolve4(M, chi0, MU, f, tt[-1], tt, cbh=cbh, with422=True, g422=g4, kscale=ks)
    s3 = fp3.evolve4(M, chi0, MU, f, tt[-1], tt, cbh=cbh, with422=False, kscale=ks)
    if s4.status != 0 or s3.status != 0:
        return {"alpha": alpha, "chi0": chi0, "f": f, "variant": variant, "true_g422": true_g, "failed": True}
    m4, m3 = s4.methods, s3.methods
    i0 = np.searchsorted(tt, te)
    i1 = np.searchsorted(tt, te * (1 + S))
    c4, c3, c3s = s4.y[4][i0], s3.y[4][i0], s3.y[4][i1]
    d = c3 - c4                                             # 422 が下げた分(>= 0 のはず)
    k = int(np.argmax(np.abs(d)))
    e2, e3 = s4.y[1][i0], s4.y[2][i0]
    grown = (e3 > e2 - 5.0) & (e2 > -60)                    # 322 が育ったあとで、422 が 322 の e^-5 倍より多い時刻
    return {"alpha": alpha, "chi0": chi0, "f": f, "variant": variant, "true_g422": true_g, "failed": False,
            "max_abs_dchi": float(np.abs(d).max()), "at_T_yr": float(te[k]), "chi4_there": float(c4[k]), "chi3_there": float(c3[k]),
            "min_chi3_minus_chi4": float(d.min()),                    # 負なら、422 が降下を遅らせた時刻がある
            "min_chi4_minus_chi3_later": float((c4 - c3s).min()),     # 負なら、(1 + S) 倍の年齢の三準位の値より下にいる時刻がある
            "max_e422_over_e322": float(np.exp((e3 - e2)[e2 > -60].max())) if np.any(e2 > -60) else 0.0,   # 322 が育ったあと
            "grew": bool(grown.any()), "final_chi4": float(c4[-1]), "final_chi3": float(c3[-1]), "solver": [m4, m3],
            "_curve": (te, c3, c4, c3s, s4.y[:4][:, i0]) if (alpha, chi0, f, variant, true_g) == FIG_CASE else None}


FIG_CASE = (0.25, 0.9, 1e16, "LO", False)


def part_ode4(nproc=6):
    cases = []
    for alpha in (0.15, 0.20, 0.25, 0.30):
        tabc = chisr.build(use_reftab=False)
        c1 = float(chisr.chic_vec(alpha, 1, tab=tabc))
        c2 = float(chisr.chic_vec(alpha, 2, tab=tabc))
        for chi0 in (round(c2 + 0.6 * (c1 - c2), 3), 0.9):
            for f in (1e15, 1e16):
                for v in ("LO", "CBH", "CBH_Qhalf"):
                    cases.append((alpha, chi0, f, v, False))
    cases += [(0.25, 0.6, 1e16, "LO", True), (0.30, 0.7, 1e16, "LO", True), (0.20, 0.6, 1e16, "CBH", True)]
    with Pool(nproc) as p:
        out = p.map(_ode_case, cases, chunksize=1)
    S = chisr.S_422
    curve = [r.pop("_curve", None) for r in out]
    res = {"S_422": S, "n_cases": len(out), "n_failed": int(sum(r["failed"] for r in out)), "cases": out}
    ok = [r for r in out if not r["failed"]]
    res["max_abs_dchi"] = {v: max(r["max_abs_dchi"] for r in ok if r["variant"] == v and not r["true_g422"]) for v in ("LO", "CBH", "CBH_Qhalf")}
    res["max_e422_over_e322"] = {v: max(r["max_e422_over_e322"] for r in ok if r["variant"] == v and not r["true_g422"]) for v in ("LO", "CBH", "CBH_Qhalf")}
    res["sandwich"] = {"min_chi3_minus_chi4": min(r["min_chi3_minus_chi4"] for r in ok),
                       "min_chi4_minus_chi3_later": min(r["min_chi4_minus_chi3_later"] for r in ok)}
    res["grew_cases"] = int(sum(r["grew"] for r in ok))
    res["grew_by_variant_alpha"] = {v: {str(a): sorted({r["grew"] for r in ok if r["variant"] == v and r["alpha"] == a and not r["true_g422"]})
                                        for a in (0.15, 0.20, 0.25, 0.30)} for v in ("LO", "CBH", "CBH_Qhalf")}
    res["true_g422_cases"] = [{k: r[k] for k in ("alpha", "chi0", "f", "variant", "max_abs_dchi", "max_e422_over_e322")} for r in ok if r["true_g422"]]
    res["cases_not_solved_by_lsoda_alone"] = int(sum(r["solver"] != [["LSODA"], ["LSODA"]] for r in ok))
    res["cases_not_solved_by_lsoda_alone"] = int(sum(r["solver"] != [["LSODA"], ["LSODA"]] for r in ok))
    res["no_growth_cases_max_abs_dchi"] = max([r["max_abs_dchi"] for r in ok if not r["grew"]] or [0.0])
    assert res["n_failed"] == 0, [r for r in out if r["failed"]]
    assert res["sandwich"]["min_chi3_minus_chi4"] > -2e-4 and res["sandwich"]["min_chi4_minus_chi3_later"] > -2e-4, res["sandwich"]
    # 図のデータ(一例)
    te, c3, c4, c3s, lne = [c for c in curve if c is not None][0]
    with open(os.path.join(FIG, "e21_ode_chi.dat"), "w") as fo:
        fo.write("t chi3 chi4 chi3later\n")
        for n in range(len(te)):
            fo.write(f"{te[n]:.5e} {c3[n]:.6f} {c4[n]:.6f} {c3s[n]:.6f}\n")
    with open(os.path.join(FIG, "e21_ode_eps.dat"), "w") as fo:
        fo.write("t e211 e322 e422 e544\n")
        for n in range(len(te)):
            fo.write(f"{te[n]:.5e} " + " ".join(f"{max(v / np.log(10), -90):.4f}" for v in lne[:, n]) + "\n")
    return res


# ---------- flags ----------
def model():
    sys.path.insert(0, PKG)
    import chi_sr
    return chi_sr.ChiSR(os.path.join(PKG, "chi_sr_tables.npz"))


def part_flags():
    m = model()
    off = lambda *x: m._eval(*x, False)                    # noqa: E731  422 の規則を外した印(E20 の版と同じ)
    res = {"note": "share of samples with regime 3 / 4 / 5 [%]; T log-uniform 1e7-1e10 yr; n = 40000 per cell; "
                   "'e20' = the same code without the 422 rule"}
    bins = [(5, 10), (10, 20), (20, 40), (40, 80), (80, 150)]

    def block(rng, draw, flo, fhi, mus):
        out = {}
        for mu12 in mus:
            row = {}
            for (lo, hi) in bins:
                n = 40000
                M = np.exp(rng.uniform(np.log(lo), np.log(hi), n))
                c0 = np.minimum(draw(n), 0.9949)
                T = 10 ** rng.uniform(7, 10, n)
                f = 10 ** rng.uniform(flo, fhi, n)
                c1, r1 = m(M, c0, T, mu12 * 1e-12, f)
                c0_, r0 = off(M, c0, T, mu12 * 1e-12, f)
                assert np.array_equal(c1, c0_)
                ch = r1 != r0
                assert np.all(np.isin(r0[ch], (1, 5)) & (r1[ch] == 4))
                row[f"M{lo}-{hi}"] = {"new": [round(100 * float(np.mean(r1 == k)), 2) for k in (3, 4, 5)],
                                      "e20": [round(100 * float(np.mean(r0 == k)), 2) for k in (3, 4, 5)],
                                      "changed_1to4": int(np.sum(ch & (r0 == 1))), "changed_5to4": int(np.sum(ch & (r0 == 5)))}
            out[f"mu{mu12}e-12"] = row
        return out

    # E19 の coverage と同じ標本(同じ乱数の列。E19 の 'new' の列が、ここの 'e20' の列と一致する)
    rng = np.random.default_rng(11)
    for sname, draw in [("gen1_beta", lambda k: rng.beta(0.69, 2.28, k)), ("fast_0.6-0.9", lambda k: rng.uniform(0.6, 0.9, k))]:
        for fname, (flo, fhi) in [("f1e11-1e13", (11, 13)), ("f1e16-1e19", (16, 19))]:
            res[f"{sname}_{fname}"] = block(rng, draw, flo, fhi, [0.1, 0.3, 0.71, 1.0, 1.5, 2.0, 3.0])
    e19 = json.load(open(os.path.join(OUT, "e19.json")))["coverage"]
    same = all(res[k][mk][bk]["e20"] == e19[k][mk][bk]["new"] for k in e19 if k != "note" for mk in e19[k] for bk in e19[k][mk])
    res["e20_column_equals_e19_coverage"] = bool(same)
    assert same
    # f = 1e13--1e16 GeV((322, 544) での降下が年齢のうちに進む範囲)
    rng = np.random.default_rng(2111)
    for sname, draw in [("gen1_beta", lambda k: rng.beta(0.69, 2.28, k)), ("fast_0.6-0.9", lambda k: rng.uniform(0.6, 0.9, k))]:
        res[f"{sname}_f1e13-1e16"] = block(rng, draw, 13, 16, [0.3, 0.71, 1.5])
    # 比べるために: 422 が育ちうる点を一律に印 4 にした場合(E15 の規則と同じ扱い)。プロジェクト版で D_422 を負にすると、
    # 422 が育ちうる点がすべて印 4 になる。増える点のうち、(322, 544) での降下(211 の臨界か初期スピンから返す値まで)が小さい点の割合も数える
    tab = chisr.build()
    rng = np.random.default_rng(2112)
    bl = {}
    for sname, draw in [("gen1_beta", lambda k: rng.beta(0.69, 2.28, k)), ("fast_0.6-0.9", lambda k: rng.uniform(0.6, 0.9, k))]:
        for (lo, hi) in [(20, 40), (40, 80)]:
            n = 40000
            M = np.exp(rng.uniform(np.log(lo), np.log(hi), n))
            c0 = np.minimum(draw(n), 0.9949)
            T = 10 ** rng.uniform(7, 10, n)
            f = 10 ** rng.uniform(13, 16, n)
            v, r_new = chisr.chi_sr(tab, M, c0, T, MU, f)
            _, r_off = chisr.chi_sr(tab, M, c0, T, MU, f, rule422=False)
            D0 = chisr.D_422
            try:
                chisr.D_422 = -1.0
                _, r_all = chisr.chi_sr(tab, M, c0, T, MU, f)
            finally:
                chisr.D_422 = D0
            c1v = chisr.chic_vec(sr.alpha_of(M, MU), 1, tab=tab)
            move = np.maximum(np.minimum(c0, c1v) - v, 0.0)
            add = (r_all == 4) & (r_off != 4)
            assert np.all((r_new == 4) <= (r_all == 4))
            bl[f"{sname}_M{lo}-{hi}"] = {"regime4_pct": {"without_rule": round(100 * float(np.mean(r_off == 4)), 2),
                                                          "rule": round(100 * float(np.mean(r_new == 4)), 2),
                                                          "blanket": round(100 * float(np.mean(r_all == 4)), 2)},
                                         "blanket_added": int(add.sum()),
                                         "blanket_added_share_moved_lt_0.015": float(np.mean(move[add] < 0.015)) if add.any() else None,
                                         "blanket_added_share_moved_lt_0.001": float(np.mean(move[add] < 0.001)) if add.any() else None}
    res["blanket_rule_f1e13-1e16_mu0.71"] = bl
    # 422 が育ちうる点の数と、そのうち印 4 にした数(alpha 0.12--0.30 の標本。export_chisr.py の確かめ 7 と同じ入力)
    n = 100000
    rng = np.random.default_rng(21)
    a = rng.uniform(0.12, 0.30, n)
    M = np.exp(rng.uniform(np.log(1), np.log(300), n))
    mu = a / sr.alpha_of(M, 1.0)
    c0 = rng.uniform(0.0, 0.995, n)
    T = np.minimum(10 ** rng.uniform(5, np.log10(1.4e10), n), 0.12 / mu)
    f = 10 ** rng.uniform(8, 19, n)
    c1, r1 = m(M, c0, T, mu, f)
    _, r0 = off(M, c0, T, mu, f)
    cl, rl = off(M, c0, T * (1 + chisr.S_422), mu, f)
    c1v = chisr.chic_vec(a, 1, tab=tab)
    below = (r0 == 1) & (c1 < c1v)                                   # 印 1 で、返す値が 211 の臨界より下((322, 544) にいる)
    width = c1 - cl
    ch = r1 != r0
    res["alpha0.12-0.30"] = {"n": n, "regime1_below_c211": int(below.sum()), "changed_to_4": int(ch.sum()),
                             "changed_from": {"1": int(np.sum(ch & (r0 == 1))), "5": int(np.sum(ch & (r0 == 5)))},
                             "width_of_changed": {"min": float(width[ch & (rl != 3)].min()), "median": float(np.median(width[ch & (rl != 3)])),
                                                  "max": float(width[ch & (rl != 3)].max())},
                             "changed_because_later_age_unavailable": int(np.sum(ch & (rl == 3))),
                             "regime1_below_c211_width_gt_D": int(np.sum(below & (width > chisr.D_422))),
                             "regime1_below_c211_width_le_D": int(np.sum(below & (width <= chisr.D_422)))}
    # 速さ(E19 の usability と同じ入力)
    rng = np.random.default_rng(3)
    n = 1000000
    M = rng.uniform(8, 50, n)
    mu = rng.uniform(0.6e-12, 1.0e-12, n)
    f = 10 ** rng.uniform(10, 19, n)
    T = 10 ** rng.uniform(5, 10, n)
    sp = {}
    for nm, c0 in (("gen1", np.minimum(rng.beta(0.69, 2.28, n), 0.9949)), ("uniform", rng.uniform(0, 0.99, n))):
        t0 = time.time()
        m(M, c0, T, mu, f)
        t1 = time.time()
        off(M, c0, T, mu, f)
        sp[nm] = {"with_rule_s": t1 - t0, "without_rule_s": time.time() - t1}
    res["time_1e6"] = sp
    return res


# ---------- validate ----------
def _validate_seed(seed):
    import e19_validate as v19
    from numpy.lib.stride_tricks import sliding_window_view
    m = model()
    off = lambda *x: m._eval(*x, False)                    # noqa: E731
    TOL = v19.TOL
    z = np.load(os.path.join(OUT, f"e19_validate_{seed}.npz"))
    M, c0, mu, f, lab, T = z["M"], z["chi0"], z["mu"], z["f"], z["kind"], z["T"]
    N, nT = len(M), len(T)

    def run(fn, Tg):
        c = np.empty((N, len(Tg)))
        r = np.empty((N, len(Tg)), np.int8)
        for i0 in range(0, N, 200):
            sl = slice(i0, min(i0 + 200, N))
            c[sl], r[sl] = fn(M[sl, None], c0[sl, None], Tg[None, :], mu[sl, None], f[sl, None])
        return c, r

    c1, r1 = run(m, T)
    c0_, r0 = run(off, T)
    res = {"seed": seed, "n_points": int(N), "n_ages": int(nT)}
    res["rule_off_reproduces_stored_e19"] = {"values_equal_float32": bool(np.array_equal(c0_.astype(np.float32), z["chisr"])),
                                             "regime_mismatch": int(np.sum(r0 != z["regime"]))}
    assert res["rule_off_reproduces_stored_e19"]["values_equal_float32"] and res["rule_off_reproduces_stored_e19"]["regime_mismatch"] == 0
    assert np.array_equal(c1, c0_)
    ch = r1 != r0
    assert np.all(np.isin(r0[ch], (1, 5)) & (r1[ch] == 4))
    tol = z["err_tol"].astype(float)
    unfl = r1 <= 2
    fin = np.isfinite(tol)
    res["coarse"] = {"cells": int(N * nT), "changed": {"1->4": int(np.sum(ch & (r0 == 1))), "5->4": int(np.sum(ch & (r0 == 5)))},
                     "points_with_change": int(np.sum(ch.any(axis=1))),
                     "regime_share": {str(k): float(np.mean(r1 == k)) for k in range(6)},
                     "regime_share_e20": {str(k): float(np.mean(r0 == k)) for k in range(6)},
                     "unflagged_cells_checked": int(np.sum(unfl & fin)), "unflagged_violations": int(np.sum(unfl & fin & (tol > TOL))),
                     "unflagged_tol_max": float(tol[unfl & fin].max()),
                     "changed_by_kind": {str(k): int(np.sum(ch[lab == k])) for k in np.unique(lab)},
                     "err_tol_of_changed_cells_max": float(np.nanmax(np.where(ch, tol, np.nan))) if ch.any() else None}
    # 細かい年齢(E19 の evaluate_dense と同じ窓)
    rpath = os.path.join(OUT, f"e19_refs_{seed}.npz")
    if os.path.exists(rpath):
        refs = np.load(rpath)["chi"].astype(float)
        half = 4
        TG, NG, TS = v19.TG, v19.NG, v19.TS
        Td = TG[half:NG - half]
        sel = (Td >= TS[0] * (1 - 1e-12)) & (Td <= TS[-1] * (1 + 1e-12))
        Td = Td[sel]
        sw = sliding_window_view(refs, 2 * half + 1, axis=1)
        lo, hi = np.nanmin(sw, axis=2)[:, sel], np.nanmax(sw, axis=2)[:, sel]
        ok = np.all(np.isfinite(sw), axis=2)[:, sel]
        cd, rd = run(m, Td)
        _, rd0 = run(off, Td)
        e = np.where(ok, np.maximum(0, np.maximum(cd - hi, lo - cd)), np.nan)
        un = rd <= 2
        fin = np.isfinite(e)
        chd = rd != rd0
        e19 = json.load(open(os.path.join(OUT, f"e19_dense_{seed}.json")))["new"]
        res["dense"] = {"n_ages": int(len(Td)), "cells": int(e.size), "changed": {"1->4": int(np.sum(chd & (rd0 == 1))), "5->4": int(np.sum(chd & (rd0 == 5)))},
                        "regime_share": {str(k): float(np.mean(rd == k)) for k in range(6)},
                        "regime_share_e20": {str(k): float(np.mean(rd0 == k)) for k in range(6)},
                        "unflagged_cells_checked": int(np.sum(un & fin)), "unflagged_violations": int(np.sum(un & fin & (e > TOL))),
                        "unflagged_tol_max": float(e[un & fin].max()),
                        "e20_unflagged_cells_checked": int(np.sum((rd0 <= 2) & fin)), "e20_unflagged_tol_max": float(e[(rd0 <= 2) & fin].max()),
                        "e19_stored": {"unflagged_cells": e19["unflagged"]["cells"], "violations": e19["unflagged"]["violations"],
                                       "tol_max": e19["unflagged"]["tol_max"]}}
        assert res["dense"]["e20_unflagged_cells_checked"] == e19["unflagged"]["cells"]
        assert abs(res["dense"]["e20_unflagged_tol_max"] - e19["unflagged"]["tol_max"]) < 1e-6
        assert res["dense"]["unflagged_violations"] == 0
    else:
        res["dense"] = "skipped (out/e19_refs_<seed>.npz is not on this machine)"
    assert res["coarse"]["unflagged_violations"] == 0
    return res


def part_validate(nproc=3):
    with Pool(nproc) as p:
        out = p.map(_validate_seed, SEEDS, chunksize=1)
    return {str(r["seed"]): r for r in out}


# ---------- survey ----------
L211, L322, L422, L544 = (2, 1, 1), (3, 2, 2), (4, 2, 2), (5, 4, 4)
SURVEY = [(3, 1, 1), (4, 1, 1), (4, 3, 3), (5, 1, 1), (5, 2, 2), (5, 3, 3)]


def linear_terms(X, pop):
    """準位 X の占有数について一次の項を作る過程: [(符号, 種類, kappa, (相手の準位 a, b))]。pop は占有されている準位。
    得: a x b -> X x (BH か inf)。失: X x a -> b x (BH か inf)。種類は "BH"(alpha^11 r+)、"res"(共鳴、alpha^7 r+)、"inf"(alpha^8)。"""
    import hyd

    def one(i, j, k):
        if hyd.resonant_n(i, j, k) is not None:
            return ("res", float(hyd.kappa_res(i, j, k)["kappa"])) if hyd.ang3(i, j, k) != 0 else None
        kap, kind = hyd.kappa(i, j, k)
        return (kind, float(kap)) if kind is not None and kap != 0 else None

    out = []
    for ia, a in enumerate(pop):
        for b in pop[ia:]:
            r = one(a, b, X)
            if r:
                out.append((+1, r[0], r[1], (a, b)))
    for a in pop:
        for b in pop:
            if a != b:
                r = one(X, a, b)
                if r:
                    out.append((-1, r[0], r[1], (a, b)))
    return out


def _leaver_checked(args):
    """準束縛状態の Im(M omega)(多倍長の連分数法)。実部が隣の n の準位に飛んでいたら NaN。"""
    import leaver
    a, chi, lv = args
    n = lv[0]
    try:
        w = leaver.bound_state_mp(a, chi, *lv, N=int(max(4000, 500 / a)), dps=50)
    except Exception:
        return float("nan")
    x = float(w.real) / a
    near = min(range(max(lv[1] + 1, n - 1), n + 2), key=lambda k: abs(x - (1 - a ** 2 / (2 * k ** 2))))
    return float(w.imag) if near == n else float("nan")


def part_survey(nproc=6):
    tab = chisr.build(use_reftab=False)
    alphas = [0.05, 0.10, 0.14, 0.18, 0.22, 0.26, 0.30]
    ws = [0.05, 0.3, 0.6, 0.9, 1.02, 1.2, 1.5, 1.9]          # w < 1: 322 と 211 の臨界のあいだ(P2)、w > 1: 211 の臨界より上(P3)
    pts = []
    for a in alphas:
        c1 = float(chisr.chic_vec(a, 1, tab=tab))
        c2 = float(chisr.chic_vec(a, 2, tab=tab))
        pts += [(a, w, float(chisr.chi_of_w(w, c1, c2))) for w in ws]
    with Pool(nproc) as p:
        im = np.array(p.map(_leaver_checked, [(a, chi, lv) for (a, w, chi) in pts for lv in SURVEY], chunksize=2)).reshape(len(pts), len(SURVEY))
    pops = {"P3": [L211, L322], "P2": [L322, L544], "P2x": [L322, L422, L544]}
    terms = {(X, fp): linear_terms(X, pop) for X in SURVEY for fp, pop in pops.items()}
    rows = []
    for ip, (a, w, chi) in enumerate(pts):
        rp = sr.rplus(chi)
        g2p = float(chisr.g_vec(a, chi, 3, 2, 2))
        jint = float(chisr.extra_growth(a, chi))                        # いまの印 4 の判定の増加率(411・433・544 の最大、C_BH 込み)
        for tag, C in (("LO", 1.0), ("CBH", cbh_of(a))):
            r = chisr._fixed(np.array([a]), np.array([chi]), c_a=C, c_c=C)
            g2, g4 = float(r["g2"][0]), float(r["g3"][0])
            if w > 1:
                occs = {"P3": {L211: float(r["e1"][0]), L322: float(r["e2"][0])}}
            else:
                occs = {"P2": {L322: float(r["f2"][0]), L544: float(r["f3"][0])}}
                e = fp3.fixed_point_422(g2, fp3.RHO_422 * g2, g4, fp3.K_C * a ** 11 * rp * C, fp3.K_D * a ** 8, fp3.K_422C * a ** 11 * rp * C,
                                        fp3.K_422D * a ** 8, fp3.K_422Q * a ** 11 * rp * C)
                if e is not None:
                    occs["P2x"] = {L322: float(e[0]), L422: float(e[1]), L544: float(e[2])}
            for il, X in enumerate(SURVEY):
                gX = 2 * im[ip, il] / a
                for fp, occ in occs.items():
                    lam = gX
                    for sgn, kind, kap, (p_, q_) in terms[(X, fp)]:
                        lam += sgn * kap * {"BH": a ** 11 * rp * C, "res": a ** 7 * rp * C, "inf": a ** 8}[kind] * occ[p_] * occ[q_]
                    rows.append({"level": name(X), "alpha": a, "w": w, "chi": chi, "variant": tag, "fp": fp, "gX_over_g322": gX / g2p,
                                 "lam_over_g322": lam / g2p, "J_integrand_over_g322": jint / g2p})
    res = {"n_points": len(pts), "alphas": alphas, "ws": ws, "n_nan_rates": int(np.sum(~np.isfinite(im))),
           "terms": {f"{name(X)}@{fp}": [{"sign": sg, "kind": kd, "kappa": kp, "with": f"{name(p_)}x{name(q_)}"} for sg, kd, kp, (p_, q_) in t]
                     for (X, fp), t in terms.items()}}
    # E15 の式(fp3.extra_rates)の項と同じであること: 411・433 の P3 と、433 の P2
    k = {(f"{name(X)}@{fp}", t["kind"], t["with"]): t["sign"] * t["kappa"] for (X, fp) in terms
         for t in res["terms"][f"{name(X)}@{fp}"]}
    for key, ref in [(("411@P3", "BH", "211x322"), -fp3.K_P1), (("411@P3", "inf", "322x211"), -fp3.K_Q1), (("433@P3", "BH", "211x322"), fp3.K_G),
                     (("433@P3", "inf", "322x211"), -fp3.K_E2), (("433@P2", "inf", "544x322"), -fp3.K_K1)]:
        assert abs(k[key] / ref - 1) < 0.03, (key, k[key], ref)
    assert len(res["terms"]["411@P3"]) == 2 and len(res["terms"]["433@P3"]) == 2 and len(res["terms"]["433@P2"]) == 1
    summ = {}
    for X in SURVEY:
        for fp in pops:
            for tag in ("LO", "CBH"):
                sel = [r for r in rows if r["level"] == name(X) and r["fp"] == fp and r["variant"] == tag and np.isfinite(r["lam_over_g322"])]
                pos = [r for r in sel if r["lam_over_g322"] > 0]
                d = {"n": len(sel), "max_lam_over_g322": max(r["lam_over_g322"] for r in sel) if sel else None,
                     "alpha_min_positive": min(r["alpha"] for r in pos) if pos else None, "n_positive": len(pos)}
                if fp == "P3":                                           # いまの判定(J)が覆っているか
                    d["n_positive_not_covered"] = int(sum(r["J_integrand_over_g322"] <= 0 for r in pos))
                    d["max_lam_over_J_integrand"] = max([r["lam_over_g322"] / r["J_integrand_over_g322"] for r in pos
                                                         if r["J_integrand_over_g322"] > 0] or [0.0])
                summ[f"{name(X)}@{fp}/{tag}"] = d
    res["summary"] = summ
    res["rows"] = rows
    # 結論を確かめる: P3 で育つ点はどれも J が覆い、増加率は J の増加率以下。P2 で正になるのは 522 だけ。四準位の固定点では何も育たない
    for X in SURVEY:
        for tag in ("LO", "CBH"):
            assert summ[f"{name(X)}@P3/{tag}"]["n_positive_not_covered"] == 0, (X, tag)
            assert summ[f"{name(X)}@P3/{tag}"]["max_lam_over_J_integrand"] <= 1.0 or tag == "LO", (X, tag)
            assert summ[f"{name(X)}@P2x/{tag}"]["n_positive"] == 0, (X, tag)
            assert summ[f"{name(X)}@P2/{tag}"]["n_positive"] == 0 or X == (5, 2, 2), (X, tag)
    return res


# ---------- x86 ----------
def part_x86():
    import hashlib
    import subprocess
    import e20_x86 as x20
    site = os.environ.get("E21_X86_SITE")
    if not site:
        return {"skipped": "E21_X86_SITE is not set"}
    site = os.path.abspath(site)
    work = os.environ.get("E21_WORK") or os.path.join(OUT, "work", "e21")
    os.makedirs(work, exist_ok=True)
    # 入力: 422 の規則が働く標本(export_chisr.py の確かめ 7 と同じ)、自己診断の入力、E19 の照合の入力(三つの種)
    sets = {}
    n = 100000
    rng = np.random.default_rng(21)
    a = rng.uniform(0.12, 0.30, n)
    M = np.exp(rng.uniform(np.log(1), np.log(300), n))
    mu = a / sr.alpha_of(M, 1.0)
    c0 = rng.uniform(0.0, 0.995, n)
    T = np.minimum(10 ** rng.uniform(5, np.log10(1.4e10), n), 0.12 / mu)
    f = 10 ** rng.uniform(8, 19, n)
    sets["rule422"] = (M, c0, T, mu, f)
    z = np.load(os.path.join(PKG, "selftest.npz"))
    sets["selftest"] = tuple(z[k] for k in ("M", "chi0", "T", "mu", "f"))
    for seed in SEEDS:
        v = np.load(os.path.join(OUT, f"e19_validate_{seed}.npz"))
        Mv, Tv = np.meshgrid(v["M"], v["T"], indexing="ij")
        b = [np.broadcast_to(v[k][:, None], Mv.shape).ravel() for k in ("chi0", "mu", "f")]
        sets[f"val_{seed}"] = (Mv.ravel(), b[0], Tv.ravel(), b[1], b[2])
    inp = os.path.join(work, "inputs.npz")
    np.savez(inp, **{f"{nm}/{k}": np.asarray(x, float) for nm, arr in sets.items() for k, x in zip(x20.ARGS, arr)})
    res = {"chi_sr_py_sha256": hashlib.sha256(open(os.path.join(PKG, "chi_sr.py"), "rb").read()).hexdigest(),
           "tables_sha256": hashlib.sha256(open(os.path.join(PKG, "chi_sr_tables.npz"), "rb").read()).hexdigest(),
           "n_inputs": {nm: int(len(arr[0])) for nm, arr in sets.items()}}
    archs = x20.arch_list(site)
    out = lambda arch: os.path.join(work, f"out_{arch}.npz")           # noqa: E731
    runs = {}
    for arch, py, st, env in archs:                                    # 警告は例外にする(strict)
        t0 = time.time()
        runs[arch] = x20.call(py, PKG, inp, out(arch), True, st, env)
        print(f"--- x86 part: {arch} ({time.time() - t0:.0f} s)", flush=True)
    res["env"] = {arch: runs[arch]["env"] for arch in runs}
    res["sets"] = {arch: runs[arch]["sets"] for arch in runs}
    names = sorted(sets)
    for x in ("x86_64", "x86_64_avx2"):
        cmp = x20.compare(out("arm64"), out(x), names)
        res[f"arm64_vs_{x}"] = cmp
        res[f"arm64_vs_{x}_total"] = x20.total(cmp, res["sets"][x], lambda nm: True)
    res["x86_64_vs_x86_64_avx2_total"] = x20.total(x20.compare(out("x86_64"), out("x86_64_avx2"), names), res["sets"]["x86_64"], lambda nm: True)
    # 規則が実際に働いている入力であること(arm64 の出力で、印 4 の数)
    za = np.load(out("arm64"))
    res["regime4_share_rule422"] = float(np.mean(za["rule422/regime"] == 4))
    # 自己診断
    res["selftest"] = {}
    for arch, py, st, env in archs:
        e = dict(env, PYTHONPATH=st) if st else env
        r = subprocess.run(list(py) + ["-s", "-B", os.path.join(PKG, "selftest.py")], capture_output=True, text=True, env=e)
        res["selftest"][arch] = {"returncode": r.returncode, "stdout": r.stdout.strip().splitlines(),
                                 "stderr_tail": r.stderr.strip().splitlines()[-3:]}
    for x in ("x86_64", "x86_64_avx2"):
        t = res[f"arm64_vs_{x}_total"]
        assert t["n_sets_raised"] == 0 and t["regime_differs"] == 0 and t["max_abs_dchi"] < 1e-9, t
        assert res["selftest"][x]["returncode"] == 0, res["selftest"][x]
    return res


def part_versions():
    import subprocess
    out = []
    for py in [x for x in os.environ.get("E21_PYTHONS", "").split(":") if x]:
        v = subprocess.run([py, "-s", "-c", "import sys, platform, numpy, scipy; print(sys.version.split()[0], "
                            "platform.machine(), numpy.__version__, scipy.__version__)"], capture_output=True, text=True,
                           check=True).stdout.split()
        r = subprocess.run([py, "-s", "-B", os.path.join(PKG, "selftest.py")], capture_output=True, text=True)
        out.append({"python": v[0], "machine": v[1], "numpy": v[2], "scipy": v[3], "returncode": r.returncode,
                    "stdout": r.stdout.strip().splitlines()})
        assert r.returncode == 0, out[-1]
    return {"selftest": out} if out else {"skipped": "E21_PYTHONS is not set"}


PARTS = {"rates": part_rates, "sr422": part_sr422, "growth": part_growth, "ode4": part_ode4, "flags": part_flags, "validate": part_validate,
         "survey": part_survey, "x86": part_x86, "versions": part_versions}


def main():
    names = sys.argv[1:] or list(PARTS)
    path = os.path.join(OUT, "e21.json")
    res = json.load(open(path)) if os.path.exists(path) else {}
    for nm in names:
        t0 = time.time()
        res[nm] = PARTS[nm]()
        print(f"--- {nm} ({time.time() - t0:.0f} s)", flush=True)
        show = res[nm]
        if nm == "rates":
            show = {k: v for k, v in show.items() if k not in ("bh", "inf")}
        elif nm == "sr422":
            show = {k: v for k, v in show.items() if k != "rows"}
        elif nm == "growth":
            show = {k: v for k, v in show.items() if k not in ("lam_table", "fp4")} | {"fp4_max": {v: res[nm]["fp4"][v]["max"] for v in res[nm]["fp4"]}}
        elif nm == "ode4":
            show = {k: v for k, v in show.items() if k != "cases"}
        elif nm == "survey":
            show = {k: v for k, v in show.items() if k not in ("rows", "terms")}
        elif nm == "x86":
            show = {k: v for k, v in show.items() if not k.startswith("arm64_vs_x86_64") or k.endswith("_total")}
        print(json.dumps(show, indent=1), flush=True)
        with open(path, "w") as fo:
            json.dump({k: res[k] for k in PARTS if k in res}, fo, indent=1)


if __name__ == "__main__":
    main()
