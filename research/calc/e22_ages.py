"""E22: BH の年齢の分布の模型(share/ages/ages.py)の確かめと、年齢について周辺化したときの窓。

  check  : 模型の部品の確かめ。宇宙の年齢(数値積分と比べる)、重みを外したときの解析的な分位点、重みつきの累積分布
           (別に書いた細かい数値積分と比べる)、標本の分布、赤方偏移についての補間。
  dist   : 既定の模型と、比べるための組(下の ONE、TWO)の年齢の分位点。図のデータ figs/data/e22_pdf.dat。
  window : 年齢について周辺化したときの、GW241011 の主星(2 世代目)が回転を保つ確率と、1 世代目が回転を失う確率を f の関数で
           出す(mu = 0.71e-12 eV、単体版の chi_SR)。年齢を一つの値に固定した E19 の数値を再現することも確かめる。
           窓が年齢の比だけで決まること(f によらない条件 T_1/T_2 >= R(M))と、年齢を k 倍・f を k^(-1/2) 倍にしても
           確率が変わらないこと(縮退)も見る。図のデータ figs/data/e22_prob_keep.dat、e22_prob_drop.dat。
出力: out/e22.json(部分ごとに書き足す)、../figs/data/e22_*.dat
実行: cd research && ~/miniforge3/envs/python313/bin/python calc/e22_ages.py [部分の名前 ...]   (名前なしで全部)
"""
import json
import os
import sys
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "..", "figs", "data")
PKG = os.path.join(HERE, "..", "..", "share", "chi_sr")
AGES = os.path.join(HERE, "..", "..", "share", "ages")
sys.path.insert(0, AGES)
sys.path.insert(0, PKG)
import ages  # noqa: E402

MU = 0.71e-12
Z1, Z2 = 0.2, 0.05            # 1 世代目の合体の赤方偏移の代表値、GW241011 の赤方偏移(D_L = 214 Mpc)
# 比べるための組(模型のパラメータ)。根拠は研究ノート E22 の表。どれも目安で、文献にある分布そのものではない
ONE = {"field": dict(),                                         # 既定: 遅れ ∝ 1/T、10 Myr 以上、星形成率の重み
       "old": dict(kappa1=0.0),                                 # 遅れが T について一様: 古い星団のように年齢が長いほうに寄る
       "flat": dict(sfr_weight=False)}                          # 星形成率の重みなし(対数一様)
TWO = {"wide": dict(),                                          # 既定: 対数一様、1 Myr から宇宙の年齢の上限まで
       "peak": dict(family2="lognormal", t2med=5e7, s2=0.7, t2lo=1e5),    # 星団の目安: 10--100 Myr に山
       "slow": dict(t2lo=3e8),                                  # 遅い側: 軽い 2 世代目の遅れの最小値の目安(3e8 yr)より上
       "agn": dict(t2lo=1e4, t2hi=1e7)}                         # AGN 円盤: 円盤の寿命(0.1--10 Myr)以内
QS = (0.05, 0.16, 0.5, 0.84, 0.95)


def chi_model():
    import chi_sr
    return chi_sr.ChiSR(os.path.join(PKG, "chi_sr_tables.npz"))


# ---------- check ----------
def _direct_cdf(am, gen, z, T, n=200001):
    """累積分布を、表を使わずに細かい格子で直接積分する(ages.py と別の書き方)。"""
    tm = float(am.cosmic_time(z))
    d = tm - am.t_start
    if gen == 1:
        lo, hi = am.tmin1, d
        w = lambda x: x ** (1 - am.kappa1) * am.sfr(tm - x)                                    # noqa: E731
    else:
        lo, hi = am.t2lo, min(am.t2hi, d - am.tmin1)

        def r1(t):                                              # 1 世代目の合体率を、その場で積分する
            if not am.parent_weight:
                return 1.0
            x = np.geomspace(am.tmin1, t - am.t_start, 20001)
            return np.trapezoid(x ** (1 - am.kappa1) * am.sfr(t - x), np.log(x))

        def w(x):
            return am._shape2(x) * np.array([r1(tm - xi) for xi in x])
    xs = np.geomspace(lo, hi, n if gen == 1 else 4001)
    ws = w(xs)
    c = np.concatenate([[0.0], np.cumsum(0.5 * (ws[1:] + ws[:-1]) * np.diff(np.log(xs)))])
    return float(np.interp(np.log(T), np.log(xs), c / c[-1]))


def part_check():
    from scipy.special import ndtri
    res = {}
    am = ages.AgeModel()
    res["describe"] = am.describe()
    # 1. 宇宙の年齢: 数値積分 int dz / ((1 + z) H(z)) と比べる
    worst = 0.0
    for z in (0.0, 0.05, 0.2, 1.0, 5.0, 20.0):
        zi = np.geomspace(1 + z, 1e8, 800001) - 1
        H = am._hub * np.sqrt(am.Om * (1 + zi) ** 3 + 1 - am.Om)
        worst = max(worst, abs(float(am.cosmic_time(z)) / float(np.trapezoid(1 / ((1 + zi) * H), zi)) - 1))
    res["cosmic_time"] = {"t0_Gyr": am.t0 / 1e9, "t_start_Gyr": am.t_start / 1e9, "t_z0.2_Gyr": float(am.cosmic_time(0.2)) / 1e9,
                          "t_z1_Gyr": float(am.cosmic_time(1.0)) / 1e9, "max_rel_diff_vs_numerical": worst,
                          "z_roundtrip_max": float(np.max(np.abs(am.z_at(am.cosmic_time(np.array([0.0, 0.3, 2.0, 15.0])))
                                                                 - np.array([0.0, 0.3, 2.0, 15.0]))))}
    assert worst < 1e-7 and res["cosmic_time"]["z_roundtrip_max"] < 1e-9
    # 星形成率の山(Madau & Dickinson: z = 1.5--2 のあいだ)
    zz = np.linspace(0, 8, 8001)
    res["sfr_peak_z"] = float(zz[np.argmax(am.sfr(am.cosmic_time(zz)))])
    assert 1.5 < res["sfr_peak_z"] < 2.0
    # 2. 重みを外したときの解析的な分位点
    p = np.linspace(0.002, 0.998, 499)
    an = {}
    for z in (0.05, 0.5, 2.0):
        a0 = ages.AgeModel(sfr_weight=False, parent_weight=False)
        Tm = float(a0.tmax(1, z))
        an[f"1G_loguniform_z{z}"] = float(np.max(np.abs(a0.quantile(1, z, p) / (a0.tmin1 ** (1 - p) * Tm ** p) - 1)))
        a1 = ages.AgeModel(kappa1=0.0, sfr_weight=False, parent_weight=False)
        an[f"1G_uniform_z{z}"] = float(np.max(np.abs(a1.quantile(1, z, p) / (a1.tmin1 + p * (Tm - a1.tmin1)) - 1)))
        b = float(a0.tmax(2, z))
        an[f"2G_loguniform_z{z}"] = float(np.max(np.abs(a0.quantile(2, z, p) / (a0.t2lo ** (1 - p) * b ** p) - 1)))
        a2 = ages.AgeModel(family2="lognormal", t2med=3e7, s2=0.5, t2lo=1e2, t2hi=1e10, sfr_weight=False, parent_weight=False)
        an[f"2G_lognormal_z{z}"] = float(np.max(np.abs(a2.quantile(2, z, p) / (3e7 * 10 ** (0.5 * ndtri(p))) - 1)))
    res["analytic_quantiles_max_rel_diff"] = an
    assert max(an.values()) < 5e-3, an
    # 3. 重みつきの累積分布: 直接の数値積分と比べる。往復(分位点 -> 累積分布)と単調性も
    cmpd = {}
    for name, kw in [("default", {}), ("old_peak", dict(kappa1=0.0, family2="lognormal", t2med=5e7, s2=0.7, t2lo=1e5)),
                     ("slow", dict(t2lo=3e8)), ("agn", dict(t2lo=1e4, t2hi=1e7))]:
        m = ages.AgeModel(**kw)
        worst = {1: 0.0, 2: 0.0}
        rt = 0.0
        for gen in (1, 2):
            for z in (0.03, 0.2, 0.47, 1.1, 2.3):
                T = m.quantile(gen, z, np.array([0.03, 0.2, 0.5, 0.8, 0.97]))
                for Ti, pi in zip(T, (0.03, 0.2, 0.5, 0.8, 0.97)):
                    worst[gen] = max(worst[gen], abs(_direct_cdf(m, gen, z, Ti) - pi))
                Tq = m.quantile(gen, z, p)
                assert np.all(np.diff(Tq) > 0)
                rt = max(rt, float(np.max(np.abs(m.cdf(gen, z, Tq) - p))))
        cmpd[name] = {"cdf_max_abs_diff_1G": worst[1], "cdf_max_abs_diff_2G": worst[2], "roundtrip_max": rt}
    res["direct_integration"] = cmpd
    assert max(max(v["cdf_max_abs_diff_1G"], v["cdf_max_abs_diff_2G"], v["roundtrip_max"]) for v in cmpd.values()) < 3e-3, cmpd
    # 4. 標本: 経験分布と累積分布の最大の差(KS 距離)。密度の積分が 1
    rng = np.random.default_rng(22)
    ks = {}
    for gen in (1, 2):
        s = np.sort(am.sample(gen, np.full(400000, 0.3), rng))
        ks[str(gen)] = float(np.max(np.abs(am.cdf(gen, 0.3, s) - (np.arange(len(s)) + 0.5) / len(s))))
        x = np.geomspace(1e5, 2e10, 400001)
        ks[f"pdf_integral_{gen}"] = float(np.trapezoid(am.pdf_lnT(gen, 0.3, x), np.log(x)))
    res["sampling_ks"] = ks
    assert ks["1"] < 4e-3 and ks["2"] < 4e-3 and abs(ks["pdf_integral_1"] - 1) < 2e-3 and abs(ks["pdf_integral_2"] - 1) < 2e-3, ks
    # 5. 節の平均が期待値になること: E[ln T] を、節(K = 32、400)と密度の積分で比べる
    nd = {}
    for gen in (1, 2):
        x = np.geomspace(1e5, 2e10, 400001)
        ex = float(np.trapezoid(am.pdf_lnT(gen, 0.3, x) * np.log(x), np.log(x)))
        nd[str(gen)] = {str(K): float(np.mean(np.log(am.nodes(gen, 0.3, K))) - ex) for K in (32, 400)}
    res["nodes_mean_lnT_minus_integral"] = nd
    assert max(abs(v) for d_ in nd.values() for v in d_.values()) < 0.02, nd
    # 6. 範囲の外: 合体が早すぎて年齢の幅がない赤方偏移では NaN
    res["nan_when_no_room"] = bool(np.isnan(am.quantile(2, 19.9, 0.5)) and np.isnan(am.quantile(1, 25.0, 0.5)))
    assert res["nan_when_no_room"]
    return res


# ---------- dist ----------
def part_dist():
    res = {"z_1G": Z1, "z_2G": Z2, "quantiles": list(QS)}
    rows = {}
    for nm, kw in ONE.items():
        m = ages.AgeModel(**kw)
        rows[nm] = {"params": kw, "by_z": {str(z): [float(x) for x in m.quantile(1, z, np.array(QS))] for z in (0.05, 0.2, 0.5, 1.0)},
                    "frac_younger": {f"{t:g}": float(m.cdf(1, Z1, t)) for t in (3e7, 1e8, 1e9)},
                    "tmax_yr": float(m.tmax(1, Z1))}
    res["1G"] = rows
    rows = {}
    for nm, kw in TWO.items():
        m = ages.AgeModel(**kw)
        rows[nm] = {"params": kw, "by_z": {str(z): [float(x) for x in m.quantile(2, z, np.array(QS))] for z in (0.05, 0.2, 0.5, 1.0)},
                    "frac_younger": {f"{t:g}": float(m.cdf(2, Z2, t)) for t in (5e7, 3e8, 1e9)},
                    "tmax_yr": float(m.tmax(2, Z2))}
    res["2G"] = rows
    # 親の合体率の重みと、星形成率の重みの効き方(既定の模型で、重みを外した場合と比べる)
    m0 = ages.AgeModel(sfr_weight=False, parent_weight=False)
    m1 = ages.AgeModel()
    res["weights_effect_median"] = {"1G_with": float(m1.quantile(1, Z1, 0.5)), "1G_without": float(m0.quantile(1, Z1, 0.5)),
                                    "2G_with": float(m1.quantile(2, Z2, 0.5)), "2G_without": float(m0.quantile(2, Z2, 0.5))}
    # 星形成の始まり z_start への感度
    ma = ages.AgeModel(z_start=10.0)
    res["z_start_10_vs_20_median_ratio"] = {"1G": float(ma.quantile(1, Z1, 0.5) / m1.quantile(1, Z1, 0.5)),
                                            "2G": float(ma.quantile(2, Z2, 0.5) / m1.quantile(2, Z2, 0.5))}
    # 図のデータ: d P / d log10 T
    x = np.geomspace(1e4, 2e10, 400)
    cols = {"T": x}
    for nm, kw in ONE.items():
        cols[f"g1_{nm}"] = ages.AgeModel(**kw).pdf_lnT(1, Z1, x) * np.log(10)
    for nm, kw in TWO.items():
        cols[f"g2_{nm}"] = ages.AgeModel(**kw).pdf_lnT(2, Z2, x) * np.log(10)
    with open(os.path.join(FIG, "e22_pdf.dat"), "w") as fo:
        fo.write(" ".join(cols) + "\n")
        for i in range(len(x)):
            fo.write(" ".join(f"{cols[k][i]:.5e}" for k in cols) + "\n")
    return res


# ---------- window ----------
LF = np.linspace(10.5, 14.0, 141)
KEEP = 0.69                    # GW241011 の主星が下回らないとするスピン(E10 から同じ)
M2G = 19.6
MASSES = (10.0, 15.0, 20.0, 25.0)   # 30 Msun(alpha = 0.16)は、(211, 322) の固定点で 411・433・544 が育ちうる印 4 になるので使わない
DELTA = 0.05                   # 1 世代目が「回転を失った」とみなす、211 の臨界からの幅(説明用。集団解析では使わない。D9)
K_AGE = 400


def _spins(model, M, chi0, T, f):
    """年齢の節 T(K,)と f(F,)の格子での chi_SR と印。形は (F, K)。"""
    c, r = model(M, chi0, T[None, :], MU, f[:, None])
    return c, r


def part_window():
    from scipy.optimize import brentq
    model = chi_model()
    f = 10 ** LF
    res = {"mu_eV": MU, "z_1G": Z1, "z_2G": Z2, "keep_threshold": KEEP, "delta": DELTA, "K_ages": K_AGE, "log10_f": [LF[0], LF[-1], len(LF)]}
    # 0. 年齢を一つの値に固定した E19 の数値を再現する(同じ関数、同じ定義)
    e19 = {(r["T2_yr"], r["chi0"]): r for r in json.load(open(os.path.join(OUT, "e19.json")))["window"]["gw241011_primary"]}
    fixed = []
    for (T2, c0), ref in e19.items():
        lf = brentq(lambda x: float(model(M2G, c0, T2, MU, 10 ** x)[0]) - KEEP, 10, 14, xtol=1e-6)
        fixed.append({"T2_yr": T2, "chi0": c0, "f_max_GeV": float(10 ** lf), "ratio_to_e19": float(10 ** lf / ref["f_max_GeV"])})
    res["fixed_age_f_max"] = fixed
    assert all(abs(r["ratio_to_e19"] - 1) < 1e-4 for r in fixed), fixed
    # 1. 年齢の比で書いた窓の条件: 1 世代目(質量 M、chi0 = 0.7)が年齢 T1 のうちに「211 の臨界 + delta」まで落ちる f の下限
    #    f_low(M) ∝ T1^(-1/2) と、GW241011 の主星が年齢 T2 のあいだ KEEP を保つ f の上限 f_max ∝ T2^(-1/2) から、
    #    窓があるのは T1 / T2 >= R(M) = (f_low sqrt(T1))^2 / (f_max sqrt(T2))^2 のとき(f によらない)
    def f_low(M, T1, delta):
        c1 = float(model.chic(model.alpha(M, MU), 1))
        return 10 ** brentq(lambda x: float(model(M, 0.7, T1, MU, 10 ** x)[0]) - (c1 + delta), 9.5, 14.5, xtol=1e-6)

    def f_max(c0, T2):
        return 10 ** brentq(lambda x: float(model(M2G, c0, T2, MU, 10 ** x)[0]) - KEEP, 10, 14, xtol=1e-6)

    ratio = {}
    for c0 in (0.72, 0.78):
        B = {T2: f_max(c0, T2) * np.sqrt(T2) for T2 in (3e7, 1e8, 3e8, 1e9)}
        row = {"f_max_sqrtT_spread": float(max(B.values()) / min(B.values()) - 1)}
        for delta in (0.02, 0.05):
            for M in (8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0):
                A = {T1: f_low(M, T1, delta) * np.sqrt(T1) for T1 in (1e9, 6e9)}
                row[f"R_delta{delta}_M{M:g}"] = float((A[6e9] / B[1e8]) ** 2)
                row[f"f_low_sqrtT_spread_delta{delta}_M{M:g}"] = float(A[1e9] / A[6e9] - 1)
        ratio[str(c0)] = row
    res["age_ratio_condition"] = ratio
    # E19 の表(T2 = 5e7 yr、chi0 = 0.78: 6 Gyr の 1 世代目の最小質量 11.9 Msun(delta = 0.05))と同じ内容であること:
    # T1 / T2 = 120 のとき R(M) = 120 になる質量
    def m_at(c0, delta, R):
        g = lambda M: (f_low(M, 6e9, delta) * np.sqrt(6e9) / (f_max(c0, 1e8) * np.sqrt(1e8))) ** 2 - R             # noqa: E731
        return float(brentq(g, 6, 30, xtol=1e-4))
    chk = {f"chi0{c0}_T2{T2:g}_delta{d_}": {"M_min_from_ratio": m_at(c0, d_, 6e9 / T2), "M_min_e19": e19[(T2, c0)][f"M_min_1G_delta{d_}"]}
           for c0 in (0.78, 0.72) for T2 in (5e7, 3e8) for d_ in (0.02, 0.05)}
    res["age_ratio_vs_e19_M_min"] = chk
    assert all(abs(v["M_min_from_ratio"] / v["M_min_e19"] - 1) < 3e-3 for v in chk.values()), chk
    # 年齢を中央値に固定して読んだ最小質量(delta = 0.05): 1 世代目と 2 世代目の年齢の中央値の比が R(M) に等しくなる質量
    med = {}
    for n1 in ("field", "old"):
        for n2 in ("wide", "peak", "slow"):
            T1 = float(ages.AgeModel(**ONE[n1]).quantile(1, Z1, 0.5))
            T2 = float(ages.AgeModel(**TWO[n2]).quantile(2, Z2, 0.5))
            row = {"T1_median_yr": T1, "T2_median_yr": T2, "ratio": T1 / T2}
            for c0 in (0.72, 0.78):
                try:
                    row[f"M_min_chi0{c0}"] = m_at(c0, DELTA, T1 / T2)
                except ValueError:                              # 6--30 Msun の外
                    row[f"M_min_chi0{c0}"] = None
            med[f"{n1}|{n2}"] = row
    res["M_min_at_median_ages"] = med
    # 2. 年齢について周辺化した確率
    keep, drop, band, flagged = {}, {}, {}, {}
    for nm, kw in TWO.items():
        T = ages.AgeModel(**kw).nodes(2, Z2, K_AGE)
        for c0 in (0.72, 0.78):
            c, r = _spins(model, M2G, c0, T, f)
            keep[f"{nm}_{c0}"] = np.mean(c >= KEEP, axis=1)
            flagged[f"2G_{nm}_{c0}"] = float(np.mean(r >= 3))
            band[f"2G_{nm}_{c0}"] = np.percentile(c, [16, 50, 84], axis=1)
    for nm, kw in ONE.items():
        T = ages.AgeModel(**kw).nodes(1, Z1, K_AGE)
        for M in MASSES:
            c1 = float(model.chic(model.alpha(M, MU), 1))
            c, r = _spins(model, M, 0.7, T, f)
            drop[f"{nm}_{M:g}"] = np.mean(c <= c1 + DELTA, axis=1)
            flagged[f"1G_{nm}_{M:g}"] = float(np.mean(r >= 3))
            band[f"1G_{nm}_{M:g}"] = np.percentile(c, [16, 50, 84], axis=1)
    res["flagged_share_max"] = max(flagged.values())
    res["flagged_share"] = flagged

    def f_at(pcurve, level, rising):
        """確率が level を横切る f(対数で線形補間)。横切らなければ None。"""
        y = pcurve - level
        idx = np.flatnonzero((y[:-1] < 0) & (y[1:] >= 0)) if rising else np.flatnonzero((y[:-1] >= 0) & (y[1:] < 0))
        if len(idx) == 0:
            return None
        i = idx[0]
        return float(10 ** (LF[i] + (LF[i + 1] - LF[i]) * (0 - y[i]) / (y[i + 1] - y[i])))

    res["f_keep"] = {k: {str(lv): f_at(v, lv, False) for lv in (0.9, 0.5, 0.1)} for k, v in keep.items()}
    res["f_drop"] = {k: {str(lv): f_at(v, lv, True) for lv in (0.1, 0.5, 0.9)} for k, v in drop.items()}
    # 3. 同じ f で両方を満たす度合い: W = max_f min(P_keep, P_drop)。年齢を固定すれば窓があるとき 1、ないとき 0
    W = {}
    for k2, pk in keep.items():
        for k1, pd in drop.items():
            m = np.minimum(pk, pd)
            i = int(np.argmax(m))
            both = np.flatnonzero((pk >= 0.5) & (pd >= 0.5))
            W[f"{k2}|{k1}"] = {"W": float(m[i]), "f_GeV": float(f[i]),
                               "f_range_both_ge_0.5": [float(f[both[0]]), float(f[both[-1]])] if len(both) else None}
    res["W"] = W
    # 4. 縮退: 2 世代目の年齢を 10 倍にして f を 10^(-1/2) 倍にすると、P_keep は同じ(親の合体率の重みは外す。上限は宇宙の年齢に
    #    かからないようにとる)
    a1 = ages.AgeModel(t2lo=1e5, t2hi=1e8, parent_weight=False)
    a2 = ages.AgeModel(t2lo=1e6, t2hi=1e9, parent_weight=False)
    deg = {}
    for c0 in (0.72, 0.78):
        p1 = np.mean(_spins(model, M2G, c0, a1.nodes(2, Z2, K_AGE), f)[0] >= KEEP, axis=1)
        p2 = np.mean(_spins(model, M2G, c0, a2.nodes(2, Z2, K_AGE), f / np.sqrt(10.0))[0] >= KEEP, axis=1)
        deg[str(c0)] = float(np.max(np.abs(p1 - p2)))
    res["degeneracy_max_abs_dP"] = deg
    assert max(deg.values()) < 0.01, deg
    # 5. 節の数: K = 32 の節で出した確率と、K = 400 の確率の差(確率は 1/K 刻みでしか出ないので、差は 1/K の程度)
    kc = {}
    for nm, kw, gen, M, c0, z in (("2G_wide_0.78", TWO["wide"], 2, M2G, 0.78, Z2), ("1G_field_20", ONE["field"], 1, 20.0, 0.7, Z1)):
        am = ages.AgeModel(**kw)
        thr = KEEP if gen == 2 else float(model.chic(model.alpha(M, MU), 1)) + DELTA
        out = {}
        for K in (32, 100):
            c = _spins(model, M, c0, am.nodes(gen, z, K), f)[0]
            pk = np.mean(c >= thr, axis=1) if gen == 2 else np.mean(c <= thr, axis=1)
            ref = keep["wide_0.78"] if gen == 2 else drop["field_20"]
            out[str(K)] = float(np.max(np.abs(pk - ref)))
        kc[nm] = out
    res["nodes_K_vs_400_max_abs_dP"] = kc
    # 図のデータ
    with open(os.path.join(FIG, "e22_prob_keep.dat"), "w") as fo:
        ks = list(keep)
        fo.write("f " + " ".join(k.replace(".", "") for k in ks) + "\n")
        for i in range(len(f)):
            fo.write(f"{f[i]:.5e} " + " ".join(f"{keep[k][i]:.5f}" for k in ks) + "\n")
    with open(os.path.join(FIG, "e22_prob_drop.dat"), "w") as fo:
        ks = list(drop)
        fo.write("f " + " ".join(k for k in ks) + "\n")
        for i in range(len(f)):
            fo.write(f"{f[i]:.5e} " + " ".join(f"{drop[k][i]:.5f}" for k in ks) + "\n")
    # f = 1e12 GeV でのスピンの分布(16%、50%、84%)。年齢の散らばりがそのままスピンの散らばりになる例
    ks = ["2G_wide_0.78", "2G_peak_0.78", "1G_field_20", "1G_old_20"]
    res["band_at_f1e12"] = {k: [float(band[k][j][int(np.argmin(np.abs(LF - 12.0)))]) for j in range(3)] for k in ks}
    return res


PARTS = {"check": part_check, "dist": part_dist, "window": part_window}


def main():
    names = sys.argv[1:] or list(PARTS)
    path = os.path.join(OUT, "e22.json")
    res = json.load(open(path)) if os.path.exists(path) else {}
    for nm in names:
        t0 = time.time()
        res[nm] = PARTS[nm]()
        print(f"--- {nm} ({time.time() - t0:.0f} s)", flush=True)
        print(json.dumps(res[nm], indent=1), flush=True)
        with open(path, "w") as fo:
            json.dump({k: res[k] for k in PARTS if k in res}, fo, indent=1)


if __name__ == "__main__":
    main()
