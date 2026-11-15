"""E23: 2 世代目の BH が生まれたときのスピンと質量(share/remnant/remnant.py)の確かめと、その分布。

  check  : 残骸の式の確かめ。スピンの式は E10 の calc/hbr16.py と同じ値、決まった入力での値(等質量・無スピンの 0.68646、
           放出エネルギー、キックの最大、スーパーキック)、merge の標本の向きが等方であること、質量の入れ替えについての対称性。
  dist   : 2 世代目の生まれたときのスピンの分布を、親の質量比・親のスピン・星団に残る条件(キック < 脱出速度)を変えて出す。
           文献の数値実験の値(球状星団の模擬計算 1906.10260)と、同じ条件で比べる。図のデータ figs/data/e23_hist.dat。
  boson  : 親(1 世代目)が合体までに超放射で回転を失っている場合。親の生まれたときのスピンと年齢(E22 の模型)から chi_SR で
           合体のときのスピンを出し、2 世代目の生まれたときのスピンの分布が f でどう変わるかを見る(mu = 0.71e-12 eV)。
           図のデータ figs/data/e23_hist_boson.dat。
  gw241011: GW241011 の主星の質量の 2 世代目が、合体のときに 0.69 以上のスピンを持つ確率を、生まれたときのスピンの分布と
           年齢の分布(E22)の両方について周辺化して f の関数で出す。生まれたスピンを一つの値に固定した場合(E22)と比べ、
           同じ f を与える「実効的な生まれたスピン」を出す。図のデータ figs/data/e23_keep.dat。
  window : E22 の「同じ f で両立する度合い」W を、gw241011 の曲線で出し直す(E22 の図のデータ e22_prob_drop.dat を使う。
           gw241011 のあとに走らせる)。
出力: out/e23.json(部分ごとに書き足す)、../figs/data/e23_*.dat
実行: cd research && ~/miniforge3/envs/python313/bin/python calc/e23_natal2g.py [部分の名前 ...]   (名前なしで全部。6 並列で約 4 分)
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
SHARE = os.path.join(HERE, "..", "..", "share")
for d in ("remnant", "ages", "chi_sr"):
    sys.path.insert(0, os.path.join(SHARE, d))
sys.path.insert(0, HERE)
import ages  # noqa: E402
import remnant  # noqa: E402

MU = 0.71e-12
M2G = 19.6                                   # GW241011 の主星の質量 [Msun]
MTOT = M2G / 0.952                           # 親の質量の和の目安(等質量・無スピンで残骸が 0.952 倍になる)
BETA = (0.69, 2.28)                          # 1 世代目の生まれたときのスピン: arXiv:2607.01317 の事後分布の中央値(E12 と同じ)
Z_PARENT = 0.3                               # 親の合体の赤方偏移の代表値(1 世代目の年齢の上限を決めるだけ)
Z2 = 0.05                                    # GW241011 の赤方偏移
KEEP = 0.69                                  # GW241011 の主星のスピンの測定の 90% 区間(0.69--0.87)の下端
THR = (0.69, 0.72, 0.78, 0.87)
N = 400000
ALPHA_PER_MSUN = 4.925490947e-6 * MU / 6.582119569e-16      # alpha = G M mu/(hbar c^3) の、1 Msun あたりの値
HIST_EDGES = np.linspace(0.40, 0.95, 111)
_MODEL = {}


def chi_model():
    """chi_SR の単体版(プロセスごとに一度だけ読む)。"""
    if "m" not in _MODEL:
        import chi_sr
        _MODEL["m"] = chi_sr.ChiSR(os.path.join(SHARE, "chi_sr", "chi_sr_tables.npz"))
    return _MODEL["m"]


def masses(rng, n, qlo=0.5, mtot=MTOT):
    """親の質量: 和が mtot、質量比が [qlo, 1] で一様(qlo = 1 なら等質量)。"""
    q = np.ones(n) if qlo >= 1 else rng.uniform(qlo, 1.0, n)
    m1 = mtot / (1 + q)
    return m1, q * m1


def spins(rng, n, kind):
    if kind == "zero":
        return np.zeros(n), np.zeros(n)
    if kind == "beta":
        return rng.beta(*BETA, n), rng.beta(*BETA, n)
    if kind == "beta25":
        return rng.beta(2, 5, n), rng.beta(2, 5, n)
    return np.full(n, float(kind)), np.full(n, float(kind))


def summ(chi, M=None, keep=None):
    """分布のまとめ。keep を渡すと、残ったものだけ。"""
    sel = np.ones(len(chi), bool) if keep is None else keep
    c = chi[sel]
    out = {"n": int(sel.sum()), "retained": float(np.mean(sel))}
    if len(c) < 200:
        return out
    out.update({"median": float(np.median(c)), "p05": float(np.percentile(c, 5)), "p16": float(np.percentile(c, 16)),
                "p84": float(np.percentile(c, 84)), "p95": float(np.percentile(c, 95)), "mean": float(np.mean(c)), "max": float(np.max(c))})
    out.update({f"P(>={t})": float(np.mean(c >= t)) for t in THR})
    if M is not None:
        out["M_over_Mtot"] = [float(np.percentile(M[sel], x)) for x in (5, 50, 95)]
    return out


# ---------- check ----------
def part_check():
    import hbr16
    res = {}
    bad = remnant.selfcheck()
    res["selfcheck_failures"] = len(bad)
    assert not bad, bad
    rng = np.random.default_rng(23)
    n = 400000
    q = rng.uniform(0.02, 1, n)
    a1, a2 = rng.uniform(0, 1, n), rng.uniform(0, 1, n)
    cb, cg, dp = rng.uniform(-1, 1, n), rng.uniform(-1, 1, n), rng.uniform(0, 2 * np.pi, n)
    ca = cb * cg + np.sqrt(1 - cb ** 2) * np.sqrt(1 - cg ** 2) * np.cos(dp)
    res["final_spin_vs_hbr16_max_abs_diff"] = float(np.max(np.abs(remnant.final_spin(q, a1, a2, cb, cg, ca) - hbr16.a_final(q, a1, a2, cb, cg, ca))))
    assert res["final_spin_vs_hbr16_max_abs_diff"] == 0.0
    af = remnant.final_spin(q, a1, a2, cb, cg, ca)
    # 当てはめは、親のスピンが 0.999 以上でほぼ揃い質量比が小さいときだけ、1 をわずかに越える(merge は 1 で切る)
    res["final_spin_range"] = [float(af.min()), float(af.max())]
    res["final_spin_gt_1"] = {"n": int(np.sum(af > 1)), "of": n, "a1_min": float(a1[af > 1].min()) if np.any(af > 1) else None}
    assert np.all(np.isfinite(af)) and af.min() >= 0 and af.max() < 1.002
    assert not np.any(af > 1) or res["final_spin_gt_1"]["a1_min"] > 0.99
    # 1 を越える大きさ: 重いほうが最大のスピンで L に揃った配置を、質量比と軽いほうのスピン(L に平行か反平行)について探す
    qg, a2g = np.meshgrid(np.linspace(0.001, 1, 4000), np.linspace(-1, 1, 401), indexing="ij")

    def top(a1):
        return remnant.final_spin(qg, a1, np.abs(a2g), 1.0, np.sign(a2g) + (a2g == 0), np.sign(a2g) + (a2g == 0))

    t1 = top(1.0)
    i = np.unravel_index(t1.argmax(), t1.shape)
    res["final_spin_max_aligned_a1=1"] = {"value": float(t1.max()), "q": float(qg[i]), "a2_signed": float(a2g[i]),
                                          "value_a1=0.999": float(top(0.999).max()), "value_a1=0.99": float(top(0.99).max())}
    assert 1.0 < res["final_spin_max_aligned_a1=1"]["value"] < 1.002 and res["final_spin_max_aligned_a1=1"]["value_a1=0.99"] <= 1.0
    assert res["final_spin_range"][1] <= res["final_spin_max_aligned_a1=1"]["value"] + 1e-4
    res["spin_values"] = {"q1_nonspin": float(remnant.final_spin(1, 0, 0, 1, 1, 1)), "q1/3_nonspin": float(remnant.final_spin(1 / 3, 0, 0, 1, 1, 1)),
                          "q1_aligned_0.5": float(remnant.final_spin(1, 0.5, 0.5, 1, 1, 1)),
                          "q->0_returns_a1": float(remnant.final_spin(1e-6, 0.37, 0.9, 0.2, -0.5, 0.1)), "k00": float(remnant.K00)}
    assert abs(res["spin_values"]["q1_nonspin"] - 0.68646) < 1e-12 and abs(res["spin_values"]["q->0_returns_a1"] - 0.37) < 1e-5
    # 角度の補正(1605.01938 の式 18)の大きさ
    d = remnant.final_spin(q, a1, a2, cb, cg, ca, correct=True) - af
    res["angle_correction"] = {"mean": float(d.mean()), "min": float(d.min()), "max": float(d.max())}
    # 質量: 等質量・無スピンで 1 - p0。等質量・揃った最大スピンで放出 9.95%(1206.3803 Sec. 2 の値)。試験粒子の極限
    res["mass_values"] = {"q1_nonspin": float(remnant.final_mass_fraction(1, 0, 0, 1, 1)),
                          "Erad_q1_aligned_a1": float(1 - remnant.final_mass_fraction(1, 1, 1, 1, 1)),
                          "Erad_over_nu_q->0": float((1 - remnant.final_mass_fraction(1e-5, 0, 0, 1, 1)) / (1e-5 / (1 + 1e-5) ** 2)),
                          "1-E_ISCO(0)": float(1 - np.sqrt(8.0 / 9.0))}
    assert abs(res["mass_values"]["q1_nonspin"] - (1 - remnant.P0)) < 1e-15 and abs(res["mass_values"]["Erad_q1_aligned_a1"] - 0.0995) < 1e-4
    assert abs(res["mass_values"]["Erad_over_nu_q->0"] / res["mass_values"]["1-E_ISCO(0)"] - 1) < 1e-3
    mf = remnant.final_mass_fraction(q, a1, a2, cb, cg)
    res["mass_fraction_range"] = [float(mf.min()), float(mf.max())]
    assert mf.min() > 0.88 and mf.max() < 1.0
    # キック: 無スピンの最大(約 175 km/s、q = 0.36。1906.10260 Sec. IV D の値)、等質量・無スピンで 0、スーパーキック V11、全体の最大。
    # Delta の符号: 重いほうが順行(L と同じ向き)ならキックは無スピンより小さく、逆行なら大きい
    qq = np.linspace(0.01, 1, 9901)
    vm = remnant.kick(qq, 0, 0, 1, 1, 0, 1)
    vmax = remnant.kick(rng.uniform(0.5, 1, 4 * n), 1.0, 1.0, rng.uniform(-1, 1, 4 * n), rng.uniform(-1, 1, 4 * n), rng.uniform(0, 2 * np.pi, 4 * n), 1.0)
    res["kick_values"] = {"nonspin_max_kms": float(vm.max()), "nonspin_max_at_q": float(qq[vm.argmax()]), "q1_nonspin": float(remnant.kick(1, 0, 0, 1, 1, 0, 1)),
                          "superkick_kms": float(remnant.kick(1, 1, 1, 0, 0, np.pi, 1)), "max_over_random_configs_kms": float(vmax.max()),
                          "q0.5_nonspin": float(remnant.kick(0.5, 0, 0, 1, 1, 0, 0)), "q0.5_heavier_0.5_prograde": float(remnant.kick(0.5, 0.5, 0, 1, 1, 0, 0)),
                          "q0.5_heavier_0.5_retrograde": float(remnant.kick(0.5, 0.5, 0, -1, 1, 0, 0))}
    assert abs(res["kick_values"]["nonspin_max_kms"] - 175.2) < 0.1 and abs(res["kick_values"]["superkick_kms"] - remnant.V11) < 1e-9
    assert 4500 < res["kick_values"]["max_over_random_configs_kms"] < 5100 and res["kick_values"]["q1_nonspin"] == 0
    assert res["kick_values"]["q0.5_heavier_0.5_prograde"] < res["kick_values"]["q0.5_nonspin"] < res["kick_values"]["q0.5_heavier_0.5_retrograde"]
    # merge: 向きが等方、質量の入れ替えについて対称(重いほうを 1 にそろえる)、速さ
    m_a, m_b = rng.uniform(5, 30, n), rng.uniform(5, 30, n)
    c_a, c_b = rng.uniform(0, 1, n), rng.uniform(0, 1, n)
    o1 = remnant.merge(m_a, m_b, c_a, c_b, np.random.default_rng(7))
    o2 = remnant.merge(m_b, m_a, c_b, c_a, np.random.default_rng(7))
    res["merge"] = {"swap_max_abs_diff_chi": float(np.max(np.abs(o1["chi"] - o2["chi"]))), "swap_max_abs_diff_M": float(np.max(np.abs(o1["M"] - o2["M"]))),
                    "mean_cosb": float(o1["cosb"].mean()), "mean_cosg": float(o1["cosg"].mean()), "q_max": float(o1["q"].max())}
    assert res["merge"]["swap_max_abs_diff_chi"] == 0 and res["merge"]["swap_max_abs_diff_M"] == 0 and res["merge"]["q_max"] <= 1
    assert abs(res["merge"]["mean_cosb"]) < 5e-3 and abs(res["merge"]["mean_cosg"]) < 5e-3
    big = 1000000
    args = (rng.uniform(5, 30, big), rng.uniform(5, 30, big), rng.uniform(0, 1, big), rng.uniform(0, 1, big))
    t0 = time.time()
    remnant.merge(*args, rng)
    res["merge"]["seconds_per_1e6"] = round(time.time() - t0, 2)
    return res


# ---------- dist ----------
VESC = (None, 50.0, 100.0, 300.0)


def kick_delta_sign(q, a1, a2, cosb, cosg, dphi, cos_theta, sgn):
    """remnant.kick と同じ式で、Delta の符号だけを sgn 倍にしたもの(sgn = +1 で remnant.kick と同じ値。part_dist で確かめる)。
    sgn = -1 は Delta = (q a2 S2 - a1 S1)/(1+q) で、1605.01067 の公開コードが v1.0.2 で直す前の符号にあたる。"""
    nu = q / (1 + q) ** 2
    sb, sg = np.sqrt(1 - cosb ** 2), np.sqrt(1 - cosg ** 2)
    d_par = sgn * (a1 * cosb - q * a2 * cosg) / (1 + q)
    d_perp = np.sqrt(np.maximum(a1 ** 2 * sb ** 2 + q ** 2 * a2 ** 2 * sg ** 2 - 2 * q * a1 * a2 * sb * sg * np.cos(dphi), 0.0)) / (1 + q)
    c_par = (a1 * cosb + q ** 2 * a2 * cosg) / (1 + q) ** 2
    c_perp = np.sqrt(np.maximum(a1 ** 2 * sb ** 2 + q ** 4 * a2 ** 2 * sg ** 2 + 2 * q ** 2 * a1 * a2 * sb * sg * np.cos(dphi), 0.0)) / (1 + q) ** 2
    vm = remnant.KICK_A * nu ** 2 * (1 - q) / (1 + q) * (1 + remnant.KICK_B * nu)
    vperp = remnant.KICK_H * nu ** 2 * d_par
    vpar = 16 * nu ** 2 * (d_perp * (remnant.V11 + 2 * remnant.VA * c_par + 4 * remnant.VB * c_par ** 2 + 8 * remnant.VC * c_par ** 3)
                           + 2 * c_perp * d_par * (remnant.C2 + 2 * remnant.C3 * c_par)) * cos_theta
    return np.sqrt(vm ** 2 + 2 * vm * vperp * np.cos(remnant.KICK_ZETA) + vperp ** 2 + vpar ** 2)


def part_dist():
    res = {"M2G": M2G, "Mtot": MTOT, "beta": list(BETA), "N": N}
    rng = np.random.default_rng(231)
    cols = {"chi": 0.5 * (HIST_EDGES[1:] + HIST_EDGES[:-1])}
    for sname in ("zero", "beta", "beta25", "0.2", "0.5"):
        for qlo in (1.0, 0.5, 0.25):
            m1, m2 = masses(rng, N, qlo)
            a1, a2 = spins(rng, N, sname)
            o = remnant.merge(m1, m2, a1, a2, rng)
            row = {}
            for v in VESC:
                keep = None if v is None else o["vkick"] < v
                row["all" if v is None else f"vesc{v:g}"] = summ(o["chi"], o["M"] / MTOT, keep)
            row["vkick_median_kms"] = float(np.median(o["vkick"]))
            row["vkick_p05_p95_kms"] = [float(np.percentile(o["vkick"], 5)), float(np.percentile(o["vkick"], 95))]
            res[f"{sname}_q{qlo:g}"] = row
            if qlo == 0.5 and sname in ("zero", "beta", "0.5"):
                cols[f"{sname.replace('.', '')}_all"] = np.histogram(o["chi"], bins=HIST_EDGES, density=True)[0]
                k = o["vkick"] < 100.0
                # 残ったものの分布は、全体に対する密度(面積 = 残る割合)で書く
                cols[f"{sname.replace('.', '')}_v100"] = np.histogram(o["chi"][k], bins=HIST_EDGES, density=True)[0] * k.mean()
    with open(os.path.join(FIG, "e23_hist.dat"), "w") as fo:
        fo.write(" ".join(cols) + "\n")
        for i in range(len(cols["chi"])):
            fo.write(" ".join(f"{cols[k][i]:.5e}" for k in cols) + "\n")
    # 角度の補正を入れた場合(基準の組)
    m1, m2 = masses(rng, N, 0.5)
    a1, a2 = spins(rng, N, "beta")
    o0 = remnant.merge(m1, m2, a1, a2, np.random.default_rng(5))
    o1 = remnant.merge(m1, m2, a1, a2, np.random.default_rng(5), correct=True)
    res["angle_correction_baseline"] = {"median_shift": float(np.median(o1["chi"]) - np.median(o0["chi"])),
                                        "P(>=0.69)": [float(np.mean(o0["chi"] >= 0.69)), float(np.mean(o1["chi"] >= 0.69))],
                                        "P(>=0.72)": [float(np.mean(o0["chi"] >= 0.72)), float(np.mean(o1["chi"] >= 0.72))],
                                        "P(>=0.78)": [float(np.mean(o0["chi"] >= 0.78)), float(np.mean(o1["chi"] >= 0.78))]}
    # 文献の数値実験と同じ条件での比較(球状星団の模擬計算 1906.10260 Sec. III A・IV D: 親が回らないとき、残骸の約 60% が残り、
    # もう一度合体する 2 世代目のスピンは 90% が 0.67--0.69。親のスピンが 0.5 のとき、残るのは 3% 未満で、中央値 0.62、90% が
    # 0.53--0.72)。質量比と脱出速度の分布は同じにできないので、質量比 [0.5, 1] の一様と、脱出速度 50・100 km/s で見る。
    # Delta の符号を逆にした場合も出す(公開コードが直す前の符号で、文献の値が説明できるかを見るため)
    rng = np.random.default_rng(2310)
    n = 2000000
    q = rng.uniform(0.5, 1.0, n)
    cb, cg, dp = rng.uniform(-1, 1, n), rng.uniform(-1, 1, n), rng.uniform(0, 2 * np.pi, n)
    ct = np.cos(rng.uniform(0, np.pi, n))
    ca = cb * cg + np.sqrt(1 - cb ** 2) * np.sqrt(1 - cg ** 2) * np.cos(dp)
    cmc = {"N": n, "literature": {"chi_birth_0": {"retained": "nearly 60%", "p05_p95": [0.67, 0.69]},
                                  "chi_birth_0.5": {"retained": "< 3%", "median": 0.62, "p05_p95": [0.53, 0.72]}}}
    for a in (0.0, 0.5):
        af = remnant.final_spin(q, a, a, cb, cg, ca)
        vk = {"as_published": remnant.kick(q, a, a, cb, cg, dp, ct), "opposite_sign": kick_delta_sign(q, a, a, cb, cg, dp, ct, -1.0)}
        assert np.allclose(vk["as_published"], kick_delta_sign(q, a, a, cb, cg, dp, ct, +1.0), rtol=1e-13, atol=1e-10)
        for tag, v in vk.items():
            for ve in (50.0, 100.0):
                cmc[f"chi{a:g}_{tag}_vesc{ve:g}"] = summ(af, keep=v < ve)
    res["compare_cmc"] = cmc
    return res


# ---------- boson ----------
LF = np.arange(11.0, 19.01, 0.5)
N_B = 200000
HIST_LF = (None, 12.0, 14.0)


def parents_at_merger(model, rng, n, f, age_kw=None, qlo=0.5, zp=Z_PARENT):
    """親の対の、合体のときのスピン。f = None なら生まれたときのまま(ボソンなし)。返り値は (m1, m2, a1, a2, info)。
    乱数の列は f によらず同じ(f = None でも年齢を引く)。info は、印 3--5 の割合、スピンが 0.01 より減った親の割合、
    親の質量の減少の目安(211 が角運動量を持ち去るとして Delta M/M = alpha Delta chi)の平均と 99 パーセンタイル。"""
    m1, m2 = masses(rng, n, qlo)
    a1, a2 = spins(rng, n, "beta")
    am = ages.AgeModel(**(age_kw or {}))
    info = {"flagged": 0.0, "spun_down": 0.0, "mass_loss_mean": 0.0, "mass_loss_p99": 0.0}
    out, dm = [], []
    for m, a in ((m1, a1), (m2, a2)):
        T = am.sample(1, np.full(n, zp), rng)
        if f is None:
            out.append(a)
            continue
        c, r = model(m, np.minimum(a, 0.9989), T, MU, f)
        info["flagged"] += 0.5 * float(np.mean(r >= 3))
        info["spun_down"] += 0.5 * float(np.mean(a - c > 0.01))
        dm.append(ALPHA_PER_MSUN * m * np.maximum(a - c, 0.0))
        out.append(c)
    if dm:
        dm = np.concatenate(dm)
        info["mass_loss_mean"], info["mass_loss_p99"] = float(dm.mean()), float(np.percentile(dm, 99))
    return m1, m2, out[0], out[1], info


def _boson_one(arg):
    tag, kw, lf, zp = arg
    rng = np.random.default_rng(232)                         # f ごとに同じ乱数(差を f の効果だけにする)
    m1, m2, a1, a2, info = parents_at_merger(None if lf is None else chi_model(), rng, N_B, None if lf is None else 10 ** lf, kw, zp=zp)
    o = remnant.merge(m1, m2, a1, a2, rng)
    row = {"log10_f": lf, "z_parent": zp, "parents": info, "parent_spin_mean": float(0.5 * (a1.mean() + a2.mean())),
           "parent_spin_p95": float(np.percentile(np.concatenate([a1, a2]), 95)),
           "all": summ(o["chi"]), "vesc100": summ(o["chi"], keep=o["vkick"] < 100.0)}
    hist = np.histogram(o["chi"], bins=HIST_EDGES, density=True)[0] if (tag == "field" and zp == Z_PARENT and lf in HIST_LF) else None
    return tag, row, hist


def part_boson(nproc=6):
    from multiprocessing import Pool
    res = {"mu_eV": MU, "z_parent": Z_PARENT, "N": N_B, "log10_f": [float(x) for x in LF],
           "m_parents": [float(MTOT / 3), float(MTOT * 2 / 3)],
           "alpha_parents": [float(ALPHA_PER_MSUN * MTOT / 3), float(ALPHA_PER_MSUN * MTOT * 2 / 3)]}
    # 211 の臨界スピンの目安(m Omega_H = mu、非相対論的な極限): 4 alpha/(1 + 4 alpha^2)
    res["chi_c211_parents"] = [float(4 * a / (1 + 4 * a * a)) for a in res["alpha_parents"]]
    res["alpha_2G"] = float(ALPHA_PER_MSUN * M2G)
    presets = (("field", {}), ("old", dict(kappa1=0.0)))
    args = [(tag, kw, lf, Z_PARENT) for tag, kw in presets for lf in [None] + [float(x) for x in LF]]
    args += [("field", {}, lf, zp) for lf in (12.0, 13.0) for zp in (0.05, 1.0)]         # 親の合体の赤方偏移を変えた場合
    with Pool(nproc) as p:
        rows = p.map(_boson_one, args, chunksize=1)
    for tag, _ in presets:
        res[tag] = [r for (t, r, h), a in zip(rows, args) if t == tag and a[3] == Z_PARENT]
    res["z_parent_variants"] = [r for (t, r, h), a in zip(rows, args) if a[3] != Z_PARENT]
    cols = {"chi": 0.5 * (HIST_EDGES[1:] + HIST_EDGES[:-1])}
    for (t, r, h), a in zip(rows, args):
        if h is not None:
            cols["none" if a[2] is None else f"lf{a[2]:g}"] = h
    with open(os.path.join(FIG, "e23_hist_boson.dat"), "w") as fo:
        fo.write(" ".join(cols) + "\n")
        for i in range(len(cols["chi"])):
            fo.write(" ".join(f"{cols[k][i]:.5e}" for k in cols) + "\n")
    return res


# ---------- gw241011 ----------
TWO = {"wide": dict(), "peak": dict(family2="lognormal", t2med=5e7, s2=0.7, t2lo=1e5), "slow": dict(t2lo=3e8)}   # E22 と同じ組
LF2 = np.linspace(10.0, 14.0, 41)
N_GW = 100000
VESC_GW = 100.0
FIXB = (0.70, 0.71, 0.72, 0.74, 0.76, 0.78)
K_AGE = 400
LEVELS = (0.9, 0.5, 0.1)
KEEP_ALT = (0.68, 0.70)


def _gw_one(lf):
    """f = 10^lf で、2 世代目が合体のときに KEEP 以上のスピンを持つ確率。親が超放射を受ける場合(sc)と受けない場合(fixed)、
    星団に残ったもの(キック < VESC_GW)だけの場合(_ret)、生まれたスピンを一つの値に固定した場合(_fix)。
    乱数は f によらず同じにする(差を f の効果だけにする)。"""
    model = chi_model()
    f = 10 ** lf
    T2 = {nm: ages.AgeModel(**kw).sample(2, np.full(N_GW, Z2), np.random.default_rng(234)) for nm, kw in TWO.items()}
    out = {}
    for tag, f_par in (("sc", f), ("fixed", None)):
        rng = np.random.default_rng(233)
        m1, m2, a1, a2, _ = parents_at_merger(model, rng, N_GW, f_par)
        o = remnant.merge(m1, m2, a1, a2, rng)
        ret = o["vkick"] < VESC_GW
        born = o["chi"] >= KEEP
        out[f"born_{tag}"], out[f"born_{tag}_ret"], out[f"retained_{tag}"] = float(born.mean()), float(born[ret].mean()), float(ret.mean())
        for nm in TWO:
            c, r = model(o["M"], np.minimum(o["chi"], 0.9989), T2[nm], MU, f)
            ok, fl = c >= KEEP, r >= 3
            out[f"{nm}_{tag}"], out[f"{nm}_{tag}_ret"] = float(ok.mean()), float(ok[ret].mean())
            # 印 3--5 の標本を「保たない」に数えた下限と、「保つ」に数えた上限(スピンは増えないので、生まれたときに KEEP 以上のものだけ)
            out[f"{nm}_{tag}_lo"], out[f"{nm}_{tag}_hi"] = float((ok & ~fl).mean()), float((ok | (fl & born)).mean())
            out[f"{nm}_{tag}_flagged"] = float(fl.mean())
            for ka in KEEP_ALT:                              # 判定の閾値を動かした場合
                out[f"{nm}_{tag}_k{ka:.2f}"] = float(np.mean(c >= ka))
        for ka in KEEP_ALT:
            out[f"born_{tag}_k{ka:.2f}"] = float(np.mean(o["chi"] >= ka))
    for nm, kw in TWO.items():                               # E22 と同じ求積(年齢の節 K_AGE 個)
        Tn = ages.AgeModel(**kw).nodes(2, np.array([Z2]), K_AGE)[0]
        for cb in FIXB:
            c, r = model(np.full(K_AGE, M2G), np.full(K_AGE, cb), Tn, MU, f)
            out[f"{nm}_fix{cb:.2f}"] = float(np.mean(c >= KEEP))
    return out


def part_gw241011(nproc=6):
    from multiprocessing import Pool
    res = {"mu_eV": MU, "keep_threshold": KEEP, "N": N_GW, "z_2G": Z2, "vesc_kms": VESC_GW, "log10_f": [float(LF2[0]), float(LF2[-1]), len(LF2)]}
    # ボソンなしの基準: 2 世代目が生まれたときに KEEP 以上である確率と、残骸の質量
    rng = np.random.default_rng(233)
    m1, m2, a1, a2, _ = parents_at_merger(None, rng, N_GW, None)
    o = remnant.merge(m1, m2, a1, a2, rng)
    ret = o["vkick"] < VESC_GW
    base, base_ret = float(np.mean(o["chi"] >= KEEP)), float(np.mean(o["chi"][ret] >= KEEP))
    res["no_boson"] = {"P_born_ge_keep": base, "P_born_ge_keep_retained": base_ret, "retained": float(ret.mean()),
                       "P_born_ge_0.87": float(np.mean(o["chi"] >= 0.87)), "M_p05_p50_p95": [float(np.percentile(o["M"], x)) for x in (5, 50, 95)],
                       "born_ge_keep_median_spin": float(np.median(o["chi"][o["chi"] >= KEEP])),
                       "born_ge_keep_median_spin_retained": float(np.median(o["chi"][ret & (o["chi"] >= KEEP)]))}
    with Pool(nproc) as p:
        rows = p.map(_gw_one, [float(x) for x in LF2], chunksize=1)
    assert all(abs(r["born_fixed"] - base) < 1e-12 and abs(r["born_fixed_ret"] - base_ret) < 1e-12 for r in rows)
    curves = {k: [r[k] for r in rows] for k in rows[0]}
    res["curves"] = curves

    def f_at(pc, level):
        y = np.array(pc) - level
        idx = np.flatnonzero((y[:-1] >= 0) & (y[1:] < 0))
        if len(idx) == 0:
            return None
        i = idx[0]
        return float(10 ** (LF2[i] + (LF2[i + 1] - LF2[i]) * (0 - y[i]) / (y[i + 1] - y[i])))

    # 確率が、ボソンなしの値の 0.9・0.5・0.1 倍になる f
    norm = {}
    for nm in TWO:
        for t in ("sc", "fixed"):
            norm[f"{nm}_{t}"], norm[f"{nm}_{t}_ret"] = base, base_ret
            norm[f"{nm}_{t}_lo"], norm[f"{nm}_{t}_hi"] = base, base
        for cb in FIXB:
            norm[f"{nm}_fix{cb:.2f}"] = 1.0
    base_alt = {ka: float(np.mean(o["chi"] >= ka)) for ka in KEEP_ALT}
    res["no_boson"]["P_born_ge_alt_threshold"] = {f"{ka:.2f}": v for ka, v in base_alt.items()}
    for nm in TWO:
        for ka in KEEP_ALT:
            norm[f"{nm}_sc_k{ka:.2f}"] = base_alt[ka]
    res["f_cross"] = {k: {f"{lv:g}": f_at(curves[k], lv * nrm) for lv in LEVELS} for k, nrm in norm.items()}
    res["P_over_base_at_smallest_f"] = {k: curves[k][0] / nrm for k, nrm in norm.items()}
    fc = res["f_cross"]
    # 生まれたスピンを固定した場合との比と、同じ f を与える「実効的な生まれたスピン」
    res["f_ratio_to_fixed_birth_spin"] = {f"{nm}{sfx}_{lv:g}": {f"{cb:.2f}": fc[f"{nm}_sc{sfx}"][f"{lv:g}"] / fc[f"{nm}_fix{cb:.2f}"][f"{lv:g}"] for cb in (0.72, 0.78)}
                                          for nm in TWO for sfx in ("", "_ret") for lv in LEVELS}
    eff = {}
    for nm in TWO:
        for sfx in ("", "_ret"):
            for lv in LEVELS:
                x = np.log([fc[f"{nm}_fix{cb:.2f}"][f"{lv:g}"] for cb in FIXB])
                assert np.all(np.diff(x) > 0)
                lx = np.log(fc[f"{nm}_sc{sfx}"][f"{lv:g}"])
                eff[f"{nm}{sfx}_{lv:g}"] = float(np.interp(lx, x, FIXB)) if x[0] <= lx <= x[-1] else None
    res["effective_birth_spin"] = eff
    # 窓に要る年齢の比 R(M)(E22)は f の上限の 2 乗に反比例する: 固定した場合の R に掛かる倍率
    res["age_ratio_factor"] = {k: {cb: 1.0 / r ** 2 for cb, r in v.items()} for k, v in res["f_ratio_to_fixed_birth_spin"].items()}
    # 判定の閾値を 0.69 から動かしたときの f の比
    res["threshold_f_ratio"] = {f"{nm}_{lv:g}": {f"{ka:.2f}": fc[f"{nm}_sc_k{ka:.2f}"][f"{lv:g}"] / fc[f"{nm}_sc"][f"{lv:g}"] for ka in KEEP_ALT}
                                for nm in TWO for lv in LEVELS}
    # 印 3--5 の標本の扱いによる幅、親の超放射を入れた場合と入れない場合の差、キックで選んだ場合
    res["flagged_max"] = {nm: max(curves[f"{nm}_sc_flagged"]) for nm in TWO}
    res["flag_bounds_max_dP"] = {nm: float(np.max(np.array(curves[f"{nm}_sc_hi"]) - np.array(curves[f"{nm}_sc_lo"]))) for nm in TWO}
    res["flag_bounds_f_cross_ratio"] = {f"{nm}_{lv:g}": [fc[f"{nm}_sc_lo"][f"{lv:g}"] / fc[f"{nm}_sc"][f"{lv:g}"], fc[f"{nm}_sc_hi"][f"{lv:g}"] / fc[f"{nm}_sc"][f"{lv:g}"]]
                                        for nm in TWO for lv in LEVELS}
    res["sc_minus_fixed_max_abs_dP"] = {nm: float(np.max(np.abs(np.array(curves[f"{nm}_sc"]) - np.array(curves[f"{nm}_fixed"])))) for nm in TWO}
    res["sc_over_fixed_f_cross"] = {f"{nm}_{lv:g}": fc[f"{nm}_sc"][f"{lv:g}"] / fc[f"{nm}_fixed"][f"{lv:g}"] for nm in TWO for lv in LEVELS}
    res["retained_over_all_f_cross"] = {f"{nm}_{lv:g}": fc[f"{nm}_sc_ret"][f"{lv:g}"] / fc[f"{nm}_sc"][f"{lv:g}"] for nm in TWO for lv in LEVELS}
    res["born_sc_range"] = [float(min(curves["born_sc"])), float(max(curves["born_sc"]))]
    with open(os.path.join(FIG, "e23_keep.dat"), "w") as fo:
        fo.write("f " + " ".join(TWO) + " " + " ".join(f"{nm}_ret" for nm in TWO) + " " + " ".join(f"{nm}_fix078" for nm in TWO)
                 + " " + " ".join(f"{nm}_fix072" for nm in TWO) + " born\n")
        for i in range(len(LF2)):
            vals = ([curves[f"{nm}_sc"][i] / base for nm in TWO] + [curves[f"{nm}_sc_ret"][i] / base_ret for nm in TWO]
                    + [curves[f"{nm}_fix0.78"][i] for nm in TWO] + [curves[f"{nm}_fix0.72"][i] for nm in TWO] + [curves["born_sc"][i] / base])
            fo.write(f"{10 ** LF2[i]:.5e} " + " ".join(f"{v:.5f}" for v in vals) + "\n")
    return res


# ---------- window ----------
def part_window():
    """E22 の「同じ f で両立する度合い」W = max_f min(P_keep, P_drop) を、生まれたスピンの分布で周辺化した P_keep/P_0 で出し直す。
    P_drop(質量 M の 1 世代目が、生まれたスピン 0.7 から「211 の臨界 + 0.05」より下に落ちている確率)は、E22 の図のデータ
    figs/data/e22_prob_drop.dat(f の刻み 0.025 dex)。こちらの曲線(gw241011 の部分、刻み 0.1 dex)は log f について線形に
    内挿する。内挿の誤差は、生まれたスピンを固定した曲線と W を、E22 の値(e22_prob_keep.dat、e22.json)と比べて見る。"""
    g = json.load(open(os.path.join(OUT, "e23.json")))["gw241011"]
    cv, nb = g["curves"], g["no_boson"]
    base, base_ret = nb["P_born_ge_keep"], nb["P_born_ge_keep_retained"]

    def load(name):
        fn = os.path.join(FIG, name)
        cols = open(fn).readline().split()
        a = np.loadtxt(fn, skiprows=1)
        return {c: a[:, i] for i, c in enumerate(cols)}

    drop, keep22 = load("e22_prob_drop.dat"), load("e22_prob_keep.dat")
    lf = np.log10(drop["f"])
    assert np.allclose(lf, np.log10(keep22["f"])) and LF2[0] <= lf[0] and lf[-1] <= LF2[-1] + 1e-9
    e22w = json.load(open(os.path.join(OUT, "e22.json")))["window"]["W"]
    res = {"log10_f": [float(lf[0]), float(lf[-1]), len(lf)], "W": {}}
    worst_curve, worst_w = 0.0, 0.0
    for nm in TWO:
        kc = {"dist": np.interp(lf, LF2, cv[f"{nm}_sc"]) / base, "dist_ret": np.interp(lf, LF2, cv[f"{nm}_sc_ret"]) / base_ret,
              "fix0.78": np.interp(lf, LF2, cv[f"{nm}_fix0.78"]), "fix0.72": np.interp(lf, LF2, cv[f"{nm}_fix0.72"])}
        for cb, col in (("0.78", f"{nm}_078"), ("0.72", f"{nm}_072")):
            worst_curve = max(worst_curve, float(np.max(np.abs(kc[f"fix{cb}"] - keep22[col]))))
        for one in drop:
            if one == "f":
                continue
            for lab, pk in kc.items():
                w = np.minimum(pk, drop[one])
                i = int(np.argmax(w))
                res["W"][f"{nm}_{lab}|{one}"] = {"W": float(w[i]), "f_GeV": float(10 ** lf[i])}
            for cb in ("0.78", "0.72"):
                worst_w = max(worst_w, abs(res["W"][f"{nm}_fix{cb}|{one}"]["W"] - e22w[f"{nm}_{cb}|{one}"]["W"]))
    res["check"] = {"fixed_curve_vs_e22_max_abs_dP": worst_curve, "fixed_W_vs_e22_max_abs_diff": worst_w}
    # 分布で周辺化した W と、固定した場合の W の差(1 世代目は既定と古い組、10--25 Msun)
    ones = [k for k in drop if k != "f" and not k.startswith("flat")]
    res["W_dist_minus_fixed"] = {cb: [float(min(res["W"][f"{nm}_dist|{o}"]["W"] - res["W"][f"{nm}_fix{cb}|{o}"]["W"] for nm in TWO for o in ones)),
                                      float(max(res["W"][f"{nm}_dist|{o}"]["W"] - res["W"][f"{nm}_fix{cb}|{o}"]["W"] for nm in TWO for o in ones))]
                                 for cb in ("0.78", "0.72")}
    res["W_ret_minus_dist"] = [float(min(res["W"][f"{nm}_dist_ret|{o}"]["W"] - res["W"][f"{nm}_dist|{o}"]["W"] for nm in TWO for o in ones)),
                               float(max(res["W"][f"{nm}_dist_ret|{o}"]["W"] - res["W"][f"{nm}_dist|{o}"]["W"] for nm in TWO for o in ones))]
    return res


PARTS = {"check": part_check, "dist": part_dist, "boson": part_boson, "gw241011": part_gw241011, "window": part_window}


def main():
    names = sys.argv[1:] or list(PARTS)
    path = os.path.join(OUT, "e23.json")
    res = json.load(open(path)) if os.path.exists(path) else {}
    for nm in names:
        t0 = time.time()
        res[nm] = PARTS[nm]()
        print(f"--- {nm} ({time.time() - t0:.0f} s)", flush=True)
        show = res[nm]
        if nm == "gw241011":
            show = {k: v for k, v in show.items() if k != "curves"}
        print(json.dumps(show, indent=1), flush=True)
        with open(path, "w") as fo:
            json.dump({k: res[k] for k in PARTS if k in res}, fo, indent=1)


if __name__ == "__main__":
    main()
