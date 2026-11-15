"""合体の残骸の質量・スピン・キック: 2 世代目の BH が生まれたときの質量とスピン(numpy だけで動く)。

使い方:
    import numpy as np, remnant
    rng = np.random.default_rng(0)
    out = remnant.merge(m_a, m_b, chi_a, chi_b, rng)     # 親の質量と、合体のときのスピンの大きさ(配列)。向きは等方に引く
    out["M"]        # 残骸の質量(親の質量と同じ単位)= 2 世代目の生まれたときの質量
    out["chi"]      # 残骸のスピン = 2 世代目の生まれたときのスピン
    out["vkick"]    # 反跳(キック)の速さ [km/s]。星団に残るのは out["vkick"] < v_esc(脱出速度)のもの
親のスピンは「合体のときの」値を渡す。超放射で回転を失う場合は、生まれたときのスピンに chi_SR を掛けたあとの値である。

式(記号は添え字 1 が重いほう、q = m2/m1 <= 1、nu = q/(1+q)^2。beta・gamma はスピン 1・2 と軌道角運動量 L のなす角、
alpha はスピンどうしのなす角。根拠と確かめは研究ノート E23、README.md):
  スピン: Hofmann, Barausse & Rezzolla 2016(arXiv:1605.01938)の式 13--16、係数は Table 1 の上段(n_M = 1、n_J = 2)。
      a_tot = (a1 cos b + a2 cos g q^2)/(1+q)^2、  a_eff = a_tot + xi nu (a1 cos b + a2 cos g)、
      |l| = |L_ISCO(a_eff) - 2 a_tot (E_ISCO(a_eff) - 1) + sum_ij k_ij nu^(1+i) a_eff^j|、
      |a_fin| = (1+q)^-2 [a1^2 + a2^2 q^4 + 2 a1 a2 q^2 cos al + 2 (a1 cos b + a2 q^2 cos g) |l| q + |l|^2 q^2]^(1/2)。
    当てはめの残差は、歳差のある一般のスピンで 0.006 程度(同論文の要旨)。
  質量: Barausse, Morozova & Rezzolla 2012(arXiv:1206.3803)の式 18。
      E_rad/M = [1 - E_ISCO(a_tot)] nu + 4 nu^2 [4 p0 + 16 p1 a_tot (a_tot + 1) + E_ISCO(a_tot) - 1]、  M_fin = M (1 - E_rad/M)。
  キック: Gerosa & Kesden 2016(arXiv:1605.01067)の式 44--47。スピンの組み合わせは同論文の式 36:
      Delta = (a1 S1 - q a2 S2)/(1+q)、  chi~ = (a1 S1 + q^2 a2 S2)/(1+q)^2   (S1・S2 はスピンの向きの単位ベクトル)。
      v_k = sqrt(v_m^2 + 2 v_m v_perp cos zeta + v_perp^2 + v_par^2)、
      v_m = A nu^2 (1-q)/(1+q) (1 + B nu)、  v_perp = H nu^2 Delta_par、
      v_par = 16 nu^2 [Delta_perp (V11 + 2 VA c + 4 VB c^2 + 8 VC c^3) + 2 chi~_perp Delta_par (C2 + 2 C3 c)] cos Theta、 c = chi~_par。
    par・perp は L に平行・垂直な成分。Theta(合体のときの落下の向きとの角)は [0, pi] で一様に引く(同論文の扱い)。
    Delta の符号に注意: 重いほうが L と同じ向きに回っているとき Delta_par > 0 で、質量比の項と一部打ち消し合う(cos zeta < 0)。
    同論文の公開コードは v1.0.2 でこの符号に直されている(コードの注記)。逆の符号で書かれた式もあるので、写すときは確かめること。
    キックの当てはめは模型による違いが大きい(数値相対論の代理模型より平均で 10 倍大きいという比較がある: arXiv:2503.21278)。
    残るかどうかの選択は、その不定性ごと使うこと(研究ノート E23 の注意)。

確かめ: python remnant.py(決まった入力での値を照らす。合えば ok と出る)。
"""
import numpy as np

K01, K02, K10, K11, K12, XI = -1.2019, -1.20764, 3.79245, 1.18385, 4.90494, 0.41616        # 1605.01938 Table 1 上段
K00 = 16.0 * (0.68646 - np.sqrt(3.0) / 2.0) - K10 / 4.0                                     # 同 式 11(等質量・無スピンで 0.68646)
EPS_ANGLE = 0.024                                                                           # 同 式 18 の角度の補正(beta、gamma)
P0, P1 = 0.04827, 0.01707                                                                   # 1206.3803 式 12
KICK_A, KICK_B, KICK_H, KICK_ZETA = 1.2e4, -0.93, 6.9e3, np.deg2rad(145.0)                  # 1605.01067 式 44--46 [km/s]
V11, VA, VB, VC, C2, C3 = 3677.76, 2481.21, 1792.45, 1506.52, 1140.0, 2481.0                # 同 式 47 [km/s]


def r_isco(a):
    a = np.asarray(a, float)
    z1 = 1 + np.cbrt(1 - a ** 2) * (np.cbrt(1 + a) + np.cbrt(1 - a))
    z2 = np.sqrt(3 * a ** 2 + z1 ** 2)
    return 3 + z2 - np.sign(a) * np.sqrt((3 - z1) * (3 + z1 + 2 * z2))


def e_isco(a):
    return np.sqrt(1 - 2.0 / (3.0 * r_isco(a)))


def l_isco(a):
    return 2.0 / (3.0 * np.sqrt(3.0)) * (1 + 2 * np.sqrt(3 * r_isco(a) - 2))


def _corrected(cosx):
    """1605.01938 の式 18: 角 Theta を Theta* = 2 arctan[(1 + eps) tan(Theta/2)] に置き換えた余弦。"""
    th = np.arccos(np.clip(cosx, -1.0, 1.0))
    return np.cos(2 * np.arctan((1 + EPS_ANGLE) * np.tan(th / 2)))


def final_spin(q, a1, a2, cosb, cosg, cosal, correct=False):
    """残骸のスピンの大きさ。correct=True で式 18 の角度の補正(当てはめの偏り -0.005 を -0.001 にする)を入れる。"""
    q = np.asarray(q, float)
    if correct:
        cosb, cosg = _corrected(cosb), _corrected(cosg)
    nu = q / (1 + q) ** 2
    s = a1 * cosb + a2 * cosg
    atot = (a1 * cosb + a2 * cosg * q ** 2) / (1 + q) ** 2
    aeff = atot + XI * nu * s
    ksum = (K00 + K01 * aeff + K02 * aeff ** 2) * nu + (K10 + K11 * aeff + K12 * aeff ** 2) * nu ** 2
    ell = np.abs(l_isco(aeff) - 2 * atot * (e_isco(aeff) - 1) + ksum)
    inner = (a1 ** 2 + a2 ** 2 * q ** 4 + 2 * a1 * a2 * q ** 2 * cosal
             + 2 * (a1 * cosb + a2 * q ** 2 * cosg) * ell * q + ell ** 2 * q ** 2)
    # inner は |a_tot + l nu|^2 (1+q)^4 で、角が互いに矛盾しなければ負にならない。角度の補正を入れたときは beta・gamma だけを
    # 動かすので、残骸のスピンがほぼ 0 の配置でわずかに負になりうる。0 で切る
    return np.sqrt(np.maximum(inner, 0.0)) / (1 + q) ** 2


def final_mass_fraction(q, a1, a2, cosb, cosg):
    """残骸の質量 / 親の質量の和。"""
    q = np.asarray(q, float)
    nu = q / (1 + q) ** 2
    atot = (a1 * cosb + a2 * cosg * q ** 2) / (1 + q) ** 2
    e = e_isco(atot)
    return 1.0 - ((1 - e) * nu + 4 * nu ** 2 * (4 * P0 + 16 * P1 * atot * (atot + 1) + e - 1))


def kick(q, a1, a2, cosb, cosg, dphi, cos_theta):
    """キックの速さ [km/s]。dphi は二つのスピンの、L のまわりの方位角の差。cos_theta は式 47 の cos Theta。"""
    q = np.asarray(q, float)
    nu = q / (1 + q) ** 2
    sb, sg = np.sqrt(1 - np.asarray(cosb, float) ** 2), np.sqrt(1 - np.asarray(cosg, float) ** 2)
    d_par = (a1 * cosb - q * a2 * cosg) / (1 + q)
    d_perp = np.sqrt(np.maximum(a1 ** 2 * sb ** 2 + q ** 2 * a2 ** 2 * sg ** 2 - 2 * q * a1 * a2 * sb * sg * np.cos(dphi), 0.0)) / (1 + q)
    c_par = (a1 * cosb + q ** 2 * a2 * cosg) / (1 + q) ** 2
    c_perp = np.sqrt(np.maximum(a1 ** 2 * sb ** 2 + q ** 4 * a2 ** 2 * sg ** 2 + 2 * q ** 2 * a1 * a2 * sb * sg * np.cos(dphi), 0.0)) / (1 + q) ** 2
    vm = KICK_A * nu ** 2 * (1 - q) / (1 + q) * (1 + KICK_B * nu)
    vperp = KICK_H * nu ** 2 * d_par
    vpar = 16 * nu ** 2 * (d_perp * (V11 + 2 * VA * c_par + 4 * VB * c_par ** 2 + 8 * VC * c_par ** 3)
                           + 2 * c_perp * d_par * (C2 + 2 * C3 * c_par)) * cos_theta
    return np.sqrt(vm ** 2 + 2 * vm * vperp * np.cos(KICK_ZETA) + vperp ** 2 + vpar ** 2)


def merge(m_a, m_b, chi_a, chi_b, rng, correct=False):
    """親の対(質量 m_a, m_b と、合体のときのスピンの大きさ chi_a, chi_b。配列)から、残骸を一つずつ作る。
    スピンの向きは等方(cos beta、cos gamma が [-1, 1] で一様、方位角の差が一様)、Theta は [0, pi] で一様に引く。
    返り値は dict: M(残骸の質量)、chi(残骸のスピン)、vkick [km/s]、q、cosb、cosg(重いほう、軽いほうの順)。"""
    m_a, m_b, chi_a, chi_b = np.broadcast_arrays(*[np.asarray(x, float) for x in (m_a, m_b, chi_a, chi_b)])
    sw = m_b > m_a
    m1, m2 = np.where(sw, m_b, m_a), np.where(sw, m_a, m_b)
    a1, a2 = np.where(sw, chi_b, chi_a), np.where(sw, chi_a, chi_b)
    q = m2 / m1
    n = m1.shape
    cosb, cosg = rng.uniform(-1, 1, n), rng.uniform(-1, 1, n)
    dphi = rng.uniform(0, 2 * np.pi, n)
    cos_theta = np.cos(rng.uniform(0, np.pi, n))
    cosal = cosb * cosg + np.sqrt(1 - cosb ** 2) * np.sqrt(1 - cosg ** 2) * np.cos(dphi)
    # スピンの当てはめは、重いほうのスピンが 0.99 を越えて L にほぼ揃い、質量比が小さい(0.03 付近)とき、1 をわずかに越えうる
    # (最大 0.0014)。1 で切る
    return {"M": (m1 + m2) * final_mass_fraction(q, a1, a2, cosb, cosg),
            "chi": np.minimum(final_spin(q, a1, a2, cosb, cosg, cosal, correct), 1.0),
            "vkick": kick(q, a1, a2, cosb, cosg, dphi, cos_theta), "q": q, "cosb": cosb, "cosg": cosg}


def selfcheck():
    """決まった入力での値を照らす(渡した先で、まず python remnant.py を走らせる)。値は研究ノート E23 の計算機でのもの。"""
    q, a1, a2, cb, cg, dp, ct = 0.6, 0.7, 0.4, 0.3, -0.5, 1.0, 0.8
    ca = cb * cg + np.sqrt(1 - cb ** 2) * np.sqrt(1 - cg ** 2) * np.cos(dp)
    qq = np.linspace(0.01, 1, 9901)
    vm = kick(qq, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0)
    rows = [("スピン: 等質量・無スピン", final_spin(1, 0, 0, 1, 1, 1), 0.68646),
            ("スピン: q = 1/3・無スピン", final_spin(1 / 3, 0, 0, 1, 1, 1), 0.540180493190051),
            ("スピン: 等質量・揃った 0.5", final_spin(1, 0.5, 0.5, 1, 1, 1), 0.8321392072441598),
            ("スピン: 一般の向き", final_spin(q, a1, a2, cb, cg, ca), 0.7495169498129935),
            ("スピン: 一般の向き、角度の補正つき", final_spin(q, a1, a2, cb, cg, ca, correct=True), 0.7457498548182114),
            ("質量: 等質量・無スピン", final_mass_fraction(1, 0, 0, 1, 1), 0.95173),
            ("質量: 等質量・揃った最大スピンの放出", 1 - final_mass_fraction(1, 1, 1, 1, 1), 0.09948),
            ("質量: 一般の向き", final_mass_fraction(q, a1, a2, cb, cg), 0.9533017474048704),
            ("キック: 無スピンの最大 [km/s]", vm.max(), 175.21082377628568),
            ("キック: 無スピンの最大になる q", qq[vm.argmax()], 0.362),
            ("キック: スーパーキック [km/s]", kick(1, 1, 1, 0, 0, np.pi, 1), 3677.76),
            ("キック: 一般の向き [km/s]", kick(q, a1, a2, cb, cg, dp, ct), 1136.6080576011552),
            ("キック: q = 0.5、重いほうだけ 0.5 で順行 [km/s]", kick(0.5, 0.5, 0.0, 1, 1, 0, 0), 91.09215906842495),
            ("キック: q = 0.5、重いほうだけ 0.5 で逆行 [km/s]", kick(0.5, 0.5, 0.0, -1, 1, 0, 0), 258.10432883191555)]
    bad = [(nm, float(v), ref) for nm, v, ref in rows if not abs(float(v) - ref) <= 1e-9 * max(1.0, abs(ref))]
    rng = np.random.default_rng(0)
    out = merge(rng.uniform(5, 30, 20000), rng.uniform(5, 30, 20000), rng.uniform(0, 1, 20000), rng.uniform(0, 1, 20000), rng)
    ok = (np.all(np.isfinite(out["M"])) and np.all(np.isfinite(out["vkick"])) and out["chi"].min() >= 0 and out["chi"].max() <= 1
          and out["q"].max() <= 1 and out["vkick"].min() >= 0)
    if not ok:
        bad.append(("merge: 値の範囲", float("nan"), float("nan")))
    return bad


if __name__ == "__main__":
    _bad = selfcheck()
    for _nm, _v, _ref in _bad:
        print(f"違う: {_nm}: {_v!r}(基準 {_ref!r})")
    print("ok" if not _bad else f"{len(_bad)} 個が基準と違う")
    raise SystemExit(1 if _bad else 0)
