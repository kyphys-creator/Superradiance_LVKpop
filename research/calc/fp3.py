"""211・322・544 の三準位の運動論と、その固定点(研究ノート E14、E15、E16)。

率は BGLS(2011.11646)の Table IV--VI(占有数の率 / mu、(M_pl/f)^4 は外に出す):
  a: 211 x 211 -> 322 x BH   4.3e-7 alpha^11 r+     (Table V)
  b: 322 x 322 -> 211 x inf  1.1e-8 alpha^8         (Table VI)
  c: 322 x 322 -> 544 x BH   1.9e-9 alpha^11 r+     (Table V)
  d: 544 x 544 -> 322 x inf  4.4e-11 alpha^8        (Table VI)
  e: 322 x 544 -> 211 x inf  6.1e-11 alpha^8        (Table VI)
  f: 544 x 544 -> 211 x inf  4.2e-13 alpha^8        (Table VI)
超放射: 211・322 は E4 の相対論的な率(表の外は解析式)、544 は解析式(1908.10370 の式 2.32)。
  de1/dt = g1 e1 - 2a e1^2 e2 + b e1 e2^2 + e e1 e2 e3 + f e1 e3^2
  de2/dt = g2 e2 + a e1^2 e2 - 2b e1 e2^2 - 2c e2^2 e3 + d e2 e3^2 - e e1 e2 e3
  de3/dt = g3 e3 + c e2^2 e3 - 2d e2 e3^2 - e e1 e2 e3 - 2f e1 e3^2
  dchi/dt = -(g1 e1 + 2 g2 e2 + 4 g3 e3)
eps = (f/M_pl)^2 x et とおくと f が消える(et の式は上で係数 (M_pl/f)^4 を落としたもの)。固定点の et は f によらず、
dchi/dt = mu (f/M_pl)^2 H(alpha, chi)、H = -(g1 et1 + 2 g2 et2 + 4 g3 et3)。
固定点の型: P3 = (211, 322)(E12)、P2 = (322, 544)。空の準位の方向の安定性で選ぶ(fixed_points)。
固定点の式と、三準位の外の準位の増加率の式は、この file の fixed_points と extra_rates だけに置き、chisr はそれを呼ぶ(E16)。
"""
import numpy as np

import si
import sr

K_A, K_B, K_C, K_D, K_E, K_F = 4.3e-7, 1.1e-8, 1.9e-9, 4.4e-11, 6.1e-11, 4.2e-13


def g_rel(alpha, chi, n, l, m):
    """超放射の率 / mu(占有数の率、スカラー)。"""
    return sr.gamma_rel_over_mu(alpha, chi, n, l, m, sr.NORM_N)


def coeffs(alpha, chi, cbh=1.0, cbh_all=False):
    """cbh: 211 x 211 -> 322 x BH に掛ける有限 alpha の補正。cbh_all=True なら BH が吸収する過程すべてに同じ補正を掛ける
    (BGLS: 補正は地平線の近くで決まるので、BH が吸収する過程どうしの比はあまり変えない)。"""
    g = np.array([g_rel(alpha, chi, 2, 1, 1), g_rel(alpha, chi, 3, 2, 2),
                  sr.gamma_over_mu(alpha, chi, 5, 4, 4, sr.NORM_N)])
    rp = sr.rplus(chi)
    k = dict(a=K_A * alpha ** 11 * rp * cbh, b=K_B * alpha ** 8, c=K_C * alpha ** 11 * rp * (cbh if cbh_all else 1.0),
             d=K_D * alpha ** 8, e=K_E * alpha ** 8, f=K_F * alpha ** 8)
    return g, k


def percap(x, g, k):
    """各準位の 1 粒子あたりの増加率(e_i で割ったもの)。x = (e1, e2, e3)。"""
    e1, e2, e3 = x
    r1 = g[0] - 2 * k["a"] * e1 * e2 + k["b"] * e2 ** 2 + k["e"] * e2 * e3 + k["f"] * e3 ** 2
    r2 = g[1] + k["a"] * e1 ** 2 - 2 * k["b"] * e1 * e2 - 2 * k["c"] * e2 * e3 + k["d"] * e3 ** 2 - k["e"] * e1 * e3
    r3 = g[2] + k["c"] * e2 ** 2 - 2 * k["d"] * e2 * e3 - k["e"] * e1 * e2 - 2 * k["f"] * e1 * e3
    return np.array([r1, r2, r3])


def fixed_points(g1, g2, g3, ka, kb, kc, kd, ke, kf):
    """三準位の固定点(配列で渡せる、f によらない et)。544 の超放射は固定点の式では落とす(alpha^20 で小さい)。
    返り値の dict:
      e1, e2: P3 = (211, 322)。u = ka e1^2 は 3u^2 + (2 g2 - 4 kb g1/ka) u - g2^2 = 0 の正の根、e2 = (u + g2)/(2 kb e1)。
      f2, f3: P2 = (322, 544)。f2 = (2/sqrt3) sqrt(g2 kd)/kc、f3 = kc f2/(2 kd)。
      lam544: P3 での 544 の方向の増加率、lam211: P2 での 211 の方向の増加率。
      use3: P3 を使う(211 の臨界より上、または 544 の方向に安定)。use2: use3 でなく、P2 が 211 の方向に安定。
      どちらでもない点(211 の臨界の下 1e-3 ほど)は、両方とも不安定な幅。
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        B = 2 * g2 - 4 * kb * g1 / ka
        D = np.sqrt(B * B + 12 * g2 * g2)
        u = np.where(B > 0, 2 * g2 * g2 / (B + D), (D - B) / 6)   # 桁落ちしない形の正の根
        e1 = np.sqrt(u / ka)
        e2 = (u + g2) / (2 * kb * e1)
        lam544 = g3 + kc * e2 ** 2 - ke * e1 * e2
        f2 = 2 / np.sqrt(3) * np.sqrt(np.maximum(g2, 0) * kd) / kc
        f3 = kc * f2 / (2 * kd)
        lam211 = g1 + kb * f2 ** 2 + ke * f2 * f3 + kf * f3 ** 2
    use3 = (g1 >= 0) | (lam544 < 0)
    use2 = (~use3) & (lam211 < 0)
    return dict(e1=e1, e2=e2, f2=f2, f3=f3, lam544=lam544, lam211=lam211, use3=use3, use2=use2)


def _fp_from(g, k):
    return fixed_points(g[0], g[1], g[2], k["a"], k["b"], k["c"], k["d"], k["e"], k["f"])


def fp_P3(g, k):
    """(211, 322) の固定点(E12 の式)。"""
    if g[1] <= 0:
        return None
    r = _fp_from(g, k)
    return np.array([float(r["e1"]), float(r["e2"]), 0.0])


def fp_P2(g, k):
    """(322, 544) の固定点(544 の超放射は落とす)。"""
    if g[1] <= 0:
        return None
    r = _fp_from(g, k)
    return np.array([0.0, float(r["f2"]), float(r["f3"])])


def H_of(x, g):
    return -(g[0] * x[0] + 2 * g[1] * x[1] + 4 * g[2] * x[2])


# ---------- 照合用の ODE(1 個の量子から、BH の質量は固定、重力波の過程は入れない) ----------
def evolve3(M, chi0, mu_eV, f_gev, T_yr, t_eval_yr, levels0=(1, 1, 1), rtol=1e-8):
    """三準位の運動論を 1 個の量子から解く。各過程の生成側の占有数に 1 個の量子(e_q = 1/(G M^2))を足して自発的な放出を入れる。
    返り値: scipy の解(y = ln e1, ln e2, ln e3, chi)。322 の臨界の 1e-5 上で止める。"""
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    mu_s = mu_eV / sr.HBAR_EVS
    a = sr.alpha_of(M, mu_eV)
    s4 = (si.MPL_GEV / f_gev) ** 4
    eq = 1.0 / (sr.N0 * M ** 2)
    c2 = brentq(lambda c: g_rel(a, c, 3, 2, 2), 1e-4, 0.9999)

    def f(t, y):
        e = np.exp(np.clip(y[:3], -700.0, 5.0))                    # 解法が大きく踏み出したときのあふれを防ぐ
        chi = min(max(y[3], c2), 0.9999)
        g, k = coeffs(a, chi)
        k = {kk: v * s4 for kk, v in k.items()}
        e1, e2, e3 = e
        q1, q2, q3 = e1 + eq, e2 + eq, e3 + eq                    # 生成される準位には 1 個の量子を足す
        Ra = k["a"] * e1 ** 2 * q2
        Rb = k["b"] * e2 ** 2 * q1
        Rc = k["c"] * e2 ** 2 * q3
        Rd = k["d"] * e3 ** 2 * q2
        Re = k["e"] * e2 * e3 * q1
        Rf = k["f"] * e3 ** 2 * q1
        sr1 = g[0] * (q1 if g[0] > 0 else e1)
        sr2 = g[1] * (q2 if g[1] > 0 else e2)
        sr3 = g[2] * (q3 if g[2] > 0 else e3)
        d1 = sr1 - 2 * Ra + Rb + Re + Rf
        d2 = sr2 + Ra - 2 * Rb - 2 * Rc + Rd - Re
        d3 = sr3 + Rc - 2 * Rd - Re - 2 * Rf
        dchi = -(sr1 + 2 * sr2 + 4 * sr3)
        return mu_s * np.array([d1 / e1, d2 / e2, d3 / e3, dchi])

    stop = lambda t, y: y[3] - (c2 + 1e-5)   # noqa: E731
    stop.terminal = True
    y0 = [np.log(n * eq) for n in levels0] + [chi0]
    return solve_ivp(f, (0.0, T_yr * sr.YR), y0, method="LSODA", rtol=rtol, atol=1e-10,
                     t_eval=np.asarray(t_eval_yr) * sr.YR, events=stop)


# ---------- 固定点の外の準位(411、433)が育つか(E15) ----------
# BGLS Table V・VI(BH が吸収する過程は alpha^11 r+、放出は alpha^8):
K_P1 = 2.5e-8   # 211 x 411 -> 322 x BH
K_Q1 = 3.8e-9   # 322 x 411 -> 211 x inf
K_G = 9.1e-8    # 211 x 322 -> 433 x BH
K_E2 = 2.6e-9   # 322 x 433 -> 211 x inf
K_H = 1.1e-9    # 211 x 433 -> 544 x BH
K_K1 = 7.8e-10  # 433 x 544 -> 322 x inf
K_K2 = 1.9e-11  # 433 x 544 -> 211 x inf


def extra_rates(alpha, rpc, x1, x2, x3, g411, g433):
    """固定点 (x1, x2, x3) = (e211, e322, e544)(f によらない et)での、411 と 433 の 1 粒子あたりの増加率 / mu(配列で渡せる)。
    rpc = r+ x C_BH(BH が吸収する過程に掛ける)。
    411: g411 - (p1 + q1) e1 e2(BGLS 式 n11evol)。433: g433 + (g - e2) e1 e2 - h e1 e3 - k1 e2 e3 - k2 e1 e3(BGLS n33 の式に 544 の項)。"""
    a11, a8 = alpha ** 11, alpha ** 8
    l411 = g411 - (K_P1 * a11 * rpc + K_Q1 * a8) * x1 * x2
    l433 = (g433 + (K_G * a11 * rpc - K_E2 * a8) * x1 * x2 - K_H * a11 * rpc * x1 * x3
            - K_K1 * a8 * x2 * x3 - K_K2 * a8 * x1 * x3)
    return l411, l433


def extra_growth(alpha, chi, x, cbh=1.0):
    """extra_rates のスカラー版。411 は Leaver の連分数法の表(E15、out/rates_411.npz)。"""
    return extra_rates(alpha, sr.rplus(chi) * cbh, x[0], x[1], x[2], g_rel(alpha, chi, 4, 1, 1), g_rel(alpha, chi, 4, 3, 3))


# ---------- 照合用: 211・322・411・433・544 の五準位の ODE(E15) ----------
LEVELS5 = [(2, 1, 1), (3, 2, 2), (4, 1, 1), (4, 3, 3), (5, 4, 4)]
# (i, j, k, kappa, BH が吸収するか): i x j -> k x (BH か inf)。BGLS Table V・VI。
PROC5 = [(0, 0, 1, K_A, True), (1, 1, 0, K_B, False), (1, 1, 4, K_C, True), (4, 4, 1, K_D, False),
         (1, 4, 0, K_E, False), (4, 4, 0, K_F, False), (0, 2, 1, K_P1, True), (2, 2, 1, 9.8e-11, True),
         (1, 2, 0, K_Q1, False), (0, 1, 3, K_G, True), (1, 3, 0, K_E2, False), (0, 3, 4, K_H, True),
         (3, 4, 1, K_K1, False), (3, 4, 0, K_K2, False), (3, 3, 0, 9.2e-11, False)]


def sr_rates5(alpha, chi):
    return np.array([g_rel(alpha, chi, 2, 1, 1), g_rel(alpha, chi, 3, 2, 2), g_rel(alpha, chi, 4, 1, 1),
                     g_rel(alpha, chi, 4, 3, 3), sr.gamma_over_mu(alpha, chi, 5, 4, 4, sr.NORM_N)])


def evolve5(M, chi0, mu_eV, f_gev, T_yr, t_eval_yr, cbh=1.0, rtol=1e-8):
    """五準位の運動論を 1 個の量子から解く(BH の質量は固定)。cbh は BH が吸収する過程すべてに掛ける。"""
    from scipy.integrate import solve_ivp
    mu_s = mu_eV / sr.HBAR_EVS
    a = sr.alpha_of(M, mu_eV)
    s4 = (si.MPL_GEV / f_gev) ** 4
    eq = 1.0 / (sr.N0 * M ** 2)
    ms = np.array([l[2] for l in LEVELS5])

    def f(t, y):
        e = np.exp(np.clip(y[:5], -700.0, 5.0))
        chi = min(max(y[5], 0.01), 0.9999)
        rp = sr.rplus(chi)
        g = sr_rates5(a, chi)
        d = np.where(g > 0, g * (e + eq), g * e)
        dchi = -np.sum(ms * d)
        for i, j, k, kap, bh in PROC5:
            R = kap * (a ** 11 * rp * cbh if bh else a ** 8) * s4 * e[i] * e[j] * (e[k] + eq)
            d[i] -= R
            d[j] -= R
            d[k] += R
        return mu_s * np.concatenate([d / e, [dchi]])

    y0 = [np.log(eq)] * 5 + [chi0]
    return solve_ivp(f, (0.0, T_yr * sr.YR), y0, method="LSODA", rtol=rtol, atol=1e-10,
                     t_eval=np.asarray(t_eval_yr) * sr.YR)
