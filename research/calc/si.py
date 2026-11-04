"""自己相互作用のあるスカラーの雲: 211・322 の二準位の準平衡と BH のスピンの減り方(研究ノート E3)。

式は Baryakhtar, Galanis, Lasenby, Simon, arXiv:2011.11646(BGLS)。抜き書きは Refs/notes/theory_rates.md §4。
  占有数 eps = N / (G M^2)、M_pl = 1/sqrt(G) = 1.22e19 GeV。lambda = mu^2/f^2(引力)。
  gamma_BH  = kappa_BH  alpha^11 (M_pl/f)^4 r+ mu   (211x211 -> 322xBH、式 21、Table V: 4.3e-7)
  gamma_inf = kappa_inf alpha^8  (M_pl/f)^4    mu   (322x322 -> 211xinf、式 18、Table VI: 1.1e-8)
  準平衡(式 55a,b): eps211 = (2/sqrt3) sqrt(g_SR g_inf) / g_BH,  eps322 = g_BH eps211 / (2 g_inf)
  スピン(式 36): da*/dt = -g_SR eps211   (322 の超放射の寄与は無視できる)
超放射の率 g_SR は、BGLS の主要項 kappa_SR alpha^8 (a* - 2 alpha r+) の代わりに、
sr.py の 1908.10370 式 2.32 の 2 倍(占有数の率)を使う。
"""
import numpy as np
import sr

MPL_GEV = 1.220890e19
KAPPA_BH = 4.3e-7
KAPPA_INF = 1.1e-8


def g_bh(alpha, chi, mu_s, f_gev):
    return KAPPA_BH * alpha ** 11 * (MPL_GEV / f_gev) ** 4 * sr.rplus(chi) * mu_s * cbh(alpha)


def g_inf(alpha, mu_s, f_gev):
    return KAPPA_INF * alpha ** 8 * (MPL_GEV / f_gev) ** 4 * mu_s


REL = False   # True で E4 の相対論的な超放射の率を使う


def g_sr(alpha, chi, mu_s, leading=False):
    """211 の超放射(占有数の率、1/s)。leading=True で BGLS の主要項(kappa_SR = 0.04)。"""
    if leading:
        return 0.04 * alpha ** 8 * (chi - 2 * alpha * sr.rplus(chi)) * mu_s
    if REL:
        return sr.gamma_rel_over_mu(alpha, chi, 2, 1, 1, sr.NORM_N) * mu_s
    return sr.gamma_over_mu(alpha, chi, 2, 1, 1, sr.NORM_N) * mu_s


def eps_eq(alpha, chi, mu_s, f_gev, leading=False):
    """準平衡の占有数 (eps211, eps322)。BGLS 式 55a,b。"""
    gs = g_sr(alpha, chi, mu_s, leading)
    gb = g_bh(alpha, chi, mu_s, f_gev)
    gi = g_inf(alpha, mu_s, f_gev)
    e211 = 2.0 / np.sqrt(3.0) * np.sqrt(gs * gi) / gb
    return e211, gb * e211 / (2.0 * gi)


def chidot(alpha, chi, mu_s, f_gev, leading=False):
    """準平衡でのスピンの変化率 da*/dt [1/s](負)。"""
    gs = g_sr(alpha, chi, mu_s, leading)
    if gs <= 0:
        return 0.0
    e211, _ = eps_eq(alpha, chi, mu_s, f_gev, leading)
    return -gs * e211


def t_spindown(M, chi_i, chi_f, mu_eV, f_gev):
    """準平衡で chi_i から chi_f まで回転を失うのにかかる時間 [yr](M は固定)。

    t = int_{chi_f}^{chi_i} dchi / |dchi/dt|。dchi/dt は f^2 に比例するので、t は f^-2 に比例する。
    chi_f が臨界スピン以下なら inf。
    """
    mu_s = mu_eV / sr.HBAR_EVS
    a = sr.alpha_of(M, mu_eV)
    if sr.gamma_over_mu(a, chi_f, 2, 1, 1) <= 0:
        return np.inf
    xs = np.linspace(chi_f, chi_i, 2001)
    rate = np.array([-chidot(a, c, mu_s, f_gev) for c in xs])
    return np.trapezoid(1.0 / rate, xs) / sr.YR


def chi_threshold_211(M, mu_eV):
    """211 の臨界スピン(質量は固定、m Omega_H = E_211 を二分法で解く)。"""
    a = sr.alpha_of(M, mu_eV)
    lo, hi = 0.0, 0.999999
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if sr.omegaH(mid) - a * sr.omega_over_mu(a, 2, 1, 1, mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# ---------- gamma_BH の相対論的な補正(感度の検定用、E7) ----------
# 2412.03655 の Fig. 3 右(Rates_BH_n4.pdf、a~ = 0.95)の gamma^{322xBH}_{211x211} の gamma_R/gamma_NR を
# 図から目で読んだ値(読み取りの誤差は 20% 程度)。BH のスピン依存性は不明なので、感度を見るためだけに使う。
CBH_ALPHA = np.array([0.03, 0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20, 0.25, 0.30, 0.35])
CBH_VALUE = np.array([1.0, 1.3, 1.7, 2.5, 3.5, 4.4, 5.2, 5.2, 3.5, 2.2, 1.4])
CBH_ON = False


def cbh(alpha):
    if not CBH_ON:
        return 1.0
    return float(np.interp(alpha, CBH_ALPHA, CBH_VALUE))
