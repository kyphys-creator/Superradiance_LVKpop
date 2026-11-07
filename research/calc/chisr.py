"""集団解析用の速い chi_SR: 自己相互作用の準平衡で、年齢 T のあとの BH のスピン(研究ノート E10)。

準平衡(E3、BGLS 式 36 と 55a)では
  dchi/dt = mu (f/M_pl)^2 H(alpha, chi),
  H = -(2/sqrt3) g^(3/2) kappa_inf^(1/2) / (kappa_BH alpha^7 r+ C_BH(alpha)) x (1 - 2 chi alpha/3 の質量の補正)
で、g = gamma_SR/mu(E4 の相対論的な率、占有数の率)。f と mu は前の因子にしか入らないので、
  I(alpha, chi) = int_chi^{CHI_TOP} dchi'/|H(alpha, chi')|
を一度表にすれば、年齢 T のあとのスピンは
  I(alpha, chi_T) = I(alpha, chi_0) + tau,   tau = mu (f/M_pl)^2 T
を chi_T について解くだけで決まる(chi_T は臨界に近づくと I が発散するので、臨界を超えない)。
表は alpha と u = (chi - chi_c)/(CHI_TOP - chi_c) の格子で ln I を持つ。
有効な範囲: 準平衡が速く成り立つ f(E5: f <~ 1e15 GeV)、二準位の alpha <~ 0.15--0.2。自由場の側は E2 の半解析を使う。
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid

import si
import sr

CHI_TOP = 0.995                       # E4 の表の上端
ALPHAS = np.geomspace(0.02, 0.30, 141)
U = np.concatenate([np.geomspace(1e-6, 1e-2, 60, endpoint=False), np.linspace(1e-2, 1.0, 400)])


def H(alpha, chi, cbh=False, mass_corr=False):
    """dchi/dt / (mu (f/M_pl)^2)(負)。"""
    g = sr.gamma_rel_over_mu(alpha, chi, 2, 1, 1, sr.NORM_N)
    if g <= 0:
        return 0.0
    c = float(np.interp(alpha, si.CBH_ALPHA, si.CBH_VALUE)) if cbh else 1.0
    h = -(2 / np.sqrt(3)) * g ** 1.5 * np.sqrt(si.KAPPA_INF) / (si.KAPPA_BH * alpha ** 7 * sr.rplus(chi) * c)
    if mass_corr:
        h *= 1 - 2 * chi * alpha / 3
    return h


def build(cbh=False, mass_corr=False, ny=600):
    """表を作る。
    lnI[alpha, u]: ln I、chi = chi_c + u (CHI_TOP - chi_c)。
    uinv[alpha, y]: ln I = y となる u(逆引き。y は共通の格子 Y)。
    """
    lnI = np.empty((len(ALPHAS), len(U)))
    chic = np.empty(len(ALPHAS))
    for i, a in enumerate(ALPHAS):
        c0 = si.chi_threshold_211(1.0, a / sr.alpha_of(1.0, 1.0))   # M = 1 Msun、mu を alpha に合わせる
        chic[i] = c0
        chis = c0 + U * (CHI_TOP - c0)
        invH = np.array([1.0 / abs(H(a, c, cbh, mass_corr)) for c in chis])
        Irev = cumulative_trapezoid(invH[::-1], chis[::-1], initial=0.0)   # 上から積分(負の値)
        I = -Irev[::-1]
        I[-1] = I[-2] * 1e-3            # u = 1(chi = CHI_TOP)の I = 0 を有限の小さい値で置き換える
        lnI[i] = np.log(I)
    Y = np.linspace(lnI[:, -2].min(), lnI[:, 0].max(), ny)
    uinv = np.empty((len(ALPHAS), ny))
    for i in range(len(ALPHAS)):
        uinv[i] = np.interp(Y, lnI[i][::-1], U[::-1], left=1.0, right=U[0])
    return {"lnI": lnI, "chic": chic, "Y": Y, "uinv": uinv}


def chi_after(tab, alpha, chi0, tau):
    """年齢のあとのスピン(配列で渡せる)。tau = mu (f/M_pl)^2 T。
    alpha の格子の外と、chi0 が臨界以下の点では chi0 をそのまま返す(回転を失わない)。
    """
    from scipy.interpolate import RegularGridInterpolator
    alpha, chi0, tau = np.broadcast_arrays(np.asarray(alpha, float), np.asarray(chi0, float),
                                           np.asarray(tau, float))
    shape = alpha.shape
    alpha, chi0, tau = alpha.ravel(), chi0.ravel(), tau.ravel()
    la = np.log(ALPHAS)
    lal = np.clip(np.log(alpha), la[0], la[-1])
    chic = np.interp(lal, la, tab["chic"])
    u0 = np.clip((chi0 - chic) / (CHI_TOP - chic), U[0], 1.0)
    f_lnI = RegularGridInterpolator((la, U), tab["lnI"], bounds_error=False, fill_value=None)
    f_u = RegularGridInterpolator((la, tab["Y"]), tab["uinv"], bounds_error=False, fill_value=None)
    lnI0 = f_lnI(np.stack([lal, u0], -1))
    yT = np.logaddexp(lnI0, np.log(np.maximum(tau, 1e-300)))
    yT = np.clip(yT, tab["Y"][0], tab["Y"][-1])
    uT = np.clip(f_u(np.stack([lal, yT], -1)), 0.0, 1.0)
    out = chic + uT * (CHI_TOP - chic)
    keep = (chi0 <= chic) | (np.log(alpha) < la[0]) | (np.log(alpha) > la[-1])
    return np.where(keep, chi0, np.minimum(out, chi0)).reshape(shape)


def tau_of(mu_eV, f_gev, T_yr):
    """tau = mu (f/M_pl)^2 T(mu は 1/s に直す)。"""
    return mu_eV / sr.HBAR_EVS * (f_gev / si.MPL_GEV) ** 2 * T_yr * sr.YR
