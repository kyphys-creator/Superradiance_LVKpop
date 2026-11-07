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


# ---------- 自由場とつないだ chi_SR(研究ノート E11) ----------
def g_vec(alpha, chi, n, l, m):
    """超放射の率 / mu(占有数の率)を配列で。連分数法の表の範囲は表、外は解析式(sr.gamma_rel_over_mu と同じ規則)。"""
    alpha, chi = np.broadcast_arrays(np.asarray(alpha, float), np.asarray(chi, float))
    lo, hi = sr._REL_RANGE[(n, l, m)]
    out = sr.gamma_over_mu(alpha, chi, n, l, m, sr.NORM_N)
    inr = (alpha >= lo) & (alpha <= hi)
    if inr.any():
        F = sr._rel_table(n, l, m)(np.stack([np.log(alpha[inr]), chi[inr]], -1))
        out = np.array(out, float)
        out[inr] = sr.NORM_N * F * alpha[inr] ** (4 * l + 5) / alpha[inr]
    return out


def chic_vec(alpha, m, n, tab=None):
    """臨界スピン(質量固定)。211 は表の範囲で表から、ほかは m Omega_H = E_n(水素原子の主要項)の解析解。"""
    alpha = np.asarray(alpha, float)
    ap = alpha * (1 - alpha ** 2 / (2 * n ** 2))          # alpha E/mu
    out = 4 * m * ap / (m ** 2 + 4 * ap ** 2)
    if tab is not None and m == 1:
        la = np.log(ALPHAS)
        inr = (alpha >= ALPHAS[0]) & (alpha <= ALPHAS[-1])
        out = np.where(inr, np.interp(np.log(np.clip(alpha, ALPHAS[0], ALPHAS[-1])), la, tab["chic"]), out)
    return out


F_FREE_322 = 1e16   # 322 から先に育つ BH を自由場として扱う f の下限 [GeV](それより小さい f は扱わない、印 3)


def _sat_state(M, chi, mu_eV, m, n, iters=5):
    """自由場で準位 (n, m=l) が飽和したあとの (M_f, chi_f)。BH の質量の減少を入れる(sr.chi_threshold を配列にしたもの)。
    角運動量 dJ を渡すと質量は dM = (E/mu) alpha_sun dJ / m 減る。臨界は m Omega_H = E_n(水素原子の主要項)。
    """
    a_sun = sr.alpha_of(1.0, mu_eV)
    J0 = chi * M ** 2
    Mf = M.copy()
    for _ in range(iters):
        a = a_sun * Mf
        ap = a * (1 - a ** 2 / (2 * n ** 2))
        cf = 4 * m * ap / (m ** 2 + 4 * ap ** 2)
        dJ = J0 - cf * Mf ** 2
        Mf = M - a_sun * (1 - a ** 2 / (2 * n ** 2)) * dJ / m
    return Mf, cf


def _free_branch(M, chi0, a, mu_eV, mu_s, T, lnGM2, c1_fixed, g1):
    """自由場: 211 が 1 個の量子から飽和(E2 の半解析と同じ数え方)、322 は t = 0 から e-fold を積む。"""
    M1, c1 = _sat_state(M, chi0, mu_eV, 1, 2)
    c1 = np.minimum(c1, chi0)
    lnN1 = np.log(np.maximum(chi0 * M ** 2 - c1 * M1 ** 2, 1e-300)) + np.log(sr.N0)
    tsat1 = lnN1 / g1
    M2, c2 = _sat_state(M1, c1, mu_eV, 2, 3)
    lnN2 = np.log(np.maximum((c1 * M1 ** 2 - c2 * M2 ** 2) / 2, 1e-300)) + np.log(sr.N0)
    g2a = np.maximum(g_vec(a, np.minimum(chi0, 0.995), 3, 2, 2), 0) * mu_s          # 211 が飽和するまで
    g2b = np.maximum(g_vec(sr.alpha_of(M1, mu_eV), c1, 3, 2, 2), 0) * mu_s          # 211 が飽和したあと
    with np.errstate(divide="ignore", invalid="ignore"):
        t2_first = np.where(g2a > 0, lnN2 / g2a, np.inf)                              # 322 が先に飽和する場合
        t2 = np.where(g2b > 0, tsat1 + np.maximum(lnN2 - g2a * tsat1, 0) / g2b, np.inf)
    res = np.where(T >= tsat1, c1, chi0)
    res = np.where((t2_first < tsat1) & (T >= t2_first), c2, res)
    res = np.where((t2_first >= tsat1) & (c1 > c2) & (T >= t2), c2, res)
    return res


def chi_sr(tab, M, chi0, T_yr, mu_eV, f_gev):
    """年齢 T のあとの BH のスピンと、どの領域で計算したかの印(配列で渡せる)。

    印: 0 = 超放射しない(chi0 が 322 の臨界以下)
        1 = 自己相互作用の準平衡(二準位、chisr の表)
        2 = 自由場(211 が飽和して臨界へ。さらに時間があれば 322 も。BH の質量の減少を入れる)
        3 = この模型の外(準平衡で alpha > 0.30、または 322 から先に育つ BH で f < F_FREE_322)。chi0 を返す
    領域の判定(BGLS の B|C の境界、式 53・56): 準平衡の占有数 eps_eq と行き過ぎの占有数 eps_th が、
    どちらも自由場で飽和する量 eps_max = chi0 - chi_c より小さければ準平衡。
    準平衡では、1 個の量子から eps_eq まで育つ時間 t_grow を年齢から引き、行き過ぎで失う eps_th を初期スピンから引く。
    """
    M, chi0, T_yr, mu_eV, f_gev = np.broadcast_arrays(*[np.asarray(x, float) for x in (M, chi0, T_yr, mu_eV, f_gev)])
    shape = M.shape
    M, chi0, T_yr, mu_eV, f_gev = (x.ravel() for x in (M, chi0, T_yr, mu_eV, f_gev))
    a = sr.alpha_of(M, mu_eV)
    mu_s = mu_eV / sr.HBAR_EVS
    T = T_yr * sr.YR
    lnGM2 = np.log(sr.N0 * M ** 2)
    c1 = chic_vec(a, 1, 2, tab)
    c2 = chic_vec(a, 2, 3)
    out = chi0.copy()
    reg = np.zeros(len(M), int)

    sr211 = chi0 > c1
    # ---- 211 が超放射する BH ----
    g1 = np.where(sr211, g_vec(a, np.minimum(chi0, 0.995), 2, 1, 1), 0.0) * mu_s       # 1/s
    ginf = si.KAPPA_INF * a ** 8 * (si.MPL_GEV / f_gev) ** 4 * mu_s
    gbh = si.KAPPA_BH * a ** 11 * (si.MPL_GEV / f_gev) ** 4 * sr.rplus(np.minimum(chi0, 0.999)) * mu_s
    emax = chi0 - c1
    with np.errstate(divide="ignore", invalid="ignore"):
        eeq = 2 / np.sqrt(3) * np.sqrt(g1 * ginf) / gbh
        eth = np.sqrt(2 * g1 * lnGM2 / gbh)
    si_reg = sr211 & (g1 > 0) & (eeq < emax) & (eth < emax)
    inval = si_reg & ((a > ALPHAS[-1]) | (a < ALPHAS[0]))
    reg[inval] = 3
    q = si_reg & ~inval
    if q.any():
        tgrow = (np.log(np.maximum(eeq[q], 1e-300)) + lnGM2[q]) / g1[q]
        Teff = np.maximum(T[q] - tgrow, 0.0)
        tau = mu_s[q] * (f_gev[q] / si.MPL_GEV) ** 2 * Teff
        c0eff = np.maximum(chi0[q] - eth[q], c1[q])
        out[q] = np.where(Teff > 0, chi_after(tab, a[q], c0eff, tau), chi0[q])
        reg[q] = 1
    fr = sr211 & ~si_reg & (g1 > 0)
    if fr.any():
        out[fr] = _free_branch(M[fr], chi0[fr], a[fr], mu_eV[fr], mu_s[fr], T[fr], lnGM2[fr], c1[fr], g1[fr])
        reg[fr] = 2
    # ---- 211 は超放射しないが 322 はする BH(重い BH) ----
    only322 = (~sr211) & (chi0 > c2)
    fr2 = only322 & (f_gev >= F_FREE_322)
    if fr2.any():
        g2 = g_vec(a[fr2], np.minimum(chi0[fr2], 0.995), 3, 2, 2) * mu_s[fr2]
        with np.errstate(divide="ignore", invalid="ignore"):
            tsat2 = (np.log(np.maximum((chi0[fr2] - c2[fr2]) / 2, 1e-300)) + lnGM2[fr2]) / g2
        out[fr2] = np.where((g2 > 0) & (T[fr2] >= tsat2), c2[fr2], chi0[fr2])
        reg[fr2] = 2
    reg[only322 & (f_gev < F_FREE_322)] = 3
    return out.reshape(shape), reg.reshape(shape)
