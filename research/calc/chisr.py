"""集団解析用の速い chi_SR: 自己相互作用のある超放射で、年齢 T のあとの BH のスピン(研究ノート E10、E11、E12)。

準平衡: BGLS(2011.11646)の式 34--36 の右辺を 0 とおいた二準位(211・322)の固定点。
  u = gamma_BH eps211^2 は 3 u^2 + (2 g322 - 4 g_inf g211 / g_BH) u - g322^2 = 0 の正の根、
  eps322 = (u + g322) / (2 g_inf eps211)、
  dchi/dt = -g211 eps211 - 2 g322 eps322 = -3 g_inf eps211 eps322^2。
g322 -> 0 で E3 の式(BGLS 式 55)に戻る。E12 で 322 の超放射を入れ、211 の臨界より下(211 は吸収される)でも
322 x 322 -> 211 x inf が 211 を補充して、322 の臨界まで下がるようにした。211 の臨界より下では
  dchi/dt ~ -(3/2) g322 sqrt(|g211| / g_inf)
で、gamma_BH によらない。どちらの側でも dchi/dt = mu (f/M_pl)^2 H(alpha, chi) で、f と mu は前の因子にしか入らないので、
  I(alpha, chi) = int_chi^{CHI_TOP} dchi' / |H(alpha, chi')|
を一度表にすれば、年齢 T のあとのスピンは I(alpha, chi_T) = I(alpha, chi_0) + tau(tau = mu (f/M_pl)^2 T)を解くだけで決まる。
表は alpha と w の格子で ln I を持つ。w は 322 の臨界(w = 0)、211 の臨界(w = 1)、CHI_TOP(w = 2)を結ぶ区分的な座標。
211 の臨界のすぐ上と下(幅 1e-8 ほど)に、回転を失う速さが最も遅い「ボトルネック」があるので、w は両側に対数で詰める。
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq

import si
import sr

CHI_TOP = 0.995                       # E4 の表の上端
ALPHAS = np.geomspace(0.01, 0.30, 161)
ALPHA_2LEVEL = 0.15                   # これより上の準平衡は印 4(433・411 が育つ: WM 2412.03655、BGLS IV C)
ALPHA_MAX = 0.5                       # これより上は印 3(211 は超放射できず、m >= 3 が効く)
F_FP_MAX = 1e17                       # 自由場の側で、211 の飽和のあと固定点の降下を使う f の上限 [GeV](E12)
KAPPA_BH_544 = 1.9e-9                 # BGLS Table V: Gamma^{544 x BH}_{322 x 322} / (alpha^11 (M_pl/f)^4 r+ mu)

_E = np.geomspace(1e-13, 1e-2, 110, endpoint=False)
W_LO = np.concatenate([_E, np.linspace(1e-2, 1 - 1e-2, 300, endpoint=False), 1 - _E[::-1]])   # (0, 1): 322 の臨界から 211 の臨界
W_HI = 1 + np.concatenate([[0.0], _E, np.linspace(1e-2, 1.0, 400)])                           # [1, 2]: 211 の臨界から CHI_TOP
W = np.concatenate([W_LO, W_HI])
U = W                                 # 旧名(E10 の e10_gwonly.py が参照する)


# ---------- 率 ----------
def g_vec(alpha, chi, n, l, m):
    """超放射の率 / mu(占有数の率)を配列で。連分数法の表の範囲は表、外は解析式(sr.gamma_rel_over_mu と同じ規則)。"""
    alpha, chi = np.broadcast_arrays(np.asarray(alpha, float), np.asarray(chi, float))
    lo, hi = sr._REL_RANGE[(n, l, m)]
    out = np.array(sr.gamma_over_mu(alpha, chi, n, l, m, sr.NORM_N), float)
    inr = (alpha >= lo) & (alpha <= hi)
    if inr.any():
        F = sr._rel_table(n, l, m)(np.stack([np.log(alpha[inr]), chi[inr]], -1))
        out[inr] = sr.NORM_N * F * alpha[inr] ** (4 * l + 5) / alpha[inr]
    return out


def H(alpha, chi, cbh=False, mass_corr=False):
    """dchi/dt / (mu (f/M_pl)^2)(負、配列で渡せる)。二準位の固定点で、322 の超放射を含む(E12)。"""
    alpha, chi = np.broadcast_arrays(np.asarray(alpha, float), np.asarray(chi, float))
    g1 = g_vec(alpha, chi, 2, 1, 1)
    g2 = g_vec(alpha, chi, 3, 2, 2)
    c = np.interp(alpha, si.CBH_ALPHA, si.CBH_VALUE) if cbh else 1.0
    kb = si.KAPPA_BH * alpha ** 11 * sr.rplus(chi) * c         # gamma_BH / (mu (M_pl/f)^4)
    ki = si.KAPPA_INF * alpha ** 8                              # gamma_inf / (mu (M_pl/f)^4)
    B = 2 * g2 - 4 * ki * g1 / kb
    D = np.sqrt(B * B + 12 * g2 * g2)
    with np.errstate(divide="ignore", invalid="ignore"):
        u = np.where(B > 0, 2 * g2 * g2 / (B + D), (D - B) / 6)   # 桁落ちしない形の正の根
        e1 = np.sqrt(u / kb)                                      # eps211 / (f/M_pl)^2
        h = -3 * (g2 + u) ** 2 / (4 * ki * e1)
    h = np.where((g2 > 0) & (u > 0), h, 0.0)
    if mass_corr:
        h = h * (1 - 2 * chi * alpha / 3)
    return h


def _zero(alpha, n, l, m):
    """準位 (n, l, m) の率が 0 になるスピン(質量固定、相対論的な率の表と同じもの)。"""
    return brentq(lambda c: float(g_vec(alpha, c, n, l, m)), 1e-4, 0.9999, xtol=1e-14)


def chi_of_w(w, c1, c2):
    return np.where(w >= 1, c1 + (w - 1) * (CHI_TOP - c1), c2 + w * (c1 - c2))


def w_of_chi(chi, c1, c2):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(chi >= c1, 1 + (chi - c1) / (CHI_TOP - c1), (chi - c2) / (c1 - c2))


# ---------- 表 ----------
def build(cbh=False, mass_corr=False, ny=1200):
    """表を作る。
    lnI[alpha, w]: ln I、chi = chi_of_w(w, chic, chic2)。
    winv[alpha, y]: ln I = y となる w(逆引き。y は共通の格子 Y)。
    chic, chic2: 211 と 322 の臨界スピン(相対論的な率の 0)。
    """
    nA = len(ALPHAS)
    lnI = np.empty((nA, len(W)))
    c1s, c2s = np.empty(nA), np.empty(nA)
    for i, a in enumerate(ALPHAS):
        c1, c2 = _zero(a, 2, 1, 1), _zero(a, 3, 2, 2)
        c1s[i], c2s[i] = c1, c2
        chis = chi_of_w(W, c1, c2)
        invH = 1.0 / np.abs(H(a, chis, cbh, mass_corr))
        I = -cumulative_trapezoid(invH[::-1], chis[::-1], initial=0.0)[::-1]   # 上から積分
        I[-1] = I[-2] * 1e-3            # w = 2(chi = CHI_TOP)の I = 0 を有限の小さい値で置き換える
        lnI[i] = np.log(I)
    Y = np.linspace(lnI[:, -2].min(), lnI[:, 0].max(), ny)
    winv = np.empty((nA, ny))
    for i in range(nA):
        winv[i] = np.interp(Y, lnI[i][::-1], W[::-1], left=W[-1], right=W[0])
    return {"lnI": lnI, "chic": c1s, "chic2": c2s, "Y": Y, "winv": winv}


def _interps(tab):
    if "_f_lnI" not in tab:
        la = np.log(ALPHAS)
        tab["_f_lnI"] = RegularGridInterpolator((la, W), tab["lnI"], bounds_error=False, fill_value=None)
        tab["_f_w"] = RegularGridInterpolator((la, tab["Y"]), tab["winv"], bounds_error=False, fill_value=None)
    return tab["_f_lnI"], tab["_f_w"]


def chi_after(tab, alpha, chi0, tau):
    """準平衡(固定点)で年齢のあとのスピン(配列で渡せる)。tau = mu (f/M_pl)^2 T。
    322 の臨界まで下がりうる。alpha の格子の外と、chi0 が 322 の臨界以下の点では chi0 をそのまま返す。
    """
    alpha, chi0, tau = np.broadcast_arrays(np.asarray(alpha, float), np.asarray(chi0, float),
                                           np.asarray(tau, float))
    shape = alpha.shape
    alpha, chi0, tau = alpha.ravel(), chi0.ravel(), tau.ravel()
    f_lnI, f_w = _interps(tab)
    la = np.log(ALPHAS)
    lal = np.clip(np.log(alpha), la[0], la[-1])
    c1 = np.interp(lal, la, tab["chic"])
    c2 = np.interp(lal, la, tab["chic2"])
    w0 = np.clip(w_of_chi(chi0, c1, c2), W[0], W[-1])
    lnI0 = f_lnI(np.stack([lal, w0], -1))
    yT = np.logaddexp(lnI0, np.log(np.maximum(tau, 1e-300)))
    yT = np.clip(yT, tab["Y"][0], tab["Y"][-1])
    wT = np.clip(f_w(np.stack([lal, yT], -1)), W[0], W[-1])
    out = chi_of_w(wT, c1, c2)
    keep = (chi0 <= c2) | (np.log(alpha) < la[0]) | (np.log(alpha) > la[-1])
    return np.where(keep, chi0, np.minimum(out, chi0)).reshape(shape)


def tau_of(mu_eV, f_gev, T_yr):
    """tau = mu (f/M_pl)^2 T(mu は 1/s に直す)。"""
    return mu_eV / sr.HBAR_EVS * (f_gev / si.MPL_GEV) ** 2 * T_yr * sr.YR


# ---------- 臨界スピン ----------
def chic_vec(alpha, m, n, tab=None):
    """臨界スピン(質量固定)。表の範囲では表(相対論的な率の 0)、外は m Omega_H = E_n(水素原子の主要項)の解析解。
    alpha_E / m >= 1/2 では、その m の準位はどのスピンでも超放射しないので 1 を返す。
    """
    alpha = np.asarray(alpha, float)
    ap = alpha * (1 - alpha ** 2 / (2 * n ** 2))          # alpha E/mu
    out = np.where(ap / m < 0.5, 4 * m * ap / (m ** 2 + 4 * ap ** 2), 1.0)
    if tab is not None and m in (1, 2):
        key = "chic" if m == 1 else "chic2"
        la = np.log(ALPHAS)
        inr = (alpha >= ALPHAS[0]) & (alpha <= ALPHAS[-1])
        out = np.where(inr, np.interp(np.log(np.clip(alpha, ALPHAS[0], ALPHAS[-1])), la, tab[key]), out)
    return out


# ---------- 自由場 ----------
def _sat_state(M, chi, mu_eV, m, n):
    """自由場で準位 (n, m=l) が飽和したあとの (M_f, chi_f)。BH の質量の減少を入れる。
    角運動量 dJ = J0 - chi_c(M_f) M_f^2 を渡すと、質量は dM = (E/mu) alpha_sun dJ / m 減る。臨界は m Omega_H = E_n
    (水素原子の主要項)。M_f を二分法で自己無撞着に解く(E12: 以前は 5 回の反復で打ち切っていた)。
    """
    a_sun = sr.alpha_of(1.0, mu_eV)
    J0 = chi * M ** 2

    def crit(Mf):
        a = a_sun * Mf
        ap = a * (1 - a ** 2 / (2 * n ** 2))
        return a, np.where(ap / m < 0.5, 4 * m * ap / (m ** 2 + 4 * ap ** 2), 1.0)

    def F(Mf):                                  # 0 になる Mf を探す(Mf = M で >= 0、Mf = M/2 で < 0)
        a, cf = crit(Mf)
        return Mf - M + a_sun * (1 - a ** 2 / (2 * n ** 2)) * (J0 - cf * Mf ** 2) / m

    lo, hi = 0.5 * M, M.copy()
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        neg = F(mid) < 0
        lo = np.where(neg, mid, lo)
        hi = np.where(neg, hi, mid)
    Mf = 0.5 * (lo + hi)
    return Mf, crit(Mf)[1]


def _free_branch(tab, M, chi0, a, mu_eV, mu_s, T, f_gev, g1):
    """自由場: 211 が 1 個の量子から飽和(E2 の半解析と同じ数え方)。そのあと
    f <= F_FP_MAX: 211 x 211 -> 322 x BH が 322 をすぐに育てるので、固定点の降下で 322 の臨界へ(E12)。
    f >  F_FP_MAX: 322 を t = 0 から 1 個の量子として e-fold を積む(E11 のまま。ODE より遅いことがある、E12)。
    """
    M1, c1 = _sat_state(M, chi0, mu_eV, 1, 2)
    c1 = np.minimum(c1, chi0)
    lnN1 = np.log(np.maximum(chi0 * M ** 2 - c1 * M1 ** 2, 1e-300)) + np.log(sr.N0)
    tsat1 = lnN1 / g1
    M2, c2 = _sat_state(M1, c1, mu_eV, 2, 3)
    lnN2 = np.log(np.maximum((c1 * M1 ** 2 - c2 * M2 ** 2) / 2, 1e-300)) + np.log(sr.N0)
    g2a = np.maximum(g_vec(a, np.minimum(chi0, 0.995), 3, 2, 2), 0) * mu_s          # 211 が飽和するまで
    a1 = sr.alpha_of(M1, mu_eV)
    g2b = np.maximum(g_vec(a1, c1, 3, 2, 2), 0) * mu_s                             # 211 が飽和したあと
    with np.errstate(divide="ignore", invalid="ignore"):
        t2_first = np.where(g2a > 0, lnN2 / g2a, np.inf)                              # 322 が先に飽和する場合
        t2 = np.where(g2b > 0, tsat1 + np.maximum(lnN2 - g2a * tsat1, 0) / g2b, np.inf)
    res = np.where(T >= tsat1, c1, chi0)
    res = np.where((t2_first < tsat1) & (T >= t2_first), c2, res)
    slow = f_gev > F_FP_MAX
    res = np.where(slow & (t2_first >= tsat1) & (c1 > c2) & (T >= t2), c2, res)
    fp = (~slow) & (t2_first >= tsat1) & (T > tsat1)
    if fp.any():
        tau1 = mu_s[fp] * (f_gev[fp] / si.MPL_GEV) ** 2 * (T[fp] - tsat1[fp])
        start = chic_vec(a1[fp], 1, 2, tab)
        drop = start - chi_after(tab, a1[fp], start + 1e-15, tau1)                   # 表の 211 の臨界からの降下
        res[fp] = np.maximum(c1[fp] - drop, c2[fp])
    return res


# ---------- 本体 ----------
def _check(M, chi0, T_yr, mu_eV, f_gev):
    for name, x in (("M", M), ("chi0", chi0), ("T_yr", T_yr), ("mu_eV", mu_eV), ("f_gev", f_gev)):
        if not np.all(np.isfinite(x)):
            raise ValueError(f"{name} に有限でない値がある")
    if np.any((chi0 < 0) | (chi0 >= 1)):
        raise ValueError("chi0 は 0 <= chi0 < 1")
    if np.any(M <= 0) or np.any(mu_eV <= 0) or np.any(f_gev <= 0) or np.any(T_yr < 0):
        raise ValueError("M, mu_eV, f_gev は正、T_yr は 0 以上")


def chi_sr(tab, M, chi0, T_yr, mu_eV, f_gev):
    """年齢 T のあとの BH のスピンと、どの領域で計算したかの印(配列で渡せる)。

    印: 0 = 超放射しない(chi0 が 211 と 322 の臨界以下。m >= 3 は入れていない)
        1 = 自己相互作用の準平衡(二準位の固定点、alpha <= 0.15)。322 の臨界まで下がりうる
        2 = 自由場(211 が飽和して臨界へ。そのあと f <= 1e17 GeV は固定点の降下、それより上は 322 を 1 個の量子から)
            または 322 から先に育つ BH で、自己相互作用が飽和までに効かない場合
        3 = この模型の外。chi0 を返す: 322 から先に育ち(chi_c322 < chi0 <= chi_c211)自己相互作用が効く BH、
            準平衡で alpha が表の外(alpha < 0.01 または > 0.30)、alpha > 0.5
        4 = 準平衡だが 0.15 < alpha <= 0.30。433・411 などが育つので二準位の近似は粗い。値は計算して返す
    領域の判定(BGLS の B|C の境界、式 53・56): 準平衡の占有数 eps_eq と行き過ぎの占有数 eps_th が、
    どちらも自由場で飽和する量 eps_max = chi0 - chi_c より小さければ準平衡。
    準平衡では、1 個の量子から eps_eq まで育つ時間 t_grow を年齢から引き、行き過ぎで失う eps_th を初期スピンから引く。
    """
    M, chi0, T_yr, mu_eV, f_gev = np.broadcast_arrays(*[np.asarray(x, float) for x in (M, chi0, T_yr, mu_eV, f_gev)])
    shape = M.shape
    M, chi0, T_yr, mu_eV, f_gev = (x.ravel() for x in (M, chi0, T_yr, mu_eV, f_gev))
    _check(M, chi0, T_yr, mu_eV, f_gev)
    a = sr.alpha_of(M, mu_eV)
    mu_s = mu_eV / sr.HBAR_EVS
    T = T_yr * sr.YR
    lnGM2 = np.log(sr.N0 * M ** 2)
    c1 = chic_vec(a, 1, 2, tab)
    c2 = chic_vec(a, 2, 3, tab)
    big = a > ALPHA_MAX
    out = chi0.copy()
    reg = np.zeros(len(M), int)
    reg[big] = 3

    sr211 = (chi0 > c1) & ~big
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
        reg[q] = np.where(a[q] > ALPHA_2LEVEL, 4, 1)
    fr = sr211 & ~si_reg & (g1 > 0)
    if fr.any():
        out[fr] = _free_branch(tab, M[fr], chi0[fr], a[fr], mu_eV[fr], mu_s[fr], T[fr], f_gev[fr], g1[fr])
        reg[fr] = 2

    # ---- 211 は超放射しないが 322 はする BH(chi_c322 < chi0 <= chi_c211) ----
    # 322 が 1 個の量子から自由に飽和するのは、自己相互作用の二つの吸い込みが、飽和の量 eps_max' = (chi0 - chi_c322)/2
    # に届くまで効かないとき: 322 x 322 -> 544 x BH の行き過ぎ eps_th' = sqrt(2 g322 ln GM^2 / g_BH') と、
    # 吸収される 211 が 322 x 322 -> 211 x inf で補充され始める占有数 sqrt(|g211| / g_inf)(E12)。
    only322 = (~sr211) & (chi0 > c2) & ~big
    if only322.any():
        i3 = np.flatnonzero(only322)
        c0 = np.minimum(chi0[i3], 0.995)
        g2 = g_vec(a[i3], c0, 3, 2, 2) * mu_s[i3]
        g1d = np.abs(g_vec(a[i3], c0, 2, 1, 1)) * mu_s[i3]
        gbh2 = KAPPA_BH_544 * a[i3] ** 11 * (si.MPL_GEV / f_gev[i3]) ** 4 * sr.rplus(c0) * mu_s[i3]
        emax2 = (chi0[i3] - c2[i3]) / 2
        with np.errstate(divide="ignore", invalid="ignore"):
            eth2 = np.sqrt(2 * g2 * lnGM2[i3] / gbh2)
            epump = np.sqrt(g1d / ginf[i3])
            tsat2 = (np.log(np.maximum(emax2, 1e-300)) + lnGM2[i3]) / g2
        free2 = (g2 > 0) & (eth2 >= emax2) & (epump >= emax2)
        out[i3] = np.where(free2 & (T[i3] >= tsat2), c2[i3], chi0[i3])
        reg[i3] = np.where(free2, 2, 3)
    return out.reshape(shape), reg.reshape(shape)
