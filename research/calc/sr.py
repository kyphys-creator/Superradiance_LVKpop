"""自由場のスカラーの超放射: 準位のエネルギー、成長率、BH のスピンの減り方(研究ノート E2)。

単位: G = c = hbar = 1 の式を、BH の質量は太陽質量、時間は秒・年で使う。
  alpha = G M mu / (hbar c)   (重力の微細構造定数)
  N0 = G Msun^2 / (hbar c)    (太陽質量の BH の角運動量 G Msun^2/c を hbar で数えた数)
雲の占有数は n = N / N0 で持つ。BH の角運動量は J~ = chi (M/Msun)^2 (単位 G Msun^2 / c)。

式の典拠(式番号は Refs/notes/theory_rates.md に原文と照合して書く):
  成長率 Gamma_nlm   Baumann, Chia, Porto, Stout, arXiv:1908.10370 (Detweiler 型の解析式)
  エネルギー omega_nlm 同上(alpha^5 まで)
"""
import numpy as np
from math import factorial
from scipy.integrate import solve_ivp

# ---------- 定数 ----------
HBAR_EVS = 6.582119569e-16      # eV s
TSUN = 4.925490947e-6           # G Msun / c^3 [s]
MSUN_KG = 1.98840987e30
MPL_KG = 2.176434e-8            # sqrt(hbar c / G)
N0 = (MSUN_KG / MPL_KG) ** 2    # G Msun^2 / (hbar c) ~ 8.35e75
YR = 3.15576e7                  # s


def alpha_of(M, mu_eV):
    """alpha = G M mu / (hbar c)。M は太陽質量、mu は eV。"""
    return TSUN * M * mu_eV / HBAR_EVS


def rplus(chi):
    """r_+ / (G M / c^2)。"""
    return 1.0 + np.sqrt(1.0 - chi ** 2)


def omegaH(chi):
    """地平線の角速度 Omega_H に G M / c^3 をかけたもの。"""
    return chi / (2.0 * rplus(chi))


def chi_crit(alpha, m):
    """m Omega_H = omega(~mu) となる臨界スピン(omega = mu の近似)。"""
    return 4.0 * m * alpha / (m ** 2 + 4.0 * alpha ** 2)


def omega_over_mu(alpha, n, l, m, chi):
    """omega_nlm / mu。alpha^5 まで(1908.10370)。"""
    f_nl = -6.0 / (2 * l + 1) + 2.0 / n
    h_l = 16.0 / (2 * l * (2 * l + 1) * (2 * l + 2)) if l > 0 else 0.0
    return (1.0 - alpha ** 2 / (2 * n ** 2) - alpha ** 4 / (8 * n ** 4)
            + f_nl * alpha ** 4 / n ** 3 + h_l * chi * m * alpha ** 5 / n ** 3)


def C_nl(n, l):
    return (2.0 ** (4 * l + 1) * factorial(n + l) / (n ** (2 * l + 4) * factorial(n - l - 1))
            * (factorial(l) / (factorial(2 * l) * factorial(2 * l + 1))) ** 2)


def gamma_over_mu(alpha, chi, n, l, m, norm=1.0):
    """超放射の率 Gamma_nlm / mu。norm=1 で 1908.10370 の式 (eqn:ScalarRate) の Gamma(振幅の率)。

    Gamma = 2 r+ C_nl g_lm (m Omega_H - omega) alpha^(4l+5)   (G M = 1 の単位)
    g_lm = prod_{k=1}^{l} [k^2 (1 - chi^2) + (chi m - 2 r+ omega)^2]
    1908.10370 は omega_nlm = E + i Gamma と書く(Bohr energy の式)ので、場の振幅が e^{Gamma t}、
    占有数 N が e^{2 Gamma t} で育つ。占有数の式には norm=2 を渡す(NORM_N)。
    小さい alpha では Gamma -> mu alpha^8 (chi - 2 alpha r+) / 48。
    """
    rp = rplus(chi)
    w = alpha * omega_over_mu(alpha, n, l, m, chi)      # G M omega
    g = 1.0
    for k in range(1, l + 1):
        g = g * (k ** 2 * (1.0 - chi ** 2) + (chi * m - 2.0 * rp * w) ** 2)
    gam_M = 2.0 * rp * C_nl(n, l) * g * (m * omegaH(chi) - w) * alpha ** (4 * l + 5)
    return norm * gam_M / alpha                          # G M Gamma / alpha = Gamma / mu


NORM_N = 2.0   # 占有数の成長率 = 2 x 振幅の率


# ---------- BH と雲の連立(自由場) ----------
def evolve_free(M0, chi0, mu_eV, levels, t_end_yr, norm=NORM_N, n_init=1.0, rtol=1e-9):
    """自由場の雲の成長と BH のスピンの減少を解く。

    levels: [(n, l, m), ...]。各準位は 1 個の量子(n_init 個)から育つ。
    返り値: solve_ivp の解(y = [ln n_1, ..., ln n_k, M, Jt]、t は秒)。
    重力波による雲の減衰は入れない(スピンの減り方には効かない)。
    """
    mu_s = mu_eV / HBAR_EVS
    k = len(levels)

    def rhs(t, y):
        lnn = y[:k]
        M, Jt = y[k], y[k + 1]
        chi = min(max(Jt / M ** 2, 0.0), 0.999999)
        a = alpha_of(M, mu_eV)
        dlnn = np.empty(k)
        dJ = 0.0
        dM = 0.0
        for i, (n, l, m) in enumerate(levels):
            G = gamma_over_mu(a, chi, n, l, m, norm) * mu_s
            dlnn[i] = G
            ni = np.exp(lnn[i])
            dJ -= m * G * ni
            dM -= (a / M) * omega_over_mu(a, n, l, m, chi) * G * ni
        return np.concatenate([dlnn, [dM, dJ]])

    y0 = np.concatenate([np.full(k, np.log(n_init / N0)), [M0, chi0 * M0 ** 2]])
    return solve_ivp(rhs, (0.0, t_end_yr * YR), y0, method="Radau",
                     rtol=rtol, atol=1e-12, dense_output=True)


def chi_at(sol, t_yr, k):
    """解 sol の時刻 t_yr(年)での BH のスピン。"""
    y = sol.sol(np.atleast_1d(t_yr) * YR)
    return y[k + 1] / y[k] ** 2


# ---------- 半解析: 準位が飽和するたびに区切る(自由場) ----------
def chi_threshold(M, chi, mu_eV, n, l, m):
    """準位 (n,l,m) が飽和したあとの (M, chi)。

    m Omega_H(chi) = G M omega(alpha, chi) を満たす chi まで回転を渡す。
    渡した角運動量 dJ と質量 dM は dM c^2 = (omega/m) dJ で結ばれる。
    """
    a_sun = alpha_of(1.0, mu_eV)            # alpha / (M/Msun)
    J0 = chi * M ** 2
    Mf, chif = M, chi
    for _ in range(4):
        a = a_sun * Mf
        lo, hi = 0.0, chif
        for _ in range(60):                 # 二分法で m Omega_H = G M omega を解く
            mid = 0.5 * (lo + hi)
            if m * omegaH(mid) - a * omega_over_mu(a, n, l, m, mid) > 0:
                hi = mid
            else:
                lo = mid
        chif = 0.5 * (lo + hi)
        # 角運動量の受け渡しで質量が減る: dM = (alpha_sun (omega/mu) / m) dJ
        w = omega_over_mu(a, n, l, m, chif)
        # J0 - Mf^2 chif = dJ,  M - Mf = a_sun w dJ / m を Mf について解く(1 回の更新)
        dJ = J0 - Mf ** 2 * chif
        Mf = M - a_sun * w * dJ / m
    return Mf, chif


def chi_final_seq(M0, chi0, mu_eV, t_yr, levels, norm=NORM_N, rel=False):
    """自由場で、時刻 t_yr(年、昇順の配列)での BH のスピン。

    各準位は 1 個の量子から e 倍の回数 A_l = int Gamma_l dt を積み、
    A_l が ln N_sat に届いたとき飽和して、BH を臨界スピンまで回転を失わせるとする。
    飽和にかかる最後の数 e-fold(ln N_sat ~ 170 に対して数個)は無視する。
    """
    t_yr = np.atleast_1d(t_yr)
    mu_s = mu_eV / HBAR_EVS
    M, chi, t = M0, chi0, 0.0
    A = np.zeros(len(levels))
    done = np.zeros(len(levels), bool)
    out = np.empty(len(t_yr))
    k = 0
    while True:
        a = alpha_of(M, mu_eV)
        gfun = gamma_rel_over_mu if rel else gamma_over_mu
        rates = np.array([gfun(a, chi, n, l, m, norm) * mu_s * YR if not d else 0.0
                          for (n, l, m), d in zip(levels, done)])
        # 次に飽和する準位
        t_next, i_next, state_next = np.inf, -1, None
        for i, ((n, l, m), d) in enumerate(zip(levels, done)):
            if d or rates[i] <= 0:
                continue
            Mf, chif = chi_threshold(M, chi, mu_eV, n, l, m)
            if chif >= chi:
                continue
            lnNsat = np.log((chi * M ** 2 - chif * Mf ** 2) * N0 / m)
            dt = max(lnNsat - A[i], 0.0) / rates[i]
            if t + dt < t_next:
                t_next, i_next, state_next = t + dt, i, (Mf, chif)
        while k < len(t_yr) and t_yr[k] < t_next:
            out[k] = chi
            k += 1
        if k == len(t_yr) or i_next < 0:
            out[k:] = chi
            return out
        A += rates * (t_next - t)
        A = np.maximum(A, 0.0)
        t = t_next
        M, chi = state_next
        done[i_next] = True


# ---------- 相対論的な率(E4 の連分数法の表) ----------
# 有効な範囲(E4 で確かめ、E9 で 433 を狭めた): 211 は alpha >= 0.03、322 は 0.1 <= alpha <= 0.52、433 は 0.30 <= alpha <= 0.40。
# 範囲の外と、表の欠けた点では解析式(1908.10370 の式 2.32)に戻す。
_REL_RANGE = {(2, 1, 1): (0.03, 0.55), (3, 2, 2): (0.1, 0.52), (4, 3, 3): (0.30, 0.40)}
_REL = {}


def _rel_table(n, l, m):
    import os
    from scipy.interpolate import RegularGridInterpolator
    key = (n, l, m)
    if key not in _REL:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", f"rates_{n}{l}{m}.npz")
        z = np.load(path)
        al, ch, im = z["alpha"], z["chi"], z["im"].copy()
        F = im / al[:, None] ** (4 * l + 5)
        for i, a in enumerate(al):                 # 欠けた点は解析式で埋める
            for j, c in enumerate(ch):
                if not np.isfinite(F[i, j]):
                    F[i, j] = gamma_over_mu(a, c, n, l, m) * a / a ** (4 * l + 5)
        _REL[key] = RegularGridInterpolator((np.log(al), ch), F, method="cubic", bounds_error=False, fill_value=None)
    return _REL[key]


def gamma_rel_over_mu(alpha, chi, n, l, m, norm=1.0):
    """相対論的な Gamma / mu(norm=1 で振幅の率、2 で占有数の率)。表の範囲外は解析式。"""
    lo, hi = _REL_RANGE.get((n, l, m), (np.inf, -np.inf))
    if not (lo <= alpha <= hi):
        return gamma_over_mu(alpha, chi, n, l, m, norm)
    F = _rel_table(n, l, m)((np.log(alpha), chi))
    return norm * float(F) * alpha ** (4 * l + 5) / alpha
