"""集団解析用の速い chi_SR: 自己相互作用のある超放射で、年齢 T のあとの BH のスピン(研究ノート E10--E16)。

固定点: BGLS(2011.11646)の式 34--36 に 544 を加えた三準位(211・322・544)の安定な固定点(fp3.fixed_points、E14)。
  211 の臨界より上は (211, 322)、下は 211 が空の (322, 544)。211 から先に育たず 322 から育つ BH も同じ固定点で降下する。
どちらの側でも dchi/dt = mu (f/M_pl)^2 H(alpha, chi) で、f と mu は前の因子にしか入らないので、
  I(alpha, chi) = int_chi^{CHI_TOP} dchi' / |H(alpha, chi')|
を一度表にすれば、年齢 T のあとのスピンは I(alpha, chi_T) = I(alpha, chi_0) + tau(tau = mu (f/M_pl)^2 T)を解くだけで決まる。
表は alpha の節(ALPHAS)と w の格子で ln I を持つ。w は 322 の臨界(w = 0)、211 の臨界(w = 1)、CHI_TOP(w = 2)を
結ぶ区分的な座標。211 の臨界のすぐ上と下(幅 1e-8 ほど)に回転を失う速さが最も遅い「ボトルネック」があるので、w は両側に
対数で詰める。chi_after は両隣の alpha の節で I を足してそれぞれ厳密に逆に引き、w を alpha について線形に混ぜる
(tau -> 0 で chi_T = chi_0 が厳密に成り立つ。E16)。
臨界スピンは、相対論的な率(g_vec)が 0 になるスピンを細かい alpha の格子(THR_ALPHAS)で表にしたものを、どこでも使う
(固定点の表、領域の判定、自由場の飽和。E16)。
"""
import os

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq

import fp3
import si
import sr

CHI_TOP = 0.995                       # E4 の表の上端
ALPHAS = np.geomspace(0.01, 0.30, 161)
# 臨界スピンの表(E16)。率は表と解析式の切り替わり(sr._REL_RANGE の端)で不連続なので、その両側に節を置く。
_THR_EDGES = sorted({b for (n, l, m), r in sr._REL_RANGE.items() if (n, l, m) in ((2, 1, 1), (3, 2, 2)) for b in r})
THR_ALPHAS = np.unique(np.concatenate([np.geomspace(0.003, 0.6, 2000)]
                                      + [[b * (1 - 1e-12), b, b * (1 + 1e-12)] for b in _THR_EDGES]))
ALPHA_MAX = 0.5                       # これより上は印 3(211 は超放射できず、m >= 3 が効く)
F_FP_MAX = 1e17                       # 自由場の側で、211 の飽和のあと固定点の降下を使う f の上限 [GeV](E12)
REFTAB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "reftab.npz")   # 値の表(E18)
# 単体版に書き出した値の表(12 MB、git にある)。基準の表がない計算機では、ここから同じ数の辞書を作る(E19)
EXPORTED_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "share", "chi_sr", "chi_sr_tables.npz")
LOW_MARGIN = 7.0                      # E19: 値の表で引けない固定点の点を chi0 のままにする条件(育つ時刻の何 e-fold 前か)
# E20: f の下限 [GeV]。f が約 1e-58 GeV より小さいと、自己相互作用の率 ∝ (M_Pl/f)^4 が倍精度をあふれ、準平衡の判定が NaN に
# なって、自由場の値を印 2 で返してしまう(正しくは chi0 のまま)。物理的に意味のない範囲なので、入力の検査で止める
F_MIN_GEV = 1e-40

_E = np.geomspace(1e-13, 1e-2, 110, endpoint=False)
W_LO = np.concatenate([_E, np.linspace(1e-2, 1 - 1e-2, 300, endpoint=False), 1 - _E[::-1]])   # (0, 1): 322 の臨界から 211 の臨界
W_HI = 1 + np.concatenate([[0.0], _E, np.linspace(1e-2, 1.0, 400)])                           # [1, 2]: 211 の臨界から CHI_TOP
W = np.concatenate([W_LO, W_HI])
U = W                                 # 旧名(E10 の e10_gwonly.py が参照する)


# ---------- 率と固定点 ----------
def g_vec(alpha, chi, n, l, m):
    """超放射の率 / mu(占有数の率)を配列で。連分数法の表の範囲は表、外は解析式(sr.gamma_rel_over_mu と同じ規則)。"""
    alpha, chi = np.broadcast_arrays(np.asarray(alpha, float), np.asarray(chi, float))
    lo, hi = sr._REL_RANGE.get((n, l, m), (np.inf, -np.inf))       # 表のない準位(433 など)は解析式
    out = np.array(sr.gamma_over_mu(alpha, chi, n, l, m, sr.NORM_N), float)
    inr = (alpha >= lo) & (alpha <= hi)
    if inr.any():
        F = sr._rel_table(n, l, m)(np.stack([np.log(alpha[inr]), chi[inr]], -1))
        out[inr] = sr.NORM_N * F * alpha[inr] ** (4 * l + 5) / alpha[inr]
    return out


def _fixed(alpha, chi, c_a=1.0, c_c=1.0):
    """三準位の率と固定点。c_a、c_c は 211x211 -> 322xBH と 322x322 -> 544xBH に掛ける有限 alpha の補正。"""
    g1 = g_vec(alpha, chi, 2, 1, 1)
    g2 = g_vec(alpha, chi, 3, 2, 2)
    g3 = np.array(sr.gamma_over_mu(alpha, chi, 5, 4, 4, sr.NORM_N), float)
    rp = sr.rplus(chi)
    a11, a8 = alpha ** 11, alpha ** 8
    r = fp3.fixed_points(g1, g2, g3, si.KAPPA_BH * a11 * rp * c_a, si.KAPPA_INF * a8, fp3.K_C * a11 * rp * c_c,
                         fp3.K_D * a8, fp3.K_E * a8, fp3.K_F * a8)
    r.update(g1=g1, g2=g2, g3=g3)
    return r


def H(alpha, chi, cbh=False, mass_corr=False):
    """dchi/dt / (mu (f/M_pl)^2)(負、配列で渡せる)。三準位の安定な固定点(E14)。
    P3 = (211, 322): 211 の臨界より上、または 544 の方向に安定なとき。P2 = (322, 544): 211 の臨界より下で、211 の方向に
    安定なとき。どちらも安定でない幅(211 の臨界の下 1e-3 ほど)では遅いほう(降下の時間への寄与は 0.2% 以下)。
    """
    alpha, chi = np.broadcast_arrays(np.asarray(alpha, float), np.asarray(chi, float))
    c = np.interp(alpha, si.CBH_ALPHA, si.CBH_VALUE) if cbh else 1.0
    r = _fixed(alpha, chi, c_a=c)
    g1, g2, g3 = r["g1"], r["g2"], r["g3"]
    with np.errstate(invalid="ignore"):
        h3 = -(g1 * r["e1"] + 2 * g2 * r["e2"])
        h2 = -(2 * g2 * r["f2"] + 4 * g3 * r["f3"])
    h = np.where(r["use3"], h3, np.where(r["use2"], h2, np.maximum(h2, h3)))
    h = np.where((g2 > 0) & np.isfinite(h), h, 0.0)
    if mass_corr:   # 放出 1 回あたり、P3 は m = 3 に対しエネルギー mu、P2 は m = 6 に対し mu を失う
        h = h * np.where(r["use3"], 1 - 2 * chi * alpha / 3, 1 - chi * alpha / 3)
    return h


def extra_growth(alpha, chi):
    """三準位の外の準位(411、433、211 の臨界より上の 544)の、固定点での 1 粒子あたりの増加率の最大 / mu(f によらない)。
    BH が吸収する過程すべてに E7 の有限 alpha の補正 C_BH を掛けて評価する(BGLS: 補正は地平線の近くで決まるので、
    BH が吸収する過程どうしの比はあまり変えない)。印 4 の判定を控えめ(多めに印をつける側)にするための選択(E15)。
    """
    alpha, chi = np.broadcast_arrays(np.asarray(alpha, float), np.asarray(chi, float))
    C = np.interp(alpha, si.CBH_ALPHA, si.CBH_VALUE)
    r = _fixed(alpha, chi, c_a=C, c_c=C)
    use3, use2 = r["use3"], r["use2"]
    x1 = np.where(use2, 0.0, r["e1"])                 # 両方とも不安定な幅では P3 の占有数
    x2 = np.where(use2, r["f2"], r["e2"])
    x3 = np.where(use2, r["f3"], 0.0)
    l411, l433 = fp3.extra_rates(alpha, sr.rplus(chi) * C, x1, x2, x3, g_vec(alpha, chi, 4, 1, 1), g_vec(alpha, chi, 4, 3, 3))
    lmax = np.maximum(np.maximum(l411, l433), np.where(use3, r["lam544"], -np.inf))
    return np.where((r["g2"] > 0) & np.isfinite(lmax), lmax, 0.0)


def switch_cell(a_nodes, alpha):
    """値の表の alpha のセルのうち、超放射の率が連分数法の表と解析式で切り替わる alpha(sr._REL_RANGE の端。211 は 0.03、
    322 は 0.1)を内側に挟むセルにいるか(E19)。基準の運動論はそこで alpha について不連続(表の率は解析式に対して
    322 で -31% から +14%、211 で -12% から +2%)なので、セルの両端の節を結ぶ補間は段差の時刻を誤る(印 5 の判定は
    reftab.compute_ks が表から測る折れで行う。この関数は報告のためだけに使う)。"""
    alpha = np.asarray(alpha, float)
    cell = np.zeros(len(a_nodes) - 1, bool)
    for b in _THR_EDGES:
        cell |= (a_nodes[:-1] < b) & (b < a_nodes[1:])
    i = np.clip(np.searchsorted(a_nodes, alpha, side="right") - 1, 0, len(a_nodes) - 2)
    return cell[i] & (alpha >= a_nodes[0]) & (alpha <= a_nodes[-1])


# ---------- 臨界スピン ----------
def _zero(alpha, n, l, m):
    """準位 (n, l, m) の率が 0 になるスピン(質量固定、相対論的な率の表と同じもの)。どのスピンでも超放射しなければ 1。
    率の表のスピンは 0.995(CHI_TOP)までなので、まず表の中で 0 を探す。表の端より上は 3 次の外挿になり、
    alpha = 0.30 の 433 のように符号を誤ることがある(E17)。臨界が表の端より上のとき(alpha ~ 0.45 以上の 211)だけ外挿で探す。"""
    g = lambda c: float(g_vec(alpha, c, n, l, m))   # noqa: E731
    if g(CHI_TOP) > 0:
        return 1e-4 if g(1e-4) > 0 else brentq(g, 1e-4, CHI_TOP, xtol=1e-14)
    if g(0.9999) > 0:
        return brentq(g, CHI_TOP, 0.9999, xtol=1e-14)
    return 1.0


def chic_vec(alpha, m, n=None, tab=None):
    """臨界スピン(質量固定)。tab があれば、相対論的な率の 0 の表(THR_ALPHAS)を log alpha について補間する(E16)。
    表の下(alpha < 0.003)は m Omega_H = E_n(水素原子の主要項)の解析解、上(alpha > 0.6)は 1(超放射しない)。
    tab がなければ解析解(alpha_E / m >= 1/2 では 1)。"""
    alpha = np.asarray(alpha, float)
    n = m + 1 if n is None else n
    ap = alpha * (1 - alpha ** 2 / (2 * n ** 2))          # alpha E/mu
    out = np.where(ap / m < 0.5, 4 * m * ap / (m ** 2 + 4 * ap ** 2), 1.0)
    if tab is not None and m in (1, 2, 3):
        ta = tab["thr_alpha"]
        thr = tab[f"thr{m}"]
        inside = np.interp(np.log(np.clip(alpha, ta[0], ta[-1])), np.log(ta), thr)
        out = np.where(alpha < ta[0], out, np.where(alpha > ta[-1], 1.0, inside))
    return out


def chi_of_w(w, c1, c2):
    return np.where(w >= 1, c1 + (w - 1) * (CHI_TOP - c1), c2 + w * (c1 - c2))


def w_of_chi(chi, c1, c2):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(chi >= c1, 1 + (chi - c1) / (CHI_TOP - c1), (chi - c2) / (c1 - c2))


# ---------- 表 ----------
def build(cbh=False, mass_corr=False, use_reftab=True, reftab_path=None):
    """表を作る。
    lnI[alpha, w]: ln I、chi = chi_of_w(w, chic, chic2)(節の alpha での臨界)。
    J[alpha, w]: 三準位の外の準位が、chi(w) から CHI_TOP までの降下のあいだに育つ e-fold の数 x (f/M_pl)^2(E15)。
    chic, chic2: 節の alpha での 211 と 322 の臨界。thr_alpha, thr1, thr2, thr3: 211・322・433 の臨界スピンの細かい表(E16、E17)。
    """
    nA = len(ALPHAS)
    lnI = np.empty((nA, len(W)))
    J = np.empty((nA, len(W)))
    c1s, c2s = np.empty(nA), np.empty(nA)
    for i, a in enumerate(ALPHAS):
        c1, c2 = _zero(a, 2, 1, 1), _zero(a, 3, 2, 2)
        c1s[i], c2s[i] = c1, c2
        chis = chi_of_w(W, c1, c2)
        invH = 1.0 / np.abs(H(a, chis, cbh, mass_corr))
        I = -cumulative_trapezoid(invH[::-1], chis[::-1], initial=0.0)[::-1]   # 上から積分
        I[-1] = I[-2] * 1e-3            # w = 2(chi = CHI_TOP)の I = 0 を有限の小さい値で置き換える
        lnI[i] = np.log(I)
        lam = np.maximum(extra_growth(a, chis), 0.0)
        J[i] = -cumulative_trapezoid((lam * invH)[::-1], chis[::-1], initial=0.0)[::-1]
    thr1 = np.array([_zero(a, 2, 1, 1) for a in THR_ALPHAS])
    thr2 = np.array([_zero(a, 3, 2, 2) for a in THR_ALPHAS])
    thr3 = np.array([_zero(a, 4, 3, 3) for a in THR_ALPHAS])
    tab = {"lnI": lnI, "J": J, "chic": c1s, "chic2": c2s, "thr_alpha": THR_ALPHAS, "thr1": thr1, "thr2": thr2,
            "thr3": thr3}
    if use_reftab:
        import reftab
        path = reftab_path or REFTAB_PATH
        if os.path.exists(path):
            tab["rt"] = reftab.load3(path)
            tab["rt_source"] = "reftab"
        elif reftab_path is None and os.path.exists(EXPORTED_PATH):
            # 基準の表(700 MB、git に入れない)がない計算機: 単体版に書き出した表から同じ数の辞書を作る(E19)。
            # プロジェクト版は量子化した値を使うので、値は基準の表から読んだときと同じ。
            tab["rt"] = reftab.load_exported(EXPORTED_PATH)
            tab["rt_source"] = "exported"
        else:                                     # 黙って半解析に戻さない
            raise FileNotFoundError(f"値の表 {path} がない。calc/e18_reftab_build.py で作るか、"
                                    "半解析だけで使うなら build(use_reftab=False)(E17 までの版、E18 の約束の外)")
    return tab


def _alpha_nodes(alpha):
    """alpha を挟む節の番号 i(i と i+1)と、log alpha での重み t。"""
    la = np.log(ALPHAS)
    x = np.clip(np.log(alpha), la[0], la[-1])
    i = np.clip(np.searchsorted(la, x, side="right") - 1, 0, len(la) - 2)
    return i, (x - la[i]) / (la[i + 1] - la[i])


def _fwd(arr, rows, w):
    """arr[row](w) を w について線形に補間(行ごと)。"""
    j = np.clip(np.searchsorted(W, w, side="right") - 1, 0, len(W) - 2)
    s = (w - W[j]) / (W[j + 1] - W[j])
    return arr[rows, j] + s * (arr[rows, j + 1] - arr[rows, j])


def _inv(lnI, rows, y):
    """lnI[row] は w について減少する。lnI[row](w) = y となる w(区分線形の補間の厳密な逆)。
    行ごとにまとめて np.searchsorted で区間を探す(行は高々 len(ALPHAS) 本)。"""
    n = lnI.shape[1]
    y = np.clip(y, lnI[rows, -1], lnI[rows, 0])
    j = np.empty(len(rows), int)
    order = np.argsort(rows, kind="stable")
    r_sorted = rows[order]
    starts = np.flatnonzero(np.r_[True, r_sorted[1:] != r_sorted[:-1]])
    ends = np.r_[starts[1:], len(rows)]
    for st, en in zip(starts, ends):
        idx = order[st:en]
        j[idx] = np.searchsorted(-lnI[r_sorted[st]], -y[idx], side="left") - 1
    j = np.clip(j, 0, n - 2)
    v0, v1 = lnI[rows, j], lnI[rows, j + 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(v0 > v1, (v0 - y) / (v0 - v1), 0.0)
    return W[j] + np.clip(s, 0.0, 1.0) * (W[j + 1] - W[j])


def chi_after(tab, alpha, chi0, tau):
    """固定点で年齢のあとのスピン(配列で渡せる)。tau = mu (f/M_pl)^2 T。
    322 の臨界まで下がりうる。alpha の節の外と、chi0 が 322 の臨界以下の点では chi0 をそのまま返す。
    """
    alpha, chi0, tau = np.broadcast_arrays(np.asarray(alpha, float), np.asarray(chi0, float),
                                           np.asarray(tau, float))
    shape = alpha.shape
    alpha, chi0, tau = alpha.ravel(), chi0.ravel(), tau.ravel()
    c1 = chic_vec(alpha, 1, tab=tab)
    c2 = chic_vec(alpha, 2, tab=tab)
    w0 = np.clip(w_of_chi(chi0, c1, c2), W[0], W[-1])
    i, t = _alpha_nodes(alpha)
    lt = np.log(np.maximum(tau, 1e-300))
    lnI = tab["lnI"]
    wa = _inv(lnI, i, np.logaddexp(_fwd(lnI, i, w0), lt))
    wb = _inv(lnI, i + 1, np.logaddexp(_fwd(lnI, i + 1, w0), lt))
    out = chi_of_w((1 - t) * wa + t * wb, c1, c2)
    keep = (chi0 <= c2) | (alpha < ALPHAS[0]) | (alpha > ALPHAS[-1])
    return np.where(keep, chi0, np.minimum(out, chi0)).reshape(shape)


def flag_fp(tab, alpha, chi_start, chi_T, s2, T_s, mu_s, lnGM2):
    """固定点の領域の印: 三準位の外の準位が年齢のうちに 1 個の量子から育ちうるなら 4、そうでなければ 1(E15)。
    降下のあいだの e-fold は (J(w_T) - J(w_start)) / (f/M_pl)^2。322 の臨界まで下がった BH は、
    そのあと 433 の超放射だけで育つ e-fold(mu g433 T、控えめに年齢全体)も見る。"""
    c1 = chic_vec(alpha, 1, tab=tab)
    c2 = chic_vec(alpha, 2, tab=tab)
    w0 = np.clip(w_of_chi(chi_start, c1, c2), W[0], W[-1])
    wT = np.clip(w_of_chi(chi_T, c1, c2), W[0], W[-1])
    i, t = _alpha_nodes(alpha)
    J = tab["J"]
    J0 = (1 - t) * _fwd(J, i, w0) + t * _fwd(J, i + 1, w0)
    JT = (1 - t) * _fwd(J, i, wT) + t * _fwd(J, i + 1, wT)
    dJ = np.where(JT > J0 * (1 + 1e-9), JT - J0, 0.0)          # 補間の丸めで 0 の区間に正の値が出ないように
    n_desc = dJ / s2
    bottom = wT <= W[0] * 1.000001
    g433 = np.maximum(g_vec(alpha, np.minimum(chi_T, 0.995), 4, 3, 3), 0.0)
    n_bot = np.where(bottom, mu_s * g433 * T_s, 0.0)
    return np.where((n_desc > lnGM2) | (n_bot > lnGM2), 4, 1)


def tau_of(mu_eV, f_gev, T_yr):
    """tau = mu (f/M_pl)^2 T(mu は 1/s に直す)。"""
    return mu_eV / sr.HBAR_EVS * (f_gev / si.MPL_GEV) ** 2 * T_yr * sr.YR


# ---------- 自由場 ----------
def _sat_state(M, chi, mu_eV, m, n, tab):
    """自由場で準位 (n, m=l) が飽和したあとの (M_f, chi_f)。BH の質量の減少を入れる。
    角運動量 dJ = J0 - chi_c(M_f) M_f^2 を渡すと、質量は dM = (E/mu) alpha_sun dJ / m 減る(E は水素原子の主要項)。
    臨界 chi_c は chic_vec(相対論的な率の 0、E16)。M_f を二分法で自己無撞着に解く。
    """
    a_sun = sr.alpha_of(1.0, mu_eV)
    J0 = chi * M ** 2

    def F(Mf):                                  # 0 になる Mf を探す(Mf = M で >= 0、Mf = M/2 で < 0)
        a = a_sun * Mf
        cf = chic_vec(a, m, n, tab)
        return Mf - M + a_sun * (1 - a ** 2 / (2 * n ** 2)) * (J0 - cf * Mf ** 2) / m

    return _bracket_root(F, 0.5 * np.asarray(M, float), np.array(M, float), a_sun, m, n, tab)


def _bracket_root(F, lo, hi, a_sun, m, n, tab, tol=1e-13, itmax=60):
    """F(lo) < 0 <= F(hi) の区間で F = 0 を Illinois 法(はさみうち法の改良)で解く(配列で渡せる)。
    二分法と同じく区間の外に出ないが、10 回ほどの評価で収束する。"""
    Flo, Fhi = F(lo), F(hi)
    side = np.zeros(lo.shape, int)            # 直前に動かした端(-1: lo、+1: hi)
    x = hi.copy()
    for _ in range(itmax):
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.where(Fhi != Flo, hi - Fhi * (hi - lo) / (Fhi - Flo), 0.5 * (lo + hi))
        x = np.clip(x, lo, hi)
        Fx = F(x)
        neg = Fx < 0
        Fhi = np.where(neg & (side == -1), 0.5 * Fhi, Fhi)    # 同じ端が続けて動いたら反対側の値を半分に
        Flo = np.where(~neg & (side == 1), 0.5 * Flo, Flo)
        lo, Flo = np.where(neg, x, lo), np.where(neg, Fx, Flo)
        hi, Fhi = np.where(neg, hi, x), np.where(neg, Fhi, Fx)
        side = np.where(neg, -1, 1)
        if np.all((hi - lo) <= tol * hi) or np.all(np.abs(Fx) <= tol * hi):
            break
    return x, chic_vec(a_sun * x, m, n, tab)


def _free_branch(tab, M, chi0, a, mu_eV, mu_s, T, f_gev, g1):
    """自由場: 211 が 1 個の量子から飽和(E2 の半解析と同じ数え方)。そのあと
    f <= F_FP_MAX: 211 x 211 -> 322 x BH が 322 をすぐに育てるので、固定点の降下で 322 の臨界へ(E12)。
    f >  F_FP_MAX: 322 を t = 0 から 1 個の量子として e-fold を積む(E11 のまま。ODE より遅いことがある、E14)。
    """
    M1, c1 = _sat_state(M, chi0, mu_eV, 1, 2, tab)
    c1 = np.minimum(c1, chi0)
    lnN1 = np.log(np.maximum(chi0 * M ** 2 - c1 * M1 ** 2, 1e-300)) + np.log(sr.N0)
    tsat1 = lnN1 / g1
    M2, c2 = _sat_state(M1, c1, mu_eV, 2, 3, tab)
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
        start = c1[fp]
        drop = start - chi_after(tab, a1[fp], start + 1e-15, tau1)                   # 211 の臨界からの降下
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
    if np.any(f_gev < F_MIN_GEV):
        raise ValueError(f"f_gev は {F_MIN_GEV:g} 以上(これより小さいと、自己相互作用の率が倍精度をあふれる)")


def chi_sr(tab, M, chi0, T_yr, mu_eV, f_gev):
    """年齢 T のあとの BH のスピンと、どの領域で計算したかの印(配列で渡せる)。

    値(E18): 超放射する点(chi0 が 322 の臨界より上)の値は、基準の運動論(fp3.evolve3m: 211・322・544 の三準位、
    BH の質量と角運動量の変化、BH へのエネルギーの戻り、重力波の過程)から直接作った表(calc/reftab.py、
    out/reftab.npz)から引く。表の範囲は alpha 0.008--0.5、f 1e8--10^20.5 GeV(それより上は上端の値)、M 1--300 Msun、
    mu T <= 0.13 eV yr。f < 1e8 GeV と alpha < 0.008 は下の半解析の値のまま(照合で確かめた)。印 1・2 の区別は
    半解析の領域の判定から付ける。
    約束(E18、照合した範囲は M 2--200 Msun、mu 1e-13--1e-11 eV、f 1e6--1e22 GeV、T 1e5--1.4e10 yr、alpha <= 0.5):
    印 0・1・2 の値は、基準の運動論と |chi_SR(T) - chi_ref(T')| <= 0.03 を満たす T' が |ln(T'/T)| <= 0.02 の中にある
    (速い段差の時刻は 2% までずれてよい)。照合の結果(越えたセルの数と最大の誤差)は E18・E19 と calc/e19_validate.py。

    印: 0 = 超放射しない(chi0 が 211 と 322 の臨界以下で、433 も年齢のうちに飽和しない)。chi0 を返す
        1 = 自己相互作用の固定点(三準位 211・322・544)。211 から育つ BH も 322 から育つ BH も、322 の臨界まで
            下がりうる。三準位の外の準位(411、433 など)は年齢のうちに育たない
        2 = 自由場(211 が飽和して臨界へ、そのあと 322 の段階)、または 322 から育つ BH で自己相互作用が飽和までに
            効かず、322 が飽和する場合
        3 = この模型の外。chi0 を返す: 固定点で alpha > 0.30(印 4 の判定の表の外)、値の表の下(alpha < 0.008。固定点では
            alpha < 0.01 で f < 1e8 GeV の点も)で雲が年齢のうちに育ちうる点(E19。固定点の占有数まで育つ時刻、自由場では
            飽和の時刻の 7 e-fold 前より若い点は、印 1・2 で chi0)、alpha > 0.5、
            211・322 は超放射しないが 433(m = 3)が年齢のうちに飽和しうる(E17)、表の範囲の中なのに値を引けない
            (M が 1--300 Msun の外、mu T > 0.13 eV yr、重みのある節がその年齢まで解けていない)、
            chi0 > 0.995 で超放射する(E18)
        4 = 値は返すが、模型の外の準位が年齢のうちに育ちうるので近似は粗い: 固定点で 411、433、211 の臨界より上の
            544 が育ちうる(有限 alpha の補正込みで判定、E15)、自由場で止まったスピンから 433 が飽和しうる(E17)、
            または自由場で 544 の超放射まで進んで 433 の臨界より下にいる(その前に 433 が飽和しているはず。E19)。
            約束の外(照合で最大 0.17)
        5 = 値は表の補間だが、段差の時刻を表の節のあいだで補間できない所。約束の外。
            (a) 速い段差の時刻が節のあいだで急に変わる所(E18): 角の節の値の幅が 0.06 を越え、しかも段差の時刻の
                曲がりが 0.05 dex を越える
            (b) 段差の時刻が alpha の節のあいだで折れている所(E19): 表から測った折れ(reftab.compute_ks)が時刻の許容 2% を
                越え、しかも時刻を折れの分だけ前後にずらすと値が 0.015 より変わる。折れの原因は、超放射の率が連分数法の
                表と解析式で切り替わる alpha(211 は 0.03、322 は 0.1。基準の運動論がそこで不連続)と、322 の率の表の
                alpha 0.10--0.125 の節の誤差(倍精度の限界、最大 9%)。E18 の版ではここが印 1・2 のまま最大 0.86 ずれた
            (c) 段のあとの平らな所で、補間が止まっていて(次の準位に届かない角の節がある)、角の節の値の幅が 0.02 を
                越える点(E19)。値は進みのいちばん遅い節の値に引かれるので、その幅までずれうる
    領域の判定(BGLS の B|C の境界、式 53・56): 準平衡の占有数 eps_eq と行き過ぎの占有数 eps_th が、
    どちらも自由場で飽和する量 eps_max = chi0 - chi_c より小さければ準平衡。
    準平衡では、1 個の量子から eps_eq まで育つ時間 t_grow を年齢から引き、行き過ぎで失う eps_th を初期スピンから引く
    (表の範囲の外で使う半解析の値)。
    """
    M, chi0, T_yr, mu_eV, f_gev = np.broadcast_arrays(*[np.asarray(x, float) for x in (M, chi0, T_yr, mu_eV, f_gev)])
    shape = M.shape
    M, chi0, T_yr, mu_eV, f_gev = (x.ravel() for x in (M, chi0, T_yr, mu_eV, f_gev))
    _check(M, chi0, T_yr, mu_eV, f_gev)
    a = sr.alpha_of(M, mu_eV)
    mu_s = mu_eV / sr.HBAR_EVS
    T = T_yr * sr.YR
    lnGM2 = np.log(sr.N0 * M ** 2)
    c1 = chic_vec(a, 1, tab=tab)
    c2 = chic_vec(a, 2, tab=tab)
    big = a > ALPHA_MAX
    out = chi0.copy()
    reg = np.zeros(len(M), int)
    reg[big] = 3

    # ---- 値の表(E18): 基準の運動論から作った表で、超放射する点の値を引く ----
    # 表の範囲の外: f が表の下端より小さい点と alpha が表の下端より小さい点は、下の半解析の値のまま(照合で確かめる)。
    # f が表の上端より大きい点は上端の値(自己相互作用が効かない)。M、mu T が表の範囲の外の点と、重みのある節が解けて
    # いない点は印 3。
    rt = tab.get("rt")
    tv = chi0.copy()
    use_t = np.zeros(len(M), bool)
    out_rt = np.zeros(len(M), bool)
    sharp_t = np.zeros(len(M), bool)
    kink_t = np.zeros(len(M), bool)              # 段差の時刻が alpha の節のあいだで折れていて、値がそれに敏感な点(E19)
    if rt is not None:
        import reftab
        cand_t = (chi0 > c2) & ~big
        if cand_t.any():
            it = np.flatnonzero(cand_t)
            lfmax = rt["LF_NODES"][-1]
            ff = np.minimum(f_gev[it], 10 ** lfmax)
            v, ok, sh, kk = reftab.lookup_rt3(rt, T_yr[it], a[it], chi0[it], ff, M[it], mu_eV[it], c1[it], c2[it])
            tv[it] = v
            use_t[it] = ok
            sharp_t[it] = sh & ok
            kink_t[it] = kk & ok
            below = (np.log10(f_gev[it]) < rt["LF_NODES"][0]) | (a[it] < rt["A_NODES"][0])
            out_rt[it] = ~ok & ~below            # 表の範囲の中なのに引けない、または M・mu T が範囲の外

    sr211 = (chi0 > c1) & ~big
    # ---- 211 が超放射する BH ----
    g1 = np.zeros(len(M))
    i1 = np.flatnonzero(sr211)
    g1[i1] = g_vec(a[i1], np.minimum(chi0[i1], 0.995), 2, 1, 1) * mu_s[i1]           # 1/s(211 が超放射する点だけ)
    ginf = si.KAPPA_INF * a ** 8 * (si.MPL_GEV / f_gev) ** 4 * mu_s
    gbh = si.KAPPA_BH * a ** 11 * (si.MPL_GEV / f_gev) ** 4 * sr.rplus(np.minimum(chi0, 0.999)) * mu_s
    emax = chi0 - c1
    with np.errstate(divide="ignore", invalid="ignore"):
        eeq = 2 / np.sqrt(3) * np.sqrt(g1 * ginf) / gbh
        eth = np.sqrt(2 * g1 * lnGM2 / gbh)
    si_reg = sr211 & (g1 > 0) & (eeq < emax) & (eth < emax)
    inval = si_reg & (a > ALPHAS[-1])                   # 固定点で alpha > 0.30: 印 4 の判定の表(J)の外
    reg[inval] = 3
    reg[sr211 & (g1 <= 0) & (chi0 > CHI_TOP)] = 3       # 臨界が率の表の端(0.995)より上: 率を引けないので模型の外(E16)
    q = si_reg & ~inval
    if q.any():
        tgrow = (np.log(np.maximum(eeq[q], 1e-300)) + lnGM2[q]) / g1[q]
        Teff = np.maximum(T[q] - tgrow, 0.0)
        tau = mu_s[q] * (f_gev[q] / si.MPL_GEV) ** 2 * Teff
        c0eff = np.maximum(chi0[q] - eth[q], c1[q])
        out[q] = np.where(Teff > 0, chi_after(tab, a[q], c0eff, tau), chi0[q])
        out[q] = np.where(use_t[q], tv[q], out[q])
        reg[q] = flag_fp(tab, a[q], np.maximum(c0eff, out[q]), out[q], (f_gev[q] / si.MPL_GEV) ** 2, Teff, mu_s[q], lnGM2[q])
        # alpha < 0.01(半解析の固定点の表の下、E19)。値の表で引けた点(alpha >= 0.008、f >= 1e8 GeV)はその値で、印は上の
        # flag_fp(J の表は alpha <= 0.038 で 0 なので、外の準位は育たない)。引けない点には半解析の値がない: 雲が固定点の
        # 占有数に届く時刻 t_grow の LOW_MARGIN e-fold 前より若ければ chi0 のまま(失うスピンは e^-7 eps_eq 以下。上で
        # Teff = 0 なので out = chi0、印 1)、そうでなければ模型の外(印 3、chi0 を返す)。
        iq = np.flatnonzero(q)
        late = (a[q] < ALPHAS[0]) & ~use_t[q] & (g1[q] * (tgrow - T[q]) < LOW_MARGIN)
        out[iq[late]] = chi0[iq[late]]
        reg[iq[late]] = 3
    fr = sr211 & ~si_reg & (g1 > 0)
    if fr.any():
        out[fr] = _free_branch(tab, M[fr], chi0[fr], a[fr], mu_eV[fr], mu_s[fr], T[fr], f_gev[fr], g1[fr])
        out[fr] = np.where(use_t[fr], tv[fr], out[fr])
        reg[fr] = 2
        if rt is not None:
            # 値の表の下(alpha < 0.008)の自由場の点(E19)。半解析の飽和の時刻は約束の 2% より粗いことがある(照合で最大 0.08)
            # ので、雲が飽和する時刻の LOW_MARGIN e-fold 前より若い点だけ chi0(印 2)にし、ほかは模型の外(印 3)。
            late = fr & (a < rt["A_NODES"][0]) & (g1 * T > lnGM2 + np.log(np.maximum(emax, 1e-300)) - LOW_MARGIN)
            out[late] = chi0[late]
            reg[late] = 3

    # ---- 211 は超放射しないが 322 はする BH(chi_c322 < chi0 <= chi_c211)、E14 ----
    # 322 が 1 個の量子から育ち、322 x 322 -> 544 x BH と 544 x 544 -> 322 x inf の (322, 544) の固定点に入る(211 は空のまま)。
    # 判定は 211 の場合と同じ形: 固定点の占有数 eps_eq' = (2/sqrt3) sqrt(g322 g_d)/g_c と行き過ぎ eps_th' = sqrt(2 g322 ln GM^2/g_c)
    # が、どちらも飽和の量 eps_max' = (chi0 - chi_c322)/2 より小さければ固定点。そのときは育つ時間を年齢から引き、
    # 行き過ぎで失う 2 eps_th' を初期スピンから引いて、同じ表で降下させる。そうでなければ自由場(322 が飽和して、
    # BH の質量の減少を入れた 322 の臨界へ。E16)。
    only322 = (~sr211) & (chi0 > c2) & ~big
    if only322.any():
        i3 = np.flatnonzero(only322)
        M3, a3, x0, T3, f3, m3, u3 = M[i3], a[i3], chi0[i3], T[i3], f_gev[i3], mu_s[i3], mu_eV[i3]
        c0 = np.minimum(x0, 0.995)
        g2 = g_vec(a3, c0, 3, 2, 2) * m3
        s4 = (si.MPL_GEV / f3) ** 4 * m3
        gc = fp3.K_C * a3 ** 11 * sr.rplus(c0) * s4
        gd = fp3.K_D * a3 ** 8 * s4
        emax2 = (x0 - c2[i3]) / 2
        with np.errstate(divide="ignore", invalid="ignore"):
            eeq2 = 2 / np.sqrt(3) * np.sqrt(np.maximum(g2, 0) * gd) / gc
            eth2 = np.sqrt(2 * np.maximum(g2, 0) * lnGM2[i3] / gc)
            tgrow2 = (np.log(np.maximum(eeq2, 1e-300)) + lnGM2[i3]) / g2
        fpq = (g2 > 0) & (eeq2 < emax2) & (eth2 < emax2)
        ok = fpq & (a3 <= ALPHAS[-1])                   # alpha > 0.30 は印 3(下の rg の初期値)
        res = x0.copy()
        rg = np.where((g2 > 0) | (x0 > CHI_TOP), 3, 0)
        if ok.any():
            Teff = np.maximum(T3[ok] - tgrow2[ok], 0.0)
            tau = m3[ok] * (f3[ok] / si.MPL_GEV) ** 2 * Teff
            c0eff = np.maximum(x0[ok] - 2 * eth2[ok], c2[i3][ok] + 1e-12)
            res[ok] = np.where(Teff > 0, chi_after(tab, a3[ok], c0eff, tau), x0[ok])
            res[ok] = np.where(use_t[i3][ok], tv[i3][ok], res[ok])
            rg[ok] = flag_fp(tab, a3[ok], np.maximum(c0eff, res[ok]), res[ok], (f3[ok] / si.MPL_GEV) ** 2, Teff, m3[ok],
                             lnGM2[i3][ok])
            # alpha < 0.01 で値の表を引けない点(E19、211 から育つ BH と同じ扱い)
            jo = np.flatnonzero(ok)
            late = (a3[ok] < ALPHAS[0]) & ~use_t[i3][ok] & (g2[ok] * (tgrow2[ok] - T3[ok]) < LOW_MARGIN)
            res[jo[late]] = x0[jo[late]]
            rg[jo[late]] = 3
        free = (g2 > 0) & ~fpq
        if free.any():
            Mf2, cf2 = _sat_state(M3[free], x0[free], u3[free], 2, 3, tab)
            cf2 = np.minimum(cf2, x0[free])
            lnN2 = np.log(np.maximum((x0[free] * M3[free] ** 2 - cf2 * Mf2 ** 2) / 2, 1e-300)) + np.log(sr.N0)
            res[free] = np.where(T3[free] >= lnN2 / g2[free], cf2, x0[free])
            res[free] = np.where(use_t[i3][free], tv[i3][free], res[free])
            rg[free] = 2
            if rt is not None:                          # 値の表の下の自由場の点(E19、211 から育つ BH と同じ扱い)
                late = free & (a3 < rt["A_NODES"][0]) & (g2 * T3 > lnGM2[i3] + np.log(np.maximum(emax2, 1e-300)) - LOW_MARGIN)
                res[late] = x0[late]
                rg[late] = 3
        out[i3] = res
        reg[i3] = rg

    # 表の範囲の中なのに値を引けない点(M・mu T が範囲の外、節が解けていない)は印 3(chi0 を返す)。E18
    bad_t = out_rt & np.isin(reg, (1, 2, 4))
    out[bad_t] = chi0[bad_t]
    reg[bad_t] = 3
    # 初期スピンが率の表の端(0.995)より上で超放射する点も印 3(E18)
    top = (chi0 > CHI_TOP) & np.isin(reg, (1, 2, 4))
    out[top] = chi0[top]
    reg[top] = 3

    # ---- m = 3(433)の超放射(模型に入っていない)、E17 ----
    # 返すスピンが 433 の臨界より上で、433 が 1 個の量子から(自由場で)年齢のうちに飽和しうる点に印をつける。
    # 印 0(超放射しない、chi0 を返す)は印 3 に、印 2(自由場、322 の臨界などで止まる)は印 4 に変える。値は変えない。
    c3 = chic_vec(a, 3, tab=tab)
    cand = ((reg == 0) | (reg == 2)) & (out > c3) & ~big
    if cand.any():
        i4 = np.flatnonzero(cand)
        g433 = g_vec(a[i4], np.minimum(out[i4], 0.995), 4, 3, 3) * mu_s[i4]
        lnN3 = np.log(np.maximum((out[i4] - c3[i4]) / 3, 1e-300)) + lnGM2[i4]
        grows = (g433 > 0) & (g433 * T[i4] > lnN3)
        reg[i4[grows]] = np.where(reg[i4[grows]] == 0, 3, 4)
    # 自由場で、返すスピンが 433 の臨界より下の点(E19)。三準位の運動論では、322 の臨界で止まったあと 544(m = 4)の超放射で
    # さらに下がる(alpha ~ 0.45 以上の古い BH)。しかし 433(m = 3)の率は 544 の率より大きく、同じスピンから先に飽和するはず
    # なので、ここまで進んだ点の値は粗い。上の判定は返すスピンしか見ないので、この点を印 2 に戻していた。印 4 にする。
    reg[(reg == 2) & (out < c3) & ~big] = 4
    # 速い段差の時刻が表の節のあいだで急に変わる所(E18): 値は表の補間のまま、印 5
    reg[np.isin(reg, (1, 2)) & sharp_t] = 5
    # 段差の時刻が alpha の節のあいだで折れていて、値がそれに敏感な所(E19): 値は表の補間のまま、印 5
    reg[np.isin(reg, (1, 2)) & kink_t] = 5
    return out.reshape(shape), reg.reshape(shape)
