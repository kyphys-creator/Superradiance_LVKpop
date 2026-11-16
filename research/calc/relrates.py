"""自己相互作用の過程の率を Kerr 時空で計算する(相対論的な補正、DESIGN.md D18 の手順 3)。

単位は G = c = hbar = M = 1(mu = alpha、a = chi)。r+- = 1 +- sqrt(1 - a^2)、Delta = (r - r+)(r - r-)、Omega_H = a/(2 r+)。
場 Phi = sum_n sqrt(N_n) (u_n e^{-i omega_n t} + c.c.)、u_n = R_n(r) S_n(theta) e^{i m_n phi}。u_n は Klein-Gordon の電荷で
1 に規格化する(N_n = eps_n G M^2):
    1 = 2 int dr |R|^2 { [(r^2 + a^2)^2 omega - 2 a r m]/Delta - a^2 omega <sin^2 theta> }、  ∮|S|^2 dOmega = 1。
(Boyer-Lindquist の切り方では地平線で対数に発散するが、束縛状態では R(r+) が小さく、切る位置 delta によらないことを確かめる。)
運動方程式 (box - mu^2) Phi = -(lambda/6) Phi^3 の、振動数 omega_X = omega_i + omega_j - omega_k、m_X = m_i + m_j - m_k の源は
    T_X = -(lambda/2) s_ij sqrt(N_i N_j N_k) u_i u_j u_k^*、 s_ij = 1(i = j)、2(i != j)
で、rho^2 = r^2 + a^2 cos^2 theta を掛けて S_ell で射影すると動径の方程式
    (Delta R')' + V R = T_ell(r) = -(lambda/2) s_ij R_i R_j R_k^* [r^2 A_ell + a^2 B_ell]、
    V = K^2/Delta + 2 a m omega - a^2 omega^2 - mu^2 r^2 - Lambda、 K = (r^2 + a^2) omega - a m、
    A_ell = ∮ S_i S_j S_k^* S_ell^*、 B_ell = ∮ cos^2 theta S_i S_j S_k^* S_ell^*
になる(Witte & Mummery 2412.03655 の式 61・64 と同じ形。式 59・60 の 1/Delta と指数の符号は s = 0 では誤り)。
地平線で内向き(R_in -> (r - r+)^{-i sigma}、sigma = 2 r+ (omega - m Omega_H)/(r+ - r-))、無限遠で減衰する(omega < mu)解
R_up から、Green 関数の解は地平線で Z R_in、Z = (1/W) int R_up T dr、W = R_in P_up - R_up P_in(P = Delta R')。
地平線を通る数の流れは N_dot = 2 (omega - m Omega_H)(r+^2 + a^2) sum_ell |Z_ell|^2。過程が 1 回起きると X が 1 個 BH に
吸われるので、lambda = 1、N = 1 での N_dot が Gamma_hat、主要項の形 gamma = kappa alpha^11 r+ (M_pl/f)^4 mu に合わせると
    kappa_eff(alpha, chi) = Gamma_hat / (alpha^8 r+)(hyd.py の kappa と同じ規約)、 補正 C = kappa_eff / kappa_hyd。
計算: x = ln(r - r+) を変数にする(dR/dx = P/(r - r-)、dP/dx = -(r - r+) V R + (r - r+) T。(r - r+) V は地平線で有限)。
角度の関数は球面調和関数(|c^2| <= a^2 alpha^4/n^2 で、扁平の補正は 1e-3 程度)、Lambda = l(l+1) - c^2 <cos^2 theta>。
放出(omega_X > mu、rate_inf): 無限遠で外向きの解 R_up ~ r^{-1 + i eta} e^{i k r}(k = sqrt(omega^2 - mu^2)、
eta = (2 omega^2 - mu^2)/k)を、漸近級数(asym_series)で k r1 >= 400 の所から始めて内向きに積分する。無限遠への数の流れは
2 |Z_inf|^2 Im(R_up^* P_up)、Z_inf = (1/W) int R_in T dr(R_up の規格化によらない)、地平線への流れは BH の過程と同じ式。
kappa_eff = Gamma_hat/alpha^5(hyd.kappa_inf と同じ規約)。部分波は l = |m_X| から l_i + l_j + l_k + 2 まで(ells を渡さないときは
l の小さい方から足し、寄与が和の ELL_TOL = 1e-9 より小さくなったら止め、a^2 cos^2 theta の項だけの l_i + l_j + l_k + 2 は飛ばす。
alpha = 0.52、chi = 0.995 で全部の部分波と 1.2e-9)。
まとめて解く(point、relcorr.py の表): 8 準位の束縛状態を共通の格子で 1 回だけ作り、過程で使い回す(過程ごとの計算と 5e-11)。
確かめ(2026-11-16): 9 個の過程で alpha -> 0 の切片が hyd.kappa_inf と 1.6e-5 以内(python3 calc/relrates.py inf)。
格子の点の数を倍、地平線の切り方 1e-8 -> 1e-11、始める k r1 を 4 倍、部分波を 4 つ増やしても、alpha <= 0.4 で率は 3e-5 以内。
"""
import math
import os
import sys

import numpy as np
from scipy.integrate import cumulative_simpson, solve_ivp, simpson
from scipy.special import gammaln, lpmv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import leaver  # noqa: E402
import sr  # noqa: E402

DELTA0 = 1e-8          # 地平線から内側の端までの r - r+(r+ を掛ける)
NX = 24000             # x の格子の点の数
ELL_TOL = 1e-9         # ells を渡さないとき、部分波は l の小さい方から足し、最後の l の寄与が和のこれより小さくなったら止める


class Kerr:
    def __init__(self, alpha, chi):
        self.mu = float(alpha)
        self.a = float(chi)
        b = np.sqrt(1.0 - self.a ** 2)
        self.rp, self.rm = 1.0 + b, 1.0 - b
        self.OmH = self.a / (2.0 * self.rp)


def Ylm(l, m, x):
    """球面調和関数の theta の部分(∮|Y|^2 dOmega = 1、e^{i m phi} は別)。x = cos theta。"""
    m_ = abs(m)
    N = np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - m_ + 1) - gammaln(l + m_ + 1)))
    return N * lpmv(m_, l, x)


_GL = np.polynomial.legendre.leggauss(96)


def ang(levels_lm, extra_cos2=False):
    """∮ prod Y (共役は呼ぶ側で m の符号に入れる。実の関数なので theta の積分だけ)。2 pi を掛ける(phi の積分)。"""
    x, w = _GL
    f = np.ones_like(x)
    for l, m in levels_lm:
        f = f * Ylm(l, m, x)
    if extra_cos2:
        f = f * x ** 2
    return 2 * np.pi * np.sum(w * f)


def cos2(l, m):
    return float(leaver.cos2_avg(l, m))


def Lam(l, m, c2):
    return l * (l + 1) - c2 * cos2(l, m)


def r_of_x(K, x):
    return K.rp + np.exp(x)


def rhs(K, omega, l, m, src=None):
    """y = (R, P) の x についての右辺。src(r) は (r - r+) を掛ける前の源(複素)。
    U = K^2/(r - r-) + (r - r+)(2 a m omega - a^2 omega^2 - mu^2 r^2 - Lambda) を、遠方で K^2 と mu^2 r^2 Delta の
    打ち消しが起きない形に展開して使う(omega ~ mu では大きな項どうしが 6 桁以上消え、解法の刻みが潰れて止まっていた):
        U = [c4 r^4 + c3 r^3 + c2 r^2 + c0]/(r - r-) + (r - r+) cst、
        c4 = omega^2 - mu^2、c3 = 2 mu^2、c2 = (2 omega^2 - mu^2) a^2 - 2 a m omega、c0 = omega^2 a^4 - 2 a^3 m omega + a^2 m^2。
    解法の時間のほとんどが右辺なので、numpy のスカラーでなく Python の float と complex で計算する。"""
    mu, a = K.mu, K.a
    omega = complex(omega)
    k2 = omega * omega - mu * mu
    L_ = Lam(l, m, a * a * k2)
    cst = 2 * a * m * omega - a * a * omega * omega - L_
    c4, c3 = k2, 2.0 * mu * mu
    c2 = (2 * omega * omega - mu * mu) * a * a - 2 * a * m * omega
    c0 = omega * omega * a ** 4 - 2 * a ** 3 * m * omega + a * a * m * m
    rp, rm = float(K.rp), float(K.rm)
    exp = math.exp

    def f(x, y):
        d = exp(x)
        r = rp + d
        r2 = r * r
        irm = 1.0 / (r - rm)
        U = (c4 * r2 * r2 + c3 * r2 * r + c2 * r2 + c0) * irm + d * cst
        y0, y1, y2, y3 = y.tolist()
        R, P = complex(y0, y1), complex(y2, y3)
        dR = P * irm
        dP = -U * R
        if src is not None:
            dP += d * src(r)
        return [dR.real, dR.imag, dP.real, dP.imag]
    return f


def _ivp(f, span, y0, t_eval=None):
    """solve_ivp(DOP853、rtol 1e-11)。途中で止まったら例外にする(止まっても返る値は途中までなので、黙って使わない)。"""
    sol = solve_ivp(f, span, y0, method="DOP853", t_eval=t_eval, rtol=1e-11, atol=1e-30)
    if sol.status != 0:
        raise RuntimeError(f"積分が止まった: {sol.message}(x = {sol.t[-1]:.4g}、目標 {span[1]:.4g})")
    return sol


def kappa_decay(K, omega):
    return np.sqrt(K.mu ** 2 - omega ** 2 + 0j)


def r_max_for(K, kap, n_eff=1.0):
    """減衰の解を始める外側の端: kappa r が 60 になる所(と、雲の大きさの 20 倍の大きい方)。"""
    return max(60.0 / kap.real, 20.0 * n_eff ** 2 / K.mu ** 2, 50.0)


def grid(K, r_max):
    return np.linspace(np.log(DELTA0 * K.rp), np.log(r_max - K.rp), NX)


def solve_up(K, omega, l, m, xg):
    """無限遠で減衰する解(内向きに積分)。R ~ r^(nu - 1) e^{-kappa r}、nu = (2 omega^2 - mu^2)/kappa(omega は複素でもよい)。"""
    kap = kappa_decay(K, omega)
    if kap.real < 0:
        kap = -kap
    nu = (2 * omega * omega - K.mu ** 2) / kap
    r1 = r_of_x(K, xg[-1])
    R1 = 1.0 + 0j
    D1 = (r1 - K.rp) * (r1 - K.rm)
    P1 = D1 * R1 * (-kap + (nu - 1) / r1)
    sol = _ivp(rhs(K, omega, l, m), (xg[-1], xg[0]), [R1.real, R1.imag, P1.real, P1.imag], t_eval=xg[::-1])
    y = sol.y[:, ::-1]
    return y[0] + 1j * y[1], y[2] + 1j * y[3]


def solve_in(K, omega, l, m, xg):
    """地平線で内向きの解(外向きに積分)。R -> (r - r+)^{-i sigma}。"""
    sig = 2 * K.rp * (omega - m * K.OmH) / (K.rp - K.rm)
    d0 = np.exp(xg[0])
    R0 = np.exp(-1j * sig * xg[0])
    P0 = (d0 + K.rp - K.rm) * (-1j * sig) * R0
    sol = _ivp(rhs(K, omega, l, m), (xg[0], xg[-1]), [R0.real, R0.imag, P0.real, P0.imag], t_eval=xg)
    return sol.y[0] + 1j * sol.y[1], sol.y[2] + 1j * sol.y[3]


def bound_omega(K, n, l, m):
    """準束縛状態の M omega(倍精度の Leaver。初期値は水素型 alpha^5 まで)。"""
    guess = K.mu * sr.omega_over_mu(K.mu, n, l, m, K.a)
    return complex(leaver.bound_state(K.mu, K.a, n, l, m, omega0=guess + 0j))


class Bound:
    """KG 電荷で規格化した準束縛状態の動径の関数(実の Re omega)。
    実の振動数で無限遠から減衰の条件で積分すると、地平線では内向きと外向きが同じ大きさの定在波になる(m = 0 では方程式が
    実なので厳密に)。本当の準束縛状態は地平線で内向きだけなので、地平線から内向きの条件で外へ、無限遠から減衰の条件で内へ
    積分し、雲の山(|R| が最大の所)でつなぐ。どちらも Leaver の複素の固有値で積分する(実の振動数では、幅の広い n00 の準位で
    外向きの積分が増える解を拾う)。つなぎ目の P の食い違い(match_err)は固有値の精度で、小さいことを確かめる。"""

    def __init__(self, K, n, l, m, xg, omega=None):
        self.n, self.l, self.m = n, l, m
        self.omega_c = bound_omega(K, n, l, m) if omega is None else omega
        w = self.omega_c.real
        self.omega = w
        Rd, Pd = solve_up(K, self.omega_c, l, m, xg)         # 複素の固有値で(幅の広い n00 でも、つなぎ目で合う)
        Rh, Ph = solve_in(K, self.omega_c, l, m, xg)
        r = r_of_x(K, xg)
        # つなぐのは雲のいちばん外の山。減衰の解は、固有値の丸めの誤差で地平線の近くに r^{-l-1} の解を拾い、小さい alpha
        # では雲の山より大きくなる(alpha = 0.01 の 322 で山が r = 6 に出た)ので、水素型の半径 n^2 a0 の 1/4 より外で探す。
        sel = np.flatnonzero(r > max(4 * K.rp, 0.25 * n * n / K.mu ** 2))
        jm = sel[np.argmax(np.abs(Rd[sel]))]
        cc = Rd[jm] / Rh[jm]
        # P の食い違いを、P の自然な大きさ Delta R / r で測る(山では P がほぼ 0 なので P 自身では割らない)
        self.match_err = float(abs(Pd[jm] - cc * Ph[jm]) / abs((r[jm] - K.rp) * (r[jm] - K.rm) * Rd[jm] / r[jm]))
        R = np.where(np.arange(len(xg)) < jm, cc * Rh, Rd)
        P = np.where(np.arange(len(xg)) < jm, cc * Ph, Pd)
        d = r - K.rp
        dens = 2 * np.abs(R) ** 2 * (((r * r + K.a ** 2) ** 2 * w - 2 * K.a * r * m) / (r - K.rm)
                                     - d * K.a ** 2 * w * (1 - cos2(l, m)))
        Nn = simpson(dens, x=xg)
        self.norm_integral = Nn
        s = 1.0 / np.sqrt(Nn)
        self.R, self.P = R * s, P * s
        # 地平線での内向きの振幅: R = A_in (r-r+)^{-i sig} + A_out (r-r+)^{+i sig}
        sig = 2 * K.rp * (w - m * K.OmH) / (K.rp - K.rm)
        dRdx = self.P[0] / (d[0] + K.rp - K.rm)
        self.A_in = (1j * sig * self.R[0] - dRdx) / (2j * sig) * np.exp(1j * sig * xg[0])
        self.A_out = (1j * sig * self.R[0] + dRdx) / (2j * sig) * np.exp(-1j * sig * xg[0])
        # 動径の節の数(雲の領域で Re(R e^{-i arg}) の符号の変化)
        ph = np.exp(-1j * np.angle(self.R[np.argmax(np.abs(self.R))]))
        sel = (r > 4 * K.rp) & (np.abs(self.R) > 1e-3 * np.abs(self.R).max())
        v = (self.R * ph).real[sel]
        self.nodes = int(np.sum(np.diff(np.sign(v)) != 0))


def decay_check(alpha, chi, n):
    """n00 の準位: 地平線を通る数の流れ 2 omega (r+^2 + a^2)|A_in|^2 と、-2 Im omega(多倍長の Leaver)の比(1 になるはず)。"""
    K = Kerr(alpha, chi)
    om = complex(leaver.bound_state_mp(alpha, chi, n, 0, 0, omega0=alpha * sr.omega_over_mu(alpha, n, 0, 0, chi),
                                       N=int(max(4000, 500 / alpha)), dps=40))
    xg = grid(K, r_max_for(K, kappa_decay(K, om.real), n))
    b = Bound(K, n, 0, 0, xg, omega=om)
    flux = 2 * b.omega * (K.rp ** 2 + K.a ** 2) * abs(b.A_in) ** 2
    return dict(alpha=alpha, chi=chi, n=n, flux=flux, minus2im=-2 * om.imag, ratio=flux / (-2 * om.imag), nodes=b.nodes,
                Aout_over_Ain=abs(b.A_out / b.A_in), match_err=b.match_err)


def _spread(v, ref, scale=None):
    """格子の中ほど(NX/10 から 9 NX/10)で、一定のはずの量 v の ref からのずれの最大(|ref| で割る)。scale(|R||P| など)を
    渡すと、scale > 1e3 |ref| の点(障壁の内側や地平線の近くの強い反射で、流れが |R||P| の打ち消しの差になる所。そこでは解法の
    相対誤差 1e-11 が |R||P|/|ref| 倍に増える)は除く。無限遠への流れは格子の外の端の値を使うので、そこが一定なら足りる。"""
    sel = slice(NX // 10, 9 * NX // 10)
    dv = np.abs(v[sel] / ref - 1)
    if scale is not None:
        dv = dv[scale[sel] < 1e3 * abs(ref)]
    return float(dv.max()) if dv.size else np.nan


def _setup(K, levs, omegas, states, xg, r_need):
    """三つの準位の (omega, 格子, 束縛状態)。xg と states を渡すと使い回す(格子の外の端が r_need より内側なら例外)。"""
    if omegas is None:
        omegas = [bound_omega(K, *lv) for lv in levs]
    r_need = max(r_need, max(r_max_for(K, kappa_decay(K, w.real), lv[0]) for w, lv in zip(omegas, levs)))
    if xg is None:
        xg = grid(K, r_need)
        if states:
            raise ValueError("states を渡すときは、その格子 xg も渡す")
    elif r_of_x(K, xg[-1]) < r_need * (1 - 1e-12):
        raise ValueError(f"格子の外の端 {r_of_x(K, xg[-1]):.4g} が要る半径 {r_need:.4g} より内側")
    B = dict(states) if states else {}
    for lv, w in zip(levs, omegas):
        if lv not in B:
            B[lv] = Bound(K, *lv, xg, omega=w)
    return omegas, xg, B


def _ells(i, j, k, mX):
    lsum = i[1] + j[1] + k[1]
    return [ell for ell in range(abs(mX), lsum + 3) if (ell + lsum) % 2 == 0]


def _forced_charge(K, wX, ell, mX, xg, Rin, Rup, T, Wm):
    """強制振動の解 R_X = [R_in int_r^inf R_up T dr' + R_up int_{r+}^r R_in T dr'] / W の KG 電荷(束縛の X、omega_X < mu)。
    int_r^inf は外の端から内向きに足す(R_in は外で e^{kappa r} で増えるので、全体から引くと桁が落ちる)。"""
    r = r_of_x(K, xg)
    d = r - K.rp
    I_in = cumulative_simpson(Rin * T * d, x=xg, initial=0)
    I_up = cumulative_simpson((Rup * T * d)[::-1], x=-xg[::-1], initial=0)[::-1]       # u = -x で外から内へ
    RX = (Rin * I_up + Rup * I_in) / Wm
    dens = 2 * np.abs(RX) ** 2 * (((r * r + K.a ** 2) ** 2 * wX - 2 * K.a * r * mX) / (r - K.rm)
                                  - d * K.a ** 2 * wX * (1 - cos2(ell, mX)))
    return float(simpson(dens, x=xg))


def rate_bh(alpha, chi, i, j, k, ells=None, states=None, xg=None, omegas=None, details=False):
    """i x j -> k x BH(omega_X < mu)の kappa_eff = Gamma_hat/(alpha^8 r+)。i, j, k は (n, l, m)。
    m_X = 0 が主要項の過程。m_X != 0(主要項で 0、閾値を越えた放出の過程)も同じ式で、地平線への流れは
    2 (omega_X - m_X Omega_H)(r+^2 + a^2)|Z|^2(超放射の m_X では負)。部分波は l = |m_X| から l_i + l_j + l_k + 2 まで。
    details の tau_X = (強制振動の X の KG 電荷)/Gamma_hat は X が雲にとどまる時間(M の単位)。BGLS の式 A14 の摂動の条件は
    (この過程が作る準位 k の 1 個あたりの率) x tau_X << 1(N_X/N_k = tau_X Gamma_k)。"""
    K = Kerr(alpha, chi)
    levs = (i, j, k)
    mX = i[2] + j[2] - k[2]
    if omegas is None:
        omegas = [bound_omega(K, *lv) for lv in levs]
    wX = omegas[0].real + omegas[1].real - omegas[2].real
    if wX >= K.mu:
        raise ValueError("omega_X >= mu(放出。rate_inf で扱う)")
    omegas, xg, B = _setup(K, levs, omegas, states, xg, r_max_for(K, kappa_decay(K, wX)))
    r = r_of_x(K, xg)
    Ri, Rj, Rk = B[i].R, B[j].R, B[k].R
    sij = 1 if i == j else 2
    adaptive = ells is None
    if adaptive:
        ells = _ells(i, j, k, mX)
    out = {}
    total = charge = 0.0
    for ell in ells:
        A = ang([(i[1], i[2]), (j[1], j[2]), (k[1], -k[2]), (ell, -mX)])
        Bc = ang([(i[1], i[2]), (j[1], j[2]), (k[1], -k[2]), (ell, -mX)], extra_cos2=True)
        if abs(A) < 1e-14 and (adaptive or abs(Bc) < 1e-14):
            continue            # 適応のときは a^2 cos^2 theta の項だけの l(l_i + l_j + l_k + 2)も飛ばす(alpha = 0.45 で 3e-11 以下)
        T = -0.5 * sij * Ri * Rj * np.conj(Rk) * (r * r * A + K.a ** 2 * Bc)
        Rin, Pin = solve_in(K, wX, ell, mX, xg)
        Rup, Pup = solve_up(K, wX, ell, mX, xg)
        W = Rin * Pup - Rup * Pin
        Wm = W[NX // 2]
        Z = simpson(Rup * T * (r - K.rp), x=xg) / Wm
        flux = 2 * (wX - mX * K.OmH) * (K.rp ** 2 + K.a ** 2) * abs(Z) ** 2
        total += flux
        q = _forced_charge(K, wX, ell, mX, xg, Rin, Rup, T, Wm) if details else np.nan
        charge += q
        out[ell] = dict(A=A, B=Bc, flux=flux, charge=q, W_spread=_spread(W, Wm))
        if adaptive and len(out) >= 2 and abs(flux) < ELL_TOL * abs(total) and not abs(q) > 1e-4 * abs(charge):
            break
    kap = total / (alpha ** 8 * K.rp)
    if details:
        return dict(kappa=kap, Gamma=total, tau_X=charge / total, ells=out, omegaX=wX, mu=K.mu,
                    nodes={lv: B[lv].nodes for lv in levs}, match_err={lv: B[lv].match_err for lv in levs},
                    norm={lv: B[lv].norm_integral for lv in levs}, r_max=float(r_of_x(K, xg[-1])))
    return kap


def _ser_div(num, den, n):
    """t のべき級数の割り算 num/den の最初の n 個の係数(den[0] != 0)。"""
    num = np.concatenate([np.asarray(num, complex), np.zeros(n)])
    out = np.zeros(n, complex)
    for i in range(n):
        acc = num[i]
        for j in range(1, min(i, len(den) - 1) + 1):
            acc -= den[j] * out[i - j]
        out[i] = acc / den[0]
    return out


def asym_series(K, omega, l, m, s, nmax=80):
    """無限遠での解 R = e^{s r} r^rho sum_n c_n r^{-n}(s^2 = mu^2 - omega^2。s = +i k が外向き、-kappa が減衰)の
    rho と c_n。R'' + p R' + q R = 0、p = Delta'/Delta、q = V/Delta を t = 1/r の級数にして、係数の漸化式で解く:
        2 s N c_N = (N-1) N c_{N-1} - sum_{j=1}^{N} P_j (N-j) c_{N-j} + sum_{j=2}^{N+1} Q_j c_{N+1-j}、
        P_1 = 2 rho + p_1、P_j = p_j、Q_j = s p_j + rho p_{j-1} + q_j + [j = 2](rho^2 - rho)、rho = -p_1/2 - q_1/(2 s)。"""
    mu, a = K.mu, K.a
    c2 = a * a * (omega * omega - mu * mu)
    c0 = 2 * a * m * omega - a * a * omega * omega - Lam(l, m, c2)
    D = np.array([1.0, -2.0, a * a])
    n = nmax + 3
    p = _ser_div([0.0, 2.0, -2.0], D, n)
    kt = np.array([omega, 0.0, omega * a * a - a * m], complex)
    num = np.convolve(kt, kt) + np.convolve(np.array([-mu * mu, 0.0, c0], complex), D)
    q = _ser_div(num, np.convolve(D, D), n)
    rho = -p[1] / 2 - q[1] / (2 * s)
    P = p.copy()
    P[1] = 2 * rho + p[1]
    Q = s * p + rho * np.concatenate([[0.0], p[:-1]]) + q
    Q[2] += rho * rho - rho
    c = np.zeros(nmax + 1, complex)
    c[0] = 1.0
    with np.errstate(over="ignore", invalid="ignore"):
        for N in range(1, nmax + 1):
            acc = (N - 1) * N * c[N - 1]
            acc -= sum(P[j] * (N - j) * c[N - j] for j in range(1, N + 1))
            acc += sum(Q[j] * c[N + 1 - j] for j in range(2, N + 2))
            c[N] = acc / (2 * s * N)
            if not np.isfinite(c[N]):              # 漸近級数は発散する。あふれたらそこで切る
                return rho, c[:N]
    return rho, c


def asym_bc(K, omega, l, m, s, r1):
    """r1 での (R, P) = (1, Delta R'/R)(漸近級数を、項が小さくなるあいだだけ足す)と、切った項の大きさ。"""
    rho, c = asym_series(K, omega, l, m, s)
    S, dS, last = 1.0 + 0j, 0j, 1.0
    for N in range(1, len(c)):
        term = c[N] * r1 ** (-N)
        if not np.isfinite(term) or abs(term) > last:
            break
        S += term
        dS += -N * c[N] * r1 ** (-N - 1)
        last = abs(term)
        if last < 1e-17:
            break
    D1 = (r1 - K.rp) * (r1 - K.rm)
    return 1.0 + 0j, D1 * (s + rho / r1 + dS / S), last


def solve_out(K, omega, l, m, xg, kr1=None):
    """無限遠で外向き(R ~ r^{-1 + i eta} e^{i k r}、omega > mu)の解。k r1 の所から漸近級数で始めて、格子の外の端まで
    積分し(出力なし)、格子の上を内向きに積分する。返り値: (R, P, 切った項の大きさ)。"""
    k = np.sqrt(omega * omega - K.mu ** 2)
    eta = (2 * omega * omega - K.mu ** 2) / k
    if kr1 is None:
        kr1 = max(400.0, 4.0 * (l * (l + 1) + eta * eta))
    r1 = max(kr1 / k, 2.0 * r_of_x(K, xg[-1]))
    R1, P1, err = asym_bc(K, omega, l, m, 1j * k, r1)
    f = rhs(K, omega, l, m)
    x1 = np.log(r1 - K.rp)
    sol = _ivp(f, (x1, xg[-1]), [R1.real, R1.imag, P1.real, P1.imag])
    y0 = sol.y[:, -1]
    sol = _ivp(f, (xg[-1], xg[0]), y0, t_eval=xg[::-1])
    y = sol.y[:, ::-1]
    return y[0] + 1j * y[1], y[2] + 1j * y[3], err


def rate_inf(alpha, chi, i, j, k, ells=None, states=None, xg=None, omegas=None, details=False, kr1=None):
    """i x j -> k x inf(omega_X > mu)の kappa_eff = Gamma_hat/alpha^5(hyd.kappa_inf と同じ規約、gamma = kappa alpha^8
    (M_pl/f)^4 mu)。Green 関数の解は、源の外で無限遠に Z_inf R_up、地平線で Z_H R_in:
        Z_inf = (1/W) int R_in T dr、 Z_H = (1/W) int R_up T dr。
    無限遠への数の流れは 2 |Z_inf|^2 Im(R_up^* P_up)(実の omega では r によらない。R_up の規格化にもよらない)、
    地平線への流れは 2 (omega - m Omega_H)(r+^2 + a^2) |Z_H|^2(超放射の条件を満たす m_X では負)。
    過程の率は二つの和(X の数の正味の生成)。返り値は kappa(和)、details では無限遠と地平線の内訳と確かめの量も。"""
    K = Kerr(alpha, chi)
    levs = (i, j, k)
    mX = i[2] + j[2] - k[2]
    if omegas is None:
        omegas = [bound_omega(K, *lv) for lv in levs]
    wX = omegas[0].real + omegas[1].real - omegas[2].real
    if wX <= K.mu:
        raise ValueError("omega_X <= mu(束縛。rate_bh で扱う)")
    omegas, xg, B = _setup(K, levs, omegas, states, xg, 0.0)
    r = r_of_x(K, xg)
    Ri, Rj, Rk = B[i].R, B[j].R, B[k].R
    sij = 1 if i == j else 2
    adaptive = ells is None
    if adaptive:
        ells = _ells(i, j, k, mX)
    out = {}
    tot_inf = tot_h = 0.0
    for ell in ells:
        A = ang([(i[1], i[2]), (j[1], j[2]), (k[1], -k[2]), (ell, -mX)])
        Bc = ang([(i[1], i[2]), (j[1], j[2]), (k[1], -k[2]), (ell, -mX)], extra_cos2=True)
        if abs(A) < 1e-14 and (adaptive or abs(Bc) < 1e-14):
            continue            # 適応のときは a^2 cos^2 theta の項だけの l(l_i + l_j + l_k + 2)も飛ばす(alpha = 0.45 で 3e-11 以下)
        T = -0.5 * sij * Ri * Rj * np.conj(Rk) * (r * r * A + K.a ** 2 * Bc)
        Rin, Pin = solve_in(K, wX, ell, mX, xg)
        Rup, Pup, aerr = solve_out(K, wX, ell, mX, xg, kr1=kr1)
        W = Rin * Pup - Rup * Pin
        Wm = W[NX // 2]
        Zinf = simpson(Rin * T * (r - K.rp), x=xg) / Wm
        Zh = simpson(Rup * T * (r - K.rp), x=xg) / Wm
        Jup = np.imag(np.conj(Rup) * Pup)
        f_inf = 2 * abs(Zinf) ** 2 * Jup[-1]
        f_h = 2 * (wX - mX * K.OmH) * (K.rp ** 2 + K.a ** 2) * abs(Zh) ** 2
        tot_inf += f_inf
        tot_h += f_h
        out[ell] = dict(A=A, B=Bc, kappa_inf=f_inf / alpha ** 5, kappa_h=f_h / alpha ** 5,
                        W_spread=_spread(W, Wm), J_spread=_spread(Jup, Jup[-1], scale=np.abs(Rup) * np.abs(Pup)),
                        asym_err=float(aerr))
        if adaptive and len(out) >= 2 and abs(f_inf + f_h) < ELL_TOL * abs(tot_inf + tot_h):
            break
    kap = (tot_inf + tot_h) / alpha ** 5
    if details:
        return dict(kappa=kap, kappa_inf=tot_inf / alpha ** 5, kappa_h=tot_h / alpha ** 5, ells=out, omegaX=wX, mu=K.mu,
                    nodes={lv: B[lv].nodes for lv in levs}, match_err={lv: B[lv].match_err for lv in levs},
                    r_max=float(r_of_x(K, xg[-1])))
    return kap


def point(alpha, chi, procs, omegas, eta_max=30.0):
    """一つの (alpha, chi) で過程をまとめて解く。束縛状態は準位ごとに 1 回だけ、共通の格子で作って使い回す(格子の外の端は
    準位の減衰で決める。束縛の X の減衰がそれより長い過程だけ、自分の格子で作り直す)。
    procs: [(i, j, k), ...]、omegas: {準位: 複素の M omega}(呼ぶ側で根を確かめたもの)。
    返り値: 過程ごとの dict(G_inf、G_h(Gamma_hat の無限遠と地平線の分)、tau_X(束縛の X だけ)、omegaX、確かめの量)と、
    準位ごとの dict(nodes、match_err)。omega_X が閾値のすぐ上(eta = (2 omega^2 - mu^2)/k > eta_max)の過程は解かない
    (外向きの解を始める半径が k^-3 で遠くなる)。"""
    K = Kerr(alpha, chi)
    levels = sorted(omegas)
    r_max = max(r_max_for(K, kappa_decay(K, omegas[lv].real), lv[0]) for lv in levels)
    xg = grid(K, r_max)
    B = {lv: Bound(K, *lv, xg, omega=omegas[lv]) for lv in levels}
    lev_out = {lv: dict(nodes=B[lv].nodes, match_err=B[lv].match_err) for lv in levels}
    out = []
    for i, j, k in procs:
        om = [omegas[i], omegas[j], omegas[k]]
        wX = om[0].real + om[1].real - om[2].real
        res = dict(omegaX=wX, G_inf=np.nan, G_h=np.nan, tau_X=np.nan, W_spread=np.nan, J_spread=np.nan, asym_err=np.nan,
                   note="")
        try:
            if wX > K.mu:
                kX = np.sqrt(wX * wX - K.mu ** 2)
                if (2 * wX * wX - K.mu ** 2) / kX > eta_max:
                    res["note"] = "near threshold"
                    out.append(res)
                    continue
                d = rate_inf(alpha, chi, i, j, k, states=B, xg=xg, omegas=om, details=True)
                res.update(G_inf=d["kappa_inf"] * alpha ** 5, G_h=d["kappa_h"] * alpha ** 5,
                           W_spread=max(v["W_spread"] for v in d["ells"].values()),
                           J_spread=max(v["J_spread"] for v in d["ells"].values()),
                           asym_err=max(v["asym_err"] for v in d["ells"].values()))
            else:
                own = r_max_for(K, kappa_decay(K, wX)) > r_max
                d = (rate_bh(alpha, chi, i, j, k, omegas=om, details=True) if own
                     else rate_bh(alpha, chi, i, j, k, states=B, xg=xg, omegas=om, details=True))
                res.update(G_inf=0.0, G_h=d["Gamma"], tau_X=d["tau_X"],
                           W_spread=max(v["W_spread"] for v in d["ells"].values()), note="own grid" if own else "")
                if own:      # 作り直した束縛状態の確かめ
                    res["match_err"] = max(d["match_err"].values())
                    res["nodes_ok"] = all(d["nodes"][lv] == lv[0] - lv[1] - 1 for lv in (i, j, k))
        except Exception as e:                       # 例外は記録して、その点だけ NaN にする
            res["note"] = "error: " + str(e)[:120]
        out.append(res)
    return out, lev_out


if __name__ == "__main__":
    import json
    what = sys.argv[1] if len(sys.argv) > 1 else "decay"
    if what == "decay":
        for a, c, n in [(0.05, 0.0, 1), (0.1, 0.5, 1), (0.2, 0.9, 1), (0.1, 0.9, 2), (0.2, 0.9, 2), (0.3, 0.5, 2)]:
            print(json.dumps(decay_check(a, c, n)))
    elif what == "bh":
        for a in [0.01, 0.02, 0.05, 0.1]:
            d = rate_bh(a, 0.9, (2, 1, 1), (2, 1, 1), (3, 2, 2), details=True)
            print(a, d["kappa"], {k: (v["flux"], v["W_spread"]) for k, v in d["ells"].items()}, d["nodes"])
    elif what == "inf":
        # 放出の過程の alpha -> 0 の極限: kappa_inf の補正 C = kappa/kappa_hyd を alpha^2 まで当てはめた切片が 1 になるはず
        import hyd
        procs = [((3, 2, 2), (3, 2, 2), (2, 1, 1)), ((4, 1, 1), (4, 1, 1), (2, 1, 1)), ((3, 2, 2), (4, 1, 1), (2, 1, 1)),
                 ((4, 1, 1), (5, 1, 1), (3, 2, 2)), ((5, 1, 1), (5, 1, 1), (3, 2, 2)), ((4, 3, 3), (5, 4, 4), (3, 2, 2)),
                 ((5, 4, 4), (5, 4, 4), (2, 1, 1)), ((4, 2, 2), (5, 4, 4), (3, 2, 2)), ((4, 1, 1), (5, 3, 3), (3, 2, 2))]
        als = np.array([0.005, 0.01, 0.02])
        for pr in procs:
            kh = float(hyd.kappa_inf(*pr)["kappa"])
            C = np.array([rate_inf(a, 0.9, *pr, details=True)["kappa_inf"] / kh for a in als])
            print(pr, "C", np.round(C, 6), "intercept", np.polyfit(als, C, 2)[-1])
