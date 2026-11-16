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
"""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.special import gammaln, lpmv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import leaver  # noqa: E402
import sr  # noqa: E402

DELTA0 = 1e-8          # 地平線から内側の端までの r - r+(r+ を掛ける)
NX = 24000             # x の格子の点の数


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
    """y = (R, P) の x についての右辺。src(r) は (r - r+) を掛ける前の源(複素)。"""
    mu, a = K.mu, K.a
    c2 = a * a * (omega * omega - mu * mu)
    L_ = Lam(l, m, c2)

    def f(x, y):
        r = K.rp + np.exp(x)
        d = r - K.rp
        Kk = (r * r + a * a) * omega - a * m
        U = Kk * Kk / (r - K.rm) + d * (2 * a * m * omega - a * a * omega * omega - mu * mu * r * r - L_)
        R, P = y[0] + 1j * y[1], y[2] + 1j * y[3]
        dR = P / (r - K.rm)
        dP = -U * R + (d * src(r) if src is not None else 0.0)
        return [dR.real, dR.imag, dP.real, dP.imag]
    return f


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
    sol = solve_ivp(rhs(K, omega, l, m), (xg[-1], xg[0]), [R1.real, R1.imag, P1.real, P1.imag], method="DOP853",
                    t_eval=xg[::-1], rtol=1e-11, atol=1e-30)
    y = sol.y[:, ::-1]
    return y[0] + 1j * y[1], y[2] + 1j * y[3]


def solve_in(K, omega, l, m, xg):
    """地平線で内向きの解(外向きに積分)。R -> (r - r+)^{-i sigma}。"""
    sig = 2 * K.rp * (omega - m * K.OmH) / (K.rp - K.rm)
    d0 = np.exp(xg[0])
    R0 = np.exp(-1j * sig * xg[0])
    P0 = (d0 + K.rp - K.rm) * (-1j * sig) * R0
    sol = solve_ivp(rhs(K, omega, l, m), (xg[0], xg[-1]), [R0.real, R0.imag, P0.real, P0.imag], method="DOP853",
                    t_eval=xg, rtol=1e-11, atol=1e-30)
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
        sel = np.flatnonzero(r > 4 * K.rp)
        jm = sel[np.argmax(np.abs(Rd[sel]))]
        cc = Rd[jm] / Rh[jm]
        self.match_err = float(abs(Pd[jm] - cc * Ph[jm]) / abs(Pd[jm]))
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


def rate_bh(alpha, chi, i, j, k, ells=(0, 2, 4), states=None, details=False):
    """i x j -> k x BH(m_X = 0、omega_X < mu)の kappa_eff。i, j, k は (n, l, m)。"""
    K = Kerr(alpha, chi)
    levs = (i, j, k)
    mX = i[2] + j[2] - k[2]
    if mX != 0:
        raise ValueError("m_X != 0")
    # 束縛状態(同じ格子で)
    omegas = [bound_omega(K, *lv) for lv in levs]
    wX = omegas[0].real + omegas[1].real - omegas[2].real
    kapX = kappa_decay(K, wX).real
    kaps = [kappa_decay(K, w.real).real for w in omegas]
    r_max = max(r_max_for(K, kappa_decay(K, wX)), max(r_max_for(K, kappa_decay(K, w.real), lv[0]) for w, lv in zip(omegas, levs)))
    xg = grid(K, r_max)
    r = r_of_x(K, xg)
    B = states or {}
    for lv, w in zip(levs, omegas):
        if lv not in B:
            B[lv] = Bound(K, *lv, xg, omega=w)
    Ri, Rj, Rk = B[i].R, B[j].R, B[k].R
    sij = 1 if i == j else 2
    parity = (i[1] + j[1] + k[1]) % 2
    out = {}
    total = 0.0
    for ell in ells:
        if ell % 2 != parity:
            continue
        A = ang([(i[1], i[2]), (j[1], j[2]), (k[1], -k[2]), (ell, 0)])
        Bc = ang([(i[1], i[2]), (j[1], j[2]), (k[1], -k[2]), (ell, 0)], extra_cos2=True)
        if abs(A) < 1e-14 and abs(Bc) < 1e-14:
            continue
        T = -0.5 * sij * Ri * Rj * np.conj(Rk) * (r * r * A + K.a ** 2 * Bc)
        Rin, Pin = solve_in(K, wX, ell, 0, xg)
        Rup, Pup = solve_up(K, wX, ell, 0, xg)
        W = Rin * Pup - Rup * Pin
        Wm = W[NX // 2]
        Z = simpson(Rup * T * (r - K.rp), x=xg) / Wm
        flux = 2 * wX * (K.rp ** 2 + K.a ** 2) * abs(Z) ** 2
        total += flux
        out[ell] = dict(A=A, B=Bc, flux=flux, W_spread=float(np.std(W[NX // 10: 9 * NX // 10]) / abs(Wm)))
    kap = total / (alpha ** 8 * K.rp)
    if details:
        return dict(kappa=kap, ells=out, omegaX=wX, mu=K.mu, nodes={lv: B[lv].nodes for lv in levs},
                    norm={lv: B[lv].norm_integral for lv in levs}, r_max=r_max)
    return kap


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
