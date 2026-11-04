"""Kerr 時空の質量のあるスカラーの準束縛状態: Leaver の連分数法(研究ノート E4)。

式は Dolan, arXiv:0705.2880(Phys. Rev. D 76, 084001)の §III(式番号は Refs/notes/theory_rates.md の
抜き書きと、この下のコメント)。単位は G = c = hbar = 1、M = 1(r, a は M、omega, mu は 1/M で測る)。
  R(r) = (r-r+)^{-i sigma} (r-r-)^{i sigma + chi - 1} e^{q r} sum a_n ((r-r+)/(r-r-))^n
  sigma = 2 r+ (omega - omega_c)/(r+ - r-),  q = -sqrt(mu^2 - omega^2)(Re q < 0 で束縛状態),
  chi = (mu^2 - 2 omega^2)/q,  b = sqrt(1 - a^2)
  alpha_n a_{n+1} + beta_n a_n + gamma_n a_{n-1} = 0
  alpha_n = n^2 + (c0+1) n + c0, beta_n = -2n^2 + (c1+2) n + c3, gamma_n = n^2 + (c2-3) n + c4
  c0..c4: Dolan の式 eqc0--eqc4。
  条件: beta_0 - alpha_0 gamma_1/(beta_1 - alpha_1 gamma_2/(beta_2 - ...)) = 0。
角度の固有値は Lambda = l(l+1) - c^2 <cos^2 theta>_lm + O(c^4)、c^2 = a^2 (omega^2 - mu^2)。
束縛状態では |c^2| ~ a^2 alpha^4 / n^2 なので、c^4 の項は 1e-5 以下(alpha <= 0.5)。
"""
import numpy as np
from scipy.optimize import newton


def cos2_avg(l, m):
    """<cos^2 theta> for Y_lm。"""
    return (2.0 * l * (l + 1) - 2.0 * m * m - 1.0) / ((2.0 * l - 1.0) * (2.0 * l + 3.0))


def cf_residual(omega, mu, a, l, m, N=4000):
    """連分数の条件の左辺(0 になる omega が準束縛状態)。"""
    b = np.sqrt(1.0 - a * a)
    q = -np.sqrt(mu * mu - omega * omega + 0j)
    if q.real > 0:
        q = -q
    w = omega
    wm = w - a * m / 2.0
    wq = w - 1j * q
    c2a = a * a * (w * w - mu * mu)
    Lam = l * (l + 1) - c2a * cos2_avg(l, m)
    c0 = 1.0 - 2j * w - 2j / b * wm
    c1 = -4.0 + 4j * (w - 1j * q * (1.0 + b)) + 4j / b * wm - 2.0 * (w * w + q * q) / q
    c2 = 3.0 - 2j * w - 2.0 * (q * q - w * w) / q - 2j / b * wm
    c3 = (2j * wq ** 3 / q + 2.0 * wq ** 2 * b + q * q * a * a + 2j * q * a * m - Lam - 1.0
          - wq ** 2 / q + 2.0 * q * b + 2j / b * (wq ** 2 / q + 1.0) * wm)
    c4 = wq ** 4 / (q * q) + 2j * w * wq ** 2 / q - 2j / b * wq ** 2 / q * wm
    n = np.arange(0, N + 1, dtype=float)
    al = n * n + (c0 + 1.0) * n + c0
    be = -2.0 * n * n + (c1 + 2.0) * n + c3
    ga = n * n + (c2 - 3.0) * n + c4
    x = 0.0 + 0j
    for k in range(N, 0, -1):
        x = al[k - 1] * ga[k] / (be[k] - x)
    return be[0] - x


def bound_state(alpha, a, n, l, m, omega0=None, N=4000, tol=1e-15):
    """準束縛状態の M omega(複素数)。n は主量子数(211 なら n=2)。

    初期値は水素原子のエネルギー alpha (1 - alpha^2/(2 n^2))。
    連分数の根は動径の節の数ごとにあるので、初期値で n を選ぶ。
    """
    if omega0 is None:
        omega0 = alpha * (1.0 - alpha ** 2 / (2.0 * n ** 2)) + 0j
    f = lambda w: cf_residual(w, alpha, a, l, m, N)   # noqa: E731
    return newton(f, omega0, x1=omega0 * (1 + 1e-7), tol=tol, maxiter=200)


def cf_residual_mp(omega, mu, a, l, m, N=4000):
    """cf_residual の多倍長版(mpmath)。322・433 のように Im omega / Re omega が 1e-15 より小さいとき用。"""
    import mpmath as mp
    w = mp.mpc(omega)
    mu = mp.mpf(mu)
    a = mp.mpf(a)
    b = mp.sqrt(1 - a * a)
    q = -mp.sqrt(mu * mu - w * w)
    if mp.re(q) > 0:
        q = -q
    j = mp.mpc(0, 1)
    wm = w - a * m / 2
    wq = w - j * q
    Lam = l * (l + 1) - a * a * (w * w - mu * mu) * mp.mpf(cos2_avg(l, m))
    c0 = 1 - 2 * j * w - 2 * j / b * wm
    c1 = -4 + 4 * j * (w - j * q * (1 + b)) + 4 * j / b * wm - 2 * (w * w + q * q) / q
    c2 = 3 - 2 * j * w - 2 * (q * q - w * w) / q - 2 * j / b * wm
    c3 = (2 * j * wq ** 3 / q + 2 * wq ** 2 * b + q * q * a * a + 2 * j * q * a * m - Lam - 1
          - wq ** 2 / q + 2 * q * b + 2 * j / b * (wq ** 2 / q + 1) * wm)
    c4 = wq ** 4 / (q * q) + 2 * j * w * wq ** 2 / q - 2 * j / b * wq ** 2 / q * wm
    x = mp.mpc(0)
    for k in range(N, 0, -1):
        al = (k - 1) ** 2 + (c0 + 1) * (k - 1) + c0
        be = -2 * k * k + (c1 + 2) * k + c3
        ga = k * k + (c2 - 3) * k + c4
        x = al * ga / (be - x)
    return c3 - x       # beta_0 = c3


def bound_state_mp(alpha, a, n, l, m, omega0=None, N=4000, dps=40):
    """bound_state の多倍長版。返り値は mpmath の mpc。"""
    import mpmath as mp
    mp.mp.dps = dps
    if omega0 is None:
        omega0 = alpha * (1.0 - alpha ** 2 / (2.0 * n ** 2))
    f = lambda w: cf_residual_mp(w, alpha, a, l, m, N)   # noqa: E731
    return mp.findroot(f, (mp.mpc(omega0), mp.mpc(omega0 * (1 + 1e-7))), solver="secant", tol=mp.mpf(10) ** (-2 * dps + 10))
