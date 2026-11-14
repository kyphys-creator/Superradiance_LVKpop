"""自己相互作用の過程の率を、水素型の波動関数(alpha の主要項)で計算する(研究ノート E21)。

    i x j -> k x BH :  gamma = kappa alpha^11 (M_pl/f)^4 r+ mu     (占有数 eps = N/(G M^2) の率)
    i x j -> k x inf:  gamma = kappa alpha^8  (M_pl/f)^4    mu
  kappa は「過程が 1 回起きる率」の係数(同じ準位が二つ入る過程は、1 回で 2 個減る。fp3.py、BGLS、Witte & Mummery と同じ規約)。

式(導出は研究ノート E21。BGLS = arXiv:2011.11646、WM = arXiv:2412.03655 の方法と同じもの):
  非相対論の場 psi は (i d/dt + nabla^2/(2 mu) + alpha/r) psi = -(lambda/(8 mu^2)) |psi|^2 psi(BGLS 式 12)、lambda = mu^2/f^2。
  psi = sum_n c_n psi_n e^{-i E_n t}、|c_n|^2 = N_n とすると、右辺は振動数 E = E_i + E_j - E_k の源
      S = -(lambda/(8 mu^2)) s_ij c_i c_j c_k^* psi_i psi_j psi_k^*、  s_ij = 1(i = j)、2(i != j)
  を含む。長さを a0 = 1/(mu alpha)、E = -(mu alpha^2/2) e、e = 1/n_i^2 + 1/n_j^2 - 1/n_k^2 と書く。
  (1) e > 0(束縛): 源の l = 0 の成分が作る強制振動は原点で有限で、BH が吸収する。吸収の率は A_H |dpsi(0)|^2
      (A_H = 8 pi r_g^2 r+。n00 の準位の減衰の率 8 alpha^5 r+ mu / n^3 を再現する)。Coulomb の Green 関数
      G(0, r) = -(mu/(2 pi r)) Gamma(1 - nu) W_{nu,1/2}(2 r/(nu a0))、nu = 1/sqrt(e) から
          kappa_BH = s_ij^2 Gamma(1 - nu)^2 A^2 J^2 / (32 pi)、
          A = int dOmega Y_i Y_j Y_k^*、  J = int_0^inf x W_{nu,1/2}(2x/nu) R_i R_j R_k dx。
      nu が整数(n00 の準位に共鳴する。n_j = n_k など)では発散する。その過程は (3) で扱う。
  (2) e < 0(放出): 波数 k = q/a0、q = sqrt(-e) の連続状態(エネルギーで規格化、引力の Coulomb 波 F_l(-1/q, q x))へ
      の遷移の率は 2 pi |<E l M|S>|^2 で、
          kappa_inf = s_ij^2 sum_l A_l^2 J_l^2 / (16 q)、
          A_l = int dOmega Y_i Y_j Y_k^* Y_{l M}^*(M = m_i + m_j - m_k)、  J_l = int_0^inf x F_l(-1/q, q x) R_i R_j R_k dx。
  (3) e = 1/n^2(共鳴): 強制振動の振動数が n00 の準位のエネルギーと主要項で一致する。ずれは微細構造 mu alpha^4 delta で、
      delta = c4(i) + c4(j) - c4(k) - c4(n00)、c4(n, l) = -1/(8 n^4) + (2/n - 6/(2l+1))/n^3(1908.10370 の式 2.28)。
      Green 関数は n00 の極 psi_n00(0) psi_n00(r)/(E - E_n00) で決まり、吸収の率は Gamma_n00 |<n00|S>|^2/(dE)^2
      (Gamma_n00 = 8 alpha^5 r+ mu/n^3 は n00 の減衰の率。幅 Gamma_n00 は dE より alpha の一乗だけ小さいので落とす)。
          gamma = kappa_res alpha^7 (M_pl/f)^4 r+ mu、  kappa_res = s_ij^2 A^2 I^2 / (32 pi n^3 delta^2)、
          I = int_0^inf x^2 R_n0 R_i R_j R_k dx。
R_nl は水素原子の動径の波動関数(a0 = 1、int R^2 x^2 dx = 1)。計算は mpmath の多倍長。
"""
import functools

import mpmath as mp
from sympy import Rational
from sympy.physics.wigner import gaunt

mp.mp.dps = 30


def radial_poly(n, l):
    """R_nl(x) = sum_p c_p x^p e^{-x/n} の係数 c_p(p = l, ..., n-1)を mpf の辞書で返す。"""
    norm = mp.mpf(2) / n ** 2 * mp.sqrt(mp.factorial(n - l - 1) / mp.factorial(n + l))
    out = {}
    k = n - l - 1                       # 一般化 Laguerre 多項式 L_k^{2l+1}(2x/n)
    for j in range(k + 1):
        cj = (-1) ** j * mp.binomial(k + 2 * l + 1, k - j) / mp.factorial(j)
        out[l + j] = norm * cj * (mp.mpf(2) / n) ** (l + j)
    return out


def _mul(p, q):
    out = {}
    for a, ca in p.items():
        for b, cb in q.items():
            out[a + b] = out.get(a + b, 0) + ca * cb
    return out


def triple(i, j, k):
    """R_i R_j R_k = P(x) e^{-s x} の (P の係数の辞書, s)。"""
    P = _mul(_mul(radial_poly(i[0], i[1]), radial_poly(j[0], j[1])), radial_poly(k[0], k[1]))
    s = mp.mpf(1) / i[0] + mp.mpf(1) / j[0] + mp.mpf(1) / k[0]
    return P, s


def _poly(P, x):
    return mp.fsum(c * x ** p for p, c in P.items())


def energy(i, j, k):
    """e = 1/n_i^2 + 1/n_j^2 - 1/n_k^2(正なら束縛 = BH、負なら放出)。"""
    return Rational(1, i[0] ** 2) + Rational(1, j[0] ** 2) - Rational(1, k[0] ** 2)


def ang3(i, j, k):
    """A = int dOmega Y_i Y_j Y_k^*(m_i + m_j = m_k でなければ 0)。"""
    if i[2] + j[2] != k[2]:
        return mp.mpf(0)
    return mp.mpf(str(gaunt(i[1], j[1], k[1], i[2], j[2], -k[2]).evalf(40))) * (-1) ** k[2]


def ang4(i, j, k, l):
    """A_l = int dOmega Y_i Y_j Y_k^* Y_{l M}^*、M = m_i + m_j - m_k。Y_i Y_j を球面調和関数で展開して計算する。"""
    mij = i[2] + j[2]
    M = mij - k[2]
    if abs(M) > l:
        return mp.mpf(0)
    tot = mp.mpf(0)
    for L in range(abs(i[1] - j[1]), i[1] + j[1] + 1):
        if abs(mij) > L:
            continue
        c1 = gaunt(i[1], j[1], L, i[2], j[2], -mij) * (-1) ** mij          # Y_i Y_j の Y_{L, mij} の成分
        if c1 == 0:
            continue
        c2 = gaunt(L, k[1], l, mij, -k[2], -M) * (-1) ** (k[2] + M)        # int Y_{L,mij} Y_k^* Y_{lM}^*
        tot += mp.mpf(str((c1 * c2).evalf(40)))
    return tot


def _breaks(s):
    xs = [mp.mpf(0)]
    for t in (0.5, 2, 5, 10, 20, 40, 80, 160):
        xs.append(mp.mpf(t) / s)
    xs.append(mp.inf)
    return xs


@functools.lru_cache(maxsize=None)
def kappa_bh(i, j, k):
    """i x j -> k x BH の係数。返り値は dict(kappa, e, nu, A, J)。束縛でない、m が合わない、角度の積分が 0 なら kappa = 0。"""
    e = energy(i, j, k)
    A = ang3(i, j, k)
    if e <= 0 or A == 0:
        return dict(kappa=mp.mpf(0), e=e, nu=None, A=A, J=None)
    nu = 1 / mp.sqrt(mp.mpf(e.p) / e.q)
    if abs(nu - mp.nint(nu)) < mp.mpf(10) ** (-20):
        raise ValueError(f"{i} x {j} -> {k} x BH は n00 の準位に共鳴する(nu = {mp.nstr(nu, 6)})。この式では扱えない")
    P, s = triple(i, j, k)
    J = mp.quad(lambda x: x * mp.whitw(nu, mp.mpf(1) / 2, 2 * x / nu) * _poly(P, x) * mp.exp(-s * x), _breaks(s + 1 / nu))
    sij = 1 if i == j else 2
    kap = sij ** 2 * mp.gamma(1 - nu) ** 2 * A ** 2 * J ** 2 / (32 * mp.pi)
    return dict(kappa=kap, e=e, nu=nu, A=A, J=J)


@functools.lru_cache(maxsize=None)
def kappa_inf(i, j, k):
    """i x j -> k x inf の係数。返り値は dict(kappa, e, q, partial = {l: 寄与})。放出でなければ kappa = 0。"""
    e = energy(i, j, k)
    if e >= 0:
        return dict(kappa=mp.mpf(0), e=e, q=None, partial={})
    q = mp.sqrt(mp.mpf(-e.p) / e.q)
    P, s = triple(i, j, k)
    M = i[2] + j[2] - k[2]
    sij = 1 if i == j else 2
    part = {}
    for l in range(abs(M), i[1] + j[1] + k[1] + 1):
        A = ang4(i, j, k, l)
        if A == 0:
            continue
        J = mp.quad(lambda x: x * mp.coulombf(l, -1 / q, q * x) * _poly(P, x) * mp.exp(-s * x), _breaks(s))
        part[l] = sij ** 2 * A ** 2 * J ** 2 / (16 * q)
    return dict(kappa=mp.fsum(part.values()), e=e, q=q, partial=part)


def c4(n, l):
    """準位のエネルギー E/mu の alpha^4 の係数(1908.10370 の式 2.28): -1/(8 n^4) + f_nl/n^3、f_nl = 2/n - 6/(2l+1)。"""
    return -mp.mpf(1) / (8 * n ** 4) + (mp.mpf(2) / n - mp.mpf(6) / (2 * l + 1)) / n ** 3


def resonant_n(i, j, k):
    """i x j -> k x BH が n00 の準位に共鳴するなら、その n。しなければ None。"""
    e = energy(i, j, k)
    if e <= 0:
        return None
    nu = mp.sqrt(mp.mpf(e.q) / e.p)
    n = int(mp.nint(nu))
    return n if abs(nu - n) < mp.mpf(10) ** (-20) else None


@functools.lru_cache(maxsize=None)
def kappa_res(i, j, k):
    """共鳴する過程 i x j -> k x BH の係数(gamma = kappa alpha^7 (M_pl/f)^4 r+ mu)。返り値は dict(kappa, n, delta, A, I)。"""
    n = resonant_n(i, j, k)
    if n is None:
        raise ValueError(f"{i} x {j} -> {k} x BH は共鳴しない")
    A = ang3(i, j, k)
    delta = c4(i[0], i[1]) + c4(j[0], j[1]) - c4(k[0], k[1]) - c4(n, 0)
    P, s = triple(i, j, k)
    Rn = radial_poly(n, 0)
    I = mp.quad(lambda x: x ** 2 * _poly(Rn, x) * mp.exp(-x / n) * _poly(P, x) * mp.exp(-s * x), _breaks(s + mp.mpf(1) / n))
    sij = 1 if i == j else 2
    return dict(kappa=sij ** 2 * A ** 2 * I ** 2 / (32 * mp.pi * n ** 3 * delta ** 2), n=n, delta=delta, A=A, I=I)


def kappa(i, j, k):
    """(係数, "BH" か "inf" か None)。束縛(e > 0)なら BH、放出(e < 0)なら inf。"""
    e = energy(i, j, k)
    if e > 0:
        return kappa_bh(i, j, k)["kappa"], "BH"
    if e < 0:
        return kappa_inf(i, j, k)["kappa"], "inf"
    return mp.mpf(0), None


def decay_n00(n):
    """確かめ用: 吸収の式 A_H |psi_n00(0)|^2 が与える n00 の準位の減衰の率 / (alpha^5 r+ mu)。8/n^3 のはず
    (1908.10370 の式 2.32 で l = 0: 2 Im omega = -2 x 2 r+ C_n0 alpha^5 mu、C_n0 = 2/n^3)。"""
    R0 = radial_poly(n, 0)[0]                     # R_n0(0)
    return 8 * mp.pi * R0 ** 2 / (4 * mp.pi)      # 8 pi r_g^2 r+ |R(0)|^2/(4 pi) を alpha^5 r+ mu の単位で
