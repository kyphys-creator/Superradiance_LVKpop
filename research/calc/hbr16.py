"""BBH 合体の残骸のスピン: Hofmann, Barausse & Rezzolla 2016(arXiv:1605.01938)の一般スピン版(研究ノート E10)。

式番号は Refs/notes/astro_inputs.md §1a の抜き書き(本文で照合済み)。G = c = 1、M1 >= M2、q = M2/M1 <= 1。
  r_ISCO, E_ISCO, L_ISCO: 式 (2)--(6)
  a_tot = (a1 cos b + a2 cos g q^2)/(1+q)^2, a_eff = a_tot + xi nu (a1 cos b + a2 cos g): 式 (14)--(15)
  |l| = |L_ISCO(a_eff) - 2 a_tot (E_ISCO(a_eff) - 1) + sum k_ij nu^(1+i) a_eff^j|: 式 (13)
  |a_fin| = (1+q)^-2 [a1^2 + a2^2 q^4 + 2 a1 a2 q^2 cos al + 2 (a1 cos b + a2 q^2 cos g) |l| q + |l|^2 q^2]^(1/2): 式 (16)
係数は Table 1 の上段(n_M = 1、n_J = 2)。k00 は式 (11) の等質量・無スピンの条件 a_fin = 0.68646 で決める。
"""
import numpy as np

K01, K02, K10, K11, K12, XI = -1.2019, -1.20764, 3.79245, 1.18385, 4.90494, 0.41616
K00 = 16.0 * (0.68646 - np.sqrt(3.0) / 2.0) - K10 / 4.0


def r_isco(a):
    a = np.asarray(a, float)
    z1 = 1 + np.cbrt(1 - a ** 2) * (np.cbrt(1 + a) + np.cbrt(1 - a))
    z2 = np.sqrt(3 * a ** 2 + z1 ** 2)
    return 3 + z2 - np.sign(a) * np.sqrt((3 - z1) * (3 + z1 + 2 * z2))


def e_isco(a):
    return np.sqrt(1 - 2.0 / (3.0 * r_isco(a)))


def l_isco(a):
    return 2.0 / (3.0 * np.sqrt(3.0)) * (1 + 2 * np.sqrt(3 * r_isco(a) - 2))


def a_final(q, a1, a2, cosb, cosg, cosal):
    """残骸のスピンの大きさ(式 16)。配列で渡してよい。"""
    q = np.asarray(q, float)
    nu = q / (1 + q) ** 2
    s = a1 * cosb + a2 * cosg
    atot = (a1 * cosb + a2 * cosg * q ** 2) / (1 + q) ** 2
    aeff = atot + XI * nu * s
    ksum = (K00 + K01 * aeff + K02 * aeff ** 2) * nu + (K10 + K11 * aeff + K12 * aeff ** 2) * nu ** 2
    ell = np.abs(l_isco(aeff) - 2 * atot * (e_isco(aeff) - 1) + ksum)
    inner = (a1 ** 2 + a2 ** 2 * q ** 4 + 2 * a1 * a2 * q ** 2 * cosal
             + 2 * (a1 * cosb + a2 * q ** 2 * cosg) * ell * q + ell ** 2 * q ** 2)
    return np.sqrt(inner) / (1 + q) ** 2
