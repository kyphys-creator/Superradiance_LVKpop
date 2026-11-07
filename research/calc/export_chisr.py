"""chi_SR の単体パッケージ(share/chi_sr/)に入れる表を書き出し、単体版がこのプロジェクトの chisr.chi_sr と一致することを確かめる。

出力: share/chi_sr/chi_sr_tables.npz
実行: cd research && ~/miniforge3/bin/python3 calc/export_chisr.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import si  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(HERE, "..", "..", "share", "chi_sr")


def rel_F(n, l, m):
    """連分数法の表 F = Im(M omega) / alpha^(4l+5)(欠けた点は解析式で埋める。sr._rel_table と同じ)。"""
    z = np.load(os.path.join(HERE, "out", f"rates_{n}{l}{m}.npz"))
    al, ch, im = z["alpha"], z["chi"], z["im"].copy()
    F = im / al[:, None] ** (4 * l + 5)
    for i, a in enumerate(al):
        for j, c in enumerate(ch):
            if not np.isfinite(F[i, j]):
                F[i, j] = sr.gamma_over_mu(a, c, n, l, m) * a / a ** (4 * l + 5)
    return al, ch, F


def main():
    tab = chisr.build()
    a211, c211, F211 = rel_F(2, 1, 1)
    a322, c322, F322 = rel_F(3, 2, 2)
    os.makedirs(PKG, exist_ok=True)
    np.savez(os.path.join(PKG, "chi_sr_tables.npz"),
             ALPHAS=chisr.ALPHAS, U=chisr.U, CHI_TOP=chisr.CHI_TOP,
             lnI=tab["lnI"], chic=tab["chic"], Y=tab["Y"], uinv=tab["uinv"],
             rel211_alpha=a211, rel211_chi=c211, rel211_F=F211, rel211_range=np.array(sr._REL_RANGE[(2, 1, 1)]),
             rel322_alpha=a322, rel322_chi=c322, rel322_F=F322, rel322_range=np.array(sr._REL_RANGE[(3, 2, 2)]),
             HBAR_EVS=sr.HBAR_EVS, TSUN=sr.TSUN, N0=sr.N0, YR=sr.YR,
             MPL_GEV=si.MPL_GEV, KAPPA_BH=si.KAPPA_BH, KAPPA_INF=si.KAPPA_INF, F_FREE_322=chisr.F_FREE_322)
    # 単体版と比べる
    sys.path.insert(0, PKG)
    import chi_sr as pkg  # noqa: E402
    model = pkg.ChiSR(os.path.join(PKG, "chi_sr_tables.npz"))
    rng = np.random.default_rng(5)
    n = 200000
    M = rng.uniform(5, 80, n)
    c0 = rng.uniform(0.0, 0.99, n)
    T = 10 ** rng.uniform(5, 10, n)
    mu = rng.uniform(0.3e-12, 1.5e-12, n)
    f = 10 ** rng.uniform(10.5, 19, n)
    c_ref, r_ref = chisr.chi_sr(tab, M, c0, T, mu, f)
    c_pkg, r_pkg = model(M, c0, T, mu, f)
    print("max |chi diff| =", np.max(np.abs(c_ref - c_pkg)), " regime mismatches =", int(np.sum(r_ref != r_pkg)))


if __name__ == "__main__":
    main()
