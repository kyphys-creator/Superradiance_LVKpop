"""chi_SR の単体パッケージ(share/chi_sr/)に入れる表を書き出し、単体版を確かめる。

確かめること(E16 で 2--5 を加えた。どれかが外れたら AssertionError で止まる):
  1. 無作為な 20 万個の標本で、単体版とプロジェクト版(chisr.chi_sr)が値も印も一致する。
  2. 臨界のすぐ上と下(alpha 0.02--0.48、211 と 322)で、両者が一致し、印 0(超放射しない)の点では
     211 と 322 の率がどちらも 0 以下である(臨界の定義の食い違いがない)。
  3. tau -> 0 で chi_after が chi0 を返す(表の引き方が厳密な逆になっている)。
  4. 322 から育つ BH が自由場で飽和したときのスピンが、BH の質量の減少を入れた _sat_state と一致する。
  5. 三準位の ODE(fp3.evolve3)と、固定点の領域の 3 例で 0.002 以内。
出力: share/chi_sr/chi_sr_tables.npz
実行: cd research && ~/miniforge3/bin/python3 calc/export_chisr.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import fp3  # noqa: E402
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


def checks(tab, model):
    rep = {}
    # 1. 無作為な標本
    rng = np.random.default_rng(5)
    n = 200000
    M = rng.uniform(5, 80, n)
    c0 = rng.uniform(0.0, 0.99, n)
    T = 10 ** rng.uniform(5, 10, n)
    mu = rng.uniform(0.3e-12, 1.5e-12, n)
    f = 10 ** rng.uniform(10.5, 19, n)
    c_ref, r_ref = chisr.chi_sr(tab, M, c0, T, mu, f)
    c_pkg, r_pkg = model(M, c0, T, mu, f)
    rep["random_max_diff"] = float(np.max(np.abs(c_ref - c_pkg)))
    rep["random_regime_mismatch"] = int(np.sum(r_ref != r_pkg))
    assert rep["random_max_diff"] == 0 and rep["random_regime_mismatch"] == 0, rep

    # 2. 臨界のすぐ上と下(臨界の表の補間の誤差より外側で)
    mu0 = 1e-12
    pts, err = [], 0.0
    for a in [0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.35, 0.42, 0.48]:
        for m in (1, 2):
            c = float(chisr.chic_vec(a, m, tab=tab))
            z = chisr._zero(a, 2, 1, 1) if m == 1 else chisr._zero(a, 3, 2, 2)
            err = max(err, abs(c - z))
            if c < 0.999:
                pts += [(a, c - 2e-5), (a, c + 2e-5)]
    rep["threshold_interp_max_err"] = err
    assert err < 1e-5, rep
    A = np.array([p[0] for p in pts])
    X = np.array([p[1] for p in pts])
    Mb = A / sr.alpha_of(1.0, mu0)
    bad = 0
    for f in (1e12, 1e16, 1e19):
        cr, rr = chisr.chi_sr(tab, Mb, X, 1e10, mu0, f)
        cp, rp = model(Mb, X, 1e10, mu0, f)
        assert np.max(np.abs(cr - cp)) == 0 and np.all(rr == rp), ("threshold band mismatch", f)
        zero = rr == 0
        bad += int(np.sum(zero & ((chisr.g_vec(A, X, 2, 1, 1) > 0) | (chisr.g_vec(A, X, 3, 2, 2) > 0))))
    rep["threshold_points"] = len(pts)
    rep["regime0_but_superradiant"] = bad
    assert bad == 0, rep

    # 3. tau -> 0
    worst = 0.0
    for a in [0.011, 0.05, 0.1, 0.2, 0.29]:
        c2 = float(chisr.chic_vec(a, 2, tab=tab))
        xs = np.linspace(c2 + 1e-6, 0.99, 300)
        worst = max(worst, float(np.max(np.abs(chisr.chi_after(tab, np.full_like(xs, a), xs, 1e-300) - xs))))
    rep["roundtrip_max"] = worst
    assert worst < 1e-12, rep

    # 4. 322 から育つ BH の自由場の飽和
    Ms, cs = np.array([20.0, 30.0, 45.0]), np.array([0.3, 0.5, 0.7])
    cfree, rfree = chisr.chi_sr(tab, Ms, cs, 1e10, 0.71e-12, 1e19)
    _, csat = chisr._sat_state(Ms, cs, np.full(3, 0.71e-12), 2, 3, tab)
    rep["free322_vs_sat_state"] = float(np.max(np.abs(cfree - csat)))
    assert np.all(rfree == 2) and rep["free322_vs_sat_state"] < 1e-12, rep

    # 5. 三準位の ODE
    worst = 0.0
    for M_, c_, f_, Ts in [(25.0, 0.8, 1e16, [1e8, 1e9]), (45.0, 0.7, 1e14, [1e9, 1e10]), (20.0, 0.3, 1e16, [1e9, 1e10])]:
        s = fp3.evolve3(M_, c_, 0.71e-12, f_, Ts[-1], Ts)
        cc, _ = model(M_, c_, np.array(Ts), 0.71e-12, f_)
        worst = max(worst, float(np.max(np.abs(cc - s.y[3]))))
    rep["ode3_max_diff"] = worst
    assert worst < 2e-3, rep
    return rep


def main():
    tab = chisr.build()
    a211, c211, F211 = rel_F(2, 1, 1)
    a322, c322, F322 = rel_F(3, 2, 2)
    a433, c433, F433 = rel_F(4, 3, 3)
    os.makedirs(PKG, exist_ok=True)
    np.savez(os.path.join(PKG, "chi_sr_tables.npz"),
             ALPHAS=chisr.ALPHAS, W=chisr.W, CHI_TOP=chisr.CHI_TOP, lnI=tab["lnI"], J=tab["J"],
             thr_alpha=tab["thr_alpha"], thr1=tab["thr1"], thr2=tab["thr2"],
             rel211_alpha=a211, rel211_chi=c211, rel211_F=F211, rel211_range=np.array(sr._REL_RANGE[(2, 1, 1)]),
             rel322_alpha=a322, rel322_chi=c322, rel322_F=F322, rel322_range=np.array(sr._REL_RANGE[(3, 2, 2)]),
             rel433_alpha=a433, rel433_chi=c433, rel433_F=F433, rel433_range=np.array(sr._REL_RANGE[(4, 3, 3)]),
             HBAR_EVS=sr.HBAR_EVS, TSUN=sr.TSUN, N0=sr.N0, YR=sr.YR,
             MPL_GEV=si.MPL_GEV, KAPPA_BH=si.KAPPA_BH, KAPPA_INF=si.KAPPA_INF, K_C=fp3.K_C, K_D=fp3.K_D,
             ALPHA_MAX=chisr.ALPHA_MAX, F_FP_MAX=chisr.F_FP_MAX)
    sys.path.insert(0, PKG)
    import chi_sr as pkg  # noqa: E402
    model = pkg.ChiSR(os.path.join(PKG, "chi_sr_tables.npz"))
    rep = checks(tab, model)
    # 入力の検査が働くこと
    for bad in [dict(chi0=1.0), dict(chi0=-0.1), dict(M=-1.0), dict(f=np.nan)]:
        args = dict(M=10.0, chi0=0.5, T=1e8, mu=1e-12, f=1e12)
        args.update(bad)
        try:
            model(args["M"], args["chi0"], args["T"], args["mu"], args["f"])
            raise AssertionError(f"input check failed for {bad}")
        except ValueError:
            pass
    rep["input_checks"] = "ok"
    for k, v in rep.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
