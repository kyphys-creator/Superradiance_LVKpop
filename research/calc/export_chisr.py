"""chi_SR の単体パッケージ(share/chi_sr/)に入れる表を書き出し、単体版を確かめる。

確かめること(E16 で 2--5 を加えた。どれかが外れたら AssertionError で止まる):
  1. 無作為な 20 万個の標本で、単体版とプロジェクト版(chisr.chi_sr)が値も印も一致する。
     1b. E18 の独立な照合(seed 1802)の入力すべてで、両者が一致する(保存してある照合の出力は E18 の版の値なので、
         E19 で印を変えたあとは、プロジェクト版をその場で計算して比べる)。
     1c. E19 で印の付け方を変えた所(率の切り替えを挟む alpha のセル、alpha < 0.012)に置いた 10 万個でも両者が一致する。
         そこで実際に E19 の規則が働いていることも見る(印 5 と、印 3 から印 1 に変わった点がある)。定数も同じ。
  2. 臨界のすぐ上と下(alpha 0.02--0.48、211・322・433)で、両者が一致し、印 0(超放射しない)の点では
     211 と 322 の率がどちらも 0 以下で、433 も年齢のうちに飽和しない(E17)。m = 3 の印の例も確かめる。
  3. tau -> 0 で chi_after が chi0 を返す(表の引き方が厳密な逆になっている)。
  4. 値の表の節の上で、基準の運動論(fp3.evolve3m)を表の分解能 0.01 以内で再現する(E18。全域の照合は
     calc/e18_validate2.py と calc/e19_validate.py)。
  5. 単体版が pickle でき(中身は表の場所だけ)、戻したものが同じ値を返す。配列で呼んでも警告が出ない(E19)。
出力: share/chi_sr/chi_sr_tables.npz(中身が変わらないときは書き直さない)、
      share/chi_sr/selftest.npz(渡した先の計算機で値と印を照らすための入力 5000 個と出力。share/chi_sr/selftest.py が使う。E19)
基準の表 out/reftab.npz(700 MB、git に入れない)がない計算機では、chisr.build() が書き出し済みの表から同じ数の辞書を作る
ので、値の表の部分は読んだものをそのまま書き戻すことになる(E19)。値の表を作り直したときは、基準の表のある計算機で走らせる。
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


def rt_arrays(rt):
    """値の表を単体版の形に: log10 tau は uint16 の符号の、準位どうしの差(届いた準位の数 rt_nreach より先は 0)。"""
    import reftab
    codes = reftab.quantize_lt(rt["LT"]).astype(np.int64)
    finite = (codes != 65535)
    nreach = finite.sum(-1)
    assert np.all(finite == (np.arange(codes.shape[-1]) < nreach[..., None])), "届いた準位が先頭に続いていない"
    c = np.where(finite, codes, 0)
    c[..., 0] = 0
    d = np.diff(c, axis=-1, prepend=0)
    d = np.where(np.arange(codes.shape[-1]) < nreach[..., None], d, 0)
    assert d.min() >= 0 and d.max() < 65536
    return {"rt_LTd": d.astype(np.uint16), "rt_nreach": nreach.astype(np.uint16), "rt_RG": rt["RG"].astype(np.float16),
            "rt_KI": rt["KI"], "rt_status": rt["status"].astype(np.int8), "rt_tau_valid": rt["tau_valid"],
            "rt_A_NODES": rt["A_NODES"], "rt_SA": rt["SA"], "rt_SB": rt["SB"], "rt_LF_NODES": rt["LF_NODES"],
            "rt_M_NODES": rt["M_NODES"], "rt_TAU_MAX": rt["TAU_MAX"], "rt_P_LEV": reftab.P_LEV2, "rt_G_LT": reftab.G_LT}


def strict(model):
    """単体版を、警告を例外にして呼ぶ(配列で呼んでも警告を出さないことを、確かめの入力すべてで見る。E19)。"""
    import warnings

    def call(*args):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            return model(*args)
    return call


def checks(tab, model_raw):
    rep = {}
    model = strict(model_raw)
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

    # 1b. 独立な照合(E18、seed 1802)の入力すべてで、プロジェクト版と同じ値・同じ印(照合の結果が単体版にそのまま移る)。
    #     保存してある出力(E18 の版)とは E19 で印が変わった所だけ違うので、プロジェクト版をその場で計算して比べる
    vpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "e18_validate2_1802.npz")
    if os.path.exists(vpath):
        v = np.load(vpath)
        Mv, Tv = np.meshgrid(v["M"], v["T"], indexing="ij")
        args = [np.broadcast_to(x[:, None], Mv.shape).ravel() for x in (v["chi0"], v["mu"], v["f"])]
        c_pkg, r_pkg = model(Mv.ravel(), args[0], Tv.ravel(), args[1], args[2])
        c_ref, r_ref = chisr.chi_sr(tab, Mv.ravel(), args[0], Tv.ravel(), args[1], args[2])
        rep["validation_inputs"] = int(Mv.size)
        rep["validation_max_diff"] = float(np.max(np.abs(c_pkg - c_ref)))
        rep["validation_regime_mismatch"] = int(np.sum(r_pkg != r_ref))
        rep["validation_changed_since_e18"] = int(np.sum((r_pkg != v["regime"].ravel()) | (c_pkg != v["chisr"].ravel())))
        assert rep["validation_max_diff"] == 0 and rep["validation_regime_mismatch"] == 0, rep

    # 1c. E19 で印の付け方を変えた所: 率の切り替えを挟む alpha のセルと、alpha < 0.012(半解析の固定点の表の下端 0.01 の両側)
    import reftab
    assert (pkg.LOW_MARGIN, pkg.KS_STEP, pkg.KS_MIN, pkg.KS_CAP, pkg.THETA_SW, pkg.THETA_P) == \
        (chisr.LOW_MARGIN, reftab.KS_STEP, reftab.KS_MIN, reftab.KS_CAP, reftab.THETA_SW, reftab.THETA_P), "E19 の定数が単体版と違う"
    an = tab["rt"]["A_NODES"]
    n = 100000
    rng = np.random.default_rng(19)
    kind = rng.integers(0, 3, n)
    edge = np.where(kind == 0, 0.03, 0.1)
    i = np.searchsorted(an, edge) - 1
    a = np.where(kind == 2, np.exp(rng.uniform(np.log(0.004), np.log(0.012), n)),
                 np.exp(rng.uniform(np.log(an[i] * 0.97), np.log(an[i + 1] * 1.03))))
    M = np.exp(rng.uniform(np.log(1), np.log(300), n))
    mu = a / sr.alpha_of(M, 1.0)
    c0 = rng.uniform(0.0, 0.995, n)
    T = np.minimum(10 ** rng.uniform(5, np.log10(1.4e10), n), 0.12 / mu)
    f = 10 ** rng.uniform(6, 21, n)
    c_ref, r_ref = chisr.chi_sr(tab, M, c0, T, mu, f)
    c_pkg, r_pkg = model(M, c0, T, mu, f)
    rep["e19_max_diff"] = float(np.max(np.abs(c_ref - c_pkg)))
    rep["e19_regime_mismatch"] = int(np.sum(r_ref != r_pkg))
    sw = chisr.switch_cell(an, a)
    rep["e19_switch_cell_regime_share"] = {str(k): float(np.mean(r_pkg[sw] == k)) for k in range(6)}
    low = a < chisr.ALPHAS[0]
    rep["e19_low_alpha_regime_share"] = {str(k): float(np.mean(r_pkg[low] == k)) for k in range(6)}
    assert rep["e19_max_diff"] == 0 and rep["e19_regime_mismatch"] == 0, rep
    assert np.any(r_pkg[sw] == 5) and np.any(r_pkg[sw] == 2), rep      # 切り替えのセルで、印 5 になる点も印 2 のままの点もある
    assert np.any(r_pkg[low] == 1) and np.any(r_pkg[low] == 3), rep    # alpha < 0.01 で、印 1(値を返す)も印 3 もある

    # 2. 臨界のすぐ上と下(臨界の表の補間の誤差より外側で)
    mu0 = 1e-12
    pts, err = [], 0.0
    for a in [0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.35, 0.42, 0.48]:
        for m in (1, 2, 3):
            c = float(chisr.chic_vec(a, m, tab=tab))
            z = chisr._zero(a, *{1: (2, 1, 1), 2: (3, 2, 2), 3: (4, 3, 3)}[m])
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
        # 印 0 の点で 433(m = 3)が年齢(1e10 年)のうちに飽和しうるものがない(E17)
        g3 = chisr.g_vec(A, np.minimum(X, 0.995), 4, 3, 3) * mu0 / sr.HBAR_EVS
        c3 = chisr.chic_vec(A, 3, tab=tab)
        lnN3 = np.log(np.maximum((X - c3) / 3, 1e-300)) + np.log(sr.N0 * Mb ** 2)
        bad += int(np.sum(zero & (X > c3) & (g3 * 1e10 * sr.YR > lnN3)))
    rep["threshold_points"] = len(pts)
    rep["regime0_but_superradiant"] = bad
    assert bad == 0, rep

    # 2b. m = 3 の印の例(E17): 50 Msun、chi0 = 0.485(433 の臨界 0.343 と 322 の臨界 0.495 のあいだ)、6 Gyr
    cm3, rm3 = model(np.array([50.0, 45.0]), np.array([0.485, 0.381]), 6e9, 0.71e-12, 1e12)
    rep["m3_flag_examples"] = [int(v) for v in rm3]
    assert list(rm3) == [3, 0], rep

    # 3. tau -> 0
    worst = 0.0
    for a in [0.011, 0.05, 0.1, 0.2, 0.29]:
        c2 = float(chisr.chic_vec(a, 2, tab=tab))
        xs = np.linspace(c2 + 1e-6, 0.99, 300)
        worst = max(worst, float(np.max(np.abs(chisr.chi_after(tab, np.full_like(xs, a), xs, 1e-300) - xs))))
    rep["roundtrip_max"] = worst
    assert worst < 1e-12, rep

    # 4. 表の節の上で、基準の運動論(evolve3m)を表の分解能 0.01 以内で再現する(E18)。進み具合の準位は段の幅の 1% 刻みで、
    #    最初の準位を越えるまでは初期のスピンを返すので、節の上でも段の幅の 1%(スピンで 0.01 以下)まではずれる
    import reftab
    rt = tab["rt"]
    worst = 0.0
    for (pc, ia, js, kf, il) in [(1, 20, 12, 30, 0), (1, 30, 6, 36, 1), (0, 25, 10, 24, 1), (1, 14, 20, 16, 0)]:
        a = rt["A_NODES"][ia]
        s_ = rt["SA"][js] if pc == 0 else 1 + rt["SB"][js]
        Mn = rt["M_NODES"][il]
        c1 = float(chisr.chic_vec(a, 1, tab=tab))
        c2 = float(chisr.chic_vec(a, 2, tab=tab))
        x0 = float(chisr.chi_of_w(s_, c1, c2))
        mu = a / sr.alpha_of(Mn, 1.0)
        Ts = np.logspace(3, 10, 15)
        s = fp3.evolve3m(Mn, x0, mu, 10 ** rt["LF_NODES"][kf], Ts[-1], Ts)
        cc, rr = model(Mn, x0, Ts, mu, 10 ** rt["LF_NODES"][kf])
        ref = s.y[4] / s.y[3] ** 2
        use = np.isin(rr, (1, 2, 4, 5))
        if use.any():
            worst = max(worst, float(np.max(np.abs(cc - ref)[use])))
    rep["on_node_vs_evolve3m"] = worst
    assert worst < 1e-2, rep

    # 5. pickle(中身は表の場所だけ。戻すと同じプロセスでは同じ表を一度だけ読む)。E19
    #    配列で呼んだときに警告が出ないことは、上の単体版の呼び出しすべてで見ている(strict)
    import pickle
    blob = pickle.dumps(model_raw)
    rep["pickle_bytes"] = len(blob)
    assert len(blob) < 10000, rep
    m2 = pickle.loads(blob)
    assert pickle.loads(blob) is m2, "同じプロセスで表を二度読んでいる"
    Mt, ct, Tt = rng.uniform(5, 80, 2000), rng.uniform(0.0, 0.99, 2000), 10 ** rng.uniform(5, 10, 2000)
    ca, ra = model(Mt, ct, Tt, 0.71e-12, 1e18)
    cb, rb = strict(m2)(Mt, ct, Tt, 0.71e-12, 1e18)
    rep["pickle_roundtrip_max_diff"] = float(np.max(np.abs(ca - cb)))
    assert rep["pickle_roundtrip_max_diff"] == 0 and np.all(ra == rb), rep
    return rep


def write_selftest(model, n=5000):
    """渡した先の計算機で照らすための入力と出力(share/chi_sr/selftest.npz)。領域を一通り含むようにとる(E19)。"""
    rng = np.random.default_rng(190)
    k = rng.integers(0, 5, n)
    M = np.where(k == 0, np.exp(rng.uniform(np.log(1), np.log(300), n)), np.where(k == 1, rng.uniform(15, 25, n), rng.uniform(5, 80, n)))
    mu = np.where(k == 0, np.exp(rng.uniform(np.log(1e-13), np.log(1e-11), n)), rng.uniform(0.3e-12, 1.5e-12, n))
    c0 = np.where(k == 2, np.minimum(rng.beta(0.69, 2.28, n), 0.9949), rng.uniform(0.0, 0.995, n))
    T = np.minimum(10 ** rng.uniform(5, np.log10(1.4e10), n), 0.12 / mu)
    f = 10 ** np.where(k == 3, rng.uniform(11, 13, n), np.where(k == 4, rng.uniform(16, 20.5, n), rng.uniform(6, 22, n)))
    chi, reg = strict(model)(M, c0, T, mu, f)
    path = os.path.join(PKG, "selftest.npz")
    new = dict(M=M, chi0=c0, T=T, mu=mu, f=f, chi=chi, regime=reg.astype(np.int8))
    unchanged = False
    if os.path.exists(path):
        with np.load(path) as old:
            unchanged = set(old.files) == set(new) and all(_same(old[x], new[x]) for x in new)
    if not unchanged:
        np.savez_compressed(path, **new)
    return {"n": n, "regime_share": [round(float(np.mean(reg == j)), 3) for j in range(6)], "file": "unchanged" if unchanged else "written"}


def _same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b, equal_nan=a.dtype.kind == "f")


def main():
    global pkg
    tab = chisr.build()
    print(f"value table source: {tab['rt_source']}"
          + ("(書き出し済みの表から読んだ。値の表の部分はそのまま書き戻す)" if tab["rt_source"] == "exported" else ""))
    a211, c211, F211 = rel_F(2, 1, 1)
    a322, c322, F322 = rel_F(3, 2, 2)
    os.makedirs(PKG, exist_ok=True)
    new = dict(ALPHAS=chisr.ALPHAS, W=chisr.W, CHI_TOP=chisr.CHI_TOP, lnI=tab["lnI"], J=tab["J"],
               thr_alpha=tab["thr_alpha"], thr1=tab["thr1"], thr2=tab["thr2"], thr3=tab["thr3"],
               rel211_alpha=a211, rel211_chi=c211, rel211_F=F211, rel211_range=np.array(sr._REL_RANGE[(2, 1, 1)]),
               rel322_alpha=a322, rel322_chi=c322, rel322_F=F322, rel322_range=np.array(sr._REL_RANGE[(3, 2, 2)]),
               HBAR_EVS=sr.HBAR_EVS, TSUN=sr.TSUN, N0=sr.N0, YR=sr.YR,
               MPL_GEV=si.MPL_GEV, KAPPA_BH=si.KAPPA_BH, KAPPA_INF=si.KAPPA_INF, K_C=fp3.K_C, K_D=fp3.K_D,
               ALPHA_MAX=chisr.ALPHA_MAX, F_FP_MAX=chisr.F_FP_MAX, **rt_arrays(tab["rt"]))
    path = os.path.join(PKG, "chi_sr_tables.npz")
    unchanged = False
    if os.path.exists(path):                     # 中身が同じなら書き直さない(zip の時刻だけが変わって git に差分が出るのを防ぐ)
        with np.load(path) as old:
            unchanged = set(old.files) == set(new) and all(_same(old[k], new[k]) for k in new)
    if not unchanged:
        np.savez_compressed(path, **new)
    print("chi_sr_tables.npz: " + ("unchanged (not rewritten)" if unchanged else "written"))
    sys.path.insert(0, PKG)
    import chi_sr as pkg  # noqa: E402
    model = pkg.ChiSR(path)
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
    rep["selftest"] = write_selftest(model)
    for k, v in rep.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
