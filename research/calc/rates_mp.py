"""超放射の率の表を、多倍長の Leaver の連分数法で作る(値の表の作り直し、DESIGN.md D18 の手順 2)。

E4 の表(倍精度)は、322 の低い alpha で Im omega / Re omega が丸めより小さく(E19)、433 は alpha ~ 0.3 で符号まで誤った
(E17)。表の範囲の外を解析式で埋めたので、切り替えの alpha で率が不連続だった(E19)。ここでは N 準位の運動論の全準位を
多倍長で解く。
  alpha の範囲: 下端は、その準位の超放射が値の表の範囲(mu T <= 0.13 eV yr)で 1 e-fold も稼げない alpha の近く
    (211 は 0.003、322・422 は 0.031、433・533 は 0.11、544 は 0.27 で 1 e-fold。411・511 は 0.004 で 1 e-fold だが、
    211 と同じ m = 1 で臨界が 211 のすぐ上なので、211 が育つあいだにしか育たず、alpha < 0.015 でスピンに効かない)。
    その下は解析式で、表の下端から 1.25 分の 1 の alpha までのあいだで表と解析式の比をなめらかに 1 へつなぐ(sr 側で行う)。
    上端は 0.52(値の表は alpha <= 0.5。alpha = 0.6 では m = 1 の準位が強く減衰し、411 と 511 の根が水素型から外れて
    同じ根に落ちる)。
  連分数の項数 N = max(4000, 500/alpha)(alpha = 0.005 で N と 2N の差 3e-5、alpha >= 0.1 で 1e-12 以下)、40 桁。
  根の取り違えを防ぐ: chi ごとに alpha の小さい方から順に解き、一つ前の alpha の根から初期値を作る。下端の根は水素型の
    エネルギー(sr.omega_over_mu、alpha^5 まで)からの初期値で解き、同じ l の隣の n の準位との間隔の 0.3 倍より水素型に
    近いことを確かめる。次の alpha からは、前の根を水素型の比で延ばした予測と、同じ間隔の 0.3 倍で比べる。
    Im omega の符号が超放射の条件 m Omega_H > Re omega と合うかも確かめる(臨界のごく近くは除く)。外れた点は NaN にして数える。
出力: out/rates_mp_<nlm>.npz(alpha, chi, re, im(M omega、M = 1)、done(chi の列ごと)、bad(確かめに外れた点))。
m = 1 の準位(211・411・511)は alpha >= 0.35 の点を 4 倍にした(refine、2026-11-17。臨界 alpha の近くで率が急に変わる)。
列ごとに途中の結果を保存し、--resume で続ける。
実行: cd research && ~/miniforge3/bin/python3 calc/rates_mp.py build <nlm,...|all> [--resume] [並列の数]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import leaver  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
DPS = 40
CHI = np.concatenate([np.linspace(0.0, 0.9, 91), np.linspace(0.905, 0.995, 19)])      # E4 と同じ 110 点
LEVELS = [(2, 1, 1), (4, 1, 1), (5, 1, 1), (3, 2, 2), (4, 2, 2), (4, 3, 3), (5, 3, 3), (5, 4, 4)]
ALPHA_LO = {(2, 1, 1): 0.004, (4, 1, 1): 0.015, (5, 1, 1): 0.015, (3, 2, 2): 0.03, (4, 2, 2): 0.03,
            (4, 3, 3): 0.10, (5, 3, 3): 0.10, (5, 4, 4): 0.22}
ALPHA_HI = 0.52


def alphas(lev):
    lo = ALPHA_LO[lev]
    n = int(round(50 * np.log10(ALPHA_HI / lo))) + 1          # 対数で 1 桁あたり 50 点
    return np.geomspace(lo, ALPHA_HI, n)


def nterms(a):
    return int(max(4000, 500.0 / a))


def spacing(a, n, l):
    """同じ l の隣の n の準位との、M omega の間隔(水素型、近い方)。"""
    gaps = [abs(1.0 / (2 * n ** 2) - 1.0 / (2 * k ** 2)) for k in (n - 1, n + 1) if k > l]
    return a ** 3 * min(gaps)


def hyd_re(a, n, l, m, c):
    return a * sr.omega_over_mu(a, n, l, m, c)


def column(args):
    """一つの chi について、alpha の小さい方から順に解く。返り値: (re, im, bad) の配列。"""
    lev, c, als = args
    n, l, m = lev
    re = np.full(len(als), np.nan)
    im = np.full(len(als), np.nan)
    bad = np.zeros(len(als), bool)
    prev = None
    t0 = time.time()
    for q, a in enumerate(als):
        h = hyd_re(a, n, l, m, c)
        gim = sr.gamma_over_mu(a, c, n, l, m, 1.0) * a                    # 解析式の Im(M omega)(初期値の虚部)
        if prev is None:
            pred = h
        else:
            ap, wp = prev
            pred = wp.real * h / hyd_re(ap, n, l, m, c)
        tol = 0.3 * spacing(a, n, l)
        w = None
        for guess in ((complex(pred, gim) if prev is None else complex(pred, prev[1].imag * (a / prev[0]) ** (4 * l + 5))),
                      complex(h, gim)):
            try:
                ww = leaver.bound_state_mp(a, c, n, l, m, omega0=guess, N=nterms(a), dps=DPS)
            except Exception:
                continue
            ww = complex(ww)
            if abs(ww.real - pred) < tol or (prev is None and abs(ww.real - h) < tol):
                w = ww
                break
        if w is None:
            bad[q] = True
            prev = None                                                     # 次は水素型からやり直す
            continue
        # Im omega の符号と超放射の条件(m Omega_H > Re omega)
        rp = 1 + np.sqrt(1 - c * c)
        margin = m * c / (2 * rp) - w.real
        if abs(margin) > 1e-3 * abs(w.real) and np.sign(margin) != np.sign(w.imag):
            bad[q] = True
        re[q], im[q] = w.real, w.imag
        prev = (a, w)
    return lev, c, re, im, bad, time.time() - t0


def build(levels, nproc, resume):
    from multiprocessing import Pool
    state = {}
    jobs = []
    for lev in levels:
        path = os.path.join(OUT, "rates_mp_%d%d%d.npz" % lev)
        al = alphas(lev)
        st = dict(path=path, al=al, re=np.full((len(al), len(CHI)), np.nan), im=np.full((len(al), len(CHI)), np.nan),
                  bad=np.zeros((len(al), len(CHI)), bool), done=np.zeros(len(CHI), bool), t0=time.time(), ndone=0)
        if resume and os.path.exists(path):
            z = np.load(path)
            for k in ("re", "im", "bad", "done"):
                st[k] = z[k]
        state[lev] = st
        jobs += [(lev, c, al) for j, c in enumerate(CHI) if not st["done"][j]]
    # 重い列(alpha の下端が小さい準位)から先に
    jobs.sort(key=lambda x: ALPHA_LO[x[0]])
    print(f"{len(jobs)} columns to solve", flush=True)
    t0 = time.time()
    with Pool(nproc) as p:
        for lev, c, re, im, bad, sec in p.imap_unordered(column, jobs):
            st = state[lev]
            j = int(np.argmin(np.abs(CHI - c)))
            st["re"][:, j], st["im"][:, j], st["bad"][:, j], st["done"][j] = re, im, bad, True
            st["ndone"] += 1
            np.savez(st["path"], alpha=st["al"], chi=CHI, re=st["re"], im=st["im"], bad=st["bad"], done=st["done"],
                     n=lev[0], l=lev[1], m=lev[2], dps=DPS)
            if st["done"].all() or st["ndone"] % 22 == 0:
                print(f"{lev}: {int(st['done'].sum())}/{len(CHI)} columns, {time.time() - t0:.0f}s; "
                      f"failed points {int(st['bad'].sum())}, NaN {int(np.sum(~np.isfinite(st['im'][:, st['done']])))}",
                      flush=True)


def refine_column(args):
    """一つの (準位, chi) の列に、alpha の点を足す(refine)。根は既存の表の点から追い、既存の点ごとにその根に合わせ直す。"""
    lev, j, c, al, re0, im0, new = args
    n, l, m = lev
    out = []
    t0 = time.time()
    for a in new:
        i = int(np.searchsorted(al, a)) - 1                         # すぐ下の既存の点
        ap, wp = al[i], complex(re0[i], im0[i])
        if out and out[-1][0] > ap:                                 # 前に解いた新しい点の方が近ければ、それから
            ap, wp = out[-1][0], out[-1][1]
        h = hyd_re(a, n, l, m, c)
        pred = wp.real * h / hyd_re(ap, n, l, m, c)
        guess = complex(pred, wp.imag * (a / ap) ** (4 * l + 5))
        try:
            w = complex(leaver.bound_state_mp(a, c, n, l, m, omega0=guess, N=nterms(a), dps=DPS))
        except Exception:
            w = complex(np.nan, np.nan)
        ok = abs(w.real - pred) < 0.3 * spacing(a, n, l)
        rp = 1 + np.sqrt(1 - c * c)
        margin = m * c / (2 * rp) - w.real
        if ok and abs(margin) > 1e-3 * abs(w.real) and np.sign(margin) != np.sign(w.imag):
            ok = False
        out.append((a, w, not ok))
    return lev, j, out, time.time() - t0


def refine(levels, a_from=0.35, factor=4, nproc=12):
    """alpha >= a_from の区間で alpha の点を factor 倍にする(2026-11-17。m = 1 の準位の臨界 alpha の近く(alpha 0.39-0.47、
    chi >= 0.98)で率が alpha に急に変わり、1 桁 50 点の格子では 3 次の補間が追いつかないため)。前の表は _v1 として残す。"""
    from multiprocessing import Pool
    for lev in levels:
        n, l, m = lev
        path = os.path.join(OUT, f"rates_mp_{n}{l}{m}.npz")
        z = {k: v for k, v in np.load(path).items()}
        al = z["alpha"]
        i0 = int(np.searchsorted(al, a_from * (1 - 1e-12), side="right")) - 1
        new = []
        for i in range(max(i0, 0), len(al) - 1):
            new.extend(np.geomspace(al[i], al[i + 1], factor + 1)[1:-1])
        new = np.array(new)
        jobs = [(lev, j, float(c), al, z["re"][:, j], z["im"][:, j], new) for j, c in enumerate(z["chi"])]
        t0 = time.time()
        with Pool(nproc) as pool:
            res = pool.map(refine_column, jobs, chunksize=1)
        al2 = np.concatenate([al, new])
        order = np.argsort(al2)
        re2 = np.full((len(al2), len(z["chi"])), np.nan)
        im2, bad2 = re2.copy(), np.zeros(re2.shape, bool)
        re2[:len(al)], im2[:len(al)], bad2[:len(al)] = z["re"], z["im"], z["bad"]
        for _, j, out, _ in res:
            for q, (a, w, b) in enumerate(out):
                re2[len(al) + q, j], im2[len(al) + q, j], bad2[len(al) + q, j] = w.real, w.imag, b
        os.replace(path, path.replace(".npz", "_v1.npz"))
        z.update(alpha=al2[order], re=re2[order], im=im2[order], bad=bad2[order], refined_from=a_from, refine_factor=factor)
        np.savez(path, **z)
        print(f"{lev}: {len(new)} alpha points added x {len(z['chi'])} chi, {time.time() - t0:.0f}s; "
              f"failed {int(bad2.sum())}, NaN {int(np.sum(~np.isfinite(im2)))}", flush=True)


def check():
    """多倍長の表と sr.gamma_mp_over_mu の確かめ(2026-11-17)。
    (1) 解けた点と、Im omega の符号と超放射の条件 m Omega_H > Re omega の一致(臨界の 1e-4 alpha より外)。
    (2) E4 の倍精度の表(sr._REL_RANGE の範囲、211・322・411)と、E4 の格子の点で比べる。
    (3) 独立に多倍長で解いた値(E19 の 322・211、E19・E21 の 422)と、sr.gamma_mp_over_mu を比べる。
    (4) 下端の切り替え: 補正 Q の幅と、下端のすぐ上と下での率の飛び(臨界から離れた所)。
    (5) 補間の誤差: alpha と chi の一つ置きの節を抜いて sr と同じ形のスプラインを作り、抜いた節で比べる(Im / alpha^(4l+5) の
        局所の最大との比)。節の間隔が倍なので、3 次の補間の誤差は実際の格子の 16 倍ほどの見積もり。
    (6) Re omega と水素型(alpha^5 まで、sr.omega_over_mu)の差。"""
    import json
    from scipy.interpolate import RectBivariateSpline
    out = {}
    for lev in LEVELS:
        n, l, m = lev
        z = np.load(os.path.join(OUT, f"rates_mp_{n}{l}{m}.npz"))
        al, ch, re, im = z["alpha"], z["chi"], z["re"], z["im"]
        p = 4 * l + 5
        s = m * sr.omegaH(ch)[None, :] - re
        far = np.abs(s) > 1e-4 * al[:, None]
        row = dict(points=int(im.size), bad=int(z["bad"].sum()), nan=int(np.sum(~np.isfinite(im))),
                   sign_mismatch=int(np.sum(far & (np.sign(im) != np.sign(s)))))
        F = im / al[:, None] ** p
        # (5) 一つ置きに抜いた節で(sr と同じ形: G = F / (m Omega_H - Re) と Re を補間して掛ける)
        with np.errstate(divide="ignore", invalid="ignore"):
            G = F / s
        for i in range(len(al)):
            g = np.isfinite(G[i]) & (np.abs(s[i]) > 1e-6 * al[i])
            G[i] = np.interp(ch, ch[g], G[i][g])
        la = np.log(al)
        loo = []
        for axis in (0, 1):
            keep = np.ones(F.shape[axis], bool)
            keep[1:-1:2] = False
            idx = np.flatnonzero(~keep)
            if axis == 0:
                pg = RectBivariateSpline(la[keep], ch, G[keep], kx=3, ky=3, s=0)(la[idx], ch)
                pr = RectBivariateSpline(la[keep], ch, re[keep], kx=3, ky=3, s=0)(la[idx], ch)
                pred = pg * (m * sr.omegaH(ch)[None, :] - pr)
                true = F[idx]
                scale = np.maximum.reduce([np.abs(F[idx - 1]), np.abs(F[idx + 1]), np.abs(true)])
            else:
                pg = RectBivariateSpline(la, ch[keep], G[:, keep], kx=3, ky=3, s=0)(la, ch[idx])
                pr = RectBivariateSpline(la, ch[keep], re[:, keep], kx=3, ky=3, s=0)(la, ch[idx])
                pred = pg * (m * sr.omegaH(ch[idx])[None, :] - pr)
                true = F[:, idx]
                scale = np.maximum.reduce([np.abs(F[:, idx - 1]), np.abs(F[:, idx + 1]), np.abs(true)])
            e = np.abs(pred - true) / scale
            loo.append(dict(median=float(np.median(e)), p99=float(np.percentile(e, 99)), max=float(e.max())))
        row["interp_every_other_node"] = dict(alpha=loo[0], chi=loo[1])
        # (6)
        hyd = np.array([[a * sr.omega_over_mu(a, n, l, m, c) for c in ch] for a in al])
        row["re_vs_alpha5"] = float(np.max(np.abs(re - hyd) / re))
        # (4) 下端の切り替え
        T = sr._mp_table(n, l, m)
        cc = np.linspace(0.0, 0.995, 400)
        q = T["Qlo"](cc)
        a0 = T["alo"]
        up = np.array([sr.gamma_mp_over_mu(a0 * (1 + 1e-9), c, n, l, m) for c in cc])
        dn = np.array([sr.gamma_mp_over_mu(a0 * (1 - 1e-9), c, n, l, m) for c in cc])
        sc = m * sr.omegaH(cc) - a0 * np.array([sr.omega_over_mu(a0, n, l, m, c) for c in cc])
        farc = np.abs(sc) > 0.05 * m * sr.omegaH(cc).max()
        row["alpha_lo"] = a0
        row["Q_at_alpha_lo"] = [float(q.min()), float(q.max())]
        row["jump_at_alpha_lo_far_from_critical"] = float(np.max(np.abs(up[farc] - dn[farc]) / np.abs(up[farc])))
        row["jump_at_alpha_lo_abs_over_max_rate"] = float(np.max(np.abs(up - dn)) / np.max(np.abs(up)))
        out[f"{n}{l}{m}"] = row
    # (2) E4 の表
    e4 = {}
    for lev, (lo, hi) in sr._REL_RANGE.items():
        n, l, m = lev
        z = np.load(os.path.join(OUT, f"rates_{n}{l}{m}.npz"))
        rows = []
        for i, a in enumerate(z["alpha"]):
            if not (lo <= a <= min(hi, ALPHA_HI)) or a < ALPHA_LO[lev]:
                continue
            for j, c in enumerate(z["chi"]):
                if not np.isfinite(z["im"][i, j]):
                    continue
                sc = m * sr.omegaH(c) - a * sr.omega_over_mu(a, n, l, m, c)
                if abs(sc) < 1e-3 * a:
                    continue
                v = sr.gamma_mp_over_mu(a, c, n, l, m) * a
                rows.append(abs(v / z["im"][i, j] - 1))
        rows = np.array(rows)
        e4[f"{n}{l}{m}"] = dict(n=int(len(rows)), median=float(np.median(rows)), p99=float(np.percentile(rows, 99)),
                                max=float(rows.max()))
    out["vs_E4_tables"] = e4
    # (3) 独立な多倍長の値
    ind = []
    e19 = json.load(open(os.path.join(OUT, "e19.json")))
    for r in e19["rates"]["rate322_analytic_vs_mp"]:
        ind.append(("322", r["alpha"], r["chi"], r["im_mp"]))
    for r in e19["rates"]["rate211_analytic_vs_mp"]:
        ind.append(("211", r["alpha"], r["chi"], r["im_mp"]))
    for r in e19["level422"]["rate422_mp"]:
        ind.append(("422", r["alpha"], r["chi"], r["im422_mp"]))
        ind.append(("322", r["alpha"], r["chi"], r["im322_mp"]))
    for r in json.load(open(os.path.join(OUT, "e21.json")))["sr422"]["rows"]:
        ind.append(("422", r["alpha"], r["chi"], r["im422"]))
        ind.append(("322", r["alpha"], r["chi"], r["im322"]))
    res = {}
    for nm, a, c, v in ind:
        lev = tuple(int(x) for x in nm)
        if a < ALPHA_LO[lev]:
            continue                                   # 表の下端より下(解析式とのつなぎの所)
        d = abs(sr.gamma_mp_over_mu(a, c, *lev) * a / v - 1)
        r = res.setdefault(nm, dict(n=0, max=0.0, where=None))
        r["n"] += 1
        if d > r["max"]:
            r["max"], r["where"] = float(d), [a, c]
    out["vs_independent_mp"] = res
    return out


if __name__ == "__main__":
    nproc = int(sys.argv[-1]) if sys.argv[-1].isdigit() else 10
    if sys.argv[1] == "build":
        levs = LEVELS if sys.argv[2] == "all" else [tuple(int(x) for x in s) for s in sys.argv[2].split(",")]
        build(levs, nproc, "--resume" in sys.argv)
    elif sys.argv[1] == "refine":
        levs = [tuple(int(x) for x in t) for t in sys.argv[2].split(",")]
        refine(levs, nproc=nproc)
    elif sys.argv[1] == "check":
        import json
        print(json.dumps(check(), indent=1))
    elif sys.argv[1] == "test":
        # 1 列ずつ、準位ごとに chi = 0.9 の列を alpha の上 8 点だけ解く(動作の確かめ)
        from multiprocessing import Pool
        with Pool(8) as p:
            for lev, c, re, im, bad, sec in p.map(column, [(lev, 0.9, alphas(lev)[::8]) for lev in LEVELS]):
                print(lev, f"{sec:.0f}s", "bad", int(bad.sum()), "NaN", int(np.sum(~np.isfinite(im))), "of", len(im))
