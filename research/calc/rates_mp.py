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


if __name__ == "__main__":
    nproc = int(sys.argv[-1]) if sys.argv[-1].isdigit() else 10
    if sys.argv[1] == "build":
        levs = LEVELS if sys.argv[2] == "all" else [tuple(int(x) for x in s) for s in sys.argv[2].split(",")]
        build(levs, nproc, "--resume" in sys.argv)
    elif sys.argv[1] == "test":
        # 1 列ずつ、準位ごとに chi = 0.9 の列を alpha の上 8 点だけ解く(動作の確かめ)
        from multiprocessing import Pool
        with Pool(8) as p:
            for lev, c, re, im, bad, sec in p.map(column, [(lev, 0.9, alphas(lev)[::8]) for lev in LEVELS]):
                print(lev, f"{sec:.0f}s", "bad", int(bad.sum()), "NaN", int(np.sum(~np.isfinite(im))), "of", len(im))
