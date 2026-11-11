"""E18: chi_SR(表の版)を基準の運動論(fp3.evolve3m)と、入力の全範囲で照合する。

約束(研究ノート E18、2026-11-10 にユーザーが選んだ許容誤差 0.03): 印 0・1・2 の出力は、基準の ODE と
|chi_SR(T) - chi_ref(T')| <= 0.03 を満たす T' が |ln(T'/T)| <= 0.02 の中にある(速い段差の時刻は 2% までずれてよい)。
印 4・5 は誤差を報告するだけ、印 3 は予測ではない。

点の取り方(層):
  uniform  : M 2--200 Msun(対数一様)、chi0 0--0.995(一様)、mu 0.1--10e-12 eV(対数一様)、f 1e6--1e22 GeV(対数一様)
  thresh   : 上と同じだが chi0 を 211 か 322 の臨界の近く(臨界との差 1e-4--0.05、対数一様)に置く
  gen1     : 1 世代目(chi0 ~ Beta(0.69, 2.28))、M 8--50、mu 0.6--1.0e-12、f 1e10--1e19
  gw241011 : M 15--25、chi0 0.6--0.95、mu 0.6--1.0e-12、f 1e10--1e19
alpha > 0.5 と chi0 > 0.995 の点は引き直す。年齢は 1e5--1.4e10 yr の対数で 40 点。
基準の ODE は年齢の密な格子で解き、年齢の窓 [T e^-0.02, T e^0.02] の中の chi_ref の幅と比べる。
最後の年齢まで解けなかった点は、評価の回数 4 倍、rtol 1e-7、rtol 1e-6 の順に解き直す。緩い rtol の曲線は、それまでの
曲線と重なる年齢で 3e-3 以内に合うときだけ先の年齢に使う。窓全体が解けていないセルは照合できないので、数を報告する
(unchecked_unflagged_cells)。

実行: cd research && ~/miniforge3/bin/python3 calc/e18_validate2.py refs <種> <層ごとの点の数> [並列の数]
      cd research && ~/miniforge3/bin/python3 calc/e18_validate2.py retry <種> [並列の数]   (refs のあとに自動で走る)
      cd research && ~/miniforge3/bin/python3 calc/e18_validate2.py eval <種> [theta_v theta_k]
出力: out/e18_refs_<種>.npz(基準の ODE)、out/e18_validate2_<種>.json(まとめ)、out/e18_validate2_<種>.npz(すべての点)
"""
import json
import os
import signal
import sys
import time
import warnings
from multiprocessing import Pool

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import fp3  # noqa: E402
import reftab  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
TS = np.logspace(5, np.log10(1.4e10), 40)
DELTA = 0.02
TOL = 0.03
NG = 2400
TG = np.logspace(np.log10(1e5 * np.exp(-2 * DELTA)), np.log10(1.4e10 * np.exp(2 * DELTA)), NG)


class _Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise _Timeout()


def ref_point(args, max_nfev=1200000, limit=600, rtol=1e-8):
    """基準の ODE(年齢の格子 TG で、limit 秒で打ち切る)。返り値: chi(TG)。解けた所まで値を持ち、その先は NaN
    (evolve3m の出力は TG の先頭からの連続した部分)。"""
    M, c0, mu, f = args
    out = np.full(NG, np.nan)
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(limit)
    try:
        s = fp3.evolve3m(M, c0, mu, f, TG[-1], TG, rtol=rtol, max_nfev=max_nfev)
        signal.alarm(0)
    except _Timeout:
        return out
    except Exception:
        signal.alarm(0)
        return out
    y = np.asarray(s.y, float)
    if y.ndim == 2 and y.shape[1] > 0:
        n = min(y.shape[1], NG)
        out[:n] = y[4, :n] / y[3, :n] ** 2
    return out


# 解き直しの段: (名前, rtol, 評価の回数, 秒)。緩い rtol の曲線は、それまでの曲線と重なる年齢で TAIL_AGREE 以内に
# 合うときだけ、解けていなかった先の年齢に使う(rtol 1e-6 は急な段差で 1e-2 ずれることがあるので、重なりで確かめる)
RETRY_STAGES = (("budget4x", 1e-8, 4800000, 5400), ("rtol1e-7", 1e-7, 4800000, 3000), ("rtol1e-6", 1e-6, 4800000, 3000))
TAIL_AGREE = 3e-3


def _retry_point(args):
    pt, rtol, nfev, limit = args
    return ref_point(pt, max_nfev=nfev, limit=limit, rtol=rtol)


def _merge(old, new, same_rtol):
    """old の解けた所はそのまま、その先を new で埋める。返り値: (曲線, 埋めたか)。"""
    fo, fn = np.isfinite(old), np.isfinite(new)
    no = NG if fo.all() else int(np.argmin(fo))
    nn = NG if fn.all() else int(np.argmin(fn))
    if nn <= no:
        return old, False
    if same_rtol:                                       # 同じ rtol で長く解けた: そのまま置き換える
        return new, True
    if no > 0 and np.max(np.abs(new[:no] - old[:no])) > TAIL_AGREE:
        return old, False
    out = old.copy()
    out[no:nn] = new[no:nn]
    return out, True
def sample(kind, n, rng, tab):
    out = []
    while len(out) < n:
        m = 4 * n
        if kind in ("uniform", "thresh"):
            M = np.exp(rng.uniform(np.log(2), np.log(200), m))
            mu = np.exp(rng.uniform(np.log(1e-13), np.log(1e-11), m))
            f = 10 ** rng.uniform(6, 22, m)
            a = sr.alpha_of(M, mu)
            if kind == "uniform":
                c0 = rng.uniform(0, chisr.CHI_TOP, m)
            else:
                which = rng.integers(1, 3, m)
                cc = np.where(which == 1, chisr.chic_vec(a, 1, tab=tab), chisr.chic_vec(a, 2, tab=tab))
                c0 = cc + rng.choice([-1, 1], m) * np.exp(rng.uniform(np.log(1e-4), np.log(0.05), m))
        elif kind == "gen1":
            M = rng.uniform(8, 50, m)
            mu = rng.uniform(0.6e-12, 1.0e-12, m)
            f = 10 ** rng.uniform(10, 19, m)
            c0 = rng.beta(0.69, 2.28, m)
        else:
            M = rng.uniform(15, 25, m)
            mu = rng.uniform(0.6e-12, 1.0e-12, m)
            f = 10 ** rng.uniform(10, 19, m)
            c0 = rng.uniform(0.6, 0.95, m)
        a = sr.alpha_of(M, mu)
        keep = (a <= chisr.ALPHA_MAX) & (c0 >= 0) & (c0 <= chisr.CHI_TOP)
        out += list(zip(M[keep], c0[keep], mu[keep], f[keep]))
    return out[:n]


KINDS = ["uniform", "thresh", "gen1", "gw241011"]


def make_refs(seed, n, nproc):
    tab = chisr.build(use_reftab=False)
    rng = np.random.default_rng(seed)
    pts, lab = [], []
    for k in KINDS:
        s = sample(k, n, rng, tab)
        pts += s
        lab += [k] * len(s)
    pts = np.array(pts)
    t0 = time.time()
    with Pool(nproc) as p:
        refs = p.map(ref_point, [tuple(x) for x in pts], chunksize=2)
    refs = np.array(refs, np.float32)
    np.savez(os.path.join(OUT, f"e18_refs_{seed}.npz"), M=pts[:, 0], chi0=pts[:, 1], mu=pts[:, 2], f=pts[:, 3],
             kind=np.array(lab), TG=TG, chi=refs, wall=time.time() - t0)
    print(f"{len(pts)} references in {time.time() - t0:.0f}s, incomplete {int(np.sum(np.isnan(refs[:, -1])))}")
    retry_refs(seed, nproc)


def retry_refs(seed, nproc=10):
    """最後の年齢まで解けなかった点を RETRY_STAGES の順に解き直す(まだ走らせていない段だけ)。
    どの段が先の年齢を埋めたかを ref_stage に残す(-1 = 最初の計算のまま)。"""
    path = os.path.join(OUT, f"e18_refs_{seed}.npz")
    z = dict(np.load(path))
    refs = z["chi"].astype(float)
    stage = z.get("ref_stage", np.full(len(refs), -1, int))
    done_st = list(z.get("stages_done", np.array([], dtype="U16")))
    if "retried" in z and "budget4x" not in done_st:      # 段の記録を入れる前に、評価の回数 4 倍の段は済んでいる
        done_st.append("budget4x")
        stage[z["retried"]] = 0
    for k, (name, rtol, nfev, limit) in enumerate(RETRY_STAGES):
        if name in done_st:
            continue
        bad = np.flatnonzero(np.isnan(refs[:, -1]))
        if len(bad) == 0:
            break
        t0 = time.time()
        jobs = [((z["M"][i], z["chi0"][i], z["mu"][i], z["f"][i]), rtol, nfev, limit) for i in bad]
        with Pool(min(nproc, len(bad))) as p:
            new = p.map(_retry_point, jobs, chunksize=1)
        for i, r in zip(bad, new):
            refs[i], used = _merge(refs[i], r, same_rtol=(rtol == 1e-8))
            if used:
                stage[i] = k
        done_st.append(name)
        nfin = np.isfinite(refs[bad]).sum(1)
        print(f"retry {name}: {len(bad)} points in {time.time() - t0:.0f}s; still incomplete "
              f"{int(np.sum(np.isnan(refs[:, -1])))}; solved ages: {nfin.tolist()} of {NG}", flush=True)
    z["chi"] = refs.astype(np.float32)
    z["ref_stage"] = stage
    z["stages_done"] = np.array(done_st, dtype="U16")
    z.pop("retried", None)
    np.savez(path, **z)


def evaluate(seed, theta=None, tab=None, quiet=False):
    z = np.load(os.path.join(OUT, f"e18_refs_{seed}.npz"))
    M, c0, mu, f, lab, refs = z["M"], z["chi0"], z["mu"], z["f"], z["kind"], z["chi"].astype(float)
    tab = tab or chisr.build()
    if theta is not None:
        orig = reftab.lookup_rt3
        reftab.lookup_rt3 = lambda *a, **k: orig(*a, theta_v=theta[0], theta_k=theta[1], **k)
    N, nT = len(M), len(TS)
    csr = np.empty((N, nT))
    reg = np.empty((N, nT), int)
    for j, T in enumerate(TS):
        csr[:, j], reg[:, j] = chisr.chi_sr(tab, M, c0, T, mu, f)
    if theta is not None:
        reftab.lookup_rt3 = orig
    ok = np.isfinite(refs[:, -1])                       # 最後の年齢まで解けた点
    lg = np.log10(TG)
    lq = np.log10(TS)
    win = np.linspace(-DELTA, DELTA, 41) / np.log(10)
    cref = np.full((N, nT), np.nan)
    tol = np.full((N, nT), np.nan)
    cell_ok = np.zeros((N, nT), bool)                   # 年齢の窓全体で基準の ODE が解けているセル
    for i in range(N):
        fin = np.isfinite(refs[i])
        n = NG if fin.all() else int(np.argmin(fin))    # 解けた所(先頭からの連続した部分)
        if n < 2:
            continue
        cell_ok[i] = lq + win[-1] <= lg[n - 1]
        if not cell_ok[i].any():
            continue
        r = refs[i, :n]
        cref[i] = np.where(cell_ok[i], np.interp(lq, lg[:n], r), np.nan)
        cw = np.interp(lq[:, None] + win[None, :], lg[:n], r)
        tol[i] = np.where(cell_ok[i], np.maximum(0, np.maximum(csr[i] - cw.max(1), cw.min(1) - csr[i])), np.nan)
    plain = np.abs(csr - cref)

    def summ(mask):
        x, y = tol[mask], plain[mask]
        g = np.isfinite(x)
        x, y = x[g], y[g]
        if len(x) == 0:
            return None
        return {"n": int(len(x)), "tol_max": float(x.max()), "tol_p99": float(np.percentile(x, 99)),
                "tol_frac_gt": float(np.mean(x > TOL)), "plain_median": float(np.median(y)),
                "plain_p99": float(np.percentile(y, 99)), "plain_max": float(y.max())}

    unfl = reg <= 2
    a = sr.alpha_of(M, mu)
    res = {"seed": seed, "n_points": int(N), "theta": theta, "ode_incomplete_points": int((~ok).sum()),
           "unchecked_cells": int(np.sum(~cell_ok)),
           "unchecked_unflagged_cells": int(np.sum(~cell_ok & unfl)),
           "contract": f"regime 0-2: |chi_SR(T) - chi_ref(T')| <= {TOL} for some |ln T'/T| <= {DELTA}",
           "regime_share": {str(r): float(np.mean(reg == r)) for r in range(6)}}
    res["by_regime"] = {str(r): summ(reg == r) for r in range(6)}
    res["unflagged"] = summ(unfl)
    res["by_kind_unflagged"] = {k: summ(unfl & (lab == k)[:, None]) for k in KINDS}
    res["regime_share_by_kind"] = {k: {str(r): float(np.mean(reg[lab == k] == r)) for r in range(6)} for k in KINDS}
    ab = a[:, None] * np.ones((1, nT))
    res["by_alpha_unflagged"] = {f"{lo}-{hi}": summ(unfl & (ab >= lo) & (ab < hi))
                                 for lo, hi in [(0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5)]}
    fb = np.log10(f)[:, None] * np.ones((1, nT))
    res["by_f_unflagged"] = {f"1e{lo}-1e{lo + 2}": summ(unfl & (fb >= lo) & (fb < lo + 2)) for lo in range(6, 22, 2)}
    bad = np.argwhere(unfl & (tol > TOL))
    res["n_violations"] = int(len(bad))
    order = np.argsort(-tol[bad[:, 0], bad[:, 1]]) if len(bad) else []
    res["worst"] = [{"kind": str(lab[i]), "M": float(M[i]), "chi0": float(c0[i]), "mu": float(mu[i]), "f": float(f[i]),
                     "alpha": float(a[i]), "T": float(TS[j]), "chisr": float(csr[i, j]), "ref": float(cref[i, j]),
                     "err_tol": float(tol[i, j]), "regime": int(reg[i, j])} for i, j in (bad[order][:40] if len(bad) else [])]
    unc = np.argwhere(~cell_ok & unfl)
    res["unchecked_unflagged"] = [{"kind": str(lab[i]), "M": float(M[i]), "chi0": float(c0[i]), "mu": float(mu[i]),
                                   "f": float(f[i]), "alpha": float(a[i]), "T": float(TS[j]), "chisr": float(csr[i, j]),
                                   "regime": int(reg[i, j])} for i, j in unc[:60]]
    tag = "" if theta is None else f"_v{theta[0]}_k{theta[1]}"
    with open(os.path.join(OUT, f"e18_validate2_{seed}{tag}.json"), "w") as fo:
        json.dump(res, fo, indent=1)
    np.savez(os.path.join(OUT, f"e18_validate2_{seed}{tag}.npz"), M=M, chi0=c0, mu=mu, f=f, kind=lab, T=TS, chisr=csr,
             regime=reg, ref=cref, err_plain=plain, err_tol=tol)
    if not quiet:
        print(json.dumps({k: v for k, v in res.items() if k not in ("worst", "by_alpha_unflagged", "by_f_unflagged")},
                         indent=1))
        for w in res["worst"][:12]:
            print(w)
    return res


if __name__ == "__main__":
    if sys.argv[1] == "refs":
        make_refs(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 10)
    elif sys.argv[1] == "retry":
        retry_refs(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 10)
    else:
        th = (float(sys.argv[3]), float(sys.argv[4])) if len(sys.argv) > 4 else None
        evaluate(int(sys.argv[2]), th)
