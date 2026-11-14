"""E19: 査読のあとの chi_SR(単体版 share/chi_sr)を、基準の運動論(fp3.evolve3m)と新しい点で照合する。

約束と、基準の曲線の作り方は E18(calc/e18_validate2.py)と同じ。違うのは次の点。
  - 新しい種。E18 の 4 つの層(uniform、thresh、gen1、gw241011)に、査読で問題の見つかった所の層を加える。
      switch211 : alpha が、211 の率の切り替え(0.03)を挟む値の表のセルと、その両側(節の 3% 外まで)
      switch322 : 同じく 322 の率の切り替え(0.1)
      lowalpha  : alpha 0.004--0.012(半解析の固定点の表の下端 0.01 と、値の表の下端 0.008 の両側)
    加えた層は M 2--200 Msun(対数一様)、chi0 0--0.995(一様)、f は半分を 1e6--1e22 GeV、半分を 1e16--10^20.5 GeV
    (対数一様。切り替えのずれは自由場の側に出る)。mu が 1e-13--1e-11 eV の外に出る点は引き直す。
  - 年齢は 1e5--1.4e10 yr の対数で 200 点(E18 は 40 点)。段差の時刻のずれは狭い年齢の幅にしか出ないので、細かくとる。
  - 評価するのは単体版。プロジェクト版(chisr.chi_sr)が同じ値・同じ印を返すことも、すべてのセルで確かめる。
  - 修正の前の単体版(git の OLD_COMMIT)も同じ点で評価し、前後を比べる。
最後の年齢まで解けなかった点は、E18 と同じ段(評価の回数 4 倍、rtol 1e-7、1e-6)で解き直す(retry)。

細かい年齢での評価(dense): 基準の曲線は年齢の 2400 点の格子 TG で持っているので、その格子の年齢のうち約束の範囲
(1e5--1.4e10 yr)に入る 2382 点すべてで単体版を呼び、年齢の窓は格子の前後 4 点(|ln(T'/T)| <= 0.0199)でとる。新しい ODE は要らない。段差の時刻の
ずれは狭い年齢の幅にしか出ないので、まれな違反の頻度はこちらで数える。

実行: cd research && ~/miniforge3/bin/python3 calc/e19_validate.py refs <種> <並列の数> [E18 の層の点の数 加えた層の点の数]
      cd research && ~/miniforge3/bin/python3 calc/e19_validate.py retry <種> <並列の数> [1 点あたりの秒の上限]
      cd research && ~/miniforge3/bin/python3 calc/e19_validate.py eval <種>
      cd research && ~/miniforge3/bin/python3 calc/e19_validate.py dense <種>
出力: out/e19_refs_<種>.npz(基準の曲線、git に入れない)、out/e19_validate_<種>.json(まとめ)、
      out/e19_validate_<種>.npz(点ごと・セルごとの値)、out/e19_dense_<種>.json(細かい年齢でのまとめ)
"""
import importlib.util
import json
import os
import subprocess
import sys
import time
import warnings
from multiprocessing import Pool

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import e18_validate2 as v18  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
PKG = os.path.join(HERE, "..", "..", "share", "chi_sr")
OLD_COMMIT = "7846eae"                                   # 査読した版(E18)
TS = np.logspace(5, np.log10(1.4e10), 200)
DELTA, TOL, NG, TG = v18.DELTA, v18.TOL, v18.NG, v18.TG
KINDS = v18.KINDS + ["switch211", "switch322", "lowalpha"]
MU_LO, MU_HI = 1e-13, 1e-11


def load_pkg(old=False):
    """単体版の ChiSR。old=True なら修正の前の版(git の OLD_COMMIT の chi_sr.py。表は同じ)。"""
    if old:
        src = subprocess.run(["git", "show", f"{OLD_COMMIT}:share/chi_sr/chi_sr.py"], cwd=HERE, capture_output=True,
                             text=True, check=True).stdout
        path = os.path.join(OUT, "_chi_sr_old.py")
        with open(path, "w") as fo:
            fo.write(src)
        spec = importlib.util.spec_from_file_location("chi_sr_old", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        os.remove(path)
    else:
        sys.path.insert(0, PKG)
        import chi_sr as mod
    return mod.ChiSR(os.path.join(PKG, "chi_sr_tables.npz"))


def sample_extra(kind, n, rng, a_nodes):
    out = []
    while len(out) < n:
        m = 4 * n
        if kind == "lowalpha":
            a = np.exp(rng.uniform(np.log(0.004), np.log(0.012), m))
        else:
            edge = 0.03 if kind == "switch211" else 0.1
            i = np.searchsorted(a_nodes, edge) - 1
            a = np.exp(rng.uniform(np.log(a_nodes[i] * 0.97), np.log(a_nodes[i + 1] * 1.03), m))
        M = np.exp(rng.uniform(np.log(2), np.log(200), m))
        mu = a / sr.alpha_of(M, 1.0)
        c0 = rng.uniform(0, chisr.CHI_TOP, m)
        f = np.where(rng.random(m) < 0.5, 10 ** rng.uniform(6, 22, m), 10 ** rng.uniform(16, 20.5, m))
        keep = (mu >= MU_LO) & (mu <= MU_HI)
        out += list(zip(M[keep], c0[keep], mu[keep], f[keep]))
    return out[:n]


def make_refs(seed, nproc, n18=750, nx=600):
    tab = chisr.build(use_reftab=False)
    a_nodes = chisr.build()["rt"]["A_NODES"]
    rng = np.random.default_rng(seed)
    pts, lab = [], []
    for k in KINDS:
        s = v18.sample(k, n18, rng, tab) if k in v18.KINDS else sample_extra(k, nx, rng, a_nodes)
        pts += s
        lab += [k] * len(s)
    pts = np.array(pts)
    t0 = time.time()
    with Pool(nproc) as p:
        refs = p.map(v18.ref_point, [tuple(x) for x in pts], chunksize=2)
    refs = np.array(refs, np.float32)
    np.savez(os.path.join(OUT, f"e19_refs_{seed}.npz"), M=pts[:, 0], chi0=pts[:, 1], mu=pts[:, 2], f=pts[:, 3],
             kind=np.array(lab), TG=TG, chi=refs, wall=time.time() - t0, ref_stage=np.full(len(pts), -1, int),
             stages_done=np.array([], dtype="U16"))
    print(f"{len(pts)} references in {time.time() - t0:.0f}s, incomplete {int(np.sum(np.isnan(refs[:, -1])))}", flush=True)


def retry_refs(seed, nproc, limit=None):
    """E18 と同じ段で解き直す(v18.RETRY_STAGES)。limit を渡すと、1 点あたりの秒の上限をそれで置き換える。"""
    path = os.path.join(OUT, f"e19_refs_{seed}.npz")
    z = dict(np.load(path))
    refs = z["chi"].astype(float)
    stage = z["ref_stage"]
    done_st = list(z["stages_done"])
    for k, (name, rtol, nfev, lim) in enumerate(v18.RETRY_STAGES):
        if name in done_st:
            continue
        bad = np.flatnonzero(np.isnan(refs[:, -1]))
        if len(bad) == 0:
            break
        t0 = time.time()
        jobs = [((z["M"][i], z["chi0"][i], z["mu"][i], z["f"][i]), rtol, nfev, int(limit or lim)) for i in bad]
        with Pool(min(nproc, len(bad))) as p:
            new = p.map(v18._retry_point, jobs, chunksize=1)
        for i, r in zip(bad, new):
            refs[i], used = v18._merge(refs[i], r, same_rtol=(rtol == 1e-8))
            if used:
                stage[i] = k
        done_st.append(name)
        print(f"retry {name}: {len(bad)} points in {time.time() - t0:.0f}s; still incomplete "
              f"{int(np.sum(np.isnan(refs[:, -1])))}", flush=True)
    z["chi"] = refs.astype(np.float32)
    z["ref_stage"] = stage
    z["stages_done"] = np.array(done_st, dtype="U16")
    np.savez(path, **z)


def tol_error(csr, refs):
    """年齢の 2% のずれを許した誤差。返り値: (tol, plain, cref, cell_ok)。E18 の evaluate と同じ計算。"""
    N, nT = csr.shape
    lg, lq = np.log10(TG), np.log10(TS)
    win = np.linspace(-DELTA, DELTA, 41) / np.log(10)
    cref = np.full((N, nT), np.nan)
    tol = np.full((N, nT), np.nan)
    cell_ok = np.zeros((N, nT), bool)
    for i in range(N):
        fin = np.isfinite(refs[i])
        n = NG if fin.all() else int(np.argmin(fin))
        if n < 2:
            continue
        cell_ok[i] = lq + win[-1] <= lg[n - 1]
        if not cell_ok[i].any():
            continue
        r = refs[i, :n]
        cref[i] = np.where(cell_ok[i], np.interp(lq, lg[:n], r), np.nan)
        cw = np.interp(lq[:, None] + win[None, :], lg[:n], r)
        tol[i] = np.where(cell_ok[i], np.maximum(0, np.maximum(csr[i] - cw.max(1), cw.min(1) - csr[i])), np.nan)
    return tol, np.abs(csr - cref), cref, cell_ok


def summ(tol, plain, mask):
    x, y = tol[mask], plain[mask]
    g = np.isfinite(x)
    x, y = x[g], y[g]
    if len(x) == 0:
        return None
    npt = int(np.sum(np.any(mask & np.isfinite(tol), axis=1)))
    vpt = int(np.sum(np.any(mask & (tol > TOL), axis=1)))
    return {"cells": int(len(x)), "points": npt, "tol_max": float(x.max()), "tol_p99": float(np.percentile(x, 99)),
            "violations": int(np.sum(x > TOL)), "points_with_violation": vpt, "plain_max": float(y.max())}


def evaluate(seed):
    z = np.load(os.path.join(OUT, f"e19_refs_{seed}.npz"))
    M, c0, mu, f, lab, refs = z["M"], z["chi0"], z["mu"], z["f"], z["kind"], z["chi"].astype(float)
    N, nT = len(M), len(TS)
    args = (M[:, None], c0[:, None], TS[None, :], mu[:, None], f[:, None])
    new, old = load_pkg(), load_pkg(old=True)
    tab = chisr.build()

    def run(fn):
        c = np.empty((N, nT))
        r = np.empty((N, nT), int)
        for i0 in range(0, N, 500):                     # メモリを抑えるために点を分けて呼ぶ(値は分け方によらない)
            sl = slice(i0, i0 + 500)
            c[sl], r[sl] = fn(*[np.broadcast_to(x[sl] if x.shape[0] == N else x, (len(M[sl]), nT)) for x in args])
        return c, r

    c_new, r_new = run(new)
    c_old, r_old = run(old)
    c_prj, r_prj = run(lambda *a: chisr.chi_sr(tab, *a))
    a = sr.alpha_of(M, mu)
    an = tab["rt"]["A_NODES"]
    res = {"seed": seed, "n_points": int(N), "n_ages": int(nT), "old_commit": OLD_COMMIT, "value_table": tab["rt_source"],
           "contract": f"regime 0-2: |chi_SR(T) - chi_ref(T')| <= {TOL} for some |ln T'/T| <= {DELTA}",
           "ode_incomplete_points": int(np.sum(~np.isfinite(refs[:, -1]))),
           "project_vs_standalone": {"max_diff": float(np.max(np.abs(c_new - c_prj))),
                                     "regime_mismatch": int(np.sum(r_new != r_prj))}}
    assert res["project_vs_standalone"]["max_diff"] == 0 and res["project_vs_standalone"]["regime_mismatch"] == 0, res

    ab = a[:, None] * np.ones((1, nT))
    fb = np.log10(f)[:, None] * np.ones((1, nT))
    sw = chisr.switch_cell(an, a)[:, None] & np.ones((1, nT), bool)
    for tag, c, r in (("new", c_new, r_new), ("old", c_old, r_old)):
        tol, plain, cref, cell_ok = tol_error(c, refs)
        unfl = r <= 2
        d = {"regime_share": {str(k): float(np.mean(r == k)) for k in range(6)},
             "unchecked_cells": int(np.sum(~cell_ok)), "unchecked_unflagged_cells": int(np.sum(~cell_ok & unfl)),
             "unflagged": summ(tol, plain, unfl),
             "by_regime": {str(k): summ(tol, plain, r == k) for k in range(6)},
             "by_kind_unflagged": {k: summ(tol, plain, unfl & (lab == k)[:, None]) for k in KINDS},
             "regime_share_by_kind": {k: {str(j): float(np.mean(r[lab == k] == j)) for j in range(6)} for k in KINDS},
             "by_alpha_unflagged": {f"{lo}-{hi}": summ(tol, plain, unfl & (ab >= lo) & (ab < hi)) for lo, hi in
                                    [(0, 0.008), (0.008, 0.01), (0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3),
                                     (0.3, 0.4), (0.4, 0.5)]},
             "by_f_unflagged": {f"1e{lo}-1e{lo + 2}": summ(tol, plain, unfl & (fb >= lo) & (fb < lo + 2))
                                for lo in range(6, 22, 2)},
             "switch_cells": {"unflagged": summ(tol, plain, unfl & sw),
                              "regime_share": {str(k): float(np.mean(r[sw] == k)) for k in range(6)},
                              "regime2": summ(tol, plain, (r == 2) & sw), "regime1": summ(tol, plain, (r == 1) & sw),
                              "regime5": summ(tol, plain, (r == 5) & sw)},
             "outside_switch_cells_unflagged": summ(tol, plain, unfl & ~sw),
             "frac_gt_tol_by_regime": {str(k): (float(np.mean(tol[(r == k) & np.isfinite(tol)] > TOL))
                                                 if np.any((r == k) & np.isfinite(tol)) else None) for k in (4, 5)}}
        # 印 3 の点で、基準の運動論が動いたか(模型の外としたのが、何も起きない点かどうか)
        m3 = (r == 3) & np.isfinite(cref)
        mv = np.abs(cref - c0[:, None])
        d["regime3"] = {"cells": int(m3.sum()), "ref_moved_gt_0.01": float(np.mean(mv[m3] > 0.01)) if m3.any() else None,
                        "by_alpha": {f"{lo}-{hi}": {"cells": int((m3 & (ab >= lo) & (ab < hi)).sum()),
                                                    "ref_moved_max": float(mv[m3 & (ab >= lo) & (ab < hi)].max())
                                                    if (m3 & (ab >= lo) & (ab < hi)).any() else None}
                                     for lo, hi in [(0, 0.008), (0.008, 0.01), (0.01, 0.3), (0.3, 0.5)]}}
        bad = np.argwhere(unfl & (tol > TOL))
        order = np.argsort(-tol[bad[:, 0], bad[:, 1]]) if len(bad) else []
        d["worst"] = [{"kind": str(lab[i]), "M": float(M[i]), "chi0": float(c0[i]), "mu": float(mu[i]), "f": float(f[i]),
                       "alpha": float(a[i]), "T": float(TS[j]), "chisr": float(c[i, j]), "ref": float(cref[i, j]),
                       "err_tol": float(tol[i, j]), "regime": int(r[i, j]), "in_switch_cell": bool(sw[i, j])}
                      for i, j in (bad[order][:40] if len(bad) else [])]
        res[tag] = d
        if tag == "new":
            keep = dict(tol=tol.astype(np.float32), ref=cref.astype(np.float32))
    ch = (r_new != r_old) | (c_new != c_old)
    res["changed_cells"] = {"total": int(ch.sum()),
                            "by_regime_old_new": {f"{i}->{j}": int(np.sum(ch & (r_old == i) & (r_new == j)))
                                                  for i in range(6) for j in range(6) if np.any(ch & (r_old == i) & (r_new == j))},
                            "value_changed": int(np.sum(c_new != c_old))}
    with open(os.path.join(OUT, f"e19_validate_{seed}.json"), "w") as fo:
        json.dump(res, fo, indent=1)
    np.savez_compressed(os.path.join(OUT, f"e19_validate_{seed}.npz"), M=M, chi0=c0, mu=mu, f=f, kind=lab, T=TS,
                        chisr=c_new.astype(np.float32), regime=r_new.astype(np.int8), regime_old=r_old.astype(np.int8),
                        ref=keep["ref"], err_tol=keep["tol"])
    show = {k: v for k, v in res.items() if k not in ("new", "old")}
    for tag in ("new", "old"):
        show[tag] = {k: res[tag][k] for k in ("regime_share", "unchecked_unflagged_cells", "unflagged", "by_kind_unflagged",
                                              "switch_cells", "outside_switch_cells_unflagged", "frac_gt_tol_by_regime",
                                              "regime3")}
    print(json.dumps(show, indent=1))
    for w in res["new"]["worst"][:15]:
        print(w)
    return res


def evaluate_dense(seed, half=4):
    """基準の曲線の格子の年齢すべてで照合する(年齢の窓は前後 half 点)。修正の前後の単体版を比べる。"""
    from numpy.lib.stride_tricks import sliding_window_view
    z = np.load(os.path.join(OUT, f"e19_refs_{seed}.npz"))
    M, c0, mu, f, lab, refs = z["M"], z["chi0"], z["mu"], z["f"], z["kind"], z["chi"].astype(float)
    N = len(M)
    Td = TG[half:NG - half]
    sel = (Td >= TS[0] * (1 - 1e-12)) & (Td <= TS[-1] * (1 + 1e-12))                # 約束の年齢の範囲(1e5--1.4e10 yr)の中だけ
    Td = Td[sel]
    nT = len(Td)
    lo = np.nanmin(sliding_window_view(refs, 2 * half + 1, axis=1), axis=2)[:, sel]      # 窓の中に NaN があるセルは下で除く
    hi = np.nanmax(sliding_window_view(refs, 2 * half + 1, axis=1), axis=2)[:, sel]
    ok = np.all(np.isfinite(sliding_window_view(refs, 2 * half + 1, axis=1)), axis=2)[:, sel]
    new, old = load_pkg(), load_pkg(old=True)
    a = sr.alpha_of(M, mu)
    sw = chisr.switch_cell(chisr.build()["rt"]["A_NODES"], a)
    res = {"seed": seed, "n_points": int(N), "n_ages": int(nT), "window_half_width_lnT": float(np.log(TG[half] / TG[0])),
           "old_commit": OLD_COMMIT, "ode_incomplete_points": int(np.sum(~np.isfinite(refs[:, -1])))}
    for tag, model in (("new", new), ("old", old)):
        tol = np.full((N, nT), np.nan, np.float32)
        reg = np.empty((N, nT), np.int8)
        for i0 in range(0, N, 100):
            sl = slice(i0, min(i0 + 100, N))
            c, r = model(M[sl, None], c0[sl, None], Td[None, :], mu[sl, None], f[sl, None])
            e = np.maximum(0, np.maximum(c - hi[sl], lo[sl] - c))
            tol[sl] = np.where(ok[sl], e, np.nan)
            reg[sl] = r
        unfl = reg <= 2
        plain = np.zeros_like(tol)                      # summ に渡す形だけ(細かい評価では使わない)
        d = {"regime_share": {str(k): float(np.mean(reg == k)) for k in range(6)},
             "unchecked_unflagged_cells": int(np.sum(unfl & ~np.isfinite(tol))),
             "unflagged": summ(tol, plain, unfl), "by_kind_unflagged": {k: summ(tol, plain, unfl & (lab == k)[:, None]) for k in KINDS},
             "switch_cells_unflagged": summ(tol, plain, unfl & sw[:, None]),
             "outside_switch_cells_unflagged": summ(tol, plain, unfl & ~sw[:, None]),
             "by_regime": {str(k): summ(tol, plain, reg == k) for k in range(6)},
             "frac_gt_tol_by_regime": {str(k): (float(np.mean(tol[(reg == k) & np.isfinite(tol)] > TOL))
                                                 if np.any((reg == k) & np.isfinite(tol)) else None) for k in (4, 5)}}
        for v in [d["unflagged"], d["switch_cells_unflagged"], d["outside_switch_cells_unflagged"]] + list(d["by_kind_unflagged"].values()) \
                + list(d["by_regime"].values()):
            if v:
                v.pop("plain_max")
        # 約束を越えた点ごとに、越えた年齢の範囲と最大の誤差
        viol = []
        for i in np.flatnonzero(np.any(unfl & (tol > TOL), axis=1)):
            b = unfl[i] & (tol[i] > TOL)
            j = int(np.nanargmax(np.where(b, tol[i], -1)))
            viol.append({"kind": str(lab[i]), "M": float(M[i]), "chi0": float(c0[i]), "mu": float(mu[i]), "f": float(f[i]),
                         "alpha": float(a[i]), "in_switch_cell": bool(sw[i]), "n_ages": int(b.sum()),
                         "T_range": [float(Td[b].min()), float(Td[b].max())], "err_tol_max": float(tol[i, j]),
                         "regime_there": sorted(set(reg[i][b].tolist()))})
        d["violating_points"] = sorted(viol, key=lambda v: -v["err_tol_max"])
        # 誤差の分布の裾(印 0・1・2)
        x = tol[unfl & np.isfinite(tol)]
        d["unflagged_tail"] = {f"gt_{t}": int(np.sum(x > t)) for t in (0.01, 0.015, 0.02, 0.025, 0.03, 0.04)}
        res[tag] = d
    with open(os.path.join(OUT, f"e19_dense_{seed}.json"), "w") as fo:
        json.dump(res, fo, indent=1)
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    cmd, seed = sys.argv[1], int(sys.argv[2])
    if cmd == "refs":
        extra = [int(x) for x in sys.argv[4:6]]
        make_refs(seed, int(sys.argv[3]), *extra)
    elif cmd == "retry":
        retry_refs(seed, int(sys.argv[3]), float(sys.argv[4]) if len(sys.argv) > 4 else None)
    elif cmd == "dense":
        evaluate_dense(seed)
    else:
        evaluate(seed)
