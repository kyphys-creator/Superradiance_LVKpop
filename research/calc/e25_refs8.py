"""E19 の照合の点の基準の曲線を、新しい基準の運動論(8 準位、kinn.evolveN_checked)で解き直し、新しい表(reftab8)と比べる
(DESIGN.md D18 の手順 6 の「照合を解き直す」。2026-11-21、ユーザー「1 と 2 を同時に」)。

今回解くのは、311・522 の確かめ(calc/e24_outside.py)で影響が見つからなかった所の点だけ:
  - 初期スピンが 211 の臨界以下(表の区分 0・-1・-2)
  - 544 の臨界より下(8 準位のどれも超放射しない)
のうち、表の alpha(0.008-0.5)と f(1e8-10^20.5 GeV)の範囲に入る点。区分 1(211 の臨界より上から始まる)は、522 が育つ範囲を
確かめてから解く。点は E19 の 3 つの種(out/e19_validate_<種>.npz の M・chi0・mu・f・層)をそのまま使う。

基準の曲線: 表の節と同じ系(reftab8._system: 8 準位・74 過程・相対論的な補正・多倍長の率)と同じ解き方
(kinn.evolveN_checked、上限は環境変数 SRLVK_CAP、既定 1e7 = 本番の表と同じ)で、年齢の格子 TG(E18 と同じ、2400 点、
yr)で解く。有効な年齢の数 n_valid(平均の確かめが通った、初めから続く年齢)と、届いた年齢までの chi を持つ。
例外を出した点は status 2 として記録して先へ進む。

比べ方(eval): 表から lookup8 で引いた値(544 の臨界より下の点は chi0 のまま)と基準の曲線を、約束と同じ物差し(時刻のずれ 2% を許した |Delta chi|、
kinn._tol_dist)で、約束の年齢の範囲(1e5-1.4e10 yr)のうち基準が有効で表も有効(lookup8 の valid)な年齢すべてで測る。
印(chi_sr の regime)の規則はまだ決め直していないので、ここで確かめるのは表と引き方の誤差だけ。

実行(research/ で):
  ~/miniforge3/bin/python3 calc/e25_refs8.py plan <場所>                点の一覧 <場所>/refs8_points.npz
  python e25_refs8.py task <番号> <本数> <場所>                         Hoffman2 のジョブ配列の 1 本(SLURM_CPUS_PER_TASK 並列)
  ~/miniforge3/bin/python3 calc/e25_refs8.py merge <場所> <出力の npz>   各本の結果を一つに
  ~/miniforge3/bin/python3 calc/e25_refs8.py eval <基準の npz> <表の npz> [出力の json]
Hoffman2 では、このファイルを checkout の外に置いて動かすときは環境変数 SRLVK_CALC に checkout の research/calc を入れる。
"""
import glob
import json
import os
import sys
import time
import warnings

import numpy as np

CALC = os.environ.get("SRLVK_CALC", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, CALC)
import e18_validate2 as v18  # noqa: E402
import reftab8 as R8  # noqa: E402
import sr  # noqa: E402

SEEDS = (1901, 1902, 1903)
TG, NG = v18.TG, v18.NG
T_LO, T_HI = 1e5, 1.4e10                        # 約束の年齢の範囲
PLAN_SEED = 25


def classify(M, chi0, mu, f):
    """点ごとの表の区分(1, 0, -1, -2。544 の臨界より下は -9、alpha が表の外は -99)と、alpha・f が表の範囲に入るか。"""
    a = mu * sr.alpha_of(M, 1.0)
    inA = (a >= R8.A_NODES[0]) & (a <= R8.A_NODES[-1])
    lf = np.log10(f)
    inF = (lf >= R8.LF_NODES[0] - 1e-12) & (lf <= R8.LF_NODES[-1] + 1e-12)
    p = np.full(len(a), -99)
    if inA.any():
        cs = R8.crit(a[inA])
        p[inA] = R8.piece_of(R8.w_of_chi(chi0[inA], cs))
    return a, p, inA & inF


def plan(out):
    """今回解く点の一覧を作る(区分 0・-1・-2 と 544 の臨界より下、alpha と f が表の範囲)。並びを混ぜて重い点を本ごとに散らす。"""
    cols = {k: [] for k in ("seed", "idx", "M", "chi0", "mu", "f", "kind", "piece", "alpha")}
    for seed in SEEDS:
        z = np.load(os.path.join(CALC, "out", f"e19_validate_{seed}.npz"))
        a, p, inr = classify(z["M"], z["chi0"], z["mu"], z["f"])
        sel = np.flatnonzero(inr & np.isin(p, [0, -1, -2, -9]))
        for k, v in (("seed", np.full(len(sel), seed)), ("idx", sel), ("M", z["M"][sel]), ("chi0", z["chi0"][sel]),
                     ("mu", z["mu"][sel]), ("f", z["f"][sel]), ("kind", z["kind"][sel]), ("piece", p[sel]),
                     ("alpha", a[sel])):
            cols[k].append(v)
    cols = {k: np.concatenate(v) for k, v in cols.items()}
    order = np.random.default_rng(PLAN_SEED).permutation(len(cols["M"]))
    cols = {k: v[order] for k, v in cols.items()}
    cols["pid"] = np.arange(len(order))         # 点の番号(結果はこの番号で持つ)
    os.makedirs(out, exist_ok=True)
    np.savez(os.path.join(out, "refs8_points.npz"), **cols)
    print(f"{len(order)} points -> {os.path.join(out, 'refs8_points.npz')}: "
          + ", ".join(f"piece {q}: {int(np.sum(cols['piece'] == q))}" for q in (0, -1, -2, -9)))
    return cols


def ref_point(args):
    """一つの点: (pid, M, chi0, mu, f, cap)。返り値の chi は TG の年齢での値(届かない年齢は NaN)。"""
    warnings.simplefilter("ignore")
    import kinn
    pid, M, c0, mu, f, cap = args
    t0 = time.time()
    out = dict(pid=int(pid), chi=np.full(NG, np.nan, np.float32), n_valid=0, n_reach=0, status=2, reason="",
               nfev=[0, 0], njump=[0, 0], secs=0.0)
    try:
        b = kinn.evolveN_checked(R8._system(), M, c0, mu, f, TG[-1], TG, max_nfev=cap)
        chi = np.asarray(kinn.chi_of(b), float)[:NG]
        out["chi"][:len(chi)] = chi
        runs = list(b.runs) + [dict(nfev=0, njump=0)] * (2 - len(b.runs))
        out.update(n_valid=int(b.n_valid), n_reach=int(len(chi)), status=int(b.status), reason=str(b.reason),
                   nfev=[int(r["nfev"]) for r in runs[:2]], njump=[int(r["njump"]) for r in runs[:2]])
    except Exception as e:                       # noqa: BLE001(どの例外でも、その点だけ失敗にする)
        out["reason"] = ("error: " + repr(e))[:200]
    out["secs"] = time.time() - t0
    return out


def _pack(res, cap):
    import platform
    return dict(pid=np.array([r["pid"] for r in res], int), chi=np.array([r["chi"] for r in res], np.float32).reshape(-1, NG),
                n_valid=np.array([r["n_valid"] for r in res], int), n_reach=np.array([r["n_reach"] for r in res], int),
                status=np.array([r["status"] for r in res], np.int8), reason=np.array([r["reason"] for r in res]),
                nfev=np.array([r["nfev"] for r in res], np.int64).reshape(-1, 2),
                njump=np.array([r["njump"] for r in res], np.int64).reshape(-1, 2),
                secs=np.array([r["secs"] for r in res]), cap=cap, host=platform.node())


def _unpack(path):
    with np.load(path) as z:
        d = {k: z[k] for k in z.files}
    return {int(p): dict(pid=int(p), chi=d["chi"][i], n_valid=int(d["n_valid"][i]), n_reach=int(d["n_reach"][i]),
                         status=int(d["status"][i]), reason=str(d["reason"][i]), nfev=[int(x) for x in d["nfev"][i]],
                         njump=[int(x) for x in d["njump"][i]], secs=float(d["secs"][i]))
            for i, p in enumerate(d["pid"])}


def task(tid, ntask, out, cap, nproc=1, every=20, every_s=900.0):
    """ジョブ配列の 1 本: refs8_points.npz のうち番号 tid(1 から)の受け持ちを nproc 並列で解き、out/refs8_<tid>.npz に
    every 点か every_s 秒ごとに保存し直す(保存したものがあれば続きから)。終了の合図(SIGTERM)で保存してやめる。
    形は reftab8.build_task と同じ。"""
    import multiprocessing as mp
    import signal
    from multiprocessing import Pool
    z = np.load(os.path.join(out, "refs8_points.npz"))
    mine = np.arange(len(z["pid"]))[tid - 1::ntask]
    path = os.path.join(out, f"refs8_{tid:05d}.npz")
    have = _unpack(path) if os.path.exists(path) else {}
    todo = [i for i in mine if int(z["pid"][i]) not in have]
    print(f"task {tid}/{ntask}: {len(mine)} points, {len(have)} already saved, {len(todo)} to solve, {nproc} workers, "
          f"cap {cap}", flush=True)

    def save():
        tmp = path + ".tmp.npz"
        np.savez(tmp, **_pack(list(have.values()), cap))
        os.replace(tmp, path)

    term = []
    signal.signal(signal.SIGTERM, lambda signum, frame: term.append(signum))
    if todo:
        last, n0, ndone = time.time(), len(have), 0
        pool = Pool(nproc, initializer=signal.signal, initargs=(signal.SIGTERM, signal.SIG_DFL))
        jobs = [(int(z["pid"][i]), float(z["M"][i]), float(z["chi0"][i]), float(z["mu"][i]), float(z["f"][i]), cap)
                for i in todo]
        it = pool.imap_unordered(ref_point, jobs, chunksize=1)
        while ndone < len(todo) and not term:
            try:
                r = it.next(timeout=5.0)
            except mp.TimeoutError:
                continue
            ndone += 1
            have[r["pid"]] = r
            if (len(have) - n0) % every == 0 or time.time() - last > every_s:
                save()
                last = time.time()
        save()
        pool.terminate()
        pool.join()
        if term:
            print(f"task {tid}: SIGTERM, saved {len(have)}/{len(mine)} points", flush=True)
            return
    print(f"task {tid}: {len(have)}/{len(mine)} points saved", flush=True)


def merge(out, dest):
    """各本の結果を点の一覧の順に一つの npz に(点の情報つき)。欠けた点は status -9(-1 は解き手が途中で止まった点で、
    n_valid までは有効)。"""
    z = dict(np.load(os.path.join(out, "refs8_points.npz")))
    n = len(z["pid"])
    res = {}
    caps, hosts = set(), set()
    for f in sorted(glob.glob(os.path.join(out, "refs8_[0-9]*.npz"))):
        if f.endswith(".tmp.npz"):
            continue
        with np.load(f) as d:
            caps.add(int(d["cap"]))
            hosts.add(str(d["host"]))
        res.update(_unpack(f))
    chi = np.full((n, NG), np.nan, np.float32)
    cols = dict(n_valid=np.zeros(n, int), n_reach=np.zeros(n, int), status=np.full(n, -9, np.int8),
                reason=np.full(n, "missing", dtype="U200"), nfev=np.zeros((n, 2), np.int64), njump=np.zeros((n, 2), np.int64),
                secs=np.zeros(n))
    for i, p in enumerate(z["pid"]):
        r = res.get(int(p))
        if r is None:
            continue
        chi[i] = r["chi"]
        for k in cols:
            cols[k][i] = r[k]
    z.update(cols, chi=chi, TG=TG, caps=np.array(sorted(caps)), n_hosts=len(hosts), provenance=json.dumps(R8.provenance()))
    np.savez(dest, **z)
    st = cols["status"]
    print(f"{n} points: status 0 {int(np.sum(st == 0))}, 1 {int(np.sum(st == 1))}, 2 {int(np.sum(st == 2))}, "
          f"-1 (valid up to n_valid) {int(np.sum(st == -1))}, missing {int(np.sum(st == -9))}; caps {sorted(caps)}; core-hours {cols['secs'].sum() / 3600:.1f} -> {dest}")


def evaluate(refs_path, tab_path, dest=None, tol=v18.TOL, delta=v18.DELTA):
    """表から引いた値と基準の曲線を、約束の物差しで比べる。"""
    import kinn
    z = np.load(refs_path)
    tab = R8.load_table(tab_path)
    inT = (TG >= T_LO) & (TG <= T_HI)
    n = len(z["pid"])
    err = np.full((n, NG), np.nan, np.float32)
    nref = np.zeros(n, int)
    ncomp = np.zeros(n, int)
    for i in range(n):
        nv = int(z["n_valid"][i])
        if z["status"][i] == -9 or nv == 0:     # 欠けた点だけ外す(status -1 も n_valid までは有効)
            continue
        T = TG
        c, valid, pc, spread, blk = R8.lookup8(tab, T, np.full(NG, z["alpha"][i]), np.full(NG, z["chi0"][i]),
                                               np.full(NG, z["f"][i]), np.full(NG, z["M"][i]), np.full(NG, z["mu"][i]))
        if z["piece"][i] == -9:                 # 544 の臨界より下: 8 準位のどれも超放射しないので、値は chi0 のまま
            c, valid = np.full(NG, float(z["chi0"][i])), np.ones(NG, bool)
        ref = z["chi"][i].astype(float)
        okref = np.zeros(NG, bool)
        okref[:nv] = True
        # 窓の中の年齢が基準の有効な範囲を出ないように、比べるのは有効な範囲から窓の幅を除いた年齢
        lt = np.log(TG)
        okref &= lt <= lt[nv - 1] - np.log1p(delta)
        d = kinn._tol_dist(TG[:nv], c[:nv], ref[:nv], delta)
        m = okref[:nv] & valid[:nv] & inT[:nv]
        err[i, :nv][m] = d[m]
        nref[i] = int((okref & inT).sum())
        ncomp[i] = int(m.sum())
    e = np.where(np.isfinite(err), err, -1.0)
    worst = e.max(axis=1)
    rows = []
    for i in np.argsort(-worst)[:15]:
        if worst[i] <= 0:
            break
        k = int(np.argmax(e[i]))
        rows.append(dict(pid=int(z["pid"][i]), seed=int(z["seed"][i]), kind=str(z["kind"][i]), piece=int(z["piece"][i]),
                         alpha=float(z["alpha"][i]), chi0=float(z["chi0"][i]), log10f=float(np.log10(z["f"][i])),
                         M=float(z["M"][i]), max_err=float(worst[i]), at_T=float(TG[k]), ref_chi=float(z["chi"][i][k])))
    by = {}
    for q in (0, -1, -2, -9):
        m = z["piece"] == q
        ages = int(np.isfinite(err[m]).sum())
        by[str(q)] = dict(points=int(m.sum()), points_compared=int(np.sum(ncomp[m] > 0)), ages=ages,
                          ages_above_tol=int(np.sum(err[m] > tol)), points_above_tol=int(np.sum(worst[m] > tol)),
                          max_err=float(worst[m].max()) if (worst[m] >= 0).any() else None,
                          p999=float(np.nanpercentile(err[m], 99.9)) if ages else None,
                          ref_ages_in_range=int(nref[m].sum()), table_coverage=float(ages / max(1, nref[m].sum())))
    summ = dict(refs=os.path.basename(refs_path), table=os.path.basename(tab_path), tol=tol, delta=delta,
                points=n, status={str(s): int(np.sum(z["status"] == s)) for s in (-9, -1, 0, 1, 2)},
                ages_compared=int(np.isfinite(err).sum()), ages_above_tol=int(np.sum(err > tol)),
                points_above_tol=int(np.sum(worst > tol)), max_err=float(worst.max()),
                median_point_max=float(np.median(worst[worst >= 0])) if (worst >= 0).any() else None,
                by_piece=by, worst=rows)
    if dest:
        json.dump(summ, open(dest, "w"), indent=1)
        np.savez(dest.replace(".json", ".npz"), pid=z["pid"], err=err)
    print(json.dumps({k: v for k, v in summ.items() if k != "worst"}, indent=1))
    return summ


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "plan":
        plan(sys.argv[2])
    elif what == "task":
        task(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], cap=int(os.environ.get("SRLVK_CAP", "10000000")),
             nproc=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    elif what == "merge":
        merge(sys.argv[2], sys.argv[3])
    elif what == "eval":
        evaluate(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else None)
    else:
        raise SystemExit(f"unknown {what}")
