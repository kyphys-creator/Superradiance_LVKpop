"""E24(予定): 8 準位の外の準位 311・522 は、8 準位の雲の上で育つか(新しい相対論的な補正で。D18 の手順 1 の
「311 は要るか確かめる」と、E21 の「522 は 422 が抑える」を、本番の表の運動論の上で確かめ直す)。

考え方: 外の準位 X が少ないあいだ、X の占有数 x = n_X / n_q(n_q は 1 個の量子)は 8 準位の雲が決める一次の式に従う:
    dx/dt = mu [ (g_X + G(t)) (x + 1) - L(t) x ]       (g_X > 0 のとき。g_X <= 0 は kinn と同じ床 g_X (x - 1e-2))
    G = sum_{a x b -> X} c_p eps_a eps_b、 L = sum_{X x a -> b} c_p eps_a (eps_b + eps_q)、 c_p = (M_pl/f)^4 kappa alpha^q [r+] C_p
8 準位の解(本番と同じ System、同じ年齢の格子)の上でこれを解き、X が雲の全体に占める割合 n_X / sum n_i の最大を見る。
X が育てば、その時刻からは一次の式が使えない(X が雲を変える)ので、そこで育つと判定する。
  rates    : 311・522 の超放射の率を rates_mp.py と同じ方法(多倍長の Leaver、40 桁、chi 110 点)で表にする。下端は
             311 が 0.03、522 が 0.05(その下は sr.gamma_mp_over_mu の解析式とのつなぎ)。311(m = 1)は alpha >= 0.35 を
             4 倍に細かくする(rates_mp.refine と同じ)。出力 out/rates_mp_311.npz、out/rates_mp_522.npz。
  procs    : X の占有数について一次の過程(得: a x b -> X、失: X x a -> b、a・b は 8 準位)を kinn と同じ規則(hyd.py)で拾う。
             出力 out/e24_procs.json。
  corr     : その過程の補正 C を relcorr.py と同じ方法(relrates.point)で、relcorr の格子の alpha >= 0.05 で解く。X の根は
             relcorr と同じく倍精度の Leaver で alpha = 0.01 から追う
             (多倍長の表とは check_roots で照らす)。出力 out/e24_corr_<nlm>.npz(relcorr.npz と同じ形)。
  runs     : 8 準位の運動論(本番の System = kinn.LEVELS8 + relcorr の補正)を、本番の節から選んだ組で mu T = 0.13 eV yr
             (値の表の範囲)まで解き、解を保存する(out/e24_runs/)。M = 1 Msun(L = ln(N0 M^2) が最小で、外の準位が
             育つのにいる e-fold がいちばん少ない側)。
  validate : 方法の確かめ: 8 準位の中の準位 Y(433・511・533・422)を外の準位とみなし、同じ一次の式で解いた x_Y が、
             Y が少ないあいだ運動論の n_Y / n_q と合うこと。
  eval     : 311・522 について一次の式を解き、節ごとの n_X / sum n_i の最大と、その年齢・印(平均の区間か)をまとめる
             (出力した年齢の点の上の後処理なので、雲が 1 目盛りに何十 e 倍も変わる所で数 e-fold の誤差がある。validate)。
  runs10   : 311・522 を運動論に加えた 10 準位(system10、加えるのは一次の過程)で同じ節を解き直す(out/e24_runs10/)。
  compare  : 8 準位と 10 準位の chi の差(時刻のずれ 0.5% を許した差、kinn._tol_dist)と、311・522 の割合の最大。
  (2026-11-21 に追加。311・522 どうしの過程も入れた 10 準位。育った 311・522 が飽和し、育ったあとの chi が出る)
  procsxx  : 311・522 の足が二つ以上の過程を kinn と同じ規則で拾う(out/e24_procsxx.json)。
  corrxx   : その過程の補正(part_corr と同じ方法、out/e24_corr_xx.npz)。
  runs10f  : 全部の過程の 10 準位(system10full)で run_nodes の節を解く(out/e24_runs10f/)。compare10f で 8 準位と比べる。
  scan     : 522 が育つ範囲を決める節(scan_nodes: 区分 1 は alpha 0.05-0.5 の全部、M = 1・300 Msun)を system10full で
             解く(out/e24_scan10f/)。scaneval で、311・522 の割合と、表の同じ節との差(包絡の chi、時刻 2% を許す)をまとめる。
出力: out/e24.json(部分ごとに書き足す)
実行: cd research && ~/miniforge3/bin/python3 calc/e24_outside.py <部分> [並列の数]
"""
import json
import os
import sys
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kinn  # noqa: E402
import si  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
RES = os.path.join(OUT, "e24.json")
RUNS = os.path.join(OUT, "e24_runs")
X_LEVELS = [(3, 1, 1), (5, 2, 2)]
ALPHA_LO_X = {(3, 1, 1): 0.03, (5, 2, 2): 0.05}
L8 = [tuple(l) for l in kinn.LEVELS8]
TAU_END = 0.13                                     # [eV yr]。値の表の範囲(reftab の TAU_MAX / 10)


def name(lv):
    return "%d%d%d" % tuple(lv)


def save_part(key, val):
    d = json.load(open(RES)) if os.path.exists(RES) else {}
    d[key] = val
    json.dump(d, open(RES, "w"), indent=1, default=float)


# ---------- rates ----------
def part_rates(nproc=12):
    import rates_mp
    rates_mp.ALPHA_LO.update(ALPHA_LO_X)
    t0 = time.time()
    rates_mp.build(X_LEVELS, nproc, resume=True)
    z = np.load(os.path.join(OUT, "rates_mp_311.npz"))
    if "refined_from" not in z.files:
        rates_mp.refine([(3, 1, 1)], nproc=nproc)
    out = {"secs": time.time() - t0}
    for lv in X_LEVELS:
        z = np.load(os.path.join(OUT, "rates_mp_%s.npz" % name(lv)))
        out[name(lv)] = dict(n_alpha=int(len(z["alpha"])), alpha_lo=float(z["alpha"][0]), n_bad=int(z["bad"].sum()),
                             n_nan=int(np.sum(~np.isfinite(z["im"]))), all_done=bool(z["done"].all()))
    save_part("rates", out)
    print(json.dumps(out), flush=True)


def _mp_check(args):
    """表の点を、項数と桁を増やした多倍長で解き直す(独立な確かめ)。"""
    import leaver
    lv, a, c, w0 = args
    w = complex(leaver.bound_state_mp(a, c, *lv, omega0=w0, N=int(2 * max(4000, 500 / a)), dps=60))
    return w


def part_ratecheck(nproc=12):
    """表の無作為な 40 点 x 2 準位を解き直し、表の値と比べる。下端で解析式との比 Q が 1 に近いことも見る。"""
    from multiprocessing import Pool
    rng = np.random.default_rng(2401)
    jobs, ref = [], []
    for lv in X_LEVELS:
        z = np.load(os.path.join(OUT, "rates_mp_%s.npz" % name(lv)))
        for _ in range(40):
            i, j = rng.integers(len(z["alpha"])), rng.integers(len(z["chi"]))
            w0 = complex(z["re"][i, j], z["im"][i, j])
            jobs.append((lv, float(z["alpha"][i]), float(z["chi"][j]), w0))
            ref.append(w0)
    with Pool(nproc) as p:
        ws = p.map(_mp_check, jobs, chunksize=1)
    rows = []
    for (lv, a, c, _), w0, w in zip(jobs, ref, ws):
        rows.append(dict(level=name(lv), alpha=a, chi=c, d_re=abs(w.real - w0.real) / abs(w0.real),
                         d_im=abs(w.imag - w0.imag) / max(abs(w0.imag), 1e-300)))
    out = {}
    for lv in X_LEVELS:
        r = [x for x in rows if x["level"] == name(lv)]
        T = sr._mp_table(*lv)
        out[name(lv)] = dict(max_rel_re=max(x["d_re"] for x in r), max_rel_im=max(x["d_im"] for x in r),
                             Q_lo_range=[float(np.min(T["Qlo"](np.linspace(0, 0.995, 50)))),
                                         float(np.max(T["Qlo"](np.linspace(0, 0.995, 50))))])
    save_part("ratecheck", out)
    print(json.dumps(out), flush=True)


# ---------- procs ----------
def lin_triples(X):
    gain = [(a, b, X) for ia, a in enumerate(L8) for b in L8[ia:]]
    loss = [(X, a, b) for a in L8 for b in L8 if a != b]
    return gain, loss


def part_procs(nproc=12):
    from multiprocessing import Pool
    todo = []
    for X in X_LEVELS:
        g, l_ = lin_triples(X)
        todo += g + l_
    with Pool(nproc) as p:
        res = p.map(kinn._proc_job, todo, chunksize=1)
    out = {}
    for X in X_LEVELS:
        g, l_ = lin_triples(X)
        keep = []
        for (i, j, k), kap, kind, n in res:
            if (i, j, k) in g or (i, j, k) in l_:
                if kind is not None and kap != 0:
                    keep.append(dict(i=name(i), j=name(j), k=name(k), kappa=kap, kind=kind, res_n=n,
                                     role="gain" if k == X else "loss"))
        out[name(X)] = keep
    json.dump(out, open(os.path.join(OUT, "e24_procs.json"), "w"), indent=1)
    summ = {x: dict(n_gain=sum(p["role"] == "gain" for p in v), n_loss=sum(p["role"] == "loss" for p in v),
                    kinds={k: sum(p["kind"] == k for p in v) for k in ("BH", "BHres", "inf")}) for x, v in out.items()}
    save_part("procs", summ)
    print(json.dumps(summ), flush=True)


def load_procs():
    d = json.load(open(os.path.join(OUT, "e24_procs.json")))
    L = lambda s: (int(s[0]), int(s[1]), int(s[2]))   # noqa: E731
    return {L(x): [dict(p, i=L(p["i"]), j=L(p["j"]), k=L(p["k"])) for p in v] for x, v in d.items()}


# ---------- corr ----------
def _roots_track(args):
    """一つの chi で、外の準位 lv の根を alpha の小さい方から倍精度の Leaver で追う(relcorr.roots_column と同じ方法: 一つ前の
    根を水素型の比で延ばした予測から始め、同じ l の隣の n との間隔の 0.3 倍より予測に近いことを確かめる。外れたら水素型から)。
    多倍長の表を待たずに解ける(表ができたら check_roots で照らす)。返り値: (chi, om, ok)。"""
    import leaver
    from scipy.optimize import newton
    warnings.simplefilter("ignore")
    c, alphas, lv = args
    n, l, m = lv
    om = np.full(len(alphas), np.nan + 0j)
    ok = np.zeros(len(alphas), bool)
    gaps = [abs(1.0 / (2 * n ** 2) - 1.0 / (2 * k ** 2)) for k in (n - 1, n + 1) if k > l]
    hre = lambda a: a * sr.omega_over_mu(a, n, l, m, c)   # noqa: E731
    prev = None
    for q, a in enumerate(alphas):
        h = hre(a)
        pred = h if prev is None else prev[1].real * h / hre(prev[0])
        tol = 0.3 * a ** 3 * min(gaps)
        guesses = ([complex(h, 0.0)] if prev is None else
                   [complex(pred, prev[1].imag * (a / prev[0]) ** (4 * l + 5)), complex(h, 0.0)])
        w = None
        for g in guesses:
            try:
                N_ = int(max(4000, 300 / a))
                ww = complex(newton(lambda w: leaver.cf_residual(w, a, c, l, m, N_), g, x1=g * (1 + 1e-7), tol=1e-20,
                                    maxiter=100, disp=False))
            except Exception:
                continue
            if abs(ww.real - pred) < tol:
                w = ww
                break
        if w is None:
            prev = None
            continue
        om[q], ok[q] = w, True
        prev = (a, w)
    return c, om, ok


def _corr_job(args):
    import relrates
    warnings.simplefilter("ignore")
    ia, ic, a, c, om8, plist, xs, omx = args
    t0 = time.time()
    omegas = {lv: complex(om8[v]) for v, lv in enumerate(L8)}
    for X, w in zip(xs, omx):
        omegas[X] = complex(w)
    res, _ = relrates.point(a, c, plist, omegas)
    return ia, ic, res, time.time() - t0


def part_corr(nproc=13, a_min=0.05, limit=None, xs=None):
    """xs: 解く外の準位(既定は両方)。準位ごとに out/e24_corr_<nlm>.npz に書く(311 の率の表ができたら、522 を待たずに解ける)。"""
    import relcorr
    from multiprocessing import Pool
    xs = X_LEVELS if xs is None else xs
    if len(xs) > 1:
        for X in xs:
            part_corr(nproc, a_min, limit, [X])
        return
    T8 = relcorr.load()
    P = load_procs()
    flat = [(X, p) for X in xs for p in P[X]]
    plist = [(p["i"], p["j"], p["k"]) for _, p in flat]
    ias = np.flatnonzero(T8["alpha"] >= a_min - 1e-12)
    A, C = T8["alpha"][ias], T8["chi"]
    shape = (len(A), len(C), len(flat))
    st = dict(C_inf=np.full(shape, np.nan), C_H=np.full(shape, np.nan), tau_X=np.full(shape, np.nan),
              omegaX=np.full(shape, np.nan), note=np.zeros(shape, np.int8), done=np.zeros(shape[:2], bool),
              root_ok=np.zeros(shape[:2] + (len(xs),), bool), sec=np.full(shape[:2], np.nan),
              omX=np.full(shape[:2], np.nan + 0j))
    path = os.path.join(OUT, "e24_corr_%s.npz" % name(xs[0]))
    if os.path.exists(path):
        z = np.load(path)
        if z["C_inf"].shape == shape:
            st = {k: z[k] for k in st}
    names = np.array([f"{name(p['i'])}x{name(p['j'])}->{name(p['k'])}" for _, p in flat])
    kinds = np.array([p["kind"] for _, p in flat])
    kap = np.array([p["kappa"] for _, p in flat])
    pw = np.array([kinn.POW[p["kind"]] for _, p in flat])
    bh = kinds != "inf"
    # 外の準位の根(chi の列ごとに alpha = 0.01 から追う)
    a_track = np.concatenate([np.geomspace(0.01, A[0], 12)[:-1], A])
    t0 = time.time()
    with Pool(min(nproc, len(C))) as pool:
        for c, om, ok in pool.imap_unordered(_roots_track, [(float(c), a_track, xs[0]) for c in C]):
            ic = int(np.argmin(np.abs(C - c)))
            st["omX"][:, ic], st["root_ok"][:, ic, 0] = om[-len(A):], ok[-len(A):]
    print(f"roots of {name(xs[0])}: {time.time() - t0:.0f}s, failed {int((~st['root_ok']).sum())} of {st['root_ok'].size}",
          flush=True)
    jobs = [(q, ic, float(A[q]), float(c), T8["om"][ia, ic], plist, xs, [st["omX"][q, ic]])
            for q, ia in enumerate(ias) for ic, c in enumerate(C)
            if not st["done"][q, ic] and st["root_ok"][q, ic].all() and T8["root_ok"][ia, ic].all()]
    if limit:
        jobs = jobs[:limit]
    print(f"corr: {len(flat)} processes, {len(jobs)} points, {nproc} workers", flush=True)

    def save():
        np.savez(path + ".tmp.npz", alpha=A, chi=C, names=names, kinds=kinds, kappa_hyd=kap,
                 X=np.array([name(X) for X, _ in flat]), role=np.array([p["role"] for _, p in flat]), **st)
        os.replace(path + ".tmp.npz", path)

    t0 = time.time()
    nd = 0
    with Pool(nproc) as pool:
        for q, ic, res, sec in pool.imap_unordered(_corr_job, jobs, chunksize=1):
            a, c = A[q], C[ic]
            st["sec"][q, ic] = sec
            if res is not None:
                rp = 1 + np.sqrt(1 - c * c)
                norm = kap * a ** (pw - 3) * np.where(bh, rp, 1.0)
                for r_i, r in enumerate(res):
                    st["C_inf"][q, ic, r_i] = r["G_inf"] / norm[r_i]
                    st["C_H"][q, ic, r_i] = r["G_h"] / norm[r_i]
                    st["tau_X"][q, ic, r_i] = r["tau_X"]
                    st["omegaX"][q, ic, r_i] = r["omegaX"]
                    note = r["note"]
                    st["note"][q, ic, r_i] = 3 if note.startswith("error") else (2 if note == "near threshold" else 0)
                st["done"][q, ic] = True
            nd += 1
            if nd % 20 == 0 or nd == len(jobs):
                save()
                el = time.time() - t0
                print(f"{int(st['done'].sum())}/{st['done'].size} points, {el:.0f}s, ~{el / nd * (len(jobs) - nd) / 60:.0f} min left",
                      flush=True)
    save()
    out = dict(n_points=int(st["done"].size), n_done=int(st["done"].sum()), n_proc=len(flat),
               n_error=int((st["note"] == 3).sum()), n_near_threshold=int((st["note"] == 2).sum()),
               median_sec=float(np.nanmedian(st["sec"])), secs=time.time() - t0)
    save_part("corr_" + name(xs[0]), out)
    print(json.dumps(out), flush=True)


def part_rootcheck():
    """補正に使った外の準位の根(倍精度で追ったもの)を、多倍長の表と比べる(relcorr の check と同じ物差し: Re の差を同じ l の
    隣の n との間隔で割ったもの)。"""
    out = {}
    for X in X_LEVELS:
        n, l, m = X
        z = np.load(os.path.join(OUT, "e24_corr_%s.npz" % name(X)))
        T = sr._mp_table(*X)
        A, C, om, ok = z["alpha"], z["chi"], z["omX"], z["root_ok"][..., 0]
        d = []
        for q, a in enumerate(A):
            gap = a ** 3 * min(abs(1.0 / (2 * n ** 2) - 1.0 / (2 * k ** 2)) for k in (n - 1, n + 1) if k > l)
            for ic, c in enumerate(C):
                if ok[q, ic] and T["alo"] <= a <= T["ahi"]:
                    d.append(abs(om[q, ic].real - float(T["splR"].ev(np.log(a), c))) / gap)
        out[name(X)] = dict(n=len(d), max_re_over_spacing=float(max(d)), n_failed=int((~ok).sum()))
    save_part("rootcheck", out)
    print(json.dumps(out), flush=True)


def load_corr():
    """X の過程の補正を relcorr.Corr で引く(relcorr と同じ整理と補間)。返り値: {X: (Corr, 過程の一覧)}。"""
    import relcorr
    P = load_procs()
    out = {}
    for X in X_LEVELS:
        z = np.load(os.path.join(OUT, "e24_corr_%s.npz" % name(X)))
        sel = z["X"] == name(X)
        T = {k: z[k][..., sel] if z[k].ndim == 3 else z[k] for k in ("C_inf", "C_H", "tau_X", "omegaX")}
        T.update(alpha=z["alpha"], chi=z["chi"], names=z["names"][sel], kinds=z["kinds"][sel])
        cr = relcorr.Corr(T)
        assert [str(s) for s in cr.names] == [f"{name(p['i'])}x{name(p['j'])}->{name(p['k'])}" for p in P[X]]
        out[X] = (cr, P[X])
    return out


# ---------- runs ----------
def run_nodes():
    """節の組: alpha 9 点 x f 5 点 x(区分 1・0 の位置 2 点ずつ + 区分 -1・-2 の位置 1 点)、M = 1 Msun。"""
    import reftab8 as R8
    A_T = [0.12, 0.16, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
    LF_T = [9.0, 11.0, 13.0, 15.0, 17.0]
    S_T = {1: [0.4, 0.95], 0: [0.4, 0.95], -1: [0.55], -2: [0.55]}
    out = []
    for at in A_T:
        ia = int(np.argmin(np.abs(np.log(R8.A_NODES / at))))
        a = R8.A_NODES[ia]
        for p in (1, 0, -1, -2):
            if a < R8.A_MIN[p] or (p == 1 and float(R8.crit(a)[0]) >= R8.CHI_TOP):
                continue
            for st in S_T[p]:
                js = int(np.argmin(np.abs(R8.S_NODES[p] - st)))
                for lt in LF_T:
                    kf = int(np.argmin(np.abs(R8.LF_NODES - lt)))
                    out.append((p, ia, js, kf, 0))
    return out


def _run_path(key):
    return os.path.join(RUNS, "run_%d_%02d_%02d_%02d_%d.npz" % key)


def _run_job(args):
    import reftab8 as R8
    warnings.simplefilter("ignore")
    key, cap = args
    path = _run_path(key)
    if os.path.exists(path):
        return key, "skip", 0.0
    a, sig, f, M, cs, chi0 = R8.node_params(*key)
    mu = a / sr.alpha_of(M, 1.0)
    tg = R8.TAU_GRID[R8.TAU_GRID <= TAU_END * (1 + 1e-12)]
    t0 = time.time()
    try:
        s = kinn.evolveN_avg(R8._system(), M, chi0, mu, f, tg[-1] / mu, tg / mu, max_nfev=cap)
        np.savez(path, y=s.y, t=s.t, avg=s.avg, status=s.status, nfev=s.nfev, njump=s.njump, mu=mu, f=f, M=M, chi0=chi0,
                 alpha0=a, sig=sig, cs=np.array(cs), secs=time.time() - t0)
        return key, "ok", time.time() - t0
    except Exception as e:                                 # noqa: BLE001
        return key, "error: " + repr(e)[:150], time.time() - t0


def part_runs(nproc=13, cap=4000000, limit=None):
    from multiprocessing import Pool
    os.makedirs(RUNS, exist_ok=True)
    nodes = run_nodes()
    if limit:
        nodes = nodes[:limit]
    print(f"runs: {len(nodes)} nodes, {nproc} workers, cap {cap}", flush=True)
    t0 = time.time()
    log = []
    with Pool(nproc) as p:
        for i, (key, msg, sec) in enumerate(p.imap_unordered(_run_job, [(k, cap) for k in nodes], chunksize=1)):
            log.append(dict(key=key, msg=msg, secs=sec))
            print(f"{i + 1}/{len(nodes)} {key} {msg} {sec:.0f}s (elapsed {time.time() - t0:.0f}s)", flush=True)
    save_part("runs", dict(n=len(nodes), n_ok=sum(x["msg"] in ("ok", "skip") for x in log),
                           errors=[x for x in log if x["msg"].startswith("error")], secs=time.time() - t0))


# ---------- 一次の式 ----------
def traj(s, meta):
    """解から年齢ごとの alpha、chi、eps_i(年齢 x 準位)、eps_q。kinn.make_rhs と同じ切り方(M は 0.5 M0 から M0、chi は 0.9999 まで)。"""
    N = len(L8)
    nq = 1.0 / sr.N0
    a_sun = sr.alpha_of(1.0, meta["mu"])
    n = np.exp(np.minimum(s.y[:N], 30.0))
    M = np.clip(s.y[N], 0.5 * meta["M"], meta["M"])
    chi = np.clip(np.maximum(s.y[N + 1], 0.0) / M ** 2, 0.0, 0.9999)
    return a_sun * M, chi, (n / M ** 2).T, nq / M ** 2


def _logmean(u0, u1):
    """区間の中で log u が時刻に比例して変わるときの、u の時間平均(対数平均)。"""
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        lr = np.log(u1) - np.log(u0)
        lm = np.where(np.abs(lr) < 1e-8, 0.5 * (u0 + u1), (u1 - u0) / lr)
    return np.where((u0 > 0) & (u1 > 0), lm, 0.5 * (u0 + u1))


def linear_x(t_s, mu, f, e, eq, cg, Ig, Jg, cl, Pl, Kl, g):
    """外の準位の x = n_X / n_q を解く。cg・cl は年齢 x 過程の係数(kappa alpha^q [r+] C)、Ig・Jg は得の過程の二つの準位、
    Pl・Kl は失の過程 X x a -> b の a と b、g は X の超放射の率 / mu(占有数の率)。返り値: ln x(年齢ごと)。
    年齢の格子の 1 区間の中で、雲の占有数の対数と係数は時刻に比例して変わるとみなす: 占有数の積の時間平均は対数平均
    (雲が 1 区間に何 e 倍も育つ所でも正しい)、係数と g は両端の平均。区間ごとに dx/dt = A x + B を厳密に進める。"""
    s4 = (si.MPL_GEV / f) ** 4
    mu_s = mu / sr.HBAR_EVS
    dt = np.diff(t_s) * mu_s
    Pg = e[:, Ig] * e[:, Jg]
    Q1 = e[:, Pl] * e[:, Kl]
    Q2 = e[:, Pl] * eq[:, None]
    mid = lambda c: 0.5 * (c[1:] + c[:-1])   # noqa: E731
    G = s4 * np.sum(mid(cg) * _logmean(Pg[:-1], Pg[1:]), axis=1)
    L = s4 * np.sum(mid(cl) * (_logmean(Q1[:-1], Q1[1:]) + _logmean(Q2[:-1], Q2[1:])), axis=1)
    gm = mid(g)
    A = gm + G - L
    B = np.where(gm > 0, gm + G, G - 1e-2 * gm)
    z = A * dt
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        # 源の項 B (e^z - 1)/A の対数を、z が大きくても溢れない形で
        lsrc = np.where(np.abs(z) < 1e-12, np.log(B * dt),
                        np.where(z > 0, np.log(B / A) + z + np.log(-np.expm1(-z)), np.log(B / -A) + np.log(-np.expm1(z))))
    lsrc = np.where(B > 0, lsrc, -np.inf)
    lnx = np.zeros(len(t_s))
    for t in range(1, len(t_s)):
        lnx[t] = np.logaddexp(lnx[t - 1] + z[t - 1], lsrc[t - 1])
    return lnx


def _sys8():
    import reftab8 as R8
    return R8._system()


def own_index(sysN, Y):
    """8 準位の中の準位 Y を外の準位とみなしたときの、得と失の過程の番号と準位(Y どうしの過程は除く)。"""
    y = sysN.levels.index(Y)
    I, J, K = sysN.I, sysN.J, sysN.K
    gq = np.flatnonzero((K == y) & (I != y) & (J != y))
    lqi = np.flatnonzero((I == y) & (J != y))        # Y x j -> k: 相手は J
    lqj = np.flatnonzero((J == y) & (I != y))        # i x Y -> k: 相手は I
    return gq, I[gq], J[gq], np.concatenate([lqi, lqj]), np.concatenate([J[lqi], I[lqj]]), np.concatenate([K[lqi], K[lqj]])


def x_index(plist):
    idx = {lv: v for v, lv in enumerate(L8)}
    gq = np.array([q for q, p in enumerate(plist) if p["role"] == "gain"], int)
    lq = np.array([q for q, p in enumerate(plist) if p["role"] == "loss"], int)
    return (gq, np.array([idx[plist[q]["i"]] for q in gq], int), np.array([idx[plist[q]["j"]] for q in gq], int),
            lq, np.array([idx[plist[q]["j"]] for q in lq], int), np.array([idx[plist[q]["k"]] for q in lq], int))


def x_coefs(cr, plist, a, chi):
    """年齢ごとの X の過程の係数 kappa alpha^q [r+] (C_inf + C_H)(年齢 x 過程)。"""
    kap = np.array([p["kappa"] for p in plist])
    pw = np.array([kinn.POW[p["kind"]] for p in plist])
    bh = np.array([p["kind"] != "inf" for p in plist])
    out = np.zeros((len(a), len(plist)))
    for t in range(len(a)):
        ci, ch = cr(a[t], chi[t])
        out[t] = kap * a[t] ** pw * np.where(bh, sr.rplus(chi[t]), 1.0) * (ci + ch)
    return out


def _load_run(path):
    z = np.load(path)

    class S:
        pass
    s = S()
    s.y, s.t, s.avg = z["y"], z["t"], z["avg"]
    return s, {k: float(z[k]) for k in ("mu", "f", "M", "chi0", "alpha0", "sig", "status", "nfev", "njump", "secs")}


def _validate_job(path):
    sysN = _sys8()
    s, meta = _load_run(path)
    a, chi, e, eq = traj(s, meta)
    nq = 1.0 / sr.N0
    N = len(L8)
    nt = len(a)
    cf = np.zeros((nt, len(sysN.procs)))
    g8 = np.zeros((nt, N))
    for t in range(nt):
        ci, ch = sysN.coef(a[t], chi[t])
        cf[t] = ci + ch
        g8[t] = sysN.sr_rate(a[t], chi[t], sysN.levels)
    lnn = s.y[:N]
    tot = np.log(np.sum(np.exp(np.minimum(lnn, 30.0)), axis=0))
    first_avg = int(np.argmax(s.avg != 0)) if np.any(s.avg != 0) else len(s.avg)
    out = []
    for Y in [(4, 3, 3), (5, 1, 1), (5, 3, 3), (4, 2, 2), (4, 1, 1)]:
        gq, Ig, Jg, lq, Pl, Kl = own_index(sysN, Y)
        y = L8.index(Y)
        lnx = linear_x(s.t, meta["mu"], meta["f"], e, eq, cf[:, gq], Ig, Jg, cf[:, lq], Pl, Kl, g8[:, y])
        lnY = lnn[y] - np.log(nq)
        frac = lnn[y] - tot                              # 本当の Y の割合
        # Y が初めて雲の 1e-4 を越えるより前で、平均の区間でない年齢だけで比べる(越えたあとは Y が雲を変え、一次の式は使えない)
        # (雲ができる前は、どの準位も 1 個の量子で割合が大きく見えるので、全体が 1e6 個を越えてから数える)
        over = np.flatnonzero((frac >= np.log(1e-4)) & (tot > np.log(nq) + np.log(1e6)))
        before = np.arange(len(frac)) < (over[0] if len(over) else len(frac))
        sel = before & (s.avg == 0) & np.isfinite(lnY)
        if sel.sum() < 5:
            continue
        d = np.abs(lnx[sel] - lnY[sel])
        # 最初の平均の区間より前(雲の振動を平均していない解だけ)と、全体(平均の区間のあとも含む)
        pre = sel & (np.arange(len(sel)) < first_avg)
        dpre = np.abs(lnx[pre] - lnY[pre])
        out.append(dict(level=name(Y), n=int(sel.sum()), max_abs_dlnx=float(d.max()), median_abs_dlnx=float(np.median(d)),
                        n_pre=int(pre.sum()), max_abs_dlnx_pre=float(dpre.max()) if pre.any() else None,
                        max_lnY=float(lnY[sel].max()), max_lnY_pre=float(lnY[pre].max()) if pre.any() else None))
    return os.path.basename(path), meta, out


def part_validate(nproc=13):
    from multiprocessing import Pool
    paths = sorted(os.path.join(RUNS, f) for f in os.listdir(RUNS) if f.endswith(".npz"))
    with Pool(nproc) as p:
        res = p.map(_validate_job, paths, chunksize=1)
    rows = [dict(run=r, alpha0=m["alpha0"], lf=float(np.log10(m["f"])), **o) for r, m, oo in res for o in oo]
    summ = {}
    for lv in sorted({r["level"] for r in rows}):
        rr = [r for r in rows if r["level"] == lv]
        pre = [r for r in rr if r["max_abs_dlnx_pre"] is not None]
        big = sorted(rr, key=lambda r: -r["max_abs_dlnx"])[:3]
        summ[lv] = dict(n_runs=len(rr), max_abs_dlnx=max(r["max_abs_dlnx"] for r in rr),
                        median_of_median=float(np.median([r["median_abs_dlnx"] for r in rr])),
                        max_lnY=max(r["max_lnY"] for r in rr), n_runs_pre=len(pre),
                        max_abs_dlnx_pre=max([r["max_abs_dlnx_pre"] for r in pre] or [None]),
                        max_lnY_pre=max([r["max_lnY_pre"] for r in pre] or [None]),
                        worst=[dict(run=r["run"], alpha0=r["alpha0"], lf=r["lf"], d=r["max_abs_dlnx"], d_pre=r["max_abs_dlnx_pre"])
                               for r in big])
    save_part("validate", dict(summary=summ, rows=rows))
    print(json.dumps(summ, indent=1), flush=True)


def _eval_job(path):
    corr = load_corr()
    s, meta = _load_run(path)
    a, chi, e, eq = traj(s, meta)
    nq = 1.0 / sr.N0
    N = len(L8)
    tot = np.log(np.sum(np.exp(np.minimum(s.y[:N], 30.0)), axis=0))
    tau = s.t / sr.YR * meta["mu"]
    res = {}
    for X in X_LEVELS:
        cr, plist = corr[X]
        cf = x_coefs(cr, plist, a, chi)
        gq, Ig, Jg, lq, Pl, Kl = x_index(plist)
        g = np.array([sr.gamma_mp_over_mu(a[t], chi[t], *X, sr.NORM_N) for t in range(len(a))])
        lnx = linear_x(s.t, meta["mu"], meta["f"], e, eq, cf[:, gq], Ig, Jg, cf[:, lq], Pl, Kl, g)
        lfrac = lnx + np.log(nq) - tot
        # 雲ができる前(全体が 1 個の量子の 1e6 倍より少ない年齢)は割合を数えない
        ok = tot > np.log(nq) + np.log(1e6)
        k = int(np.argmax(np.where(ok, lfrac, -np.inf))) if ok.any() else 0
        res[name(X)] = dict(max_log10_frac=float(lfrac[k] / np.log(10)) if ok.any() else None,
                            tau_at_max=float(tau[k]), avg_at_max=int(s.avg[k]), max_lnx=float(lnx.max()),
                            final_log10_frac=float(lfrac[-1] / np.log(10)), n_avg=int(np.sum(s.avg != 0)))
    return os.path.basename(path), meta, res


def part_eval(nproc=13):
    from multiprocessing import Pool
    paths = sorted(os.path.join(RUNS, f) for f in os.listdir(RUNS) if f.endswith(".npz"))
    with Pool(nproc) as p:
        res = p.map(_eval_job, paths, chunksize=1)
    rows = [dict(run=r, alpha0=m["alpha0"], lf=float(np.log10(m["f"])), chi0=m["chi0"], status=m["status"],
                 tau_end=None, **{f"{x}_{k}": v for x, d in o.items() for k, v in d.items()}) for r, m, o in res]
    summ = {}
    for X in X_LEVELS:
        v = [r[f"{name(X)}_max_log10_frac"] for r in rows if r[f"{name(X)}_max_log10_frac"] is not None]
        top = sorted(rows, key=lambda r: -(-1e9 if r[f"{name(X)}_max_log10_frac"] is None else r[f"{name(X)}_max_log10_frac"]))[:8]
        summ[name(X)] = dict(n=len(v), max_log10_frac=max(v), n_above_1e_3=int(sum(x > -3 for x in v)),
                             n_above_1e_6=int(sum(x > -6 for x in v)),
                             top=[dict(run=r["run"], alpha0=r["alpha0"], lf=r["lf"], chi0=r["chi0"],
                                       log10_frac=r[f"{name(X)}_max_log10_frac"], tau=r[f"{name(X)}_tau_at_max"],
                                       avg=r[f"{name(X)}_avg_at_max"]) for r in top])
    save_part("eval", dict(summary=summ, rows=rows))
    print(json.dumps(summ, indent=1), flush=True)


# ---------- runs10・compare: 311・522 を運動論に加えて解き直す ----------
RUNS10 = os.path.join(OUT, "e24_runs10")


class CombinedCorr:
    """複数の relcorr.Corr を一つにまとめる(kinn.System は names で過程を対応させる)。"""

    def __init__(self, parts):
        self.parts = parts
        self.names = [str(nm) for c in parts for nm in c.names]

    def __call__(self, alpha, chi):
        out = [c(alpha, chi) for c in self.parts]
        return np.concatenate([o[0] for o in out]), np.concatenate([o[1] for o in out])

    def tau(self, alpha, chi):
        return np.concatenate([c.tau(alpha, chi) for c in self.parts])


_SYS10 = None


def system10():
    """本番の 8 準位の運動論に 311・522 を加えたもの。加える過程は、311・522 の占有数について一次の過程(得: 8 準位の
    二つから作られる、失: 8 準位の一つと散乱する)だけ。311・522 が少ないあいだは一次の式を運動論の刻みで正確に解くのと同じで、
    育てば 8 準位への反作用も入る(311・522 どうしの過程は入らない)。"""
    global _SYS10
    if _SYS10 is None:
        import relcorr
        levels = L8 + X_LEVELS
        procs = [dict(p) for p in kinn.process_table()]
        corr = load_corr()
        for X in X_LEVELS:
            for p in corr[X][1]:
                procs.append(dict(i=levels.index(p["i"]), j=levels.index(p["j"]), k=levels.index(p["k"]), kappa=p["kappa"],
                                  kind=p["kind"], res_n=p["res_n"]))
        cc = CombinedCorr([relcorr.Corr(relcorr.load())] + [corr[X][0] for X in X_LEVELS])
        _SYS10 = kinn.System(levels, procs, corr=cc)
    return _SYS10


def _run10_job(args):
    import reftab8 as R8
    warnings.simplefilter("ignore")
    key, cap = args
    path = os.path.join(RUNS10, os.path.basename(_run_path(key)))
    if os.path.exists(path):
        return key, "skip", 0.0
    a, sig, f, M, cs, chi0 = R8.node_params(*key)
    mu = a / sr.alpha_of(M, 1.0)
    tg = R8.TAU_GRID[R8.TAU_GRID <= TAU_END * (1 + 1e-12)]
    t0 = time.time()
    try:
        s = kinn.evolveN_avg(system10(), M, chi0, mu, f, tg[-1] / mu, tg / mu, max_nfev=cap)
        np.savez(path, y=s.y, t=s.t, avg=s.avg, status=s.status, nfev=s.nfev, njump=s.njump, mu=mu, f=f, M=M, chi0=chi0,
                 alpha0=a, sig=sig, cs=np.array(cs), secs=time.time() - t0)
        return key, "ok", time.time() - t0
    except Exception as e:                                 # noqa: BLE001
        return key, "error: " + repr(e)[:150], time.time() - t0


def part_runs10(nproc=13, cap=4000000, limit=None):
    from multiprocessing import Pool
    os.makedirs(RUNS10, exist_ok=True)
    nodes = run_nodes()
    if limit:
        nodes = nodes[:limit]
    print(f"runs10: {len(nodes)} nodes, {nproc} workers, cap {cap}", flush=True)
    t0 = time.time()
    log = []
    with Pool(nproc) as p:
        for i, (key, msg, sec) in enumerate(p.imap_unordered(_run10_job, [(k, cap) for k in nodes], chunksize=1)):
            log.append(dict(key=key, msg=msg, secs=sec))
            print(f"{i + 1}/{len(nodes)} {key} {msg} {sec:.0f}s (elapsed {time.time() - t0:.0f}s)", flush=True)
    save_part("runs10", dict(n=len(nodes), n_ok=sum(x["msg"] in ("ok", "skip") for x in log),
                             errors=[x for x in log if x["msg"].startswith("error")], secs=time.time() - t0))


def _compare_job(fn):
    s8, m8 = _load_run(os.path.join(RUNS, fn))
    z = np.load(os.path.join(RUNS10, fn))
    nq = 1.0 / sr.N0
    n = min(s8.y.shape[1], z["y"].shape[1])
    N8, N10 = 8, 10
    c8 = s8.y[N8 + 1, :n] / s8.y[N8, :n] ** 2
    c10 = z["y"][N10 + 1, :n] / z["y"][N10, :n] ** 2
    t = s8.t[:n]
    d_tol = kinn._tol_dist(t, c8, c10, 5e-3)
    d_tol2 = kinn._tol_dist(t, c8, c10, 2e-2)         # 表の約束と同じ、時刻の 2% のずれを許す
    d_raw = np.abs(c8 - c10)
    k2 = int(np.nanargmax(d_tol2))
    lnn = z["y"][:N10, :n]
    tot = np.log(np.sum(np.exp(np.minimum(lnn, 30.0)), axis=0))
    formed = tot > np.log(nq) + np.log(1e6)
    out = dict(run=fn, alpha0=m8["alpha0"], lf=float(np.log10(m8["f"])), chi0=m8["chi0"], n=int(n),
               status8=int(m8["status"]), status10=int(z["status"]), n8=int(s8.y.shape[1]), n10=int(z["y"].shape[1]),
               njump8=int(m8["njump"]), njump10=int(z["njump"]),
               max_dchi_tol=float(np.nanmax(d_tol)), max_dchi_raw=float(np.max(d_raw)),
               max_dchi_tol2=float(d_tol2[k2]), avg_at_max2=bool(s8.avg[k2] or z["avg"][k2]),
               chi_end8=float(c8[-1]), chi_end10=float(c10[-1]))
    tau = t / sr.YR * m8["mu"]
    out["tau_at_max2"] = float(tau[k2])
    for v, X in zip((8, 9), X_LEVELS):
        lf = np.where(formed, lnn[v] - tot, -np.inf)
        k = int(np.argmax(lf))
        out[f"{name(X)}_max_log10_frac"] = float(lf[k] / np.log(10)) if formed.any() else None
        out[f"{name(X)}_tau_at_max"] = float(tau[k])
    return out


def part_compare(nproc=13):
    from multiprocessing import Pool
    fns = sorted(f for f in os.listdir(RUNS10) if f.endswith(".npz") and os.path.exists(os.path.join(RUNS, f)))
    with Pool(nproc) as p:
        rows = p.map(_compare_job, fns, chunksize=1)
    summ = dict(n=len(rows), n_status_mismatch=int(sum(r["status8"] != r["status10"] for r in rows)),
                max_dchi_tol=max(r["max_dchi_tol"] for r in rows), max_dchi_raw=max(r["max_dchi_raw"] for r in rows),
                n_dchi_tol_above_5e_3=int(sum(r["max_dchi_tol"] > 5e-3 for r in rows)),
                max_dchi_tol2=max(r["max_dchi_tol2"] for r in rows),
                n_dchi_tol2_above_5e_3=int(sum(r["max_dchi_tol2"] > 5e-3 for r in rows)),
                n_dchi_tol2_above_0_03=int(sum(r["max_dchi_tol2"] > 0.03 for r in rows)),
                n_dchi_tol2_above_5e_3_outside_avg=int(sum(r["max_dchi_tol2"] > 5e-3 and not r["avg_at_max2"]
                                                           for r in rows)),
                worst_dchi=sorted(rows, key=lambda r: -r["max_dchi_tol2"])[:8])
    for X in X_LEVELS:
        v = [r[f"{name(X)}_max_log10_frac"] for r in rows if r[f"{name(X)}_max_log10_frac"] is not None]
        summ[name(X)] = dict(max_log10_frac=max(v), n_above_1e_3=int(sum(x > -3 for x in v)),
                             n_above_1e_6=int(sum(x > -6 for x in v)),
                             top=sorted(rows, key=lambda r: -(-1e9 if r[f"{name(X)}_max_log10_frac"] is None else r[f"{name(X)}_max_log10_frac"]))[:5])
    save_part("compare", dict(summary=summ, rows=rows))
    print(json.dumps({k: v for k, v in summ.items() if k not in ("worst_dchi",)}, indent=1, default=float)[:3000], flush=True)


# ---------- 311・522 どうしの過程も入れた 10 準位(2026-11-21、ユーザー「1 と 2 を同時に」) ----------
# system10 は 311・522 について一次の過程だけなので、育った 311・522 が飽和せず、育ったあとの chi は出ない。ここでは、
# 10 準位のあいだの過程を残らず入れる(311・522 の足が二つ以上の過程を足す)。
#   procsxx  : その過程を kinn と同じ規則(kinn._proc_job)で拾う。出力 out/e24_procsxx.json
#   corrxx   : その過程の補正 C を part_corr と同じ方法で(311・522 の根は e24_corr_<nlm>.npz に保存したものを使う)。
#              出力 out/e24_corr_xx.npz
#   runs10f  : system10full(10 準位・全部の過程)で run_nodes の節を解く(out/e24_runs10f/)。scan は節の一覧を広げたもの
#   compare10f: 8 準位・一次だけの 10 準位(runs10)と比べる
L10 = L8 + X_LEVELS
RUNS10F = os.path.join(OUT, "e24_runs10f")
PROCXX = os.path.join(OUT, "e24_procsxx.json")


def xx_triples():
    """10 準位の過程の候補のうち、311・522 の足が二つ以上のもの(kinn.process_table と同じ並び: i <= j(10 準位の番号)、
    k は i・j と違う)。足が 0 本(8 準位の過程)と 1 本(lin_triples)と合わせて、候補の全部になることを確かめる。"""
    keys = [(a, b, c) for a in range(len(L10)) for b in range(a, len(L10)) for c in range(len(L10)) if c not in (a, b)]
    nx = [sum(x >= len(L8) for x in k) for k in keys]
    n_lin = sum(len(g) + len(l_) for g, l_ in (lin_triples(X) for X in X_LEVELS))
    n8 = sum(1 for a in range(len(L8)) for b in range(a, len(L8)) for c in range(len(L8)) if c not in (a, b))
    assert sum(v == 0 for v in nx) == n8 and sum(v == 1 for v in nx) == n_lin
    return [(L10[a], L10[b], L10[c]) for (a, b, c), v in zip(keys, nx) if v >= 2]


def part_procsxx(nproc=12):
    from multiprocessing import Pool
    todo = xx_triples()
    with Pool(min(nproc, len(todo))) as p:
        res = p.map(kinn._proc_job, todo, chunksize=1)
    keep = [dict(i=name(i), j=name(j), k=name(k), kappa=kap, kind=kind, res_n=n)
            for (i, j, k), kap, kind, n in res if kind is not None and kap != 0]
    json.dump(keep, open(PROCXX, "w"), indent=1)
    summ = dict(n_candidates=len(todo), n_nonzero=len(keep), kinds={k: sum(p["kind"] == k for p in keep) for k in ("BH", "BHres", "inf")},
                names=[f"{p['i']}x{p['j']}->{p['k']}" for p in keep])
    save_part("procsxx", summ)
    print(json.dumps(summ), flush=True)


def load_procsxx():
    L = lambda s: (int(s[0]), int(s[1]), int(s[2]))   # noqa: E731
    return [dict(p, i=L(p["i"]), j=L(p["j"]), k=L(p["k"])) for p in json.load(open(PROCXX))]


def part_corrxx(nproc=13, limit=None):
    """311・522 どうしの過程の補正(part_corr と同じ格子・同じ正規化)。根は e24_corr_311/522.npz の omX を使う。"""
    import relcorr
    from multiprocessing import Pool
    T8 = relcorr.load()
    P = load_procsxx()
    plist = [(p["i"], p["j"], p["k"]) for p in P]
    zx = {X: np.load(os.path.join(OUT, "e24_corr_%s.npz" % name(X))) for X in X_LEVELS}
    A, C = zx[X_LEVELS[0]]["alpha"], zx[X_LEVELS[0]]["chi"]
    assert all(np.array_equal(zx[X]["alpha"], A) and np.array_equal(zx[X]["chi"], C) for X in X_LEVELS)
    ias = np.array([int(np.argmin(np.abs(T8["alpha"] - a))) for a in A])
    assert np.allclose(T8["alpha"][ias], A) and np.array_equal(T8["chi"], C)
    shape = (len(A), len(C), len(P))
    st = dict(C_inf=np.full(shape, np.nan), C_H=np.full(shape, np.nan), tau_X=np.full(shape, np.nan),
              omegaX=np.full(shape, np.nan), note=np.zeros(shape, np.int8), done=np.zeros(shape[:2], bool),
              sec=np.full(shape[:2], np.nan))
    path = os.path.join(OUT, "e24_corr_xx.npz")
    if os.path.exists(path):
        z = np.load(path)
        if z["C_inf"].shape == shape:
            st = {k: z[k] for k in st}
    names = np.array([f"{name(p['i'])}x{name(p['j'])}->{name(p['k'])}" for p in P])
    kinds = np.array([p["kind"] for p in P])
    kap = np.array([p["kappa"] for p in P])
    pw = np.array([kinn.POW[p["kind"]] for p in P])
    bh = kinds != "inf"
    rok = np.all([zx[X]["root_ok"][..., 0] for X in X_LEVELS], axis=0)
    jobs = [(q, ic, float(A[q]), float(c), T8["om"][ia, ic], plist, X_LEVELS, [zx[X]["omX"][q, ic] for X in X_LEVELS])
            for q, ia in enumerate(ias) for ic, c in enumerate(C)
            if not st["done"][q, ic] and rok[q, ic] and T8["root_ok"][ia, ic].all()]
    if limit:
        jobs = jobs[:limit]
    print(f"corrxx: {len(P)} processes, {len(jobs)} points, {nproc} workers", flush=True)

    def save():
        np.savez(path + ".tmp.npz", alpha=A, chi=C, names=names, kinds=kinds, kappa_hyd=kap, **st)
        os.replace(path + ".tmp.npz", path)

    t0 = time.time()
    nd = 0
    with Pool(nproc) as pool:
        for q, ic, res, sec in pool.imap_unordered(_corr_job, jobs, chunksize=1):
            a, c = A[q], C[ic]
            st["sec"][q, ic] = sec
            if res is not None:
                rp = 1 + np.sqrt(1 - c * c)
                norm = kap * a ** (pw - 3) * np.where(bh, rp, 1.0)
                for r_i, r in enumerate(res):
                    st["C_inf"][q, ic, r_i] = r["G_inf"] / norm[r_i]
                    st["C_H"][q, ic, r_i] = r["G_h"] / norm[r_i]
                    st["tau_X"][q, ic, r_i] = r["tau_X"]
                    st["omegaX"][q, ic, r_i] = r["omegaX"]
                    note = r["note"]
                    st["note"][q, ic, r_i] = 3 if note.startswith("error") else (2 if note == "near threshold" else 0)
                st["done"][q, ic] = True
            nd += 1
            if nd % 20 == 0 or nd == len(jobs):
                save()
                el = time.time() - t0
                print(f"{int(st['done'].sum())}/{st['done'].size} points, {el:.0f}s, ~{el / nd * (len(jobs) - nd) / 60:.0f} min left",
                      flush=True)
    save()
    out = dict(n_points=int(st["done"].size), n_done=int(st["done"].sum()), n_proc=len(P),
               n_error=int((st["note"] == 3).sum()), n_near_threshold=int((st["note"] == 2).sum()),
               median_sec=float(np.nanmedian(st["sec"])), secs=time.time() - t0)
    save_part("corrxx", out)
    print(json.dumps(out), flush=True)


def load_corrxx():
    import relcorr
    P = load_procsxx()
    z = np.load(os.path.join(OUT, "e24_corr_xx.npz"))
    T = {k: z[k] for k in ("C_inf", "C_H", "tau_X", "omegaX", "alpha", "chi", "names", "kinds")}
    cr = relcorr.Corr(T)
    assert [str(s) for s in cr.names] == [f"{name(p['i'])}x{name(p['j'])}->{name(p['k'])}" for p in P]
    return cr, P


_SYS10F = None


def system10full():
    """10 準位(8 準位 + 311・522)のあいだの過程を残らず入れた運動論。補正は 8 準位の relcorr、311・522 の一次の過程
    (e24_corr_311/522)、311・522 どうしの過程(e24_corr_xx)。"""
    global _SYS10F
    if _SYS10F is None:
        import relcorr
        procs = [dict(p) for p in kinn.process_table()]
        corr = load_corr()
        cxx, pxx = load_corrxx()
        for plist in [corr[X][1] for X in X_LEVELS] + [pxx]:
            for p in plist:
                procs.append(dict(i=L10.index(p["i"]), j=L10.index(p["j"]), k=L10.index(p["k"]), kappa=p["kappa"],
                                  kind=p["kind"], res_n=p["res_n"]))
        cc = CombinedCorr([relcorr.Corr(relcorr.load())] + [corr[X][0] for X in X_LEVELS] + [cxx])
        _SYS10F = kinn.System(L10, procs, corr=cc)
    return _SYS10F


def _run10f_job(args):
    import reftab8 as R8
    warnings.simplefilter("ignore")
    key, cap, outdir = args
    path = os.path.join(outdir, os.path.basename(_run_path(key)))
    if os.path.exists(path):
        return key, "skip", 0.0
    a, sig, f, M, cs, chi0 = R8.node_params(*key)
    mu = a / sr.alpha_of(M, 1.0)
    tg = R8.TAU_GRID[R8.TAU_GRID <= TAU_END * (1 + 1e-12)]
    t0 = time.time()
    try:
        s = kinn.evolveN_avg(system10full(), M, chi0, mu, f, tg[-1] / mu, tg / mu, max_nfev=cap)
        np.savez(path, y=s.y, t=s.t, avg=s.avg, status=s.status, nfev=s.nfev, njump=s.njump, mu=mu, f=f, M=M, chi0=chi0,
                 alpha0=a, sig=sig, cs=np.array(cs), secs=time.time() - t0)
        return key, "ok", time.time() - t0
    except Exception as e:                                 # noqa: BLE001
        return key, "error: " + repr(e)[:150], time.time() - t0


def part_runs10f(nproc=13, cap=4000000, nodes=None, outdir=RUNS10F, tag="runs10f"):
    from multiprocessing import Pool
    os.makedirs(outdir, exist_ok=True)
    nodes = run_nodes() if nodes is None else nodes
    print(f"{tag}: {len(nodes)} nodes, {nproc} workers, cap {cap}", flush=True)
    t0 = time.time()
    log = []
    with Pool(nproc) as p:
        for i, (key, msg, sec) in enumerate(p.imap_unordered(_run10f_job, [(k, cap, outdir) for k in nodes], chunksize=1)):
            log.append(dict(key=key, msg=msg, secs=sec))
            print(f"{i + 1}/{len(nodes)} {key} {msg} {sec:.0f}s (elapsed {time.time() - t0:.0f}s)", flush=True)
    save_part(tag, dict(n=len(nodes), n_ok=sum(x["msg"] in ("ok", "skip") for x in log),
                        errors=[x for x in log if x["msg"].startswith("error")], secs=time.time() - t0))


def _chi_frac(z, N):
    """保存した解から chi と、準位ごとの雲の全体に占める割合(log10、雲ができる前は -inf)。"""
    nq = 1.0 / sr.N0
    c = z["y"][N + 1] / z["y"][N] ** 2
    lnn = z["y"][:N]
    tot = np.log(np.sum(np.exp(np.minimum(lnn, 30.0)), axis=0))
    formed = tot > np.log(nq) + np.log(1e6)
    frac = np.where(formed[None, :], (lnn - tot[None, :]) / np.log(10), -np.inf)
    return c, frac, formed


def _compare10f_job(args):
    fn, d10f = args
    z8 = np.load(os.path.join(RUNS, fn))
    zf = np.load(os.path.join(d10f, fn))
    c8, f8, _ = _chi_frac(z8, 8)
    cf, ff, formed = _chi_frac(zf, 10)
    n = min(len(c8), len(cf))
    t = zf["t"][:n]
    d2 = kinn._tol_dist(t, c8[:n], cf[:n], 2e-2)
    k2 = int(np.nanargmax(d2))
    tau = zf["t"] / sr.YR * float(zf["mu"])
    out = dict(run=fn, alpha0=float(zf["alpha0"]), lf=float(np.log10(zf["f"])), chi0=float(zf["chi0"]), M=float(zf["M"]),
               n8=len(c8), n10f=len(cf), status10f=int(zf["status"]), njump10f=int(zf["njump"]), secs10f=float(zf["secs"]),
               max_dchi_8_vs_10f=float(d2[k2]), tau_at_max=float(tau[k2]), chi_end8=float(c8[-1]), chi_end10f=float(cf[-1]),
               tau_end10f=float(tau[len(cf) - 1]))
    p10 = os.path.join(RUNS10, fn)
    if os.path.exists(p10):
        z10 = np.load(p10)
        c10, _, _ = _chi_frac(z10, 10)
        out["chi_end10_linear"] = float(c10[-1])
    for v, X in zip((8, 9), X_LEVELS):
        k = int(np.argmax(ff[v]))
        out[f"{name(X)}_max_log10_frac"] = float(ff[v][k]) if formed.any() else None
        out[f"{name(X)}_tau_at_max"] = float(tau[k])
        out[f"{name(X)}_end_log10_frac"] = float(ff[v][len(cf) - 1]) if formed[len(cf) - 1] else None
    # どの準位が最後に雲の大部分を持つか
    if formed[len(cf) - 1]:
        out["dominant_end"] = [name(L10[i]) for i in np.argsort(-ff[:, len(cf) - 1])[:3]]
    return out


def part_compare10f(nproc=13, d10f=RUNS10F, tag="compare10f"):
    from multiprocessing import Pool
    fns = sorted(f for f in os.listdir(d10f) if f.endswith(".npz") and os.path.exists(os.path.join(RUNS, f)))
    with Pool(nproc) as p:
        rows = p.map(_compare10f_job, [(f, d10f) for f in fns], chunksize=1)
    grow = lambda r, X: r[f"{name(X)}_max_log10_frac"] is not None and r[f"{name(X)}_max_log10_frac"] > -2   # noqa: E731
    summ = dict(n=len(rows), n_522_above_1pct=sum(grow(r, (5, 2, 2)) for r in rows),
                n_311_above_1pct=sum(grow(r, (3, 1, 1)) for r in rows),
                n_dchi_above_0_03=sum(r["max_dchi_8_vs_10f"] > 0.03 for r in rows),
                n_dchi_above_5e_3=sum(r["max_dchi_8_vs_10f"] > 5e-3 for r in rows),
                n_dchi_above_0_03_without_X=sum(r["max_dchi_8_vs_10f"] > 0.03 and not grow(r, (5, 2, 2)) and not grow(r, (3, 1, 1))
                                                for r in rows),
                n_status_not0=sum(r["status10f"] != 0 for r in rows))
    save_part(tag, dict(summary=summ, rows=rows))
    print(json.dumps(summ, indent=1), flush=True)
    return rows


SCAN10F = os.path.join(OUT, "e24_scan10f")


def scan_nodes():
    """522 が育つ範囲を決めるための節(run_nodes で抜けていた所を埋める): 区分 1 は alpha 0.05-0.5 の全部の節(26 点)、
    f 5 点、位置 2 点、M = 1・300 Msun。区分 0・-1・-2 は run_nodes と同じ組を M = 300 Msun で。"""
    import reftab8 as R8
    LF_T = [9.0, 11.0, 13.0, 15.0, 17.0]
    S_T = {1: [0.4, 0.95], 0: [0.4, 0.95], -1: [0.55], -2: [0.55]}
    out = []
    for ia, a in enumerate(R8.A_NODES):
        if a < 0.05 or float(R8.crit(a)[0]) >= R8.CHI_TOP:
            continue
        for st in S_T[1]:
            js = int(np.argmin(np.abs(R8.S_NODES[1] - st)))
            for lt in LF_T:
                kf = int(np.argmin(np.abs(R8.LF_NODES - lt)))
                for il in (0, 1):
                    out.append((1, ia, js, kf, il))
    out += [(p, ia, js, kf, 1) for (p, ia, js, kf, il) in run_nodes() if p != 1]
    return out


def _table_chi(LT, p, sig, cs, lt10):
    """表の節の準位(LT: 包絡が P8 を越える log10 tau)から、log10 tau の格子での chi(包絡)。"""
    import reftab8 as R8
    LT = np.where(np.isneginf(LT), -1e9, LT)       # 準位 0 は初めから(levels8 は -inf を入れる)
    m = np.isfinite(LT)
    if not m.any():
        return np.full(len(lt10), np.nan)
    # 最後に越えた準位より後は、包絡がその準位にとどまる(準位の刻み 0.01 の分だけ粗い)。比べるのは有効な年齢までに限る
    b = np.interp(lt10, LT[m], R8.P8[m], left=0.0, right=R8.P8[m][-1])
    return np.asarray(R8.chi_of_w(R8.w_of_b(b, p, sig), cs), float)


def _scaneval_job(args):
    import reftab8 as R8
    warnings.simplefilter("ignore")
    fn, tabpath = args
    key = tuple(int(x) for x in fn[4:-4].split("_"))
    p, ia, js, kf, il = key
    zf = np.load(os.path.join(SCAN10F, fn))
    cf, ff, formed = _chi_frac(zf, 10)
    tau = zf["t"][:len(cf)] / sr.YR * float(zf["mu"])
    a, sig, f, M, cs, chi0 = R8.node_params(*key)
    # 10 準位の曲線から表と同じ形の準位を作り、表の節と同じ物差し(包絡の chi、時刻 2% を許す)で比べる
    LTf, _, _ = R8.levels8(tau, cf, p, sig, cs)
    with np.load(tabpath) as z:
        tg = R8.PIECE_TAG[p]
        LT8 = z[f"LT_{tg}"][key[1:]].astype(float)
        st8, tv8 = int(z[f"status_{tg}"][key[1:]]), float(z[f"tau_valid_{tg}"][key[1:]])
    lt10 = np.log10(tau)
    ok = tau <= (tv8 if st8 == 1 else tau[-1]) * (1 + 1e-12)
    c8e = _table_chi(LT8, p, sig, cs, lt10)
    cfe = _table_chi(LTf, p, sig, cs, lt10)
    good = ok & np.isfinite(c8e) & np.isfinite(cfe)
    d = np.full(len(tau), np.nan)
    if good.sum() > 2:
        d[good] = kinn._tol_dist(tau[good], c8e[good], cfe[good], 2e-2)
    k = int(np.nanargmax(d)) if np.isfinite(d).any() else 0
    out = dict(key=list(key), alpha=a, sigma=sig, lf=float(np.log10(f)), M=M, chi0=chi0, status8=st8, status10f=int(zf["status"]),
               n10f=len(cf), max_dchi=float(np.nanmax(d)) if np.isfinite(d).any() else None, tau_at_max=float(tau[k]),
               chi_end_table=float(c8e[good][-1]) if good.any() else None, chi_end10f=float(cf[-1]), secs=float(zf["secs"]))
    for v, X in zip((8, 9), X_LEVELS):
        kk = int(np.argmax(ff[v]))
        out[f"{name(X)}_max_log10_frac"] = float(ff[v][kk]) if formed.any() else None
        out[f"{name(X)}_tau_at_max"] = float(tau[kk])
        above = np.flatnonzero(ff[v] > -2)
        out[f"{name(X)}_tau_1pct"] = float(tau[above[0]]) if len(above) else None   # 雲の 1% を越えた最初の年齢
    if formed[len(cf) - 1]:
        out["dominant_end"] = [name(L10[i]) for i in np.argsort(-ff[:, len(cf) - 1])[:3]]
    return out


def part_scaneval(nproc=13, tabpath=os.path.join(OUT, "reftab8.npz")):
    from multiprocessing import Pool
    fns = sorted(f for f in os.listdir(SCAN10F) if f.endswith(".npz"))
    with Pool(nproc) as p:
        rows = p.map(_scaneval_job, [(f, tabpath) for f in fns], chunksize=1)
    grow = lambda r, X: r[f"{name(X)}_max_log10_frac"] is not None and r[f"{name(X)}_max_log10_frac"] > -2   # noqa: E731
    summ = dict(n=len(rows), n_522_above_1pct=sum(grow(r, (5, 2, 2)) for r in rows),
                n_311_above_1pct=sum(grow(r, (3, 1, 1)) for r in rows),
                n_dchi_above_0_03=sum((r["max_dchi"] or 0) > 0.03 for r in rows),
                n_dchi_above_0_03_without_X=sum((r["max_dchi"] or 0) > 0.03 and not grow(r, (5, 2, 2)) and not grow(r, (3, 1, 1))
                                                for r in rows))
    save_part("scaneval", dict(summary=summ, rows=rows))
    print(json.dumps(summ, indent=1), flush=True)
    return rows


# ---------- 本番の表と同じ解き方の 10 準位(Hoffman2 のジョブ配列) ----------
# 表の節と同じ解き方(evolveN_checked、上限 1e7、年齢 TAU_GRID を TAU_MAX まで)で system10full を解き、表と同じ形の準位
# (levels8)を作る。表の同じ節と、包絡の chi を時刻 2% を許して比べられる(chk10eval)。
def chk10_nodes():
    """scan_nodes に、run_nodes の区分 0・-1・-2 の M = 1 Msun の節を足したもの。"""
    return scan_nodes() + [k for k in run_nodes() if k[0] != 1]


def _chk10_point(args):
    import reftab8 as R8
    warnings.simplefilter("ignore")
    key, cap = args
    p, ia, js, kf, il = key
    a, sig, f, M, cs, chi0 = R8.node_params(*key)
    mu = a / sr.alpha_of(M, 1.0)
    t0 = time.time()
    NT = len(R8.TAU_GRID)
    out = dict(key=tuple(key), LT=np.full(len(R8.P8), np.nan, np.float32), RZ=np.full(len(R8.G_LT), np.nan, np.float32),
               chi=np.full(NT, np.nan, np.float32), frac=np.full((2, NT), np.nan, np.float32), status=2, n_valid=0,
               tau_valid=0.0, reason="", nfev=[0, 0], njump=[0, 0], secs=0.0)
    try:
        b = kinn.evolveN_checked(system10full(), M, chi0, mu, f, R8.TAU_MAX / mu, R8.TAU_GRID / mu, max_nfev=cap)
        nv = b.n_valid
        chi = np.asarray(kinn.chi_of(b), float)
        LT, RZ, _ = R8.levels8(R8.TAU_GRID[:nv], chi[:nv], p, sig, cs)
        N = len(L10)
        lnn = np.asarray(b.y[:N], float)
        tot = np.log(np.sum(np.exp(np.minimum(lnn, 30.0)), axis=0))
        formed = tot > np.log(1.0 / sr.N0) + np.log(1e6)
        fr = np.where(formed[None, :], (lnn[len(L8):] - tot[None, :]) / np.log(10), -np.inf)
        runs = list(b.runs) + [dict(nfev=0, njump=0)] * (2 - len(b.runs))
        out.update(LT=LT.astype(np.float32), RZ=RZ.astype(np.float32), status=0 if nv == NT else 1, n_valid=int(nv),
                   tau_valid=float(R8.TAU_GRID[nv - 1]) if nv else 0.0, reason=str(b.reason),
                   nfev=[int(r["nfev"]) for r in runs[:2]], njump=[int(r["njump"]) for r in runs[:2]])
        out["chi"][:len(chi)] = chi
        out["frac"][:, :fr.shape[1]] = fr
    except Exception as e:                                 # noqa: BLE001
        out["reason"] = ("error: " + repr(e))[:200]
    out["secs"] = time.time() - t0
    return out


def _chk10_pack(res, cap):
    import platform
    return dict(key=np.array([r["key"] for r in res], int).reshape(-1, 5), LT=np.array([r["LT"] for r in res], np.float32),
                RZ=np.array([r["RZ"] for r in res], np.float32), chi=np.array([r["chi"] for r in res], np.float32),
                frac=np.array([r["frac"] for r in res], np.float32), status=np.array([r["status"] for r in res], np.int8),
                n_valid=np.array([r["n_valid"] for r in res], int), tau_valid=np.array([r["tau_valid"] for r in res]),
                reason=np.array([r["reason"] for r in res]), nfev=np.array([r["nfev"] for r in res], np.int64).reshape(-1, 2),
                njump=np.array([r["njump"] for r in res], np.int64).reshape(-1, 2), secs=np.array([r["secs"] for r in res]),
                cap=cap, host=platform.node())


def _chk10_unpack(path):
    with np.load(path) as z:
        d = {k: z[k] for k in z.files}
    keys = ("LT", "RZ", "chi", "frac", "status", "n_valid", "tau_valid", "reason", "nfev", "njump", "secs")
    return {tuple(int(x) for x in k): dict({kk: d[kk][i] for kk in keys}, key=tuple(int(x) for x in k))
            for i, k in enumerate(d["key"])}


def chk10_task(tid, ntask, out, cap, nproc=1, every=10, every_s=900.0):
    """ジョブ配列の 1 本(reftab8.build_task と同じ形: 受け持ちを nproc 並列で解き、out/chk10_<tid>.npz に保存し直す。
    続きから解ける。終了の合図で保存してやめる)。節の一覧は chk10_nodes を重い節が散るように混ぜたもの。"""
    import multiprocessing as mp
    import signal
    from multiprocessing import Pool
    nodes = chk10_nodes()
    order = np.random.default_rng(24).permutation(len(nodes))
    mine = [nodes[i] for i in order[tid - 1::ntask]]
    path = os.path.join(out, f"chk10_{tid:05d}.npz")
    have = _chk10_unpack(path) if os.path.exists(path) else {}
    todo = [k for k in mine if tuple(k) not in have]
    print(f"task {tid}/{ntask}: {len(mine)} nodes, {len(have)} already saved, {len(todo)} to solve, {nproc} workers, cap {cap}",
          flush=True)

    def save():
        tmp = path + ".tmp.npz"
        np.savez(tmp, **_chk10_pack(list(have.values()), cap))
        os.replace(tmp, path)

    term = []
    signal.signal(signal.SIGTERM, lambda signum, frame: term.append(signum))
    if todo:
        last, n0, ndone = time.time(), len(have), 0
        pool = Pool(nproc, initializer=signal.signal, initargs=(signal.SIGTERM, signal.SIG_DFL))
        it = pool.imap_unordered(_chk10_point, [(k, cap) for k in todo], chunksize=1)
        while ndone < len(todo) and not term:
            try:
                r = it.next(timeout=5.0)
            except mp.TimeoutError:
                continue
            ndone += 1
            have[r["key"]] = r
            if (len(have) - n0) % every == 0 or time.time() - last > every_s:
                save()
                last = time.time()
        save()
        pool.terminate()
        pool.join()
        if term:
            print(f"task {tid}: SIGTERM, saved {len(have)}/{len(mine)} nodes", flush=True)
            return
    print(f"task {tid}: {len(have)}/{len(mine)} nodes saved", flush=True)


def part_chk10eval(dirs, tabpath=os.path.join(OUT, "reftab8.npz"), margin=10.0):
    """chk10 の結果を表の同じ節と比べる: 包絡の chi(どちらも levels8 の準位から)を、両方が有効で表を引く範囲
    (tau <= TAU_MAX / margin)の年齢で、時刻 2% を許して比べる。311・522 の割合の最大と、雲の 1% を越えた年齢も。"""
    import glob
    import reftab8 as R8
    res = {}
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, "chk10_[0-9]*.npz"))):
            if not f.endswith(".tmp.npz"):
                res.update(_chk10_unpack(f))
    z = np.load(tabpath)
    tg = R8.TAU_GRID
    lt = np.log10(tg)
    rows = []
    for key, r in sorted(res.items()):
        p = key[0]
        a, sig, f, M, cs, chi0 = R8.node_params(*key)
        T = R8.PIECE_TAG[p]
        LT8 = z[f"LT_{T}"][key[1:]].astype(float)
        st8, tv8 = int(z[f"status_{T}"][key[1:]]), float(z[f"tau_valid_{T}"][key[1:]])
        tv10 = float(r["tau_valid"]) if int(r["status"]) == 1 else (tg[-1] if int(r["status"]) == 0 else 0.0)
        tmax = min(tv8 if st8 == 1 else tg[-1], tv10, R8.TAU_MAX / margin)
        m = tg <= tmax * (1 + 1e-12)
        row = dict(key=list(key), alpha=a, sigma=sig, lf=float(np.log10(f)), M=M, chi0=chi0, status8=st8, status10=int(r["status"]),
                   tau_valid10=float(r["tau_valid"]), tau_compared_to=float(tmax), reason10=str(r["reason"]),
                   secs10=float(r["secs"]), nfev10=[int(x) for x in r["nfev"]])
        if m.sum() > 2 and int(r["status"]) >= 0 and int(r["status"]) != 2:
            c8 = _table_chi(LT8, p, sig, cs, lt[m])
            c10 = _table_chi(r["LT"].astype(float), p, sig, cs, lt[m])
            ok = np.isfinite(c8) & np.isfinite(c10)
            if ok.sum() > 2:
                dd = kinn._tol_dist(tg[m][ok], c8[ok], c10[ok], 2e-2)
                k = int(np.argmax(dd))
                row.update(max_dchi=float(dd[k]), tau_at_max=float(tg[m][ok][k]), chi8_at_max=float(c8[ok][k]),
                           chi10_at_max=float(c10[ok][k]))
        for v, X in enumerate(X_LEVELS):
            fr = r["frac"][v]
            mm = np.isfinite(fr) & (tg <= tmax * (1 + 1e-12))
            row[f"{name(X)}_max_log10_frac"] = float(fr[mm].max()) if mm.any() else None
            above = np.flatnonzero(mm & (fr > -2))
            row[f"{name(X)}_tau_1pct"] = float(tg[above[0]]) if len(above) else None
        rows.append(row)
    grow = lambda r, X: (r[f"{name(X)}_max_log10_frac"] or -99) > -2   # noqa: E731
    dmax = lambda r: r.get("max_dchi", 0.0)   # noqa: E731
    summ = dict(n=len(rows), status10={str(s): sum(r["status10"] == s for r in rows) for s in (0, 1, 2)},
                n_compared=sum("max_dchi" in r for r in rows),
                n_522_above_1pct=sum(grow(r, (5, 2, 2)) for r in rows), n_311_above_1pct=sum(grow(r, (3, 1, 1)) for r in rows),
                n_dchi_above_0_03=sum(dmax(r) > 0.03 for r in rows), n_dchi_above_0_01=sum(dmax(r) > 0.01 for r in rows),
                n_dchi_above_0_03_with_522=sum(dmax(r) > 0.03 and grow(r, (5, 2, 2)) for r in rows),
                max_dchi=max(dmax(r) for r in rows), core_hours10=sum(r["secs10"] for r in rows) / 3600,
                worst=sorted(rows, key=lambda r: -dmax(r))[:12])
    save_part("chk10eval", dict(summary=summ, rows=rows))
    print(json.dumps({k: v for k, v in summ.items() if k != "worst"}, indent=1), flush=True)
    return rows


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "chk10task":                 # Hoffman2: e24_outside.py chk10task <番号> <本数> <場所>
        chk10_task(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], cap=int(os.environ.get("SRLVK_CAP", "10000000")),
                   nproc=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
        sys.exit(0)
    npr = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    extra = sys.argv[3:]
    if what == "rates":
        part_rates(npr)
    elif what == "ratecheck":
        part_ratecheck(npr)
    elif what == "procs":
        part_procs(npr)
    elif what == "corr":
        xs = [tuple(int(c) for c in extra[0])] if extra else None
        part_corr(npr, xs=xs, limit=int(extra[1]) if len(extra) > 1 else None)
    elif what == "runs":
        part_runs(npr, limit=int(extra[0]) if extra else None)
    elif what == "rootcheck":
        part_rootcheck()
    elif what == "validate":
        part_validate(npr)
    elif what == "eval":
        part_eval(npr)
    elif what == "runs10":
        part_runs10(npr, limit=int(extra[0]) if extra else None)
    elif what == "compare":
        part_compare(npr)
    elif what == "procsxx":
        part_procsxx(npr)
    elif what == "corrxx":
        part_corrxx(npr, limit=int(extra[0]) if extra else None)
    elif what == "runs10f":
        part_runs10f(npr)
    elif what == "compare10f":
        part_compare10f(npr)
    elif what == "scan":
        part_runs10f(npr, nodes=scan_nodes(), outdir=SCAN10F, tag="scan10f")
    elif what == "scaneval":
        part_scaneval(npr)
    elif what == "chk10eval":               # e24_outside.py chk10eval 1 <場所> [...]
        part_chk10eval(extra)
    else:
        raise SystemExit(f"unknown part {what}")
