"""N 準位の基準の運動論(値の表の作り直し、DESIGN.md D18 の手順 1)。

準位と過程をデータで渡す。考え方と式は fp3.evolve3m(E18)と同じで、それを N 準位に一般化したもの:
  変数は ln n_i(n = N/N0、N0 = G Msun^2/(hbar c))、BH の質量 M [Msun]、角運動量 J [G Msun^2/c]。
  超放射: dn_i = g_i mu (n_i + n_q) dt(g_i > 0)、減衰する準位(g_i <= 0)は g_i (n_i - 1e-2 n_q)(evolve3m と同じ床)。
          BH は dM = -(omega_i/mu) alpha_sun dn_i、dJ = -m_i dn_i。
  自己相互作用 p = (i, j, k): i x j -> k x X。eps = n/M^2 で
          R_p = kappa_p alpha^q_p [r+] C_p(alpha, chi) (M_pl/f)^4 eps_i eps_j (eps_k + eps_q)、
          q = 11(BH が吸収、r+ を掛ける)、7(BH が吸収、n00 に共鳴、r+ を掛ける)、8(無限遠への放出)。
          1 回で i と j を 1 個ずつ(i = j なら 2 個)減らし、k を 1 個増やす。BH が吸収する過程では、BH がエネルギー
          omega_i + omega_j - omega_k を受け取る(l = 0 の強制振動なので m = 0、角運動量は変わらない)。
          C_p は相対論的な補正の係数(既定は 1。手順 3 で入れる)。
  重力波の過程(BGLS Table IV): evolve3m と同じく 211・322 の消滅と 322 -> 211 の遷移だけ(準位にあれば)。
過程の係数 kappa は hyd.py(水素型、主要項)で総当たりして計算し、out/kinn_procs.json に置く(process_table)。
確かめ(python3 calc/kinn.py check): 三準位で evolve3m と同じ率を渡すと evolve3m を再現する、BH の質量を固定すると
evolve5・evolve4 と合う、無限遠への過程と重力波を切ると角運動量とエネルギーが保存する。
"""
import json
import os
import sys
import time
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fp3  # noqa: E402
import si  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PROC_CACHE = os.path.join(HERE, "out", "kinn_procs.json")
LEVELS8 = [(2, 1, 1), (3, 2, 2), (4, 1, 1), (4, 2, 2), (4, 3, 3), (5, 1, 1), (5, 3, 3), (5, 4, 4)]
POW = {"BH": 11, "BHres": 7, "inf": 8}


def name(lev):
    return "%d%d%d" % tuple(lev)


def _proc_job(ijk):
    import hyd
    i, j, k = ijk
    if hyd.energy(i, j, k) > 0 and hyd.resonant_n(i, j, k) is not None:
        r = hyd.kappa_res(i, j, k)
        return ijk, float(r["kappa"]), "BHres", int(r["n"])
    kap, kind = hyd.kappa(i, j, k)
    return ijk, float(kap), kind, None


def process_table(levels=LEVELS8, nproc=10, rebuild=False):
    """準位のあいだの、主要項で 0 でない過程をすべて返す: dict(i, j, k(levels の番号), kappa, kind, res_n)。
    i <= j、k は i, j と違う。値は out/kinn_procs.json に名前で保存し、あれば読む。"""
    cache = {}
    if os.path.exists(PROC_CACHE) and not rebuild:
        cache = json.load(open(PROC_CACHE))
    keys = [(a, b, c) for a in range(len(levels)) for b in range(a, len(levels)) for c in range(len(levels))
            if c not in (a, b)]
    todo = [(levels[a], levels[b], levels[c]) for a, b, c in keys
            if f"{name(levels[a])}x{name(levels[b])}->{name(levels[c])}" not in cache]
    if todo:
        from multiprocessing import Pool
        with Pool(min(nproc, len(todo))) as p:
            for (i, j, k), kap, kind, n in p.map(_proc_job, todo, chunksize=1):
                cache[f"{name(i)}x{name(j)}->{name(k)}"] = {"kappa": kap, "kind": kind, "res_n": n}
        os.makedirs(os.path.dirname(PROC_CACHE), exist_ok=True)
        json.dump(cache, open(PROC_CACHE, "w"), indent=0, sort_keys=True)
    out = []
    for a, b, c in keys:
        v = cache[f"{name(levels[a])}x{name(levels[b])}->{name(levels[c])}"]
        if v["kind"] is not None and v["kappa"] != 0:
            out.append(dict(i=a, j=b, k=c, kappa=v["kappa"], kind=v["kind"], res_n=v["res_n"]))
    return out


def procs_from_fp3(table, levels):
    """fp3 の (i, j, k, kappa, BH か) の表(PROC5、PROC4)を、process_table と同じ形にする(共鳴する過程は入っていない)。"""
    return [dict(i=i, j=j, k=k, kappa=kap, kind="BH" if bh else "inf", res_n=None) for i, j, k, kap, bh in table]


def procs_evolve3m():
    """evolve3m が使う 6 個の過程と係数(fp3.coeffs の K_A--K_F)。準位は 211・322・544 の順。"""
    return [dict(i=0, j=0, k=1, kappa=fp3.K_A, kind="BH", res_n=None), dict(i=1, j=1, k=0, kappa=fp3.K_B, kind="inf", res_n=None),
            dict(i=1, j=1, k=2, kappa=fp3.K_C, kind="BH", res_n=None), dict(i=2, j=2, k=1, kappa=fp3.K_D, kind="inf", res_n=None),
            dict(i=1, j=2, k=0, kappa=fp3.K_E, kind="inf", res_n=None), dict(i=2, j=2, k=0, kappa=fp3.K_F, kind="inf", res_n=None)]


def sr_default(alpha, chi, levels):
    """超放射の率 / mu(占有数の率)。evolve3m と同じ関数(表の範囲は表、外は解析式)。手順 2 で多倍長の表に替える。"""
    return np.array([fp3.g_rel(alpha, chi, *lev) for lev in levels])


class System:
    """準位と過程から、右辺に要る配列を前もって作る。"""

    def __init__(self, levels, procs, sr_rate=sr_default, corr=None):
        self.levels = [tuple(l) for l in levels]
        self.nl = len(levels)
        self.ms = np.array([l[2] for l in levels], float)
        self.procs = procs
        P = len(procs)
        self.I = np.array([p["i"] for p in procs], int)
        self.J = np.array([p["j"] for p in procs], int)
        self.K = np.array([p["k"] for p in procs], int)
        self.kap = np.array([p["kappa"] for p in procs], float)
        self.pw = np.array([POW[p["kind"]] for p in procs], float)
        self.bh = np.array([p["kind"] != "inf" for p in procs])
        # 化学量論の行列 S[level, process]: i と j を 1 個ずつ(i = j なら 2 個)減らし、k を 1 個増やす
        S = np.zeros((self.nl, P))
        for q, p in enumerate(procs):
            S[p["i"], q] -= 1
            S[p["j"], q] -= 1
            S[p["k"], q] += 1
        self.S = S
        self.sr_rate = sr_rate
        self.corr = corr                                   # corr(alpha, chi) -> 長さ P の配列、または None
        idx = {l: n for n, l in enumerate(self.levels)}
        self.i211, self.i322 = idx.get((2, 1, 1)), idx.get((3, 2, 2))

    def coef(self, alpha, chi):
        c = self.kap * alpha ** self.pw * np.where(self.bh, sr.rplus(chi), 1.0)
        if self.corr is not None:
            c = c * self.corr(alpha, chi)
        return c


def _continue(fun, y0, te, T_end, nvar, rtol, atol, max_nfev):
    """evolve3m と同じ続け方: LSODA で止まったら最後の出力の時刻から Radau、BDF と切り替える(1 区間の評価は
    max_nfev/4 まで、全体で max_nfev まで)。"""
    from scipy.integrate import solve_ivp
    seg_nfev = max(max_nfev // 4, 1000)
    count = [0]
    t_stop = [np.inf]

    def budget(t, y):
        return 1.0 if not np.isfinite(t_stop[0]) else t_stop[0] - t
    budget.terminal = True

    def fc(t, y):
        count[0] += 1
        if count[0] == seg_nfev:
            t_stop[0] = t
        return fun(t, y)

    methods = ("LSODA", "Radau", "BDF")
    t_cur, y_cur = 0.0, np.asarray(y0, float)
    ts, ys = [], []
    done, total, k, stall = 0, 0, 0, 0
    while done < len(te) and total < max_nfev and stall < len(methods):
        count[0] = 0
        t_stop[0] = np.inf
        sol = solve_ivp(fc, (t_cur, T_end), y_cur, method=methods[k % len(methods)], rtol=rtol, atol=atol,
                        t_eval=te[done:], events=budget)
        total += count[0]
        sy = np.asarray(sol.y, float)
        if sy.ndim != 2 or sy.shape[1] == 0:
            sy = np.zeros((nvar, 0))
        good = np.all(np.isfinite(sy), axis=0)
        nnew = int(np.argmin(good)) if not np.all(good) else len(good)
        if nnew > 0:
            ts.append(np.asarray(sol.t)[:nnew])
            ys.append(sy[:, :nnew])
            done += nnew
            t_cur, y_cur = ts[-1][-1], sy[:, nnew - 1]
            stall = 0
        else:
            stall += 1
        if done < len(te):
            k += 1
    ok = done == len(te)
    y = np.concatenate(ys, axis=1) if ys else np.zeros((nvar, 0))
    t = np.concatenate(ts) if ts else np.zeros(0)
    return SimpleNamespace(status=0 if ok else -1, y=y, t=t, nfev=total, nseg=k + 1,
                           message="ok" if ok else f"stopped at t = {t_cur / sr.YR:.3e} yr")


def evolveN(sysN, M0, chi0, mu_eV, f_gev, T_yr, t_eval_yr, rtol=1e-8, gw=True, max_nfev=1200000, fixed_mass=False,
            inf_off=False):
    """N 準位の運動論を 1 個の量子から解く。返り値: y = (ln n_1..ln n_N, M, J)、t [s]、status(0 = 最後まで解けた)。
    fixed_mass: BH の質量を M0 に固定する(evolve5・evolve4 との照合用。J は変える)。inf_off: 無限遠への放出を切る
    (保存則の確かめ用)。"""
    mu_s = mu_eV / sr.HBAR_EVS
    a_sun = sr.alpha_of(1.0, mu_eV)
    s4 = (si.MPL_GEV / f_gev) ** 4
    nq = 1.0 / sr.N0
    N = sysN.nl
    lev = sysN.levels
    I, J_, K, S, bh = sysN.I, sysN.J, sysN.K, sysN.S, sysN.bh
    i1, i2 = sysN.i211, sysN.i322
    keep = sysN.bh.astype(float) if inf_off else np.ones(len(sysN.procs))

    def f(t, y):
        n = np.exp(np.clip(y[:N], -800.0, 30.0))
        M = M0 if fixed_mass else min(max(y[N], 0.5 * M0), M0)
        Jb = max(y[N + 1], 0.0)
        chi = min(max(Jb / M ** 2, 0.0), 0.9999)
        a = a_sun * M
        g = sysN.sr_rate(a, chi, lev)
        w = np.array([sr.omega_over_mu(a, nn, l, m, chi) for (nn, l, m) in lev])
        e = n / M ** 2
        eq = nq / M ** 2
        R = sysN.coef(a, chi) * s4 * e[I] * e[J_] * (e[K] + eq) * keep
        sr_ = np.where(g > 0, g * (n + nq), g * (n - 1e-2 * nq))
        M2 = M ** 2
        d = sr_ + M2 * (S @ R)
        if gw:
            if i1 is not None:
                d[i1] -= M2 * 2 * 1e-2 * a ** 14 * e[i1] ** 2
            if i2 is not None:
                d[i2] -= M2 * 2 * 3e-8 * a ** 18 * e[i2] ** 2
            if i1 is not None and i2 is not None:
                Gt = 5e-6 * a ** 10 * e[i2] * (e[i1] + eq)
                d[i1] += M2 * Gt
                d[i2] -= M2 * Gt
        dM = 0.0 if fixed_mass else (-a_sun * np.dot(w, sr_) + a_sun * M2 * np.sum((R * bh) * (w[I] + w[J_] - w[K])))
        dJ = -np.dot(sysN.ms, sr_)
        out = mu_s * np.concatenate([d / n, [dM, dJ]])
        return np.where(np.isfinite(out), out, 0.0)

    y0 = [np.log(nq)] * N + [M0, chi0 * M0 ** 2]
    te = np.clip(np.asarray(t_eval_yr, float) * sr.YR, 0.0, T_yr * sr.YR)
    return _continue(f, y0, te, T_yr * sr.YR, N + 2, rtol, 1e-12, max_nfev)


def chi_of(sol, M0=None):
    """解からスピン J/M^2。"""
    N = sol.y.shape[0] - 2
    return sol.y[N + 1] / sol.y[N] ** 2


# ---------- 確かめ ----------
def check_evolve3m():
    """三準位・evolve3m と同じ係数で、evolve3m を再現するか(12 例 x 60 年齢)。"""
    sys3 = System([(2, 1, 1), (3, 2, 2), (5, 4, 4)], procs_evolve3m())
    cases = [(19.6, 0.78, 0.71e-12, 1e12), (19.6, 0.78, 0.71e-12, 1e16), (19.6, 0.78, 0.71e-12, 1e18),
             (45.0, 0.70, 0.71e-12, 1e14), (30.0, 0.50, 0.71e-12, 1e15), (8.0, 0.95, 1e-12, 1e13),
             (60.0, 0.29, 0.71e-12, 5e18), (25.0, 0.86, 1.0e-12, 3e17), (12.0, 0.6, 3e-12, 1e11),
             (5.0, 0.9, 5e-12, 1e20), (100.0, 0.8, 0.3e-12, 1e16), (15.0, 0.45, 2e-12, 1e19)]
    Ts = np.logspace(3, 10, 60)
    worst, rows = 0.0, []
    for M, c0, mu, f in cases:
        t0 = time.time()
        a = fp3.evolve3m(M, c0, mu, f, Ts[-1], Ts)
        t1 = time.time()
        b = evolveN(sys3, M, c0, mu, f, Ts[-1], Ts)
        t2 = time.time()
        ca, cb = a.y[4] / a.y[3] ** 2, chi_of(b)
        n = min(len(ca), len(cb))
        d = float(np.max(np.abs(ca[:n] - cb[:n]))) if n else np.nan
        worst = max(worst, d)
        rows.append(dict(M=M, chi0=c0, mu=mu, f=f, max_diff=d, n_ages=[len(ca), len(cb)], secs=[t1 - t0, t2 - t1]))
    return dict(worst=worst, rows=rows)


def check_fixed_mass():
    """BH の質量を固定して、evolve5(211・322・411・433・544)と evolve4(211・322・422・544)と比べる。"""
    out = {}
    for nm, levels, table, fn in (("evolve5", fp3.LEVELS5, fp3.PROC5, fp3.evolve5), ("evolve4", fp3.LEVELS4, fp3.PROC4, fp3.evolve4)):
        if nm == "evolve5":
            rate = lambda a, c, lv: fp3.sr_rates5(a, c)                       # noqa: E731
        else:                                                                    # evolve4: 422 は RHO_422 x 322 の率
            rate = lambda a, c, lv: np.array([fp3.g_rel(a, c, 2, 1, 1), fp3.g_rel(a, c, 3, 2, 2),   # noqa: E731
                                              fp3.RHO_422 * fp3.g_rel(a, c, 3, 2, 2), fp3.g_rel(a, c, 5, 4, 4)])
        sysN = System(levels, procs_from_fp3(table, levels), sr_rate=rate)
        worst = 0.0
        for M, c0, mu, f in [(30.0, 0.9, 0.71e-12, 1e15), (40.0, 0.85, 0.71e-12, 1e16), (25.0, 0.95, 1.0e-12, 1e15)]:
            Ts = np.logspace(4, 10, 40)
            ref = fn(M, c0, mu, f, Ts[-1], Ts)
            cref = ref.y[len(levels)]
            b = evolveN(sysN, M, c0, mu, f, Ts[-1], Ts, fixed_mass=True, gw=False)
            cb = chi_of(b)
            n = min(len(cref), len(cb))
            worst = max(worst, float(np.max(np.abs(cref[:n] - cb[:n]))))
        out[nm] = worst
    return out


def check_conservation():
    """8 準位、無限遠への過程と重力波を切ると、BH と雲の角運動量とエネルギーの和が保存するか。"""
    procs = process_table()
    sys8 = System(LEVELS8, procs)
    res = []
    for M, c0, mu, f in [(30.0, 0.9, 0.71e-12, 1e15), (20.0, 0.8, 1.0e-12, 1e16)]:
        Ts = np.logspace(4, 10, 30)
        b = evolveN(sys8, M, c0, mu, f, Ts[-1], Ts, gw=False, inf_off=True)
        N = len(LEVELS8)
        n = np.exp(b.y[:N])
        Mb, Jb = b.y[N], b.y[N + 1]
        a_sun = sr.alpha_of(1.0, mu)
        Jtot = Jb + np.array([lv[2] for lv in LEVELS8]) @ n
        w = np.array([[sr.omega_over_mu(a_sun * Mb[t], *lv, Jb[t] / Mb[t] ** 2) for t in range(len(Mb))] for lv in LEVELS8])
        Etot = Mb + a_sun * np.sum(w * n, axis=0)
        res.append(dict(M=M, chi0=c0, f=f, dJ_rel=float(np.max(np.abs(Jtot / Jtot[0] - 1))),
                        dE_rel=float(np.max(np.abs(Etot / Etot[0] - 1))), chi_end=float(Jb[-1] / Mb[-1] ** 2)))
    return res


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "check"
    if what == "procs":
        p = process_table(rebuild="--rebuild" in sys.argv)
        print(len(p), "processes")
    else:
        rep = {"evolve3m": check_evolve3m(), "fixed_mass": check_fixed_mass(), "conservation": check_conservation()}
        print(json.dumps(rep, indent=1))
