"""N 準位の基準の運動論(値の表の作り直し、DESIGN.md D18 の手順 1)。

準位と過程をデータで渡す。考え方と式は fp3.evolve3m(E18)と同じで、それを N 準位に一般化したもの:
  変数は ln n_i(n = N/N0、N0 = G Msun^2/(hbar c))、BH の質量 M [Msun]、角運動量 J [G Msun^2/c]。
  超放射: dn_i = g_i mu (n_i + n_q) dt(g_i > 0)、減衰する準位(g_i <= 0)は g_i (n_i - 1e-2 n_q)(evolve3m と同じ床)。
          BH は dM = -(omega_i/mu) alpha_sun dn_i、dJ = -m_i dn_i。
  自己相互作用 p = (i, j, k): i x j -> k x X。eps = n/M^2 で
          R_p = kappa_p alpha^q_p [r+] C_p(alpha, chi) (M_pl/f)^4 eps_i eps_j (eps_k + eps_q)、
          q = 11(BH が吸収、r+ を掛ける)、7(BH が吸収、n00 に共鳴、r+ を掛ける)、8(無限遠への放出)。
          1 回で i と j を 1 個ずつ(i = j なら 2 個)減らし、k を 1 個増やす。
          C_p = C_inf + C_H は相対論的な補正の係数(relcorr.Corr。既定は補正なし: 放出の過程は C_inf = 1、BH が吸収する
          過程は C_H = 1)。C_inf の分は X が無限遠へ出る。C_H の分は X が地平線に吸われ、BH がエネルギー
          omega_X = omega_i + omega_j - omega_k と角運動量 m_X = m_i + m_j - m_k を受け取る(放出の過程の超放射の m_X では
          C_H < 0 で、BH が失う)。BH が吸収する過程はすべて m_X = 0。
          BGLS の式 A14 の摂動の条件(作られる準位の 1 個あたりの率 x X が雲にとどまる時間 tau_X << 1)は a14 で確かめる。
  重力波の過程(BGLS Table IV): evolve3m と同じく 211・322 の消滅と 322 -> 211 の遷移だけ(準位にあれば)。
過程の係数 kappa は hyd.py(水素型、主要項)で総当たりして計算し、out/kinn_procs.json に置く(process_table)。
確かめ(python3 calc/kinn.py check): 三準位で evolve3m と同じ率を渡すと evolve3m を再現する、BH の質量を固定すると
evolve5・evolve4 と合う、無限遠への過程と重力波を切ると角運動量とエネルギーが保存する。
補正の確かめ(python3 calc/kinn.py corrcheck、2026-11-16): 恒等の補正で corr なしと同じ解(差 0)、右辺の保存則
3e-16(地平線の経路)、M0 = 30・chi0 = 0.95 の alpha0 0.03-0.45 x f 1e13-1e17 の 55 例で A14 の x は 4e-4 以下。
雲が大きく振動する(バースト、周期 0.02-5e4 年)区間では解法が止まる(右辺の形や積分法では直らない)ので、振動を平均して
解く(evolveN_avg、DESIGN.md D20)。表には evolveN_checked を使う: 計算に上限を置き、平均を使った点は細かい設定でもう一度
解いて比べ、両方が届き差が 5e-3 以内の年齢までを有効とする(確かめ: python3 calc/kinn.py avgcheck [--refs])。
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
        self.mX = self.ms[self.I] + self.ms[self.J] - self.ms[self.K]
        self.names = [f"{name(self.levels[p['i']])}x{name(self.levels[p['j']])}->{name(self.levels[p['k']])}"
                      for p in procs]
        self.sr_rate = sr_rate
        # corr(alpha, chi) -> (C_inf, C_H)。corr.names があれば過程を名前で対応させる(表にない過程は誤り)
        self.corr = corr
        self.cidx = None
        if corr is not None and hasattr(corr, "names"):
            pos = {str(nm): q for q, nm in enumerate(corr.names)}
            miss = [nm for nm in self.names if nm not in pos]
            if miss:
                raise KeyError(f"補正の表にない過程: {miss}")
            self.cidx = np.array([pos[nm] for nm in self.names], int)
        idx = {l: n for n, l in enumerate(self.levels)}
        self.i211, self.i322 = idx.get((2, 1, 1)), idx.get((3, 2, 2))

    def coef(self, alpha, chi):
        """(無限遠の経路の係数, 地平線の経路の係数)。R_p = 係数 x (M_pl/f)^4 eps_i eps_j (eps_k + eps_q)。"""
        c = self.kap * alpha ** self.pw * np.where(self.bh, sr.rplus(chi), 1.0)
        if self.corr is None:
            return np.where(self.bh, 0.0, c), np.where(self.bh, c, 0.0)
        ci, ch = self.corr(alpha, chi)
        if self.cidx is not None:
            ci, ch = ci[self.cidx], ch[self.cidx]
        return c * ci, c * ch

    def tau(self, alpha, chi):
        """BH が吸収する過程と共鳴する過程の tau_X(M の単位)。それ以外と、corr に tau がなければ nan。"""
        if self.corr is None or not hasattr(self.corr, "tau"):
            return np.full(len(self.procs), np.nan)
        t = self.corr.tau(alpha, chi)
        if self.cidx is not None:
            t = t[self.cidx]
        return np.where(self.bh, t, np.nan)


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
        try:
            sol = solve_ivp(fc, (t_cur, T_end), y_cur, method=methods[k % len(methods)], rtol=rtol, atol=atol,
                            t_eval=te[done:], events=budget)
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):   # 陰解法の線形代数が溢れた: 止まったと同じに扱う
            sol = SimpleNamespace(t=np.zeros(0), y=np.zeros((nvar, 0)))
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


def make_rhs(sysN, M0, mu_eV, f_gev, gw=True, fixed_mass=False, inf_off=False):
    """evolveN の右辺 f(t, y)(y = (ln n_1..ln n_N, M, J)、返り値は d y/dt [1/s])。引数は evolveN と同じ。
    d ln n_i/dt は 1 個あたりの率で直接書く(n で割らない): 準位 i が過程 p で消える率は c_p eps_j (eps_k + eps_q)、
    k が作られる率は c_p eps_i eps_j (1 + n_q/n_k)。n が倍精度の 0 になっても(ln n < -745)微分が 0 にならない。"""
    mu_s = mu_eV / sr.HBAR_EVS
    a_sun = sr.alpha_of(1.0, mu_eV)
    s4 = (si.MPL_GEV / f_gev) ** 4
    nq = 1.0 / sr.N0
    lnq = np.log(nq)
    N = sysN.nl
    lev = sysN.levels
    I, J_, K, mX = sysN.I, sysN.J, sysN.K, sysN.mX
    P = len(I)
    Si, Sj, Sk = np.zeros((N, P)), np.zeros((N, P)), np.zeros((N, P))
    Si[I, np.arange(P)] = 1.0
    Sj[J_, np.arange(P)] = 1.0
    Sk[K, np.arange(P)] = 1.0
    i1, i2 = sysN.i211, sysN.i322
    keep_inf = 0.0 if inf_off else 1.0

    def f(t, y):
        ln = y[:N]
        n = np.exp(np.minimum(ln, 30.0))
        r = np.exp(np.minimum(lnq - ln, 300.0))            # n_q / n(大きすぎると陰解法の線形代数で溢れる)
        M = M0 if fixed_mass else min(max(y[N], 0.5 * M0), M0)
        Jb = max(y[N + 1], 0.0)
        chi = min(max(Jb / M ** 2, 0.0), 0.9999)
        a = a_sun * M
        g = sysN.sr_rate(a, chi, lev)
        w = np.array([sr.omega_over_mu(a, nn, l, m, chi) for (nn, l, m) in lev])
        M2 = M ** 2
        e = n / M2
        eq = nq / M2
        ci, ch = sysN.coef(a, chi)
        ct = (ci * keep_inf + ch) * s4
        ek = e[K] + eq
        d = np.where(g > 0, g * (1.0 + r), g * (1.0 - 1e-2 * r))
        d += -Si @ (ct * e[J_] * ek) - Sj @ (ct * e[I] * ek) + Sk @ (ct * e[I] * e[J_] * (1.0 + r[K]))
        if gw:
            if i1 is not None:
                d[i1] -= 2e-2 * a ** 14 * e[i1]
            if i2 is not None:
                d[i2] -= 6e-8 * a ** 18 * e[i2]
            if i1 is not None and i2 is not None:
                d[i1] += 5e-6 * a ** 10 * e[i2] * (1.0 + r[i1])
                d[i2] -= 5e-6 * a ** 10 * (e[i1] + eq)
        sr_ = np.where(g > 0, g * (n + nq), g * (n - 1e-2 * nq))
        Rh = ch * s4 * e[I] * e[J_] * ek                   # 地平線の経路(BH が omega_X と m_X を受け取る)
        dM = 0.0 if fixed_mass else (-a_sun * np.dot(w, sr_) + a_sun * M2 * np.dot(Rh, w[I] + w[J_] - w[K]))
        dJ = -np.dot(sysN.ms, sr_) + M2 * np.dot(Rh, mX)
        out = mu_s * np.concatenate([d, [dM, dJ]])
        return np.where(np.isfinite(out), out, 0.0)
    return f


def evolveN(sysN, M0, chi0, mu_eV, f_gev, T_yr, t_eval_yr, rtol=1e-8, gw=True, max_nfev=1200000, fixed_mass=False,
            inf_off=False):
    """N 準位の運動論を 1 個の量子から解く。返り値: y = (ln n_1..ln n_N, M, J)、t [s]、status(0 = 最後まで解けた)。
    fixed_mass: BH の質量を M0 に固定する(evolve5・evolve4 との照合用。J は変える)。inf_off: 無限遠への経路を切る
    (地平線の経路は残す。保存則の確かめ用)。"""
    f = make_rhs(sysN, M0, mu_eV, f_gev, gw=gw, fixed_mass=fixed_mass, inf_off=inf_off)
    N = sysN.nl
    nq = 1.0 / sr.N0
    y0 = [np.log(nq)] * N + [M0, chi0 * M0 ** 2]
    te = np.clip(np.asarray(t_eval_yr, float) * sr.YR, 0.0, T_yr * sr.YR)
    return _continue(f, y0, te, T_yr * sr.YR, N + 2, rtol, 1e-12, max_nfev)


def _hann_fit(t, Y):
    """Y の各行を t に Hann の重み(と台形の幅)つきで直線に当てはめる。返り値: (中心の時刻, 平均, 傾き)。
    窓の中に K 周期ある振動の漏れは (1/2 pi K)^3 ほどに抑えられる。"""
    dt = np.gradient(t)
    w = 0.5 * (1.0 - np.cos(2.0 * np.pi * (t - t[0]) / (t[-1] - t[0]))) * dt
    W = w.sum()
    tc = np.dot(w, t) / W
    mean = Y @ w / W
    x = t - tc
    slope = (Y - mean[:, None]) @ (w * x) / np.dot(w, x * x)
    return tc, mean, slope


def evolveN_avg(sysN, M0, chi0, mu_eV, f_gev, T_yr, t_eval_yr, rtol=1e-8, atol=1e-12, gw=True,
                stall_steps=3000, stall_dlnt=0.05, burst_calls=5000, burst_max=160000, relax=0.25,
                dchi=1e-3, dlnM=5e-4, dln=0.25, dln_small=5.0, jump_frac=0.25, max_nfev=50000000, burst_rtol=1e-6, trace=None):
    """evolveN と同じ運動論を、雲が速く振動する(バースト)区間は平均して解く(DESIGN.md D20)。
    8 準位の雲は、自己相互作用の準平衡の近くで大きく振動する(少ない準位で数十 e 倍、周期 0.02 年から 5e4 年)ことがあり、
    1e10 年までそのまま解くと最大 1e15 歩かかる。BH にはほとんど効かない(J の振れ 1e-8 から 2e-2)ので、振動を平均する。
    ふつうは LSODA で解き、stall_steps 歩で ln t が stall_dlnt しか進まなくなったら平均の区間に入る:
      (1) バースト: 今の状態から burst_calls 回(足りなければ倍に、burst_max まで)解く。初めの relax の割合(前の跳びの
          あとのなじみ)を除き、各変数を Hann の重みで直線に当てはめて平均と傾きを出す。前半と後半の傾きが M・J で
          3 割以内に合うまで延ばす。
      (2) 跳び: M と J、それに傾きが前半と後半で合いはっきり動く準位を、傾きで D だけ進める。D は chi の変化 dchi、
          M の相対変化 dlnM、進める準位の ln n の変化 dln(最大の準位の 1e-10 倍より少なく育つ準位は dln_small まで、ただし
          最大の e^-10 倍の手前まで。そういう準位はほかの準位に効かず、増え方は多い準位で決まる)、今の時刻の jump_frac 倍の
          うち小さい方に、(3) の係数を掛ける。
          前のバーストとの傾きの差から曲がりを出して 2 次で進める(曲がりの項は 1 次の項を越えないように抑える)。
          減っていく準位は 1 個の量子の 1e-3 倍より下へは進めない(そこから下は 1 個からの種で決まり、深さが効かない)。
          傾きが合わない準位(振動が窓で平均しきれない)は進めない(次のバーストでなじむ)。
      (3) 跳びの検査: 跳んだ先のバーストで測った M・J の傾きが、跳びに使った傾きと(揺れを除いて)3 割より違えば、
          跳びを半分にしてやり直し、係数も半分にする(雲がなじまない中間の速さの成分があると、傾きで進めるのは前進
          Euler と同じで振れが育つ。曲がりの強い所も同じ)。1 割以内なら係数を 1.5 倍(1 まで)。
      バーストが jump_frac t を楽に越えるか、D が窓より短ければ、ふつうの解き方に戻る。解法が失敗したら LSODA、
      小さい初めの歩の LSODA、Radau、BDF の順に替える。
    出力の年齢: ふつうの区間は解そのもの。平均の区間は当てはめた直線(雲の振動を平均した値)。
    返り値: evolveN と同じ y、t、status に、avg(年齢ごとの印: 0 ふつう、1 平均、2 平均で M・J の傾きが前半と後半で
    合わないまま)、chi_amp(平均した区間の chi の振れの半分)、njump、nrej(やり直した跳び)、nfev を加える。"""
    from scipy.integrate import BDF, LSODA, Radau
    fr = make_rhs(sysN, M0, mu_eV, f_gev, gw=gw)
    N = sysN.nl
    nq = 1.0 / sr.N0
    lnq = np.log(nq)
    T_end = T_yr * sr.YR
    te = np.clip(np.asarray(t_eval_yr, float) * sr.YR, 0.0, T_end)
    out = np.full((N + 2, len(te)), np.nan)
    flag = np.zeros(len(te), int)
    camp = np.zeros(len(te))
    t, y = 0.0, np.array([lnq] * N + [M0, chi0 * M0 ** 2])
    total, njump, k = 0, 0, 0
    mode = "normal"

    rtol_b = rtol if burst_rtol is None else burst_rtol

    def new_solver(t0, y0, attempt, rt=rtol):
        """失敗したら次の方法でやり直す(evolveN の _continue と同じ考え): LSODA、小さい初めの歩の LSODA、Radau、BDF。"""
        m = attempt % 4
        if m == 0:
            return LSODA(fr, t0, y0, T_end, rtol=rt, atol=atol)
        if m == 1:
            return LSODA(fr, t0, y0, T_end, rtol=rt, atol=atol, first_step=max(1e-9 * t0, 1.0))
        return (Radau if m == 2 else BDF)(fr, t0, y0, T_end, rtol=rt, atol=atol)

    def step(sol):
        try:
            sol.step()
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            sol.status = "failed"
        if sol.status != "failed" and not np.all(np.isfinite(sol.y)):
            sol.status = "failed"
        return sol.status != "failed"

    def collect(sol, buf):
        nonlocal k
        if k < len(te) and te[k] <= sol.t:
            dense = sol.dense_output()
            while k < len(te) and te[k] <= sol.t:
                buf[k] = dense(te[k])
                k += 1

    def jump(p, D):
        """バーストの当てはめ p で D だけ跳ぶ。年齢の出力は、バーストの始めから跳びの先まで当てはめた直線、その先は消す。"""
        nonlocal k
        a_ = np.clip(p["a"], -np.abs(p["s"]) / D, np.abs(p["s"]) / D)   # 曲がりの項は 1 次の項を越えない
        dyv = p["s"] * D + a_ * ((p["t_end"] - p["tc"]) * D + 0.5 * D * D)
        y_new = p["y_end"].copy()
        y_new[N:] += dyv[N:]
        for i in np.flatnonzero(p["proj"]):
            dy = dyv[i]
            y_new[i] = max(p["y_end"][i] + dy, min(p["y_end"][i], lnq - 7.0)) if dy < 0 else p["y_end"][i] + dy
        t_new = p["t_end"] + D
        kk = p["k0"]
        while kk < len(te) and te[kk] <= t_new:
            x = te[kk] - p["tc"]
            v = p["m"] + p["s"] * x + 0.5 * a_ * x * x
            v[:N] = np.where(p["proj"], v[:N], p["m"][:N])
            out[:, kk] = v
            flag[kk] = p["flag"]
            camp[kk] = p["amp"]
            kk += 1
        out[:, kk:] = np.nan
        flag[kk:] = 0
        camp[kk:] = 0.0
        k = kk
        return t_new, y_new

    attempt, prev, dscale, nrej = 0, None, 1.0, 0
    while t < T_end and total < max_nfev and k < len(te) and attempt < 4:
        if mode == "normal":
            sol = new_solver(t, y, attempt)
            hist = []
            buf = {}
            t_ok, y_ok = t, y
            while sol.status == "running" and total + sol.nfev < max_nfev:
                if not step(sol):
                    break
                collect(sol, buf)
                hist.append(sol.t)
                t_ok, y_ok = sol.t, sol.y.copy()
                if len(hist) > stall_steps and hist[-stall_steps - 1] > 0 and \
                        np.log(hist[-1] / hist[-stall_steps - 1]) < stall_dlnt:
                    mode = "avg"
                    break
            for kk, v in buf.items():
                out[:, kk] = v
            prev = None
            total += sol.nfev
            t, y = t_ok, y_ok                               # 失敗した解法の y は使わない
            attempt = (0 if len(hist) > 100 else attempt + 1) if sol.status == "failed" else 0
            continue
        # ---- 平均の区間 ----
        span_max = jump_frac * t if t > 0 else T_end
        ts, ys, buf = [t], [y.copy()], {}
        k0 = k
        lim_calls = burst_calls
        easy = False
        used, batt = 0, 0
        bsol = new_solver(t, y, 0, rtol_b)
        while True:
            while used + bsol.nfev < lim_calls and ts[-1] - t < span_max and bsol.status != "finished" and batt < 4:
                if not step(bsol):                          # 失敗: 最後の状態から次の方法で続ける
                    used += bsol.nfev
                    batt += 1
                    bsol = new_solver(ts[-1], ys[-1], batt, rtol_b)
                    continue
                ts.append(bsol.t)
                ys.append(bsol.y.copy())
                collect(bsol, buf)
            if bsol.status == "finished" or ts[-1] - t >= span_max or batt >= 4:
                easy = True
                break
            T_, Y_ = np.array(ts), np.array(ys).T
            sel = T_ >= T_[0] + relax * (T_[-1] - T_[0])
            Ts, Ys = T_[sel], Y_[:, sel]
            if not np.all(np.isfinite(Ys)) or Ts[-1] <= Ts[0] * (1 + 1e-12):   # 進まないバースト(歩が時刻の刻みより小さい)
                easy = True
                break
            if len(Ts) < 50:
                lim_calls *= 2
                if lim_calls > burst_max:
                    easy = True
                    break
                continue
            tc, m, s = _hann_fit(Ts, Ys)
            h = len(Ts) // 2
            _, _, sa = _hann_fit(Ts[:h], Ys[:, :h])
            _, _, sb = _hann_fit(Ts[h:], Ys[:, h:])
            if not (np.all(np.isfinite(s)) and np.all(np.isfinite(sa)) and np.all(np.isfinite(sb))):
                easy = True
                break
            ok_slow = np.all(np.abs(sa[N:] - sb[N:]) <= 0.3 * np.abs(s[N:]) + 1e-30)
            if ok_slow or lim_calls >= burst_max:
                break
            lim_calls *= 2
        used += bsol.nfev
        total += used
        if easy:
            for kk, v in buf.items():
                out[:, kk] = v
            if len(ts) == 1:                                # バーストが 1 歩も進めない
                attempt += 1
            t, y = ts[-1], ys[-1].copy()
            mode, prev = "normal", None
            continue
        if prev is not None:
            # 跳びの検査: 跳んだ先のバーストの M・J の傾きが、跳びに使った傾きと(このバーストの傾きの揺れを除いて)3 割より
            # 違えば、跳びが大きすぎた(雲がなじまない中間の速さの成分、または傾きの曲がり)。跳びを半分にしてやり直す。
            chk = np.zeros(N + 2, bool)                      # M・J と、前の跳びで進めた多い準位
            chk[N:] = True
            chk[:N] = prev["proj"] & (prev["y_end"][:N] > np.max(prev["y_end"][:N]) - 23.0)
            q = float(np.max((np.abs(s[chk] - prev["s"][chk]) - np.abs(sa[chk] - sb[chk]))
                             / (np.maximum(np.abs(s[chk]), np.abs(prev["s"][chk])) + 1e-300)))
            if q > 0.3:
                nrej += 1
                dscale *= 0.5
                D_new = 0.5 * prev["D"]
                if D_new >= prev["W"]:
                    prev["D"] = D_new
                    t, y = jump(prev, D_new)
                    continue
                for kk, v in prev["buf"].items():          # 跳ぶ得がない: 前のバーストの終わりからふつうに解く
                    out[:, kk] = v
                k = prev["k_end"]
                out[:, k:] = np.nan
                flag[k:] = 0
                camp[k:] = 0.0
                t, y = prev["t_end"], prev["y_end"].copy()
                mode, prev = "normal", None
                continue
            if q < 0.1:
                dscale = min(1.0, 1.5 * dscale)
        # 傾きの変化(前のバーストとの差)から曲がりを出し、跳びを 2 次にする(1 次だと、加速する変化で遅れる)
        acc = (s - prev["s"]) / (tc - prev["tc"]) if prev is not None and tc > prev["tc"] else np.zeros_like(s)
        y_end, t_end = Y_[:, -1], T_[-1]
        W = t_end - ts[0]
        lim = [span_max, T_end - t_end]
        if s[N + 1] != 0:
            lim.append(dchi * m[N] ** 2 / abs(s[N + 1]))
        if s[N] != 0:
            lim.append(dlnM * m[N] / abs(s[N]))
        lnmax = np.max(y_end[:N])
        proj = (np.abs(s[:N]) * span_max > 0.1) & (np.abs(sa[:N] - sb[:N]) <= 0.3 * np.abs(s[:N]))
        for i in np.flatnonzero(proj):
            if y_end[i] > lnmax - 23.0:                     # 最大の 1e-10 倍より多い準位: 1 回で dln まで
                lim.append(dln / abs(s[i]))
            elif s[i] > 0:                                  # それより少なく育つ準位: 1 回で dln_small まで、最大の e^-10 倍の手前まで
                lim.append(min(dln_small, max(dln, lnmax - 10.0 - y_end[i])) / s[i])
        D = dscale * max(min(lim), 0.0)
        if D < W:                                           # 跳ぶ得がない: ふつうに続ける
            for kk, v in buf.items():
                out[:, kk] = v
            t, y = t_end, y_end.copy()
            mode, prev = "normal", None
            continue
        chi_res = Ys[N + 1] / Ys[N] ** 2                   # chi の振れ: 直線を引いた残りの幅の半分
        xc = (Ts - Ts.mean()) / (Ts[-1] - Ts[0])
        bc = np.dot(xc, chi_res - chi_res.mean()) / max(np.dot(xc, xc), 1e-300)
        amp = 0.5 * float(np.ptp(chi_res - chi_res.mean() - bc * xc))
        prev = dict(t_end=t_end, y_end=y_end.copy(), tc=tc, m=m, s=s, a=acc, proj=proj, flag=1 if ok_slow else 2, amp=amp,
                    W=W, D=D, k0=k0, k_end=k, buf=buf)
        t, y = jump(prev, D)
        if trace is not None:
            trace.append(dict(t=t_end, W=W, D=D, dscale=dscale, lim=[float(x) for x in lim], chi0=float(y_end[N + 1] / y_end[N] ** 2),
                              chi1=float(y[N + 1] / y[N] ** 2), proj=proj.tolist(), s=s.tolist(), ok=bool(ok_slow), y_end=y_end.tolist()))
        njump += 1
    done = int(np.sum(np.all(np.isfinite(out), axis=0)))
    ok = done == len(te)
    return SimpleNamespace(status=0 if ok else -1, y=out[:, :done], t=te[:done], avg=flag[:done], chi_amp=camp[:done],
                           njump=njump, nrej=nrej, nfev=total, message="ok" if ok else f"stopped at t = {t / sr.YR:.3e} yr")


# 平均の確かめに使う細かい設定: 跳びを半分、少なく育つ準位も dln で縛る(緩めた規則を使わない)、バーストの許容誤差を 1e-7
AVG_FINE = dict(dchi=5e-4, dlnM=2.5e-4, dln=0.125, dln_small=0.125, burst_calls=10000, burst_rtol=1e-7)


def _tol_dist(t, ca, cb, dt_rel):
    """年齢ごとの |ca - cb| を、時刻の相対的なずれ dt_rel まで許して測る(片方を log t で直線補間し、窓の中の最小。
    両方向の大きい方)。"""
    lt = np.log(np.maximum(t, 1e-300))
    q = lt[:, None] + np.linspace(-np.log1p(dt_rel), np.log1p(dt_rel), 11)[None, :]

    def one(cx, cy):
        return np.min(np.abs(cx[:, None] - np.interp(q, lt, cy)), axis=1)
    return np.maximum(one(ca, cb), one(cb, ca))


def evolveN_checked(sysN, M0, chi0, mu_eV, f_gev, T_yr, t_eval_yr, tol=5e-3, dt_rel=5e-3, max_nfev=20000000, fine=None,
                    **kw):
    """evolveN_avg に計算の上限と、平均の確かめの印をつけたもの(DESIGN.md D20、2026-11-17 のユーザー決定「Cap and flag」)。
    既定の設定で解き(呼び出し max_nfev 回まで)、平均を使ったら細かい設定(AVG_FINE)でもう一度解いて chi を比べる。
    年齢ごとの差は時刻のずれ dt_rel を許して測る(_tol_dist)。n_valid は、両方が届き差が tol 以内の、初めから続く年齢の数で、
    表ではこの年齢までを使う(reftab の tau_valid)。平均を使わなければ解そのもので、届いた年齢まで有効。
    tol = 5e-3 と dt_rel = 0.5% は、表の約束(|Δχ| ≤ 0.03、段差の時刻 2%)の 6 分の 1 と 4 分の 1(表の補間の誤差 0.015-0.018
    に足しても約束の中)。
    返り値: y、t は細かい方の解(届かない年齢は既定の解。平均を使わなければ既定の解)、status(0 = 最後の年齢まで有効)、
    n_valid、dchi(年齢ごとの差、比べられない年齢は NaN)、avg(値を出した方の年齢ごとの平均の印)、reason
    (有効な範囲を決めたもの: "" / "cap"(上限で止まった)/ "fail"(解法が止まった)/ "disagree")、runs(二つの解の status・njump・nrej・nfev)。"""
    a = evolveN_avg(sysN, M0, chi0, mu_eV, f_gev, T_yr, t_eval_yr, max_nfev=max_nfev, **kw)
    na = a.y.shape[1]
    runs = [dict(status=a.status, njump=a.njump, nrej=a.nrej, nfev=a.nfev)]
    if not np.any(a.avg > 0):
        return SimpleNamespace(status=a.status, y=a.y, t=a.t, n_valid=na, dchi=np.zeros(na), avg=a.avg, chi_amp=a.chi_amp,
                               reason="" if a.status == 0 else ("cap" if a.nfev >= max_nfev else "fail"), runs=runs,
                               message=a.message)
    b = evolveN_avg(sysN, M0, chi0, mu_eV, f_gev, T_yr, t_eval_yr, max_nfev=max_nfev, **{**kw, **(AVG_FINE if fine is None else fine)})
    nb = b.y.shape[1]
    runs.append(dict(status=b.status, njump=b.njump, nrej=b.nrej, nfev=b.nfev))
    n = min(na, nb)
    d = np.full(max(na, nb), np.nan)
    if n:
        d[:n] = _tol_dist(a.t[:n], chi_of(a)[:n], chi_of(b)[:n], dt_rel)
    bad = np.flatnonzero(~(d[:n] <= tol))
    n_valid = int(bad[0]) if len(bad) else n
    if nb >= na:
        y, t, avg, amp = b.y, b.t, b.avg, b.chi_amp
    else:
        y, t = np.concatenate([b.y, a.y[:, nb:]], axis=1), a.t
        avg, amp = np.concatenate([b.avg, a.avg[nb:]]), np.concatenate([b.chi_amp, a.chi_amp[nb:]])
    nt = len(np.atleast_1d(t_eval_yr))
    reason = "disagree" if len(bad) else ("" if n_valid == nt else ("cap" if max(a.nfev, b.nfev) >= max_nfev else "fail"))
    msg = "ok" if n_valid == nt else (f"valid to t = {t[n_valid - 1] / sr.YR:.3e} yr ({reason})" if n_valid else f"no valid age ({reason})")
    return SimpleNamespace(status=0 if n_valid == nt else -1, y=y, t=t, n_valid=n_valid, dchi=d, avg=avg, chi_amp=amp,
                           reason=reason, runs=runs, message=msg)


def chi_of(sol, M0=None):
    """解からスピン J/M^2。"""
    N = sol.y.shape[0] - 2
    return sol.y[N + 1] / sol.y[N] ** 2


def a14(sysN, sol, mu_eV, f_gev):
    """BGLS の式 A14 の摂動の条件を、解の各時刻で調べる: x_p = Gamma_k tau_X(= N_X / N_k)。
    Gamma_k = (過程 p が準位 k を作る率)/(n_k + n_q) は誘導放出の 1 個あたりの率、tau_X は X が雲にとどまる時間
    (BH が吸収する過程と共鳴する過程だけ。放出の過程は X が出ていく)。率の時間の単位は 1/mu、tau_X は M で、mu M = alpha
    なので x_p = alpha (c_inf + c_H) (M_pl/f)^4 eps_i eps_j tau_X。返り値: (時刻ごとの x の最大値, その過程の番号)。"""
    s4 = (si.MPL_GEV / f_gev) ** 4
    a_sun = sr.alpha_of(1.0, mu_eV)
    N = sysN.nl
    xmax, arg = np.zeros(sol.y.shape[1]), np.zeros(sol.y.shape[1], int)
    for t in range(sol.y.shape[1]):
        n = np.exp(sol.y[:N, t])
        M = sol.y[N, t]
        chi = min(max(sol.y[N + 1, t] / M ** 2, 0.0), 0.9999)
        a = a_sun * M
        ci, ch = sysN.coef(a, chi)
        e = n / M ** 2
        x = a * (ci + ch) * s4 * e[sysN.I] * e[sysN.J] * sysN.tau(a, chi)
        x = np.where(np.isfinite(x), x, -np.inf)
        arg[t] = int(np.argmax(x))
        xmax[t] = x[arg[t]]
    return xmax, arg


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


def corr_grid():
    """補正の確かめの例: M0 = 30、chi0 = 0.95 で、alpha0 を 11 点(0.03-0.45)、f を 5 点(1e13-1e17 GeV)。"""
    mu1 = sr.alpha_of(30.0, 1e-12) / 1e-12                 # alpha = mu1 x mu [eV](M = 30)
    return [(30.0, 0.95, float(a / mu1), f) for a in (0.03, 0.05, 0.08, 0.12, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
            for f in (1e13, 1e14, 1e15, 1e16, 1e17)]


def _corr_job(case):
    import relcorr
    M, c0, mu, f = case
    procs = process_table()
    s0, s1 = System(LEVELS8, procs), System(LEVELS8, procs, corr=relcorr.Corr(relcorr.load()))
    Ts = np.logspace(3, 10, 120)
    t0 = time.time()
    b0 = evolveN_checked(s0, M, c0, mu, f, Ts[-1], Ts)
    t1 = time.time()
    b1 = evolveN_checked(s1, M, c0, mu, f, Ts[-1], Ts)
    t2 = time.time()
    c0s, c1s = chi_of(b0), chi_of(b1)
    n = min(b0.n_valid, b1.n_valid)
    x, arg = a14(s1, b1, mu, f)
    x, arg = x[:b1.n_valid], arg[:b1.n_valid]                # 有効な年齢だけで
    k = int(np.argmax(x)) if len(x) else 0
    return dict(M=M, chi0=c0, mu=mu, f=f, alpha0=sr.alpha_of(M, mu), status=[b0.status, b1.status],
                reason=[b0.reason, b1.reason], t_valid_yr=[float(Ts[b.n_valid - 1]) if b.n_valid else 0.0 for b in (b0, b1)],
                chi_end=[float(c0s[-1]) if len(c0s) else np.nan, float(c1s[-1]) if len(c1s) else np.nan],
                max_dchi=float(np.max(np.abs(c1s[:n] - c0s[:n]))) if n else np.nan,
                check_dchi=[float(np.nanmax(b.dchi[:b.n_valid])) if b.n_valid else np.nan for b in (b0, b1)],
                a14_max=float(x[k]) if len(x) else np.nan, a14_proc=s1.names[arg[k]] if len(x) else "",
                a14_t_yr=float(b1.t[k] / sr.YR) if len(x) else np.nan, secs=[t1 - t0, t2 - t1],
                n_avg=[int((b0.avg > 0).sum()), int((b1.avg > 0).sum())], runs=[b0.runs, b1.runs],
                chi_amp_max=[float(b0.chi_amp.max()) if len(b0.chi_amp) else 0.0, float(b1.chi_amp.max()) if len(b1.chi_amp) else 0.0])


def check_corr(nproc=8):
    """補正を入れた運動論の確かめ。
    (1) 恒等の補正(放出は C_inf = 1、BH は C_H = 1。名前の対応を確かめるため逆順で渡す)が corr なしと同じ解になる。
    (2) 右辺の保存則: 補正を入れ、無限遠の経路と重力波を切ると、どの状態でも dM + alpha_sun sum omega_i dn_i = 0、
        dJ + sum m_i dn_i = 0(地平線の経路の omega_X、m_X の勘定。放出の過程の C_H < 0 も含む)。乱数の 400 状態で。
    (3) 解に沿った角運動量の保存(無限遠と重力波を切る)。誤差は積分の精度(rtol)の程度。エネルギーは M + alpha_sun sum omega n
        が、BH が変わると雲の束縛のエネルギーが変わる分(運動論の dM = -omega dN が入れない、補正なしでも同じ)だけずれる。
    (4) corr_grid の 55 例で(evolveN_checked で解く)、補正のあり・なしのスピンの差(両方が有効な年齢で)と、BGLS A14 の
        x = Gamma_k tau_X の最大値、有効な年齢の範囲。"""
    import relcorr
    procs = process_table()
    C = relcorr.Corr(relcorr.load())
    bh = np.array([p["kind"] != "inf" for p in procs])

    class Ident:
        names = C.names[::-1]                              # 逆順に並べて、名前の対応も確かめる

        def __call__(self, a, c):
            return np.where(bh, 0.0, 1.0)[::-1], np.where(bh, 1.0, 0.0)[::-1]
    s0, si_, s1 = System(LEVELS8, procs), System(LEVELS8, procs, corr=Ident()), System(LEVELS8, procs, corr=C)
    Ts = np.logspace(4, 10, 30)
    M, c0, mu, f = 30.0, 0.9, 0.71e-12, 1e15
    b0, bi = evolveN(s0, M, c0, mu, f, Ts[-1], Ts), evolveN(si_, M, c0, mu, f, Ts[-1], Ts)
    n = min(b0.y.shape[1], bi.y.shape[1])
    ident_diff = float(np.max(np.abs(b0.y[:, :n] - bi.y[:, :n])))
    # (2) 右辺の保存則
    N = len(LEVELS8)
    ms = np.array([lv[2] for lv in LEVELS8])
    rng = np.random.default_rng(1)
    worst = [0.0, 0.0]
    for _ in range(400):
        M, a, chi = 30.0, rng.uniform(0.02, 0.5), rng.uniform(0.0, 0.99)
        mu = a / sr.alpha_of(M, 1.0)
        f = 10 ** rng.uniform(12.5, 17)
        ln = np.log(1.0 / sr.N0) + rng.uniform(0, 1, N) * (np.log(0.05 * M ** 2 / a) - np.log(1.0 / sr.N0))
        y = np.concatenate([ln, [M, chi * M ** 2]])
        out = make_rhs(s1, M, mu, f, gw=False, inf_off=True)(0.0, y)
        dn = np.exp(ln) * out[:N]
        w = np.array([sr.omega_over_mu(a, *lv, chi) for lv in LEVELS8])
        a_sun = sr.alpha_of(1.0, mu)
        rE = abs(out[N] + a_sun * np.dot(w, dn)) / (abs(out[N]) + a_sun * np.sum(np.abs(w * dn)))
        rJ = abs(out[N + 1] + np.dot(ms, dn)) / (abs(out[N + 1]) + np.sum(np.abs(ms * dn)))
        worst = [max(worst[0], rE), max(worst[1], rJ)]
    # (3) 解に沿った保存
    cons = []
    for M, c0, mu, f in [(30.0, 0.9, 0.71e-12, 1e15), (40.0, 0.95, 1.0e-12, 1e15), (50.0, 0.95, 1.2e-12, 1e16)]:
        b = evolveN(s1, M, c0, mu, f, Ts[-1], Ts, gw=False, inf_off=True)
        nn = np.exp(b.y[:N])
        Mb, Jb = b.y[N], b.y[N + 1]
        a_sun = sr.alpha_of(1.0, mu)
        Jtot = Jb + ms @ nn
        w = np.array([[sr.omega_over_mu(a_sun * Mb[t], *lv, Jb[t] / Mb[t] ** 2) for t in range(len(Mb))] for lv in LEVELS8])
        Etot = Mb + a_sun * np.sum(w * nn, axis=0)
        cons.append(dict(M=M, chi0=c0, mu=mu, f=f, status=b.status, dJ_rel=float(np.max(np.abs(Jtot / Jtot[0] - 1))),
                         dE_rel=float(np.max(np.abs(Etot / Etot[0] - 1))), chi_end=float(Jb[-1] / Mb[-1] ** 2)))
    from multiprocessing import Pool
    cases = corr_grid()
    with Pool(min(nproc, len(cases))) as p:
        rows = p.map(_corr_job, cases, chunksize=1)
    return dict(identity_max_diff=ident_diff, rhs_conservation=dict(E=worst[0], J=worst[1]), trajectory=cons, cases=rows)


AVG_FORCED = [(a, f, c) for c in ("corr", "none") for (a, f) in ((0.03, 1e13), (0.05, 1e15), (0.08, 1e13), (0.08, 1e17),
              (0.12, 1e13), (0.12, 1e15), (0.16, 1e16), (0.25, 1e13), (0.25, 1e15), (0.30, 1e14), (0.40, 1e13), (0.40, 1e15))]
AVG_STALLED = [(0.35, 1e13, "none"), (0.20, 1e15, "corr"), (0.16, 1e14, "corr"), (0.45, 1e16, "none")]
AVG_REF_T2 = {(0.20, 1e15, "corr"): 5e6, (0.45, 1e16, "none"): 6e7, (0.16, 1e14, "corr"): 1.5e8}
FORCE = dict(stall_steps=300, stall_dlnt=2.0, burst_calls=400, jump_frac=1.0)


def _avg_sys(cflag):
    import relcorr
    return System(LEVELS8, process_table(), corr=relcorr.Corr(relcorr.load()) if cflag == "corr" else None)


def _avg_forced_job(case):
    """平均の区間を無理に使い(ふつうに解ける例で)、evolveN と全年齢の chi を比べる。"""
    a0, f, cflag = case
    s_ = _avg_sys(cflag)
    mu = a0 / (sr.alpha_of(30.0, 1e-12) / 1e-12)
    Ts = np.logspace(3, 10, 120)
    cr = chi_of(evolveN(s_, 30.0, 0.95, mu, f, Ts[-1], Ts))
    b = evolveN_avg(s_, 30.0, 0.95, mu, f, Ts[-1], Ts, **FORCE)
    cb = chi_of(b)
    n = min(len(cr), len(cb))
    return dict(case=case, status=b.status, njump=b.njump, n_avg=int((b.avg > 0).sum()),
                max_dchi=float(np.max(np.abs(cb[:n] - cr[:n]))) if n == 120 else None)


def _avg_stall_job(case):
    """止まる例を evolveN_checked で解く(既定の設定と細かい設定を比べ、有効な年齢の範囲を決める)。"""
    a0, f, cflag = case
    s_ = _avg_sys(cflag)
    mu = a0 / (sr.alpha_of(30.0, 1e-12) / 1e-12)
    Ts = np.logspace(3, 10, 120)
    t0 = time.time()
    b = evolveN_checked(s_, 30.0, 0.95, mu, f, Ts[-1], Ts)
    return dict(case=case, status=b.status, reason=b.reason, n_valid=b.n_valid,
                t_valid_yr=float(Ts[b.n_valid - 1]) if b.n_valid else 0.0, runs=b.runs, secs=time.time() - t0,
                chi=chi_of(b).tolist(), dchi=[None if not np.isfinite(x) else float(x) for x in b.dchi],
                n_avg=int((b.avg > 0).sum()), chi_amp=float(b.chi_amp.max()) if len(b.chi_amp) else 0.0)


def _avg_ref_job(args):
    """止まる例の参照: 振動まで解く(ふつうの LSODA)。wall 秒で届いた年齢まで。"""
    from scipy.integrate import LSODA
    case, t2, wall = args
    a0, f, cflag = case
    s_ = _avg_sys(cflag)
    mu = a0 / (sr.alpha_of(30.0, 1e-12) / 1e-12)
    Ts = np.logspace(3, 10, 120)
    te = Ts[Ts <= t2] * sr.YR
    sol = LSODA(make_rhs(s_, 30.0, mu, f), 0.0, np.array([np.log(1.0 / sr.N0)] * 8 + [30.0, 0.95 * 900.0]), te[-1],
                rtol=1e-8, atol=1e-12)
    out, k, t0 = [], 0, time.time()
    while sol.status == "running" and time.time() - t0 < wall:
        sol.step()
        if k < len(te) and te[k] <= sol.t:
            d = sol.dense_output()
            while k < len(te) and te[k] <= sol.t:
                yy = d(te[k])
                out.append(float(yy[9] / yy[8] ** 2))
                k += 1
    return dict(case=case, reached_yr=float(sol.t / sr.YR), chi=out)


def check_avg(nproc=8, refs=False, wall=2400):
    """evolveN_avg(バーストの平均)と evolveN_checked(上限と印)の確かめ(DESIGN.md D20)。
    (1) ふつうに解ける 24 例で平均の区間を無理に使い(FORCE)、evolveN と全年齢の chi を比べる。
    (2) 止まる 4 例を evolveN_checked で解く: 有効な年齢の範囲、それを決めたもの(上限・設定で値が動く)、既定と細かい
        設定の差。
    (3) refs=True: 止まる 3 例で、振動まで解いた参照(wall 秒で届いた年齢まで)と、出した値を比べる(有効な年齢と、全体)。"""
    from multiprocessing import Pool
    with Pool(nproc) as p:
        rs = p.map_async(_avg_stall_job, AVG_STALLED, chunksize=1)
        rr = p.map_async(_avg_ref_job, [(c, t2, wall) for c, t2 in AVG_REF_T2.items()], chunksize=1) if refs else None
        rf = p.map(_avg_forced_job, AVG_FORCED, chunksize=1)
        stall = rs.get()
        ref = rr.get() if refs else []
    out = dict(forced=rf, forced_max=max(r["max_dchi"] for r in rf if r["max_dchi"] is not None),
               forced_failed=[r["case"] for r in rf if r["max_dchi"] is None])
    st = []
    for d in stall:
        nv = d["n_valid"]
        row = dict(case=d["case"], status=d["status"], reason=d["reason"], t_valid_yr=d["t_valid_yr"], runs=d["runs"],
                   secs=round(d["secs"]), n_avg=d["n_avg"], chi_amp=d["chi_amp"], chi_end=d["chi"][-1] if d["chi"] else None,
                   check_dchi_valid=max(d["dchi"][:nv]) if nv else None,
                   check_dchi_all=max((x for x in d["dchi"] if x is not None), default=None))
        for r in ref:
            if r["case"] == d["case"]:
                m = min(len(d["chi"]), len(r["chi"]))
                dd = np.abs(np.array(d["chi"][:m]) - np.array(r["chi"][:m]))
                row["ref_reached_yr"] = r["reached_yr"]
                row["vs_ref_valid"] = float(dd[:min(m, nv)].max()) if min(m, nv) else None
                row["vs_ref_all"] = float(dd.max()) if m else None
        st.append(row)
    out["stalled"] = st
    return out


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "check"
    if what == "procs":
        p = process_table(rebuild="--rebuild" in sys.argv)
        print(len(p), "processes")
    elif what == "avgcheck":
        nums = [int(a) for a in sys.argv[2:] if a.isdigit()]
        print(json.dumps(check_avg(nums[0] if nums else 8, refs="--refs" in sys.argv), indent=1))
    elif what == "corrcheck":
        print(json.dumps(check_corr(int(sys.argv[2]) if len(sys.argv) > 2 else 8), indent=1))
    else:
        rep = {"evolve3m": check_evolve3m(), "fixed_mass": check_fixed_mass(), "conservation": check_conservation()}
        print(json.dumps(rep, indent=1))
