"""自己相互作用の過程の相対論的な補正の表 C(alpha, chi)(値の表の作り直し、DESIGN.md D18 の手順 3)。

kinn.process_table の 74 個の過程(BH が吸収 18、n00 に共鳴 13、無限遠への放出 43)を、relrates.point で Kerr 時空の
Green 関数から解く。補正は主要項(hyd.py)との比で、過程の率を
    Gamma_hat = kappa_hyd alpha^(q-3) [r+] (C_inf + C_H)   (q = 8 放出、11 BH、7 共鳴。[r+] は BH と共鳴だけ)
と書いたときの C_inf(無限遠へ出る X)と C_H(地平線に吸われる X)。BH が吸収する過程は C_inf = 0。放出の過程の C_H は
m_X = 0 で正(BH が X を吸う)、m_X > 0 で負(超放射で BH から X が出る)。運動論では C_H の分を BH の経路として扱う
(BH が omega_X のエネルギーと m_X の角運動量を受け取る)。
tau_X(束縛の X だけ)は強制振動の X が雲にとどまる時間(M の単位)。BGLS の式 A14 の摂動の条件は
(この過程が作る準位 k の 1 個あたりの率) x tau_X << 1 で、運動論の中で確かめる。
格子: chi は CHI の 11 点、alpha は 0.004-0.1 を対数で 50 点、0.105-0.52 を 0.005 おき(n00 の極をまたぐ所の幅は
alpha で 0.03 以上(polescan、2026-11-16)なので、0.005 おきで解ける。l = 2 の極(820・720・620)は幅が 1e-7 以下で
解けないが、その部分波の寄与は極のごく近くでしか効かない)。
準位の根: chi ごとに alpha の小さい方から倍精度の Leaver で追う(一つ前の根を水素型の比で延ばした予測から始め、
同じ l の隣の n との間隔の 0.3 倍より予測に近いことを確かめる。外れたら水素型からやり直す)。動径の節の数も確かめる。
出力: out/relcorr.npz(途中も保存し、--resume で続ける)。記録: out/relcorr_build.log。
運動論からは Corr(load()) で引く(clean で閾値の先・解けなかった点・とがりを整理してから、C をそのまま 3 次で補間)。
確かめ(2026-11-16、python3 calc/relcorr.py check): 1474 点すべて解け、根は多倍長の表(rates_mp)と準位の間隔の 1.5% 以内
(8 準位すべて)、つなぎ目 6e-8、Wronskian 9e-6、無限遠への流れ 6e-6。
alpha <= 0.02 に 3 次を当てはめた切片が 1 と 1.2e-3。格子の外の 130 点(relcorr.py holdout、9530 個の値)で補間の誤差は
中央値 2.6e-6、99% で 4.4e-3、幅の狭い l = 2 の共鳴のある 411x411->322(alpha >~ 0.39)と 511x511->422(alpha >~ 0.30)を
除いて最大 0.9%(閾値の 0.03 以内の放出の側では 5%)。節の数は、強く減衰する準位(alpha >= 0.425 で超放射しない chi の
211・411・511)で波動関数の位相が回って数え方が崩れる(59 点)が、根は多倍長の表と同じ枝。
限界: 二つの過程の幅の狭い共鳴は解けず背景の値を使う(そこでは A14 も破れうる)。閾値を越えた 411x511->322(alpha 0.365-0.42)・
411x533->322(0.455-0.485)・411x544->322(0.465-0.495)の束縛の側は、極で平均した率。
2026-11-17: 独立な監査で relrates の誤り 2 件(tau_X、減衰する準位の規格化)を直し、表を作り直す(上の数値は 2026-11-16 の版、
古い表は out/relcorr_v1_20261116.npz)。閾値の先の埋め方は、X の部分波が l >= 2 で極が孤立しているので前提が成り立たず、
扱いは未決(DESIGN.md D19 の追記)。
実行: cd research && ~/miniforge3/bin/python3 calc/relcorr.py build [--resume] [並列の数]
"""
import os
import sys
import time
import warnings

import numpy as np
from scipy.optimize import newton

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kinn  # noqa: E402
import leaver  # noqa: E402
import relrates  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "relcorr.npz")
LEVELS = kinn.LEVELS8
CHI = np.array([0.0, 0.2, 0.4, 0.55, 0.7, 0.8, 0.87, 0.92, 0.955, 0.98, 0.995])
ALPHA = np.concatenate([np.geomspace(0.004, 0.1, 50), np.round(np.arange(0.105, 0.5201, 0.005), 3)])
# 点ごとの記録の符号
NOTE = {"": 0, "own grid": 1, "near threshold": 2, "error": 3, "bound X": 4, "unbound X": 5}


def procs():
    P = kinn.process_table()
    return P, [(LEVELS[p["i"]], LEVELS[p["j"]], LEVELS[p["k"]]) for p in P]


def spacing(a, n, l):
    gaps = [abs(1.0 / (2 * n ** 2) - 1.0 / (2 * k ** 2)) for k in (n - 1, n + 1) if k > l]
    return a ** 3 * min(gaps)


def hyd_re(a, lv, c):
    return a * sr.omega_over_mu(a, *lv, c)


def roots_column(c, alphas=None):
    """一つの chi について、8 準位の根を alpha の小さい方から追う。返り値: (chi, om (Na, 8) 複素, ok (Na, 8))。"""
    warnings.simplefilter("ignore")
    ALPHA_ = ALPHA if alphas is None else np.asarray(alphas)
    om = np.full((len(ALPHA_), len(LEVELS)), np.nan + 0j)
    ok = np.zeros((len(ALPHA_), len(LEVELS)), bool)
    for v, lv in enumerate(LEVELS):
        n, l, m = lv
        prev = None
        for q, a in enumerate(ALPHA_):
            h = hyd_re(a, lv, c)
            pred = h if prev is None else prev[1].real * h / hyd_re(prev[0], lv, c)
            tol = 0.3 * spacing(a, n, l)
            guesses = [complex(h, 0.0)] if prev is None else [complex(pred, prev[1].imag * (a / prev[0]) ** (4 * l + 5)),
                                                              complex(h, 0.0)]
            w = None
            for g in guesses:
                try:
                    # Newton を丸めの所まで回す(tol 1e-15 では 511 の Re omega が alpha = 0.004 で 4e-16 ずれ、共鳴する過程の
                    # 離調 omega_X - omega_n00(~1e-14)が数 % 動いた。1e-20 では多倍長の根と 1e-18 で合う)
                    N_ = int(max(4000, 300 / a))
                    ww = complex(newton(lambda w: leaver.cf_residual(w, a, c, l, m, N_), g, x1=g * (1 + 1e-7),
                                        tol=1e-20, maxiter=100, disp=False))
                except Exception:
                    continue
                if abs(ww.real - pred) < tol:
                    w = ww
                    break
            if w is None:
                prev = None
                continue
            om[q, v], ok[q, v] = w, True
            prev = (a, w)
    return c, om, ok


def job(args):
    """一つの (alpha, chi)。返り値: 添字と、過程ごと・準位ごとの配列。"""
    warnings.simplefilter("ignore")
    ia, ic, a, c, om_row = args
    t0 = time.time()
    _, plist = procs()
    omegas = {lv: complex(om_row[v]) for v, lv in enumerate(LEVELS)}
    res, lev = relrates.point(a, c, plist, omegas)
    return ia, ic, res, lev, time.time() - t0


def new_state(P):
    Na, Nc, Np, Nl = len(ALPHA), len(CHI), len(P), len(LEVELS)
    f = lambda: np.full((Na, Nc, Np), np.nan)   # noqa: E731
    return dict(C_inf=f(), C_H=f(), tau_X=f(), omegaX=f(), W_spread=f(), J_spread=f(), asym_err=f(),
                note=np.zeros((Na, Nc, Np), np.int8), om=np.full((Na, Nc, Nl), np.nan + 0j),
                root_ok=np.zeros((Na, Nc, Nl), bool), nodes=np.full((Na, Nc, Nl), -1, np.int16),
                match_err=np.full((Na, Nc, Nl), np.nan), done=np.zeros((Na, Nc), bool), sec=np.full((Na, Nc), np.nan))


def save(st, P):
    tmp = OUT + ".tmp.npz"
    np.savez(tmp, alpha=ALPHA, chi=CHI, levels=np.array(LEVELS),
             names=np.array([f"{kinn.name(LEVELS[p['i']])}x{kinn.name(LEVELS[p['j']])}->{kinn.name(LEVELS[p['k']])}"
                             for p in P]),
             kinds=np.array([p["kind"] for p in P]), kappa_hyd=np.array([p["kappa"] for p in P]),
             note_codes=np.array(list(NOTE.keys())), **st)
    os.replace(tmp, OUT)


def build(nproc, resume):
    from multiprocessing import Pool
    P, plist = procs()
    st = new_state(P)
    if resume and os.path.exists(OUT):
        z = np.load(OUT)
        assert np.array_equal(z["alpha"], ALPHA) and np.array_equal(z["chi"], CHI)
        for key in st:
            st[key] = z[key]
    # 準位の根(chi ごとに alpha で追う)
    t0 = time.time()
    with Pool(min(nproc, len(CHI))) as pool:
        for c, om, ok in pool.imap_unordered(roots_column, list(CHI)):
            ic = int(np.argmin(np.abs(CHI - c)))
            st["om"][:, ic], st["root_ok"][:, ic] = om, ok
    print(f"roots: {time.time() - t0:.0f}s, failed {int((~st['root_ok']).sum())} of {st['root_ok'].size}", flush=True)
    jobs = [(ia, ic, float(a), float(c), st["om"][ia, ic]) for ia, a in enumerate(ALPHA) for ic, c in enumerate(CHI)
            if not st["done"][ia, ic] and st["root_ok"][ia, ic].all()]
    print(f"{len(jobs)} points to solve with {nproc} workers", flush=True)
    pw = np.array([kinn.POW[p["kind"]] for p in P])
    kap = np.array([p["kappa"] for p in P])
    t0 = time.time()
    ndone = 0
    with Pool(nproc) as pool:
        for ia, ic, res, lev, sec in pool.imap_unordered(job, jobs, chunksize=1):
            a, c = ALPHA[ia], CHI[ic]
            rp = 1 + np.sqrt(1 - c * c)
            norm = kap * a ** (pw - 3) * np.where(np.array([p["kind"] != "inf" for p in P]), rp, 1.0)
            for q, (p, r) in enumerate(zip(P, res)):
                st["C_inf"][ia, ic, q] = r["G_inf"] / norm[q]
                st["C_H"][ia, ic, q] = r["G_h"] / norm[q]
                for key in ("tau_X", "omegaX", "W_spread", "J_spread", "asym_err"):
                    st[key][ia, ic, q] = r[key]
                note = r["note"]
                code = NOTE["error"] if note.startswith("error") else NOTE.get(note, 0)
                # X の束縛・非束縛が主要項と入れ替わった点(自分の格子で解いた点も含めて)。Corr は omegaX から直接判定する
                if code in (0, 1) and (r["omegaX"] > a) != (p["kind"] == "inf"):
                    code = NOTE["bound X"] if p["kind"] == "inf" else NOTE["unbound X"]
                st["note"][ia, ic, q] = code
                if code == NOTE["error"]:
                    print(f"  error at alpha {a:.4f}, chi {c}: {p} {note}", flush=True)
            for v, lv in enumerate(LEVELS):
                st["nodes"][ia, ic, v] = lev[lv]["nodes"]
                st["match_err"][ia, ic, v] = lev[lv]["match_err"]
            st["done"][ia, ic], st["sec"][ia, ic] = True, sec
            ndone += 1
            if ndone % 10 == 0 or ndone == len(jobs):
                save(st, P)
            if ndone % 50 == 0 or ndone == len(jobs):
                el = time.time() - t0
                print(f"{int(st['done'].sum())}/{st['done'].size} points, {el:.0f}s, "
                      f"~{el / ndone * (len(jobs) - ndone) / 3600:.1f} h left; errors {int((st['note'] == 3).sum())}",
                      flush=True)
    save(st, P)


HOLD = os.path.join(HERE, "out", "relcorr_holdout.npz")
HOLD_CHI = np.array([0.3, 0.75, 0.9, 0.95, 0.99])
# 格子の alpha の 6 区間ごとの中点(幾何平均)と、以前の確かめ(chi = 0.95、alpha = 0.05、0.1、0.2、0.3)
HOLD_ALPHA = np.unique(np.concatenate([np.sqrt(ALPHA[:-1] * ALPHA[1:])[3::6], [0.05, 0.1, 0.2, 0.3]]))


def _roots_hold(c):
    return roots_column(c, HOLD_ALPHA)


def holdout(nproc):
    """格子の外の点(HOLD_ALPHA x HOLD_CHI)を解いて保存する(補間の誤差の確かめ用)。根は細かい alpha で追ってから拾う。"""
    from multiprocessing import Pool
    P, plist = procs()
    fine = np.unique(np.concatenate([ALPHA, HOLD_ALPHA]))
    with Pool(min(nproc, len(HOLD_CHI))) as pool:
        cols = pool.map(_roots_fine, [(c, fine) for c in HOLD_CHI])
    idx = np.searchsorted(fine, HOLD_ALPHA)
    jobs = []
    om_all = np.full((len(HOLD_ALPHA), len(HOLD_CHI), len(LEVELS)), np.nan + 0j)
    for ic, (c, om, ok) in enumerate(cols):
        for ia, q in enumerate(idx):
            om_all[ia, ic] = om[q]
            if ok[q].all():
                jobs.append((ia, ic, float(HOLD_ALPHA[ia]), float(c), om[q]))
    print(f"holdout: {len(jobs)} points", flush=True)
    kap = np.array([p["kappa"] for p in P])
    pw = np.array([kinn.POW[p["kind"]] for p in P])
    bh = np.array([p["kind"] != "inf" for p in P])
    Cinf = np.full((len(HOLD_ALPHA), len(HOLD_CHI), len(P)), np.nan)
    CH, tau, note = Cinf.copy(), Cinf.copy(), np.zeros(Cinf.shape, np.int8)
    with Pool(nproc) as pool:
        for ia, ic, res, lev, sec in pool.imap_unordered(job, jobs, chunksize=1):
            a, c = HOLD_ALPHA[ia], HOLD_CHI[ic]
            norm = kap * a ** (pw - 3) * np.where(bh, 1 + np.sqrt(1 - c * c), 1.0)
            for q, r in enumerate(res):
                Cinf[ia, ic, q], CH[ia, ic, q], tau[ia, ic, q] = r["G_inf"] / norm[q], r["G_h"] / norm[q], r["tau_X"]
                note[ia, ic, q] = NOTE["error"] if r["note"].startswith("error") else NOTE.get(r["note"], 0)
    np.savez(HOLD, alpha=HOLD_ALPHA, chi=HOLD_CHI, C_inf=Cinf, C_H=CH, tau_X=tau, note=note, om=om_all)


def _roots_fine(args):
    return roots_column(*args)


def m_of_X(name):
    """'411x533->322' -> m_X = m_i + m_j - m_k。"""
    return int(name[2]) + int(name[6]) - int(name[11])


def clean(T, spike=1.5):
    """補間に渡す前の表の整理。返り値: (C の全体 tot, 地平線の割合 fH, 埋めた点の印 filled, とがりの一覧)。
    (1) 放出の過程で X が束縛になった点(omegaX <= alpha。note は使わない: 2026-11-16 の版の build は自分の格子で解いた点に
        bound X をつけていなかった): Rydberg の極が密に並ぶので、極で平均した率(閾値の上の放出の率の続き)を使い、
        束縛の区間を閾値の上の最後の点の値で埋める。X の行き先は、m_X = 0 なら n00 の準位を経て地平線(fH = 1)、
        m_X != 0 なら長く生きる準位に残って BH に入らない(fH = 0、BH の勘定は無限遠へ出るのと同じ)。
    (2) 解けなかった点(閾値のすぐ上、例外)は、同じ chi の列で alpha の両隣から C を線形に埋める。
    (3) とがり: 7 点の中央値の spike 倍を越える山と、両隣の spike 倍を越える 1 点の山は、格子で解けない幅の狭い共鳴(l >= 2 の部分波)とみなし、両隣の
        とがりでない点から線形に埋めて一覧に残す(411x411->322 の alpha >~ 0.39 と 511x511->422 の alpha >~ 0.30)。
        n00 の幅の広い極(幅が 5 点以上)は 7 点の中央値の 1.2 倍以下なので触らない。率が 0 を通る谷も触らない。"""
    A = T["alpha"]
    tot = T["C_inf"] + T["C_H"]
    with np.errstate(invalid="ignore", divide="ignore"):
        fH = T["C_H"] / tot
    names, kinds = T["names"], T["kinds"]
    xbound = (T["omegaX"] <= A[:, None, None]) & (kinds == "inf")[None, None, :]
    filled = np.zeros(tot.shape, bool)
    spikes = []
    u = np.log(A)
    for q in range(tot.shape[2]):
        for ic in range(tot.shape[1]):
            t, f, b = tot[:, ic, q], fH[:, ic, q], xbound[:, ic, q]
            # (1) 束縛の区間ごとに
            ia = 0
            while ia < len(A):
                if b[ia]:
                    ib = ia
                    while ib < len(A) and b[ib]:
                        ib += 1
                    good = np.flatnonzero(np.isfinite(t[:ia]) & (t[:ia] > 0) & ~b[:ia])
                    t[ia:ib], f[ia:ib] = t[good[-1]], (1.0 if m_of_X(names[q]) == 0 else 0.0)
                    filled[ia:ib, ic, q] = True
                    ia = ib
                else:
                    ia += 1
            # (2) 解けなかった点
            bad = ~np.isfinite(t)
            if bad.any():
                t[bad] = np.interp(u[bad], u[~bad], t[~bad])
                f[bad] = np.interp(u[bad], u[~bad], f[~bad])
                filled[bad, ic, q] = True
            # (3) とがり
            med = np.full(len(A), np.inf)            # 端の 3 点は窓が片側になる(上り坂を山と取り違える)ので調べない
            med[3:-3] = [np.median(t[k - 3:k + 4]) for k in range(3, len(A) - 3)]
            peak = np.zeros(len(A), bool)               # 両隣の spike 倍を越える 1 点の山(端の近くも)
            peak[1:-1] = t[1:-1] > spike * np.maximum(t[:-2], t[2:])
            sp = (((t > spike * med) & (med > 0)) | peak) & ~filled[:, ic, q]
            if sp.any():
                for k in np.flatnonzero(sp):
                    spikes.append((str(names[q]), float(A[k]), float(T["chi"][ic]), float(t[k]), float(med[k])))
                t[sp] = np.interp(u[sp], u[~sp], t[~sp])
                f[sp] = np.interp(u[sp], u[~sp], f[~sp])
                filled[sp, ic, q] = True
    return tot, fH, filled, spikes


class Corr:
    """表の補間: corr(alpha, chi) -> (C_inf, C_H)(過程の数の長さの配列)。
    clean で整理した C(全体)と f_H(地平線の割合)を、(log alpha, chi) の 3 次で補間する。C は log でなくそのまま補間する
    (率が 0 を通る谷(振幅の符号が変わる所)は C が alpha の 2 次なのでそのままの 3 次で正しく、log では谷の両側が崩れる)。
    負になったら 0 にする。表の外: alpha < ALPHA[0] では C - 1 を alpha に比例して 0 へ、alpha > ALPHA[-1] と chi の外は
    端の値。tau_X は別に log で補間する(tau(alpha, chi)、束縛の X の過程だけ)。"""

    def __init__(self, T):
        from scipy.interpolate import RegularGridInterpolator
        A, X = T["alpha"], T["chi"]
        tot, fH, self.filled, self.spikes = clean(T)
        self.A, self.X, self.P = A, X, tot.shape[2]
        self.names = [str(s) for s in T["names"]]          # kinn.System はこれで過程を対応させる
        self.ip = RegularGridInterpolator((np.log(A), X), np.concatenate([tot, fH], axis=2), method="cubic")
        with np.errstate(invalid="ignore", divide="ignore"):
            lt = np.log(np.abs(T["tau_X"]))
        lt = np.where(np.isfinite(lt), lt, 0.0)
        self.ip_tau = RegularGridInterpolator((np.log(A), X), lt, method="linear")
        self.has_tau = np.isfinite(T["tau_X"]).any(axis=(0, 1))

    def _clip(self, alpha, chi):
        return [[np.log(min(max(alpha, self.A[0]), self.A[-1])), min(max(chi, self.X[0]), self.X[-1])]]

    def __call__(self, alpha, chi):
        v = self.ip(self._clip(alpha, chi))[0]
        tot, fH = np.maximum(v[:self.P], 0.0), np.clip(v[self.P:], -1.0, 1.0)
        if alpha < self.A[0]:
            tot = 1 + (tot - 1) * alpha / self.A[0]
        return tot * (1 - fH), tot * fH

    def tau(self, alpha, chi):
        return np.where(self.has_tau, np.exp(self.ip_tau(self._clip(alpha, chi))[0]), np.nan)


def load(path=OUT):
    z = np.load(path)
    return {k: z[k] for k in z.files}


def _loo(y, x):
    """一つ抜いて両隣の 4 点の 3 次で予測した値との差(端の 2 点は NaN)。y は x の方向が最初の軸。"""
    out = np.full_like(y, np.nan)
    for q in range(2, len(x) - 2):
        idx = [q - 2, q - 1, q + 1, q + 2]
        xs = x[idx]
        # ラグランジュの 3 次
        pred = 0.0
        for a_ in range(4):
            w = 1.0
            for b_ in range(4):
                if b_ != a_:
                    w = w * (x[q] - xs[b_]) / (xs[a_] - xs[b_])
            pred = pred + w * y[idx[a_]]
        out[q] = y[q] - pred
    return out


def check(path=OUT):
    """表の確かめ: 根と束縛状態、Wronskian と流れの一定、alpha -> 0 で C -> 1、一つ抜きの補間の誤差、点ごとの記録。"""
    T = load(path)
    A, X, names, kinds = T["alpha"], T["chi"], T["names"], T["kinds"]
    done = T["done"]
    print(f"points done {int(done.sum())}/{done.size}")
    exp_nodes = np.array([n - l - 1 for n, l, m in T["levels"]])
    d = done[..., None]
    print("roots ok:", bool(T["root_ok"].all()), "| nodes wrong:", int(((T["nodes"] != exp_nodes) & d).sum()),
          f"| max match_err {np.nanmax(np.where(d, T['match_err'], np.nan)):.1e}")
    # 根を多倍長の表(rates_mp、あれば)と比べる: 水素型からのずれを rates_mp の alpha の節から log alpha で線形に補間
    for v, lv in enumerate(LEVELS):
        f = os.path.join(HERE, "out", "rates_mp_%d%d%d.npz" % lv)
        if not os.path.exists(f):
            continue
        m = np.load(f)
        worst, nchk = 0.0, 0
        for ic, c in enumerate(X):
            jc = int(np.argmin(np.abs(m["chi"] - c)))
            if abs(m["chi"][jc] - c) > 1e-9 or not m["done"][jc]:
                continue
            ok = np.isfinite(m["re"][:, jc])
            am = m["alpha"][ok]
            res = m["re"][ok, jc] - np.array([hyd_re(a, lv, c) for a in am])
            sel = (A >= am[0]) & (A <= am[-1])
            pred = np.interp(np.log(A[sel]), np.log(am), res) + np.array([hyd_re(a, lv, c) for a in A[sel]])
            gap = np.array([spacing(a, lv[0], lv[1]) for a in A[sel]])
            if sel.any():
                worst = max(worst, float(np.max(np.abs(T["om"][sel, ic, v].real - pred) / gap)))
                nchk += int(sel.sum())
        print(f"  roots of {kinn.name(lv)} vs multiprecision table: {nchk} compared, max |d Re omega| / spacing {worst:.1e}")
    for key in ("W_spread", "J_spread", "asym_err"):
        v = np.where(done[..., None], T[key], np.nan)
        print(f"max {key} {np.nanmax(v):.1e}" if np.isfinite(v).any() else f"{key}: none")
    codes = T["note_codes"]
    for c_, nm in enumerate(codes):
        if c_ == 0:
            continue
        m = (T["note"] == c_) & done[..., None]
        if m.any():
            per = {names[q]: int(m[..., q].sum()) for q in range(len(names)) if m[..., q].any()}
            print(f"note '{nm}': {int(m.sum())} points:", per)
    C = T["C_inf"] + T["C_H"]
    if not done.all():
        return T
    # alpha -> 0: alpha <= 0.02 の点に alpha の 3 次を当てはめた切片が 1(共鳴する過程は 1 次の項がある)
    lo = A <= 0.02
    c0 = np.polynomial.polynomial.polyfit(A[lo], C[lo].reshape(lo.sum(), -1), 3)[0].reshape(C.shape[1:])
    print(f"alpha = {A[0]}: max |C - 1| = {np.nanmax(np.abs(C[0] - 1)):.2e}; intercept of cubic fit over alpha <= 0.02:"
          f" max |c0 - 1| = {np.nanmax(np.abs(c0 - 1)):.2e}")
    # 補間の誤差: 整理した表(clean)で、一つ抜きの 3 次の予測との差を、近くの 5 点の |C| の最大で割る(間隔が倍なので、
    # なめらかな所では本当の誤差の約 16 倍)。埋めた点とその両隣は数えない
    tot, fH, filled, spikes = clean(T)
    print(f"cleaned: filled {int(filled.sum())} values, spikes removed {len(spikes)} "
          f"({sorted(set(sp[0] for sp in spikes))})")
    near = filled.copy()
    near[1:] |= filled[:-1]
    near[:-1] |= filled[1:]

    def local_max(y, axis):
        y = np.moveaxis(np.abs(y), axis, 0)
        out = np.full_like(y, np.nan)
        for k in range(2, y.shape[0] - 2):
            out[k] = y[k - 2:k + 3].max(axis=0)
        return np.moveaxis(out, 0, axis)

    for axis, label, x in ((0, "alpha", np.log(A)), (1, "chi", X)):
        d = np.moveaxis(_loo(np.moveaxis(tot, axis, 0), x), 0, axis)
        e = np.abs(d) / local_max(tot, axis)
        e[near] = np.nan
        v = e[np.isfinite(e)]
        print(f"leave-one-out in {label} (error / local max of C): median {np.median(v):.1e}, 99% {np.percentile(v, 99):.1e},"
              f" 99.9% {np.percentile(v, 99.9):.1e}, max {v.max():.2e}")
        for q in np.argsort(-np.nanmax(e, axis=(0, 1)))[:6]:
            ia, ic_ = np.unravel_index(np.nanargmax(e[..., q]), e[..., q].shape)
            print(f"    {names[q]:16s} {kinds[q]:5s} max {np.nanmax(e[..., q]):.2e} at alpha {A[ia]:.4f}, chi {X[ic_]}")
    # 格子の外の点(holdout)と補間
    if os.path.exists(HOLD):
        H = np.load(HOLD)
        cr = Corr(T)
        Ct = H["C_inf"] + H["C_H"]
        err = np.full(Ct.shape, np.nan)
        for ia, a in enumerate(H["alpha"]):
            ja = np.clip(np.searchsorted(A, a), 1, len(A) - 1)
            for ic, c in enumerate(H["chi"]):
                jc = np.clip(np.searchsorted(X, c), 1, len(X) - 1)
                ci, ch = cr(a, c)
                scale = np.abs(tot[ja - 1:ja + 1, jc - 1:jc + 1]).max(axis=(0, 1))
                err[ia, ic] = np.abs(ci + ch - Ct[ia, ic]) / np.maximum(np.abs(Ct[ia, ic]), scale)
        xb = (H["om"][..., 0].real < 0)       # 使わない(形を合わせるだけ)
        skip = np.zeros(Ct.shape, bool)
        for q in range(Ct.shape[2]):           # 閾値の周りと、とがりのある過程の大きい alpha は別に見る
            if names[q] in ("411x511->322", "411x533->322", "411x544->322"):
                skip[:, :, q] = H["alpha"][:, None] > 0.33
        ok = np.isfinite(err) & ~skip
        v = err[ok]
        print(f"holdout ({int(ok.sum())} values at {len(H['alpha'])} alpha x {len(H['chi'])} chi off the grid; error /"
              f" max(|C|, nearby nodes)): median {np.median(v):.1e}, 99% {np.percentile(v, 99):.1e}, max {v.max():.2e}")
        for q in np.argsort(-np.nanmax(np.where(ok, err, np.nan), axis=(0, 1)))[:8]:
            e = np.where(ok[..., q], err[..., q], np.nan)
            ia, ic_ = np.unravel_index(np.nanargmax(e), e.shape)
            print(f"    {names[q]:16s} {kinds[q]:5s} max {np.nanmax(e):.2e} at alpha {H['alpha'][ia]:.4f}, chi {H['chi'][ic_]}"
                  f" (true {Ct[ia, ic_, q]:.4g})")
        del xb
    return T


if __name__ == "__main__":
    nproc = int(sys.argv[-1]) if sys.argv[-1].isdigit() else 4
    if sys.argv[1] == "build":
        build(nproc, "--resume" in sys.argv)
    elif sys.argv[1] == "check":
        check()
    elif sys.argv[1] == "holdout":
        holdout(nproc)
