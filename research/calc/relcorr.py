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


class Corr:
    """表の補間: corr(alpha, chi) -> (C_inf, C_H)(過程の数の長さの配列)。
    補間するのは log(C_inf + C_H)(全体)と f_H = C_H/(C_inf + C_H)(地平線の割合)を (log alpha, chi) の 3 次で。
    表の外: alpha < ALPHA[0] では C - 1 を alpha に比例して 0 へ、alpha > ALPHA[-1] と chi の外は端の値。
    解けなかった点(閾値のすぐ上、例外)は、同じ chi の列で alpha の両隣から log C を線形に埋める(filled に記録)。
    放出の過程が閾値を越えて X が束縛になった点(note = bound X)は、X が雲から逃げられず地平線に吸われる。Rydberg の極が
    密に並ぶので、極で平均した率(閾値の上の放出の率の続き)を使う: 全体は閾値の上の最後の点の値で一定、f_H = 1。"""

    def __init__(self, T):
        from scipy.interpolate import RegularGridInterpolator
        A, X = T["alpha"], T["chi"]
        tot = T["C_inf"] + T["C_H"]
        with np.errstate(invalid="ignore", divide="ignore"):
            fH = T["C_H"] / tot
        # 放出の過程で X が束縛になった点(note でなく omegaX から判定する。2026-11-16 の版の build は、自分の格子で解いた
        # 点に bound X の記録をつけていなかった)
        xbound = (T["omegaX"] <= A[:, None, None]) & (T["kinds"] == "inf")[None, None, :]
        self.filled = np.zeros(tot.shape, bool)
        for q in range(tot.shape[2]):
            for ic in range(len(X)):
                t, f = tot[:, ic, q], fH[:, ic, q]
                bound = xbound[:, ic, q]
                if bound.any():
                    ia0 = np.flatnonzero(bound)[0]
                    last = np.flatnonzero(np.isfinite(t[:ia0]) & (t[:ia0] > 0))[-1]
                    t[ia0:], f[ia0:] = t[last], 1.0
                    self.filled[ia0:, ic, q] = True
                bad = ~np.isfinite(t) | (t <= 0)
                if bad.any():
                    good = ~bad
                    t[bad] = np.exp(np.interp(np.log(A[bad]), np.log(A[good]), np.log(t[good])))
                    f[bad] = np.interp(np.log(A[bad]), np.log(A[good]), f[good])
                    self.filled[bad, ic, q] = True
        self.A, self.X = A, X
        self.y = np.log(tot)
        self.fH = fH
        vals = np.concatenate([self.y, fH], axis=2)
        self.ip = RegularGridInterpolator((np.log(A), X), vals, method="cubic")
        self.P = tot.shape[2]

    def __call__(self, alpha, chi):
        a = min(max(alpha, self.A[0]), self.A[-1])
        c = min(max(chi, self.X[0]), self.X[-1])
        v = self.ip([[np.log(a), c]])[0]
        tot, fH = np.exp(v[:self.P]), v[self.P:]
        if alpha < self.A[0]:
            tot = 1 + (tot - 1) * alpha / self.A[0]
        return tot * (1 - fH), tot * fH


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
    # alpha -> 0
    lo = slice(0, 4)
    if done[lo].all():
        dev = np.abs(C[0] - 1)
        q_, c_ = np.unravel_index(np.nanargmax(dev.T), dev.T.shape)
        print(f"alpha = {A[0]}: max |C - 1| = {np.nanmax(dev):.2e} ({names[q_]}, chi {X[c_]})")
        ic = np.polynomial.polynomial.polyfit(A[lo], C[lo].reshape(4, -1), 2)[0]
        print(f"intercept of quadratic fit over alpha <= {A[3]:.4f}: max |c0 - 1| = {np.nanmax(np.abs(ic - 1)):.2e}")
    # 一つ抜きの補間(alpha は log alpha、chi はそのまま。log C で)
    with np.errstate(invalid="ignore", divide="ignore"):
        y = np.log(np.abs(C))
    full_a = done.all(axis=1)
    if full_a.sum() >= 5:
        ya = np.where(done[..., None], y, np.nan)
        la = _loo(ya, np.log(A))
        print("leave-one-out in alpha (|d log C|, spacing doubled, so a bound ~16x the real error):")
        for q in np.argsort(-np.nanmax(np.abs(la), axis=(0, 1)))[:12]:
            v = np.abs(la[..., q])
            ia, ic_ = np.unravel_index(np.nanargmax(v), v.shape)
            print(f"  {names[q]:16s} {kinds[q]:5s} max {np.nanmax(v):.2e} at alpha {A[ia]:.4f}, chi {X[ic_]};"
                  f" median {np.nanmedian(v):.1e}")
    if done.all(axis=0).sum() >= 5 or done.all():
        lc = _loo(np.moveaxis(np.where(done[..., None], y, np.nan), 1, 0), X)
        print("leave-one-out in chi:")
        for q in np.argsort(-np.nanmax(np.abs(lc), axis=(0, 1)))[:12]:
            v = np.abs(lc[..., q])
            ic_, ia = np.unravel_index(np.nanargmax(v), v.shape)
            print(f"  {names[q]:16s} {kinds[q]:5s} max {np.nanmax(v):.2e} at alpha {A[ia]:.4f}, chi {X[ic_]};"
                  f" median {np.nanmedian(v):.1e}")
    return T


if __name__ == "__main__":
    nproc = int(sys.argv[-1]) if sys.argv[-1].isdigit() else 4
    if sys.argv[1] == "build":
        build(nproc, "--resume" in sys.argv)
    elif sys.argv[1] == "check":
        check()
    elif sys.argv[1] == "holdout":
        holdout(nproc)
