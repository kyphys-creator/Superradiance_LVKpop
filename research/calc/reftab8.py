"""N 準位の基準の運動論から作る chi_SR の新しい表(DESIGN.md D18 の手順 4・5、2026-11-17)。

reftab.py(E18、三準位、今の表)の第 3 の形を、8 準位の運動論(kinn.evolveN_checked)に広げたもの。今の表(reftab.py、
chisr.py)は作り直しが終わるまでそのまま使うので、別のモジュールにした。違い:
- スピンの座標 w は 0.995・211・322・433・544 の臨界を節に持つ区分的な一次式(w = 2, 1, 0, -1, -2)。臨界は多倍長の率
  (sr.gamma_mp_over_mu)が 0 になるスピン(質量固定)を細かい alpha の格子で表にしたもの(crit)。433 と 533 の臨界は
  差が 0.0014 以下なので m = 3 は 433 で代表する。544 の臨界より下は w = -2 - (c4 - chi)/c4(chi = 0 で -3)。
- 初期スピン chi0 の位置 w0 で区分 p を決める: w0 が (p, p+1] にあれば p(1: 211 から育つ、0: 322 から、-1: 433 から、
  -2: 544 から)。区分の中の位置 sigma = w0 - p(0, 1]。
- 最初の段階の正規化: w >= p で b = (w0 - w)/sigma、w < p で b = 1 + (p - w)。どの区分も b = 0 から始まり、最初の段階の
  終わり(すぐ下の臨界)で b = 1、その下の臨界を一つ越えるごとに +1。進み具合は reftab の第 3 の形と同じく、最後に下へ
  通る包絡で持ち(P8: 0-5 を 0.01 刻み)、包絡からの一時的な落ち込みは G_LT の時刻で持つ。
- 節: alpha 0.008-0.5 の対数で 47 点、sigma 25 点(区分ごと)、f 1e8-1e20.5 GeV を 0.25 dex で 51 点、M 1・300 Msun
  (M は L = ln(N0 M^2) を通してだけ効く)。区分 1 は 211 の臨界が 0.995 より下の alpha だけ、区分 -1 は alpha >= 0.08、
  -2 は alpha >= 0.18(それより小さい alpha では 433・544 が値の表の年齢 mu T <= 0.13 eV yr のうちに育たない)。
- 節ごとの運動論は kinn.evolveN_checked(8 準位、相対論的な補正、多倍長の率、計算の上限と平均の確かめ)。有効な年齢までの
  曲線から準位を作り、status(0: 全年齢、1: 途中まで)と tau_valid を持つ。
- 引き方は reftab._core2・lookup_rt3 と同じ骨組み(角の 8 節で準位の時刻を混ぜる、M は tau を線形に混ぜる)。違うのは
  補間が止まった所(重みのある角の節のどれかが次の準位に届かない): 全部の角の節が届いている最後の準位で止め、それを自分の
  時刻ですでに越えた角の節の超過分を重みで足す(reftab の版は進みのいちばん遅い節の値に引かれた。E19)。
試し(2026-11-17、trial・trial2、1 本の上限 4e6 回): 16 セル・256 節を解き、セルの中の 48 点を直接解いた値と約束の物差し
(時刻のずれ 2% を許した |Delta chi|)で比べて、節が有効な年齢で最大 0.0039(433 から育つ区分の降りる所でも 0.0035)。
544 から育つ区分の変化する所と、印 5 の判定(曲がり KI・折れ KS、全体の表が要る)はまだ確かめていない(DESIGN.md D21)。
実行: cd research && ~/miniforge3/bin/python3 calc/reftab8.py crit | trial | trial2 | trial_eval | trial2_eval [並列の数]
本番(Hoffman2、h2/run.py の build8): reftab8.py plan <場所> <見積もりの場所> で解く節の一覧を作り、ジョブ配列の各本が
build_task で受け持ちを解く(1 本に複数のコア、節ごとに保存し直す)。最後に reftab8.py merge <場所> ... <表のファイル> で
区分ごとの配列と曲がり KI・折れ KS に束ねる(load_table で lookup8 の形に読む)。
"""
import json
import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reftab as R  # noqa: E402
import sr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
CHI_TOP = 0.995
LEV_CRIT = {1: (2, 1, 1), 2: (3, 2, 2), 3: (4, 3, 3), 4: (5, 4, 4)}
CRIT_PATH = os.path.join(OUT, "reftab8_crit.npz")
CRIT_ALPHAS = np.geomspace(0.003, 0.52, 3000)
P8 = np.round(np.arange(0.0, 5.0 + 1e-9, 0.01), 6)
G_LT = R.G_LT
TAU_MIN, TAU_MAX, N_TAU, TAU_GRID = R.TAU_MIN, R.TAU_MAX, R.N_TAU, R.TAU_GRID


def _fine(coarse):
    c = np.asarray(coarse, float)
    return np.sort(np.concatenate([c, np.sqrt(c[1:] * c[:-1])]))


A_NODES = np.geomspace(0.008, 0.5, 47)
_S_COARSE = [1e-3, 3e-3, 0.01, 0.03, 0.07, 0.15, 0.25, 0.4, 0.55, 0.7, 0.85, 0.95]
S_NODES = {1: _fine(_S_COARSE + [1.0]), 0: _fine(_S_COARSE + [0.999]), -1: _fine(_S_COARSE + [0.999]),
           -2: _fine(_S_COARSE + [0.999])}
LF_NODES = np.round(np.arange(8.0, 20.51, 0.25), 6)
M_NODES = np.array([1.0, 300.0])
PIECES = (1, 0, -1, -2)
A_MIN = {1: 0.0, 0: 0.0, -1: 0.08, -2: 0.18}


# ---------- 臨界スピン ----------
def _analytic_crit(alpha, n, m):
    ap = alpha * (1 - alpha ** 2 / (2 * n ** 2))
    return np.where(ap / m < 0.5, 4 * m * ap / (m ** 2 + 4 * ap ** 2), 1.0)


def _root(alpha, lev):
    """準位 lev の多倍長の率が 0 になるスピン(質量固定)。0.995 まで超放射しなければ 0.9999 まで外挿で探し、それでも
    なければ解析式(1 以下)。"""
    from scipy.optimize import brentq
    g = lambda c: sr.gamma_mp_over_mu(alpha, c, *lev)   # noqa: E731
    cs = np.concatenate([np.geomspace(1e-5, 0.05, 60, endpoint=False), np.linspace(0.05, CHI_TOP, 200)])
    v = np.array([g(c) for c in cs])
    k = np.flatnonzero((v[:-1] <= 0) & (v[1:] > 0))
    if len(k):
        return brentq(g, cs[k[0]], cs[k[0] + 1], xtol=1e-13)
    if v[0] > 0:
        return float(_analytic_crit(alpha, lev[0], lev[2]))
    if g(0.9999) > 0:
        return brentq(g, CHI_TOP, 0.9999, xtol=1e-13)
    return max(float(_analytic_crit(alpha, lev[0], lev[2])), 0.9999)


def build_crit():
    """臨界の表(CRIT_ALPHAS x m = 1-4)を作って保存する。"""
    tab = np.array([[_root(a, LEV_CRIT[m]) for m in (1, 2, 3, 4)] for a in CRIT_ALPHAS])
    if not np.all(np.diff(tab[:, ::-1], axis=1) > 0):
        raise ValueError("臨界が c4 < c3 < c2 < c1 の順に並ばない alpha がある")
    np.savez(CRIT_PATH, alpha=CRIT_ALPHAS, c=tab)
    return tab


_CRIT = None


def crit(alpha):
    """臨界スピン (c1, c2, c3, c4)(211・322・433・544、質量固定)。表の範囲で log alpha について線形に補間する。"""
    global _CRIT
    if _CRIT is None:
        z = np.load(CRIT_PATH)
        _CRIT = (np.log(z["alpha"]), z["c"])
    la, c = _CRIT
    x = np.log(np.asarray(alpha, float))
    if np.any(x < la[0] - 1e-12) or np.any(x > la[-1] + 1e-12):
        raise ValueError("alpha が臨界の表(0.003-0.52)の外")
    return tuple(np.interp(x, la, c[:, m]) for m in range(4))


# ---------- 座標 ----------
def w_of_chi(chi, cs):
    c1, c2, c3, c4 = (np.asarray(x, float) for x in cs)
    chi = np.asarray(chi, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.select([chi >= c1, chi >= c2, chi >= c3, chi >= c4],
                         [1 + (chi - c1) / (CHI_TOP - c1), (chi - c2) / (c1 - c2), -1 + (chi - c3) / (c2 - c3),
                          -2 + (chi - c4) / (c3 - c4)], default=-2 - (c4 - chi) / c4)


def chi_of_w(w, cs):
    c1, c2, c3, c4 = (np.asarray(x, float) for x in cs)
    w = np.asarray(w, float)
    return np.select([w >= 1, w >= 0, w >= -1, w >= -2],
                     [c1 + (w - 1) * (CHI_TOP - c1), c2 + w * (c1 - c2), c3 + (w + 1) * (c2 - c3),
                      c4 + (w + 2) * (c3 - c4)], default=c4 * (1 + (w + 2)))


def piece_of(w0):
    """区分(w0 が (p, p+1] なら p)。-2 以下(544 の臨界以下、8 準位はどれも超放射しない)と 2 を越える所は -9。"""
    w0 = np.asarray(w0, float)
    p = np.ceil(w0) - 1
    return np.where((w0 > -2) & (w0 <= 2), p, -9).astype(int)


def b_of_w(w, p, sig):
    w = np.asarray(w, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(w >= p, (p + sig - w) / sig, 1 + (p - w))


def w_of_b(b, p, sig):
    b = np.asarray(b, float)
    return np.where(b <= 1, p + sig - b * sig, p - (b - 1))


# ---------- 節の曲線から準位へ ----------
def levels8(tau, chi, p, sig, cs):
    """reftab.levels3 と同じ形(最後に下へ通る包絡 b_lp(t) = min_{t' >= t} b(t') の各準位 P8 を越える log10 tau と、
    G_LT の時刻での包絡からの落ち込み D の符号を変えたもの)。返り値: (LT, RZ, 包絡の最大値)。"""
    out_t = np.full(len(P8), np.nan)
    rz = np.full(len(G_LT), np.nan)
    ok = np.isfinite(chi)
    tau, chi = tau[ok], chi[ok]
    if len(tau) == 0:
        return out_t, rz, 0.0
    b = np.concatenate([[0.0], b_of_w(w_of_chi(chi, cs), p, sig)])
    lt = np.concatenate([[-np.inf], np.log10(tau)])
    env = np.minimum.accumulate(b[::-1])[::-1]
    env = np.maximum.accumulate(np.maximum(env, 0.0))
    k = P8 <= env[-1]
    idx = np.clip(np.searchsorted(env, P8[k], side="left"), 1, len(env) - 1)
    e0, e1 = env[idx - 1], env[idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        fr = np.where(e1 > e0, (P8[k] - e0) / (e1 - e0), 1.0)
    l0, l1 = lt[idx - 1], lt[idx]
    out_t[k] = np.where(np.isfinite(l0), l0 + fr * (l1 - l0), l1)
    out_t[0] = -np.inf
    g = G_LT <= lt[-1]
    rz[g] = np.interp(G_LT[g], lt[1:], (env - b)[1:], left=0.0)
    return out_t, rz, float(env[-1])


# ---------- 節を解く ----------
_SYS = None


def _system():
    global _SYS
    if _SYS is None:
        import kinn
        import relcorr
        _SYS = kinn.System(kinn.LEVELS8, kinn.process_table(), corr=relcorr.Corr(relcorr.load()))
    return _SYS


def node_params(p, ia, js, kf, il):
    a = float(A_NODES[ia])
    sig = float(S_NODES[p][js])
    cs = tuple(float(x) for x in crit(a))
    return a, sig, 10.0 ** LF_NODES[kf], float(M_NODES[il]), cs, float(chi_of_w(p + sig, cs))


def solve_node(args):
    """一つの節: (p, ia, js, kf, il, max_nfev)。返り値は準位と確かめの量。"""
    warnings.simplefilter("ignore")
    import kinn
    p, ia, js, kf, il, cap = args
    a, sig, f, M, cs, chi0 = node_params(p, ia, js, kf, il)
    mu = a / sr.alpha_of(M, 1.0)
    t0 = time.time()
    b = kinn.evolveN_checked(_system(), M, chi0, mu, f, TAU_MAX / mu, TAU_GRID / mu, max_nfev=cap)
    nv = b.n_valid
    LT, RZ, env = levels8(TAU_GRID[:nv], kinn.chi_of(b)[:nv], p, sig, cs)
    return dict(key=(p, ia, js, kf, il), LT=LT.astype(np.float32), RZ=RZ.astype(np.float32), env=env,
                status=0 if nv == N_TAU else 1, tau_valid=float(TAU_GRID[nv - 1]) if nv else 0.0, reason=b.reason,
                runs=b.runs, chi0=chi0, secs=time.time() - t0)


def solve_node_safe(args):
    """solve_node が例外を出したら、失敗した節(status 2、準位は NaN)として返す(本番の 1 本が一つの節で止まって、
    投げ直しても先へ進めなくなるのを防ぐ。引き方は status 0・1 の節しか使わない)。"""
    t0 = time.time()
    try:
        return solve_node(args)
    except Exception as e:                       # noqa: BLE001(どの例外でも、その節だけ失敗にする)
        return dict(key=tuple(args[:5]), LT=np.full(len(P8), np.nan, np.float32), RZ=np.full(len(G_LT), np.nan, np.float32),
                    env=0.0, status=2, tau_valid=0.0, reason=("error: " + repr(e))[:200], runs=[], chi0=float("nan"),
                    secs=time.time() - t0)


# ---------- 引き方 ----------
def _lt_corner(LT, a_, s_, f_, il, tl, k):
    v0, v1 = LT[a_, s_, f_, il, k], LT[a_, s_, f_, il + 1, k]
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        return np.where(tl <= 0, v0, np.where(tl >= 1, v1, np.log10((1 - tl) * 10.0 ** v0 + tl * 10.0 ** v1)))


def _corner_b(LT, RZ, a_, s_, f_, il, tl, lq):
    """一つの角の節(M は混ぜる)の、時刻 lq での b と準位の番号(reftab._corner_b を P8 で)。"""
    n = len(lq)
    nK = LT.shape[-1]
    lo = np.zeros(n, int)
    hi = np.full(n, nK)
    while np.any(hi - lo > 1):
        mid = (lo + hi) // 2
        v = _lt_corner(LT, a_, s_, f_, il, tl, np.minimum(mid, nK - 1))
        r = (hi - lo > 1) & (np.nan_to_num(v, nan=np.inf) <= lq)
        lo = np.where(r, mid, lo)
        hi = np.where((hi - lo > 1) & ~r, mid, hi)
    nxt = np.minimum(lo + 1, nK - 1)
    l0 = np.where(lo == 0, -np.inf, _lt_corner(LT, a_, s_, f_, il, tl, lo))
    l1 = _lt_corner(LT, a_, s_, f_, il, tl, nxt)
    ok1 = np.isfinite(l1) & (nxt > lo) & np.isfinite(l0) & (l1 > l0)
    with np.errstate(invalid="ignore", divide="ignore"):
        fr = np.clip(np.where(ok1, (lq - l0) / (l1 - l0), 0.0), 0, 1)
    benv = P8[lo] + fr * (P8[nxt] - P8[lo])
    ilr = np.where(tl < 0.5, il, il + 1)
    j = np.nan_to_num(np.clip((lq - G_LT[0]) / (G_LT[1] - G_LT[0]), 0, len(G_LT) - 1.000001), nan=0.0)
    j0 = np.floor(j).astype(int)
    u = j - j0
    rv = (1 - u) * RZ[a_, s_, f_, ilr, j0] + u * RZ[a_, s_, f_, ilr, np.minimum(j0 + 1, len(G_LT) - 1)]
    return benv - np.nan_to_num(rv, nan=0.0), lo


def core8(T, alpha, sig, p, f_gev, M, mu_eV, cs, LT, RZ, a_nodes, s_nodes, lf_nodes, L_nodes):
    """reftab._core2 の 4 区分版。返り値: (chi, blocked, 角の節の b の重みつき平均)。blocked の点の chi は重みつき平均の値。"""
    n = len(T)
    lq = np.log10(np.maximum(mu_eV * T, 1e-300))
    lf_c = np.clip(np.log10(f_gev), lf_nodes[0], lf_nodes[-1])
    ia, ta = R._bracket(np.log(alpha), np.log(a_nodes))
    js, ts = R._bracket(np.log(np.maximum(sig, 1e-300)), np.log(s_nodes))
    kf, tf = R._bracket(lf_c, lf_nodes)
    il, tl = R._bracket(R.L_of_M(M), L_nodes)
    corners = []
    for da in (0, 1):
        for ds in (0, 1):
            for df in (0, 1):
                wc = (ta if da else 1 - ta) * (ts if ds else 1 - ts) * (tf if df else 1 - tf)
                corners.append((ia + da, js + ds, kf + df, wc))

    def at(k):
        tot = np.zeros(n)
        bad = np.zeros(n, bool)
        for (a_, s_, f_, wc) in corners:
            v = _lt_corner(LT, a_, s_, f_, il, tl, k)
            pos = wc > 0
            bad |= pos & np.isnan(v)
            tot += np.where(pos, wc * np.nan_to_num(v, nan=0.0, neginf=-1e30), 0.0)
        return np.where(bad, np.nan, tot)

    nK = LT.shape[-1]
    lo = np.zeros(n, int)
    hi = np.full(n, nK)
    while np.any(hi - lo > 1):
        mid = (lo + hi) // 2
        v = at(np.minimum(mid, nK - 1))
        r = (hi - lo > 1) & (v <= lq)
        lo = np.where(r, mid, lo)
        hi = np.where((hi - lo > 1) & ~r, mid, hi)
    kk = lo
    nxt = np.minimum(kk + 1, nK - 1)
    l0 = np.where(kk == 0, -np.inf, at(kk))
    l1 = at(nxt)
    ok1 = np.isfinite(l1) & (nxt > kk) & np.isfinite(l0) & (l1 > l0)
    with np.errstate(invalid="ignore", divide="ignore"):
        fr = np.clip(np.where(ok1, (lq - l0) / (l1 - l0), 0.0), 0, 1)
    benv = P8[kk] + fr * (P8[nxt] - P8[kk])
    d = np.where(np.isfinite(l0), lq - l0, 0.0)
    rz = np.zeros(n)
    bw = np.zeros(n)
    ex = np.zeros(n)
    for (a_, s_, f_, wc) in corners:
        pos = wc > 0
        bc, _ = _corner_b(LT, RZ, a_, s_, f_, il, tl, lq)
        bw += np.where(pos, wc * bc, 0.0)
        ex += np.where(pos, wc * np.maximum(bc - P8[kk], 0.0), 0.0)
        pk = pos & (kk > 0)
        if not pk.any():
            continue
        for il_, wl in ((il, 1 - tl), (il + 1, tl)):
            tk = LT[a_, s_, f_, il_, kk]
            tn = LT[a_, s_, f_, il_, nxt]
            tc = tk + d
            tc = np.where(np.isfinite(tn) & (nxt > kk), np.minimum(tc, tn), tc)
            j = np.nan_to_num(np.clip((tc - G_LT[0]) / (G_LT[1] - G_LT[0]), 0, len(G_LT) - 1.000001), nan=0.0)
            j0 = np.floor(j).astype(int)
            u = j - j0
            r0 = RZ[a_, s_, f_, il_, j0]
            r1 = RZ[a_, s_, f_, il_, np.minimum(j0 + 1, len(G_LT) - 1)]
            rz += np.where(pk, wc * wl * np.nan_to_num((1 - u) * r0 + u * r1, nan=0.0), 0.0)
    # 補間が止まった所(重みのある角の節のどれかが次の準位に届かない。最後に止まる位置が節で違うときなど): 全部の角の節が
    # 届いている最後の準位 kk で止め、kk を自分の時刻ですでに越えた角の節の超過分だけを重みで足す(段差の途中の角の節は
    # kk に置く)。値は戻らず、全部の角の節が越えたあとは角の節の b の重みつき平均につながる。角の節の b をそのまま平均すると、
    # 段差の前の節と後の節を混ぜて中途半端な値になる(試しで 0.12 ずれた、2026-11-17)
    blocked = np.isnan(l1) & (nxt > kk)
    b = np.where(blocked, P8[kk] + ex, benv - rz)
    return np.asarray(chi_of_w(w_of_b(b, p, sig), cs), float), blocked, bw


def lookup8(tab, T, alpha, chi0, f_gev, M, mu_eV, margin=10.0, theta_v=0.06, theta_p=R.THETA_P):
    """表からスピンを引く。tab: dict(p -> dict(LT, RZ, status, tau_valid, a_nodes, s_nodes, lf_nodes))。
    返り値: (chi, valid, piece, spread, blocked)。valid は、表の範囲の中で重みのある節がどれも tau_q の margin 倍まで
    解けている点。spread は角の節の値の幅(印 5 の判定に使う。曲がり KI・折れ KS の印は全体の表で足す)。"""
    T, alpha, chi0, f_gev, M, mu_eV = np.broadcast_arrays(*[np.asarray(x, float) for x in (T, alpha, chi0, f_gev, M, mu_eV)])
    cs = crit(alpha)
    w0 = w_of_chi(chi0, cs)
    pp = piece_of(w0)
    sig = w0 - pp
    tau = mu_eV * T
    chi = chi0.copy()
    valid = np.zeros(len(T), bool)
    spread = np.zeros(len(T))
    blocked = np.zeros(len(T), bool)
    Ln = R.L_of_M(M_NODES)
    for p, t in tab.items():
        an, sn, lfn = t["a_nodes"], t["s_nodes"], t["lf_nodes"]
        lf = np.log10(f_gev)
        q = np.flatnonzero((pp == p) & (alpha >= an[0]) & (alpha <= an[-1]) & (sig >= sn[0]) & (sig <= sn[-1])
                           & (M >= M_NODES[0]) & (M <= M_NODES[-1]) & (lf >= lfn[0] - 1e-12) & (lf <= lfn[-1] + 1e-12)
                           & (tau <= TAU_MAX / margin))
        if len(q) == 0:
            continue
        csq = tuple(c[q] for c in cs)
        c, blk, _ = core8(T[q], alpha[q], sig[q], p, f_gev[q], M[q], mu_eV[q], csq, t["LT"], t["RZ"], an, sn, lfn, Ln)
        ia, ta = R._bracket(np.log(alpha[q]), np.log(an))
        js, ts = R._bracket(np.log(sig[q]), np.log(sn))
        kf, tf = R._bracket(np.clip(lf[q], lfn[0], lfn[-1]), lfn)
        il, tl = R._bracket(R.L_of_M(M[q]), Ln)
        lq = np.log10(np.maximum(tau[q], 1e-300))
        ok = np.ones(len(q), bool)
        vmin = np.full(len(q), np.inf)
        vmax = np.full(len(q), -np.inf)
        for da in (0, 1):
            for ds in (0, 1):
                for df in (0, 1):
                    wc = (ta if da else 1 - ta) * (ts if ds else 1 - ts) * (tf if df else 1 - tf)
                    a_, s_, f_ = ia + da, js + ds, kf + df
                    for il_, wl in ((il, 1 - tl), (il + 1, tl)):
                        st, tv = t["status"][a_, s_, f_, il_], t["tau_valid"][a_, s_, f_, il_]
                        good = (st == 0) | ((st == 1) & (tv >= margin * tau[q]))
                        ok &= ~((wc * wl > 0) & ~good)
                    bc, _ = _corner_b(t["LT"], t["RZ"], a_, s_, f_, il, tl, lq)
                    xc = np.asarray(chi_of_w(w_of_b(bc, p, sig[q]), csq), float)
                    pos = wc > 0
                    vmin = np.where(pos, np.minimum(vmin, xc), vmin)
                    vmax = np.where(pos, np.maximum(vmax, xc), vmax)
        chi[q] = c
        valid[q] = ok & np.isfinite(c)
        spread[q] = vmax - vmin
        blocked[q] = blk
    return chi, valid, pp, spread, blocked


# ---------- 全節と、費用の見積もり(Hoffman2 のジョブ配列) ----------
def all_nodes():
    """作る節の一覧 (p, ia, js, kf, il)。区分 1 は 211 の臨界が 0.995 より下の alpha、区分 -1・-2 は A_MIN 以上の alpha。"""
    out = []
    for p in PIECES:
        for ia, a in enumerate(A_NODES):
            if a < A_MIN[p] or (p == 1 and float(crit(a)[0]) >= CHI_TOP):
                continue
            for js in range(len(S_NODES[p])):
                for kf in range(len(LF_NODES)):
                    for il in range(len(M_NODES)):
                        out.append((p, ia, js, kf, il))
    return out


def sample_nodes(n, seed=2026):
    """全節から n 個を無作為に(重複なし)。"""
    allv = all_nodes()
    idx = np.sort(np.random.default_rng(seed).choice(len(allv), n, replace=False))
    return [allv[i] for i in idx]


def cost_task(tid, ntask, out, n=2000, cap=20000000, seed=2026):
    """費用の見積もりの 1 本: sample_nodes(n) のうち番号 tid(1 から)の受け持ちを解き、out/cost8_<tid>.npz に保存する
    (準位も保存するので、上限を変えなければ本番で使い回せる)。"""
    import platform
    path = os.path.join(out, f"cost8_{tid:05d}.npz")
    if os.path.exists(path):
        print("already done:", path)
        return
    nodes = sample_nodes(n, seed)[tid - 1::ntask]
    res = [solve_node(k + (cap,)) for k in nodes]
    tmp = path + ".tmp.npz"
    np.savez(tmp, key=np.array([r["key"] for r in res]), LT=np.array([r["LT"] for r in res]), RZ=np.array([r["RZ"] for r in res]),
             env=np.array([r["env"] for r in res]), status=np.array([r["status"] for r in res]),
             tau_valid=np.array([r["tau_valid"] for r in res]), secs=np.array([r["secs"] for r in res]),
             nfev=np.array([[x["nfev"] for x in r["runs"]] + [0] * (2 - len(r["runs"])) for r in res]),
             njump=np.array([[x["njump"] for x in r["runs"]] + [0] * (2 - len(r["runs"])) for r in res]),
             reason=np.array([r["reason"] for r in res]), chi0=np.array([r["chi0"] for r in res]), cap=cap,
             host=platform.node(), n_all=len(all_nodes()))
    os.replace(tmp, path)
    print(f"{len(res)} nodes -> {path}, {sum(r['secs'] for r in res):.0f}s")


# ---------- 本番の作り直し(Hoffman2 のジョブ配列。1 本に複数のコア、節ごとに保存し直す) ----------
BUILD_SEED = 7
PIECE_TAG = {1: "p1", 0: "p0", -1: "m1", -2: "m2"}       # 保存するときの区分の名前


def _pack(res, cap):
    """節の結果のリストを保存用の配列に。"""
    import platform
    return dict(key=np.array([r["key"] for r in res], int).reshape(-1, 5), LT=np.array([r["LT"] for r in res], np.float32),
                RZ=np.array([r["RZ"] for r in res], np.float32), env=np.array([r["env"] for r in res]),
                status=np.array([r["status"] for r in res], np.int8), tau_valid=np.array([r["tau_valid"] for r in res]),
                secs=np.array([r["secs"] for r in res]),
                nfev=np.array([[x["nfev"] for x in r["runs"]] + [0] * (2 - len(r["runs"])) for r in res], np.int64).reshape(-1, 2),
                njump=np.array([[x["njump"] for x in r["runs"]] + [0] * (2 - len(r["runs"])) for r in res], np.int64).reshape(-1, 2),
                reason=np.array([r["reason"] for r in res]), chi0=np.array([r["chi0"] for r in res]), cap=cap,
                host=platform.node())


def _unpack(path):
    """保存した結果を {節: 結果の dict} に。"""
    z = np.load(path)
    out = {}
    for i, k in enumerate(z["key"]):
        out[tuple(int(x) for x in k)] = dict(key=tuple(int(x) for x in k), LT=z["LT"][i], RZ=z["RZ"][i], env=float(z["env"][i]),
                                             status=int(z["status"][i]), tau_valid=float(z["tau_valid"][i]), secs=float(z["secs"][i]),
                                             runs=[dict(nfev=int(n), njump=int(j)) for n, j in zip(z["nfev"][i], z["njump"][i]) if n > 0],
                                             reason=str(z["reason"][i]), chi0=float(z["chi0"][i]), cap=int(z["cap"]))
    return out


def _result_files(dirs):
    import glob
    fs = []
    for d in dirs:
        fs += glob.glob(os.path.join(d, "cost8_*.npz")) + glob.glob(os.path.join(d, "build8_*.npz"))
    return sorted(f for f in fs if not f.endswith(".tmp.npz"))


def plan_build(out, skip_dirs=()):
    """本番で解く節の一覧を作って out/build8_nodes.npy に保存する(投げる前に一度だけ。全部の本が同じ一覧を読む)。
    skip_dirs の結果(費用の見積もりなど)で解いた節は除き、BUILD_SEED で並びを混ぜる(重い節を本ごとに散らす)。
    skip_dirs のジョブが全部終わってから作る(途中で作ると、終わっていない本の節をもう一度解く。結果は壊れない)。"""
    done = set()
    for f in _result_files(skip_dirs):
        done.update(tuple(int(x) for x in k) for k in np.load(f)["key"])
    todo = np.array([k for k in all_nodes() if k not in done], int)
    np.random.default_rng(BUILD_SEED).shuffle(todo)
    np.save(os.path.join(out, "build8_nodes.npy"), todo)
    print(f"{len(todo)} nodes to solve ({len(done)} already solved) -> {os.path.join(out, 'build8_nodes.npy')}")
    return todo


def build_task(tid, ntask, out, cap, nproc=1, every=50, every_s=900.0):
    """本番の 1 本: build8_nodes.npy のうち番号 tid(1 から)の受け持ちを nproc 並列で解き、out/build8_<tid>.npz に every 節か
    every_s 秒ごとに保存し直す。保存したものがあれば、その続きから解く(時間の上限で切られても投げ直せば続く)。
    終了の合図(SIGTERM。Slurm が時間の上限の前に送る)を受けたら、その場で保存してやめる。
    投げ直すときは本数 ntask を変えない(受け持ちが本数で決まる。変えても束ねるときに重複は除かれ、壊れはしない)。
    節が例外を出したら、その節だけ失敗(status 2)として記録して先へ進む(solve_node_safe)。"""
    import signal
    from multiprocessing import Pool
    nodes = [tuple(int(x) for x in k) for k in np.load(os.path.join(out, "build8_nodes.npy"))[tid - 1::ntask]]
    path = os.path.join(out, f"build8_{tid:05d}.npz")
    have = _unpack(path) if os.path.exists(path) else {}
    todo = [k for k in nodes if k not in have]
    print(f"task {tid}/{ntask}: {len(nodes)} nodes, {len(have)} already saved, {len(todo)} to solve, {nproc} workers, cap {cap}",
          flush=True)

    def save():
        tmp = path + ".tmp.npz"
        np.savez(tmp, **_pack(list(have.values()), cap))
        os.replace(tmp, path)

    # 終了の合図は印を立てるだけにし、結果を 5 秒ずつ区切って待って印を見る(結果を待つあいだ親は止まっていて、
    # macOS ではその待ちが合図で中断されず、例外を投げる形では合図が効かなかった。2026-11-17 の試し)
    term = []
    signal.signal(signal.SIGTERM, lambda signum, frame: term.append(signum))
    if todo:
        import multiprocessing as mp
        last, n0, ndone = time.time(), len(have), 0
        pool = Pool(nproc, initializer=signal.signal, initargs=(signal.SIGTERM, signal.SIG_DFL))
        it = pool.imap_unordered(solve_node_safe, [k + (cap,) for k in todo], chunksize=1)
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
            print(f"task {tid}: SIGTERM, saved {len(have)}/{len(nodes)} nodes", flush=True)
            return
    print(f"task {tid}: {len(have)}/{len(nodes)} nodes saved", flush=True)


def merge_build(dirs, out_path):
    """結果(cost8_*.npz と build8_*.npz)を区分ごとの配列に束ね、曲がり KI と折れ KS(reftab の関数)を足して out_path に保存する。
    配列: LT_<区分>[alpha, sigma, f, M, 準位](float32、NaN は届かない)、RZ_<区分>、status_<区分>(-1 は解いていない節)、
    tau_valid_<区分>、KI_<区分>、KS_<区分>。区分の名前は PIECE_TAG(p1、p0、m1、m2)。status は 0(全年齢)、1(途中まで)、
    2(例外で失敗)、-1(解いていない)。全節の一覧にない結果は数えて飛ばす。出どころ(provenance)も保存する。"""
    res = {}
    caps = set()
    for f in _result_files(dirs):
        z = np.load(f)
        caps.add(int(z["cap"]))
        for i, k in enumerate(z["key"]):
            res[tuple(int(x) for x in k)] = (z["LT"][i], z["RZ"][i], int(z["status"][i]), float(z["tau_valid"][i]),
                                             float(z["secs"][i]), int(z["nfev"][i].sum()))
    want = all_nodes()
    wset = set(want)
    extra = [k for k in res if k not in wset]          # 全節の一覧にない結果(数えて飛ばす)
    for k in extra:
        del res[k]
    missing = [k for k in want if k not in res]
    nfail = sum(1 for v in res.values() if v[2] == 2)
    arrs = dict(A_NODES=A_NODES, LF_NODES=LF_NODES, M_NODES=M_NODES, P8=P8, G_LT=G_LT, TAU_MAX=TAU_MAX,
                caps=np.array(sorted(caps)), n_nodes=len(want), n_missing=len(missing), n_failed=nfail, n_extra=len(extra))
    z = np.load(CRIT_PATH)
    arrs.update(crit_alpha=z["alpha"], crit_c=z["c"], provenance=json.dumps(provenance(), ensure_ascii=False))
    for p in PIECES:
        tg = PIECE_TAG[p]
        shp = (len(A_NODES), len(S_NODES[p]), len(LF_NODES), len(M_NODES))
        LT = np.full(shp + (len(P8),), np.nan, np.float32)
        RZ = np.full(shp + (len(G_LT),), np.nan, np.float32)
        st = np.full(shp, -1, np.int8)
        tv = np.zeros(shp, np.float32)
        secs = np.zeros(shp, np.float32)
        nfev = np.zeros(shp, np.int64)
        for k, (lt, rz, s_, t_, sc, nf) in res.items():
            if k[0] != p:
                continue
            idx = k[1:]
            LT[idx], RZ[idx], st[idx], tv[idx], secs[idx], nfev[idx] = lt, rz, s_, t_, sc, nf
        arrs.update({f"S_NODES_{tg}": S_NODES[p], f"LT_{tg}": LT, f"RZ_{tg}": RZ, f"status_{tg}": st, f"tau_valid_{tg}": tv,
                     f"secs_{tg}": secs, f"nfev_{tg}": nfev,
                     f"KI_{tg}": R.compute_ki(LT, (np.log(A_NODES), np.log(S_NODES[p]), LF_NODES)),
                     f"KS_{tg}": R.compute_ks(LT, np.log(A_NODES))})
    tmp = out_path + ".tmp.npz"
    np.savez(tmp, **arrs)
    os.replace(tmp, out_path)
    print(f"merged {len(res)} nodes ({len(missing)} of {len(want)} missing, {nfail} failed, {len(extra)} not in the node list "
          f"skipped, caps {sorted(caps)}) -> {out_path}")
    return len(missing)


def provenance():
    """表の出どころ: コードのコミット(と手元の変更の有無)、入力のチェックサム、運動論の設定。"""
    import hashlib
    import inspect
    import subprocess
    import kinn
    def sha(path):
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for blk in iter(lambda: fh.read(1 << 20), b""):
                h.update(blk)
        return h.hexdigest()
    try:
        commit = subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", HERE, "status", "--porcelain", "--untracked-files=no", "--", "."],
                                    capture_output=True, text=True).stdout.strip())
    except Exception:                            # noqa: BLE001
        commit, dirty = "unknown", None
    inputs = [os.path.join(OUT, f) for f in ["relcorr.npz", "kinn_procs.json", "reftab8_crit.npz"]
              + [f"rates_mp_{n}{l}{m}.npz" for (n, l, m) in kinn.LEVELS8]]
    sig = {k: v.default for k, v in inspect.signature(kinn.evolveN_checked).parameters.items() if v.default is not inspect._empty}
    siga = {k: v.default for k, v in inspect.signature(kinn.evolveN_avg).parameters.items() if v.default is not inspect._empty}
    import scipy
    return dict(commit=commit, dirty=dirty, inputs={os.path.basename(f): sha(f) for f in inputs if os.path.exists(f)},
                evolveN_checked=sig, evolveN_avg={k: v for k, v in siga.items() if k != "trace"}, AVG_FINE=kinn.AVG_FINE,
                numpy=np.__version__, scipy=scipy.__version__, levels=[list(l) for l in kinn.LEVELS8])


def load_table(path):
    """merge_build の表を lookup8 の形 {区分: dict(LT, RZ, status, tau_valid, a_nodes, s_nodes, lf_nodes, KI, KS)} に。"""
    z = np.load(path)
    return {p: dict(LT=z[f"LT_{PIECE_TAG[p]}"], RZ=z[f"RZ_{PIECE_TAG[p]}"], status=z[f"status_{PIECE_TAG[p]}"],
                    tau_valid=z[f"tau_valid_{PIECE_TAG[p]}"], a_nodes=z["A_NODES"], s_nodes=z[f"S_NODES_{PIECE_TAG[p]}"],
                    lf_nodes=z["LF_NODES"], KI=z[f"KI_{PIECE_TAG[p]}"], KS=z[f"KS_{PIECE_TAG[p]}"]) for p in PIECES}


# ---------- 試し: 少数のセルで引き方を確かめる ----------
# セル = 区分 p と、alpha・sigma・f の隣り合う 2 節(M は 2 節)。16 節を解き、セルの中の点を直接解いた値と比べる。
TRIAL_CELLS = [(1, 0.05, 0.5, 12.0), (1, 0.20, 0.3, 15.0), (1, 0.10, 0.5, 20.0), (0, 0.10, 0.5, 14.0), (0, 0.35, 0.5, 16.0),
               (0, 0.15, 0.002, 13.0), (-1, 0.30, 0.5, 14.0), (-1, 0.40, 0.5, 16.0), (-1, 0.45, 0.3, 13.0),
               (-2, 0.35, 0.5, 15.0), (-2, 0.45, 0.5, 17.0), (0, 0.25, 0.1, 18.0)]


# 433・544 から育つ区分で、飽和が速く節が最後まで解ける所(f が大きめ)。最初の試しの区分 -1・-2 のセルは、試しの上限で
# 節が早く切れたか 544 が育たず、スピンが変わる年齢まで有効にならなかった
TRIAL_CELLS2 = [(-1, 0.30, 0.5, 16.0), (-1, 0.35, 0.3, 17.5), (-1, 0.42, 0.6, 19.0), (-2, 0.45, 0.5, 18.5)]


def trial_cells(which=1):
    out = []
    for p, a, s, lf in (TRIAL_CELLS if which == 1 else TRIAL_CELLS2):
        ia = int(np.clip(np.searchsorted(A_NODES, a) - 1, 0, len(A_NODES) - 2))
        js = int(np.clip(np.searchsorted(S_NODES[p], s) - 1, 0, len(S_NODES[p]) - 2))
        kf = int(np.clip(np.searchsorted(LF_NODES, lf) - 1, 0, len(LF_NODES) - 2))
        out.append((p, ia, js, kf))
    return out


def _trial_ref(args):
    """セルの中の点を直接解く(evolveN_checked)。"""
    warnings.simplefilter("ignore")
    import kinn
    p, a, sig, f, M, cap = args
    cs = tuple(float(x) for x in crit(a))
    chi0 = float(chi_of_w(p + sig, cs))
    mu = a / sr.alpha_of(M, 1.0)
    taus = np.geomspace(1e-14, TAU_MAX / 10, 400)
    b = kinn.evolveN_checked(_system(), M, chi0, mu, f, taus[-1] / mu, taus / mu, max_nfev=cap)
    return dict(args=args[:5], chi0=chi0, taus=taus[:b.n_valid].tolist(), chi=kinn.chi_of(b)[:b.n_valid].tolist(),
                reason=b.reason)


def trial(nproc=12, cap=4000000, npts=3, seed=17, which=1):
    from multiprocessing import Pool
    cells = trial_cells(which)
    jobs = sorted({(p, ia + da, js + ds, kf + df, il, cap) for (p, ia, js, kf) in cells
                   for da in (0, 1) for ds in (0, 1) for df in (0, 1) for il in (0, 1)})
    rng = np.random.default_rng(seed)
    refs = []
    for (p, ia, js, kf) in cells:
        for _ in range(npts):
            a = float(np.exp(rng.uniform(np.log(A_NODES[ia]), np.log(A_NODES[ia + 1]))))
            sig = float(np.exp(rng.uniform(np.log(S_NODES[p][js]), np.log(S_NODES[p][js + 1]))))
            lf = float(rng.uniform(LF_NODES[kf], LF_NODES[kf + 1]))
            M = float(np.exp(rng.uniform(np.log(2.0), np.log(200.0))))
            refs.append((p, a, sig, 10.0 ** lf, M, cap))
    t0 = time.time()
    with Pool(nproc) as pool:
        rr = pool.map_async(_trial_ref, refs, chunksize=1)
        nodes = pool.map(solve_node, jobs, chunksize=1)
        refres = rr.get()
    print(f"trial: {len(jobs)} nodes and {len(refs)} reference points in {time.time() - t0:.0f}s", flush=True)
    import pickle
    sfx = "" if which == 1 else str(which)
    pickle.dump(dict(nodes=nodes, refs=refres), open(os.path.join(OUT, f"reftab8_trial{sfx}_nodes.pkl"), "wb"))
    return trial_eval(nodes, refres, which)


def trial_eval(nodes=None, refres=None, which=1):
    """試しの節と参照の点(保存したもの)で、引き方を確かめ直す(節を解き直さない)。"""
    import pickle
    sfx = "" if which == 1 else str(which)
    if nodes is None:
        z = pickle.load(open(os.path.join(OUT, f"reftab8_trial{sfx}_nodes.pkl"), "rb"))
        nodes, refres = z["nodes"], z["refs"]
    cells = trial_cells(which)
    nd = {r["key"]: r for r in nodes}
    rows = []
    for (p, ia, js, kf) in cells:
        sub = dict(a_nodes=A_NODES[ia:ia + 2], s_nodes=S_NODES[p][js:js + 2], lf_nodes=LF_NODES[kf:kf + 2])
        LT = np.full((2, 2, 2, 2, len(P8)), np.nan, np.float32)
        RZ = np.full((2, 2, 2, 2, len(G_LT)), np.nan, np.float32)
        st = np.zeros((2, 2, 2, 2), np.int8)
        tv = np.zeros((2, 2, 2, 2))
        for da in (0, 1):
            for ds in (0, 1):
                for df in (0, 1):
                    for il in (0, 1):
                        r = nd[(p, ia + da, js + ds, kf + df, il)]
                        LT[da, ds, df, il], RZ[da, ds, df, il] = r["LT"], r["RZ"]
                        st[da, ds, df, il], tv[da, ds, df, il] = r["status"], r["tau_valid"]
        sub.update(LT=LT, RZ=RZ, status=st, tau_valid=tv)
        for ref in refres:
            q, a, sig, f, M = ref["args"]
            if q != p or not (sub["a_nodes"][0] <= a <= sub["a_nodes"][1] and sub["s_nodes"][0] <= sig <= sub["s_nodes"][1]
                              and sub["lf_nodes"][0] <= np.log10(f) <= sub["lf_nodes"][1]):
                continue
            taus = np.array(ref["taus"])
            if len(taus) < 2:
                rows.append(dict(cell=(p, ia, js, kf), point=ref["args"], n=0))
                continue
            mu = a / sr.alpha_of(M, 1.0)
            T = taus / mu
            c, valid, pc, spread, blk = lookup8({p: sub}, T, np.full(len(T), a), np.full(len(T), ref["chi0"]),
                                                np.full(len(T), f), np.full(len(T), M), np.full(len(T), mu))
            import kinn
            d = kinn._tol_dist(taus, c, np.array(ref["chi"]), 0.02)
            m = valid
            rows.append(dict(cell=(p, ia, js, kf), point=[p, a, sig, f, M], n=int(m.sum()), n_ages=len(T),
                             max_err=float(d[m].max()) if m.any() else None, worst_tau=float(taus[m][np.argmax(d[m])]) if m.any() else None,
                             blocked=int((blk & m).sum()), spread_max=float(spread[m].max()) if m.any() else None,
                             chi0=ref["chi0"], chi_end=float(ref["chi"][-1]), reason=ref["reason"]))
    secs = [r["secs"] for r in nodes]
    summ = dict(nodes=len(nodes), node_secs=dict(median=float(np.median(secs)), p90=float(np.percentile(secs, 90)),
                                                 max=float(np.max(secs)), total_h=float(np.sum(secs) / 3600)),
                nodes_partial=int(sum(r["status"] for r in nodes)), rows=rows)
    json.dump(summ, open(os.path.join(OUT, f"reftab8_trial{sfx}.json"), "w"), indent=1, default=float)
    return summ


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "crit"
    nproc = int(sys.argv[-1]) if sys.argv[-1].isdigit() else 12
    if what == "crit":
        t0 = time.time()
        tab = build_crit()
        for a in (0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5):
            print(a, [round(float(x), 4) for x in crit(a)])
        print(f"{len(CRIT_ALPHAS)} alpha, {time.time() - t0:.0f}s")
    elif what == "plan":                 # reftab8.py plan <出力の場所> [解いた結果の場所 ...]
        plan_build(sys.argv[2], sys.argv[3:])
    elif what == "merge":                # reftab8.py merge <結果の場所> [...] <出力のファイル>
        merge_build(sys.argv[2:-1], sys.argv[-1])
    elif what in ("trial", "trial_eval", "trial2", "trial2_eval"):
        which = 2 if what.startswith("trial2") else 1
        s = trial(nproc, which=which) if not what.endswith("eval") else trial_eval(which=which)
        print(json.dumps({k: v for k, v in s.items() if k != "rows"}, indent=1))
        for r in s["rows"]:
            print(r)
