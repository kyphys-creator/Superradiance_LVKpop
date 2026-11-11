"""基準の運動論(fp3.evolve3m)から作る chi_SR の表(E18)。

節 (alpha, s, f, L) ごとに基準の ODE を 1 本解き、スピンの道のりを「進み具合」p の格子 P_K で持つ。
- 時間は tau = mu [eV] x t [yr](率はどれも mu に比例するので、ほかの依存は L = ln(N0 M^2) を通してだけ。N0 M^2 は
  BH の角運動量 G M^2/hbar で、1 個の量子から育つ e-fold の数に入る)。
- スピンの位置 w は臨界に対する座標(chisr.w_of_chi: 322 の臨界が 0、211 の臨界が 1、0.995 が 2)。初期スピンの位置 s。
- 最初の段階の正規化 b(w): 211 から育つ節(s > 1)では w >= 1 で b = (s - w)/(s - 1)、w < 1 で b = 1 + (1 - w)。
  322 から育つ節(s < 1)では w >= 0 で b = (s - w)/s、w < 0 で b = 1 - w。どの節も b = 0 から始まり、最初の段階の
  終わり(臨界)で b = 1 になる。
- 進み具合 p = b の全変動(単調に増える)。スピンが戻る(b が減る)ところも p は増える。R = (p - b)/2 は戻りの累計。
- 各 P_K を初めて越える時刻 log10 tau と、そのときの R を持つ。届かない準位は NaN。
引き方: 両隣の節で log10 tau を(L については tau を線形に)混ぜ、log10 tau_q を挟む準位の間で p と R を線形に補間して、
b = p - 2R から w、chi に戻す。
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import chisr  # noqa: E402
import fp3  # noqa: E402
import sr  # noqa: E402

TAU_MIN, TAU_MAX = 1e-16, 1.3                 # [eV yr]。1.4e10 yr x 3e-12 eV の 30 倍まで
N_TAU = 6500
TAU_GRID = np.logspace(np.log10(TAU_MIN), np.log10(TAU_MAX), N_TAU)
P_LEV = np.round(np.arange(0.0, 4.0 + 1e-9, 0.01), 6)


def L_of_M(M):
    return np.log(sr.N0 * np.asarray(M, float) ** 2)


def b_of_w(w, s):
    """最初の段階を正規化した座標 b(配列可)。"""
    w = np.asarray(w, float)
    s = np.asarray(s, float)
    hi = s > 1
    with np.errstate(divide="ignore", invalid="ignore"):
        bh = np.where(w >= 1, (s - w) / (s - 1), 1 + (1 - w))
        bl = np.where(w >= 0, (s - w) / s, 1 - w)
    return np.where(hi, bh, bl)


def w_of_b(b, s):
    b = np.asarray(b, float)
    s = np.asarray(s, float)
    hi = s > 1
    wh = np.where(b <= 1, s - b * (s - 1), 1 - (b - 1))
    wl = np.where(b <= 1, s - b * s, -(b - 1))
    return np.where(hi, wh, wl)


def node_curve(alpha, s, f_gev, M, tab, max_nfev=1200000):
    """1 つの節の ODE。返り値: (tau の格子, chi(tau), 成功したか, 届いた最後の tau)。"""
    c1 = float(chisr.chic_vec(alpha, 1, tab=tab))
    c2 = float(chisr.chic_vec(alpha, 2, tab=tab))
    chi0 = float(chisr.chi_of_w(s, c1, c2))
    mu = alpha / sr.alpha_of(M, 1.0)
    sol = fp3.evolve3m(M, chi0, mu, f_gev, TAU_MAX / mu, TAU_GRID / mu, max_nfev=max_nfev)
    n = sol.y.shape[1]
    chi = sol.y[4] / sol.y[3] ** 2 if n else np.zeros(0)
    return TAU_GRID[:n], chi, sol.status == 0, (TAU_GRID[n - 1] if n else 0.0), c1, c2, chi0


def levels(tau, chi, s, c1, c2):
    """道のりを準位 P_LEV に直す。返り値: (log10 tau at P_K, R at P_K)。届かない準位は NaN。"""
    out_t = np.full(len(P_LEV), np.nan)
    out_r = np.full(len(P_LEV), np.nan)
    if len(tau) == 0:
        return out_t, out_r
    w = np.asarray(chisr.w_of_chi(chi, c1, c2), float)
    b = b_of_w(w, s)
    b = np.concatenate([[0.0], b])                       # tau = 0 で b = 0
    lt = np.concatenate([[-np.inf], np.log10(tau)])
    p = np.concatenate([[0.0], np.cumsum(np.abs(np.diff(b)))])
    R = (p - b) / 2
    pmax = p[-1]
    k = P_LEV <= pmax
    idx = np.searchsorted(p, P_LEV[k], side="left")     # p[idx-1] < P <= p[idx]
    idx = np.clip(idx, 1, len(p) - 1)
    p0, p1 = p[idx - 1], p[idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        fr = np.where(p1 > p0, (P_LEV[k] - p0) / (p1 - p0), 1.0)
    l0, l1 = lt[idx - 1], lt[idx]
    lt_k = np.where(np.isfinite(l0), l0 + fr * (l1 - l0), l1)
    lt_k[P_LEV[k] <= 0] = -np.inf                        # P = 0 は始めから
    out_t[k] = lt_k
    out_r[k] = R[idx - 1] + fr * (R[idx] - R[idx - 1])
    out_r[P_LEV <= 0] = 0.0
    return out_t, out_r


# ---------- 引き方 ----------
def _bracket(x, nodes):
    """格子 nodes(昇順)で x を挟む添字 i と重み t(x = (1-t) nodes[i] + t nodes[i+1])。端の外は端に止める。"""
    x = np.clip(x, nodes[0], nodes[-1])
    i = np.clip(np.searchsorted(nodes, x, side="right") - 1, 0, len(nodes) - 2)
    t = (x - nodes[i]) / (nodes[i + 1] - nodes[i])
    return i, t


def _core(T, alpha, s, f_gev, M, mu_eV, c1, c2, LT, RR, a_nodes, sig_nodes, lf_nodes, L_nodes, piece_hi,
          ST=None, TV=None, margin=10.0):
    """表 LT[ia, js, kf, iL, K](log10 tau)、RR[ia, js, kf, iL, K](R)からスピンを引く。piece_hi: s > 1 の表か。
    準位は二分法で探す(メモリは点の数に比例)。ST、TV(節の状態と、解けた最後の tau)を渡すと、重みのある節がどれも
    tau_q の margin 倍まで解けているかも返す。返り値: (chi, valid)。"""
    T, alpha, s, f_gev, M, mu_eV, c1, c2 = (np.asarray(x, float) for x in (T, alpha, s, f_gev, M, mu_eV, c1, c2))
    n = len(T)
    tau = mu_eV * T
    lf = np.log10(f_gev)
    lf_c = np.clip(lf, lf_nodes[0], lf_nodes[-1])
    sig = np.log(np.maximum(s - 1, 1e-300)) if piece_hi else np.log(np.maximum(s, 1e-300))
    ia, ta = _bracket(np.log(alpha), np.log(a_nodes))
    js, ts = _bracket(sig, sig_nodes)
    kf, tf = _bracket(lf_c, lf_nodes)
    il, tl = _bracket(L_of_M(M), L_nodes)
    corners = []
    for da in (0, 1):
        for ds in (0, 1):
            for df in (0, 1):
                wc = (ta if da else 1 - ta) * (ts if ds else 1 - ts) * (tf if df else 1 - tf)
                corners.append((ia + da, js + ds, kf + df, wc))
    valid = np.ones(n, bool)
    if ST is not None:
        for (a_, s_, f_, wc) in corners:
            for il_, wl in ((il, 1 - tl), (il + 1, tl)):
                st, tv = ST[a_, s_, f_, il_], TV[a_, s_, f_, il_]
                good = (st == 0) | ((st == 1) & (tv >= margin * tau))
                valid &= ~((wc * wl > 0) & ~good)
    rows = np.arange(n)

    def at(k, arr, logt):
        tot = np.zeros(n)
        bad = np.zeros(n, bool)
        for (a_, s_, f_, wc) in corners:
            v0 = arr[a_, s_, f_, il, k]
            v1 = arr[a_, s_, f_, il + 1, k]
            if logt:
                with np.errstate(over="ignore", invalid="ignore"):
                    v = np.where(tl <= 0, v0, np.where(tl >= 1, v1, np.log10((1 - tl) * 10.0 ** v0 + tl * 10.0 ** v1)))
            else:
                v = np.where(tl <= 0, v0, np.where(tl >= 1, v1, (1 - tl) * v0 + tl * v1))
            pos = wc > 0
            bad |= pos & np.isnan(v)
            tot += np.where(pos, wc * np.nan_to_num(v, nan=0.0, neginf=-1e30), 0.0)
        return np.where(bad, np.nan, tot)

    nK = LT.shape[-1]
    lq = np.log10(np.maximum(tau, 1e-300))
    lo = np.zeros(n, int)                                # 準位 0(P = 0)は始めから越えている
    hi = np.full(n, nK)
    while np.any(hi - lo > 1):
        mid = (lo + hi) // 2
        m = np.minimum(mid, nK - 1)
        v = at(m, LT, True)
        r = (hi - lo > 1) & (v <= lq)                   # NaN は届いていない
        lo = np.where(r, mid, lo)
        hi = np.where((hi - lo > 1) & ~r, mid, hi)
    kk = lo
    nxt = np.minimum(kk + 1, nK - 1)
    l0 = np.where(kk == 0, -np.inf, at(kk, LT, True))
    l1 = at(nxt, LT, True)
    ok1 = np.isfinite(l1) & (nxt > kk) & np.isfinite(l0) & (l1 > l0)
    with np.errstate(invalid="ignore", divide="ignore"):
        fr = np.where(ok1, (lq - l0) / (l1 - l0), 0.0)
    fr = np.clip(fr, 0, 1)
    p = P_LEV[kk] + fr * (P_LEV[nxt] - P_LEV[kk])
    r_lo = np.where(kk == 0, 0.0, at(kk, RR, False))
    r_hi = at(nxt, RR, False)
    Rq = np.where(fr > 0, r_lo + fr * (r_hi - r_lo), r_lo)
    b = p - 2 * Rq
    w = w_of_b(b, s)
    chi = np.asarray(chisr.chi_of_w(w, c1, c2), float)
    return chi, valid & np.isfinite(chi)


def lookup(T, alpha, s, f_gev, M, mu_eV, c1, c2, LT, RR, a_nodes, sig_nodes, lf_nodes, L_nodes, piece_hi):
    """試作との互換: スピンだけ返す。"""
    return _core(T, alpha, s, f_gev, M, mu_eV, c1, c2, LT, RR, a_nodes, sig_nodes, lf_nodes, L_nodes, piece_hi)[0]


def lookup_rt(rt, T, alpha, chi0, f_gev, M, mu_eV, c1, c2):
    """表(load の辞書)からスピンを引く。返り値: (chi, valid)。valid = 表の範囲の中で、重みのある節がどれも解けている。
    範囲: alpha が節の範囲、chi0 が 322 の臨界より上で CHI_TOP 以下、M が M_NODES の範囲、f が LF_NODES の範囲、
    tau = mu T <= TAU_MAX / 10。"""
    T, alpha, chi0, f_gev, M, mu_eV, c1, c2 = np.broadcast_arrays(*[np.asarray(x, float) for x in
                                                                    (T, alpha, chi0, f_gev, M, mu_eV, c1, c2)])
    s = np.asarray(chisr.w_of_chi(chi0, c1, c2), float)
    an, lfn, Mn = rt["A_NODES"], rt["LF_NODES"], rt["M_NODES"]
    lf = np.log10(f_gev)
    inside = ((alpha >= an[0]) & (alpha <= an[-1]) & (s > 0) & (chi0 <= chisr.CHI_TOP) & (M >= Mn[0]) & (M <= Mn[-1])
              & (lf >= lfn[0] - 1e-12) & (lf <= lfn[-1] + 1e-12) & (mu_eV * T <= rt["TAU_MAX"] / 10))
    chi = chi0.copy()
    valid = np.zeros(len(T), bool)
    Ln = L_of_M(Mn)
    for piece, nodes in ((0, np.log(rt["SA"])), (1, np.log(rt["SB"]))):
        q = inside & ((s > 1) if piece else (s <= 1))
        if not q.any():
            continue
        c, v = _core(T[q], alpha[q], s[q], f_gev[q], M[q], mu_eV[q], c1[q], c2[q], rt["LT"][piece], rt["RR"][piece],
                     an, nodes, lfn, Ln, bool(piece), rt["status"][piece], rt["tau_valid"][piece])
        chi[q] = c
        valid[q] = v
    return chi, valid


# ---------- 保存の形(単体版と同じ数にするため、プロジェクト版も量子化した値を使う) ----------
LT_LO, LT_HI = -17.0, 1.0                      # log10 tau の範囲
LT_STEP = (LT_HI - LT_LO) / 65533              # 2.7e-4 dex(tau で 0.06%)


def quantize_lt(LT):
    """log10 tau を uint16 に: 0 = -inf(始めの準位)、65535 = NaN(届かない)、ほかは LT_LO から LT_STEP 刻み。"""
    q = np.full(LT.shape, 65535, np.uint16)
    fin = np.isfinite(LT)
    q[fin] = (1 + np.round((np.clip(LT[fin], LT_LO, LT_HI) - LT_LO) / LT_STEP)).astype(np.uint16)
    q[np.isneginf(LT)] = 0
    return q


def dequantize_lt(q):
    q = np.asarray(q)
    v = LT_LO + (q.astype(float) - 1) * LT_STEP
    v[q == 0] = -np.inf
    v[q == 65535] = np.nan
    return v


def load(path=None, quantized=True):
    """表を読む。quantized: log10 tau を uint16、R を float16 に丸めた値を使う(単体版と同じ)。"""
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "reftab.npz")
    z = np.load(path)
    rt = {k: z[k] for k in z.files}
    if quantized:
        rt["LT"] = dequantize_lt(quantize_lt(rt["LT"]))
        rt["RR"] = rt["RR"].astype(np.float16).astype(float)
    return rt


# ---------- 第 2 の形(E18): 進み具合は b の包絡(それまでの最大値)、戻りは別に持つ ----------
# 包絡 b_env はスピンの最小値の座標なので、戻りの大きさが節で違っても準位がずれない。戻り Rz = b_env - b(>= 0)は
# log10 tau の粗い格子 G_LT で持ち、引くときは「最後に越えた準位からの経過時間」をそろえて節のあいだで混ぜる。
P_LEV2 = np.round(np.arange(0.0, 4.0 + 1e-9, 0.01), 6)
G_LT = np.round(np.arange(-16.0, 0.2 + 1e-9, 0.05), 6)


def levels2(tau, chi, s, c1, c2):
    """返り値: (log10 tau at P_LEV2(包絡が越える時刻、届かなければ NaN)、G_LT での戻り Rz、包絡の最大値)。"""
    out_t = np.full(len(P_LEV2), np.nan)
    rz = np.full(len(G_LT), np.nan)
    ok = np.isfinite(chi)
    tau, chi = tau[ok], chi[ok]
    if len(tau) == 0:
        return out_t, rz, 0.0
    w = np.asarray(chisr.w_of_chi(chi, c1, c2), float)
    b = np.concatenate([[0.0], b_of_w(w, s)])
    lt = np.concatenate([[-np.inf], np.log10(tau)])
    env = np.maximum.accumulate(b)
    k = P_LEV2 <= env[-1]
    idx = np.clip(np.searchsorted(env, P_LEV2[k], side="left"), 1, len(env) - 1)   # env[idx-1] < P <= env[idx]
    e0, e1 = env[idx - 1], env[idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        fr = np.where(e1 > e0, (P_LEV2[k] - e0) / (e1 - e0), 1.0)
    l0, l1 = lt[idx - 1], lt[idx]
    out_t[k] = np.where(np.isfinite(l0), l0 + fr * (l1 - l0), l1)
    out_t[0] = -np.inf
    g = G_LT <= lt[-1]
    rz[g] = np.interp(G_LT[g], lt[1:], (env - b)[1:], left=0.0)
    return out_t, rz, float(env[-1])


def _core2(T, alpha, s, f_gev, M, mu_eV, c1, c2, LT, RG, a_nodes, sig_nodes, lf_nodes, L_nodes, piece_hi):
    """第 2 の形の引き方。返り値: (chi, b_env, Rz)。"""
    T, alpha, s, f_gev, M, mu_eV, c1, c2 = (np.asarray(x, float) for x in (T, alpha, s, f_gev, M, mu_eV, c1, c2))
    n = len(T)
    lq = np.log10(np.maximum(mu_eV * T, 1e-300))
    lf_c = np.clip(np.log10(f_gev), lf_nodes[0], lf_nodes[-1])
    sig = np.log(np.maximum(s - 1, 1e-300)) if piece_hi else np.log(np.maximum(s, 1e-300))
    ia, ta = _bracket(np.log(alpha), np.log(a_nodes))
    js, ts = _bracket(sig, sig_nodes)
    kf, tf = _bracket(lf_c, lf_nodes)
    il, tl = _bracket(L_of_M(M), L_nodes)
    corners = []
    for da in (0, 1):
        for ds in (0, 1):
            for df in (0, 1):
                wc = (ta if da else 1 - ta) * (ts if ds else 1 - ts) * (tf if df else 1 - tf)
                corners.append((ia + da, js + ds, kf + df, wc))

    def lt_corner(a_, s_, f_, k):
        """節の log10 tau(L について tau を線形に混ぜる)。"""
        v0, v1 = LT[a_, s_, f_, il, k], LT[a_, s_, f_, il + 1, k]
        with np.errstate(over="ignore", invalid="ignore"):
            return np.where(tl <= 0, v0, np.where(tl >= 1, v1, np.log10((1 - tl) * 10.0 ** v0 + tl * 10.0 ** v1)))

    def at(k):
        tot = np.zeros(n)
        bad = np.zeros(n, bool)
        for (a_, s_, f_, wc) in corners:
            v = lt_corner(a_, s_, f_, k)
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
    benv = P_LEV2[kk] + fr * (P_LEV2[nxt] - P_LEV2[kk])
    # 戻り: 最後に越えた準位 kk からの経過時間 d をそろえ、各節の Rz をその時刻で引く(次の準位を越える時刻で止める)
    d = np.where(np.isfinite(l0), lq - l0, 0.0)
    rz = np.zeros(n)
    for (a_, s_, f_, wc) in corners:
        pos = (wc > 0) & (kk > 0)
        if not pos.any():
            continue
        for il_, wl in ((il, 1 - tl), (il + 1, tl)):
            tk = LT[a_, s_, f_, il_, kk]
            tn = LT[a_, s_, f_, il_, nxt]
            tc = tk + d
            tc = np.where(np.isfinite(tn) & (nxt > kk), np.minimum(tc, tn), tc)
            j = np.clip((tc - G_LT[0]) / (G_LT[1] - G_LT[0]), 0, len(G_LT) - 1.000001)
            j0 = np.floor(j).astype(int)
            u = j - j0
            r0 = RG[a_, s_, f_, il_, j0]
            r1 = RG[a_, s_, f_, il_, np.minimum(j0 + 1, len(G_LT) - 1)]
            rv = np.nan_to_num((1 - u) * r0 + u * r1, nan=0.0)
            rz += np.where(pos, wc * wl * rv, 0.0)
    b = benv - rz
    w = w_of_b(b, s)
    return np.asarray(chisr.chi_of_w(w, c1, c2), float), benv, rz


def levels3(tau, chi, s, c1, c2):
    """第 3 の形: 進み具合は「最後に下へ通る」包絡 b_lp(t) = min_{t' >= t} b(t')(スピンでは未来の最大値)。
    一時的な行き過ぎ(戻ってくる落ち込み)は包絡に入らないので、落ち込みの深さが節で違っても準位がそろう。
    落ち込み D = b - b_lp(>= 0)は G_LT で持ち、_core2 には RG = -D として渡す(b = b_lp + D)。
    返り値: (log10 tau at P_LEV2、G_LT での -D、包絡の最大値)。"""
    out_t = np.full(len(P_LEV2), np.nan)
    rz = np.full(len(G_LT), np.nan)
    ok = np.isfinite(chi)
    tau, chi = tau[ok], chi[ok]
    if len(tau) == 0:
        return out_t, rz, 0.0
    w = np.asarray(chisr.w_of_chi(chi, c1, c2), float)
    b = np.concatenate([[0.0], b_of_w(w, s)])
    lt = np.concatenate([[-np.inf], np.log10(tau)])
    env = np.minimum.accumulate(b[::-1])[::-1]
    env = np.maximum.accumulate(np.maximum(env, 0.0))         # 数値の揺れを除いて単調に
    k = P_LEV2 <= env[-1]
    idx = np.clip(np.searchsorted(env, P_LEV2[k], side="left"), 1, len(env) - 1)
    e0, e1 = env[idx - 1], env[idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        fr = np.where(e1 > e0, (P_LEV2[k] - e0) / (e1 - e0), 1.0)
    l0, l1 = lt[idx - 1], lt[idx]
    out_t[k] = np.where(np.isfinite(l0), l0 + fr * (l1 - l0), l1)
    out_t[0] = -np.inf
    g = G_LT <= lt[-1]
    rz[g] = np.interp(G_LT[g], lt[1:], (env - b)[1:], left=0.0)
    return out_t, rz, float(env[-1])


# ---------- 第 3 の形の表を使う(E18) ----------
KI_BLOCK = 4                                   # 曲がりの指標は準位 4 つずつの最大で持つ
KI_STEP = 0.01                                 # [dex]、uint8(255 = 2.55 dex 以上)


def compute_ki(LT, axes):
    """曲がりの指標 KI[ia, js, kf, iL, block]: 各軸で、節の log10 tau が両隣を結ぶ直線からずれる量 [dex] の最大
    (端の節は内側の隣の値)。準位 KI_BLOCK 個ずつの最大を uint8 で返す。axes = (x_alpha, x_s, x_f)(補間に使う座標)。"""
    lt = np.where(np.isfinite(LT), LT, np.nan).astype(float)
    ki = np.zeros(LT.shape, float)
    for ax, x in enumerate(axes):
        n = lt.shape[ax]
        if n < 3:
            continue
        sl = lambda a, b: tuple([slice(None)] * ax + [slice(a, b)])
        t1, t2, t3 = lt[sl(0, n - 2)], lt[sl(1, n - 1)], lt[sl(2, n)]
        shape = [1] * lt.ndim
        shape[ax] = n - 2
        x1, x2, x3 = (np.asarray(x[a:b]).reshape(shape) for a, b in ((0, n - 2), (1, n - 1), (2, n)))
        with np.errstate(invalid="ignore"):
            dev = np.abs(t2 - (t1 + (t3 - t1) * (x2 - x1) / (x3 - x1)))
        dev = np.nan_to_num(dev, nan=0.0)
        full = np.concatenate([dev[sl(0, 1)], dev, dev[sl(n - 3, n - 2)]], axis=ax)
        # 真ん中の節のずれは、その節を角に持つ両側の区間に効くので、隣どうしの最大もとる
        nb = np.maximum(full, np.concatenate([full[sl(1, n)], full[sl(n - 1, n)]], axis=ax))
        nb = np.maximum(nb, np.concatenate([full[sl(0, 1)], full[sl(0, n - 1)]], axis=ax))
        ki = np.maximum(ki, nb)
    K = ki.shape[-1]
    nb_ = -(-K // KI_BLOCK)
    pad = np.zeros(ki.shape[:-1] + (nb_ * KI_BLOCK - K,))
    kb = np.concatenate([ki, pad], axis=-1).reshape(ki.shape[:-1] + (nb_, KI_BLOCK)).max(-1)
    return np.minimum(np.round(kb / KI_STEP), 255).astype(np.uint8)


def load3(path=None, quantized=True):
    """第 3 の形の表を読む(log10 tau は uint16、戻りは float16 に丸めた値。単体版と同じ数になる)。曲がりの指標 KI を足す。"""
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "reftab.npz")
    z = np.load(path)
    rt = {k: z[k] for k in z.files}
    LT = rt["LT"].astype(float)
    RG = rt["RG"].astype(float)
    if quantized:
        LT = dequantize_lt(quantize_lt(LT))
        RG = RG.astype(np.float16).astype(float)
    rt["LT"], rt["RG"] = LT, RG
    if "KI" not in rt:
        rt["KI"] = np.stack([compute_ki(LT[pc], (np.log(rt["A_NODES"]), np.log(rt["SA"] if pc == 0 else rt["SB"]),
                                                rt["LF_NODES"])) for pc in (0, 1)])
    return rt


def _lt_corner(LT, a_, s_, f_, il, tl, k):
    """節 (a_, s_, f_) の準位 k の log10 tau(M について tau を線形に混ぜる)。"""
    v0, v1 = LT[a_, s_, f_, il, k], LT[a_, s_, f_, il + 1, k]
    with np.errstate(over="ignore", invalid="ignore"):
        return np.where(tl <= 0, v0, np.where(tl >= 1, v1, np.log10((1 - tl) * 10.0 ** v0 + tl * 10.0 ** v1)))


def _corner_b(LT, RG, a_, s_, f_, il, tl, lq):
    """一つの節(M は混ぜる)の、時刻 lq での b と準位の番号。"""
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
    benv = P_LEV2[lo] + fr * (P_LEV2[nxt] - P_LEV2[lo])
    ilr = np.where(tl < 0.5, il, il + 1)
    j = np.clip((lq - G_LT[0]) / (G_LT[1] - G_LT[0]), 0, len(G_LT) - 1.000001)
    j0 = np.floor(j).astype(int)
    u = j - j0
    rv = (1 - u) * RG[a_, s_, f_, ilr, j0] + u * RG[a_, s_, f_, ilr, np.minimum(j0 + 1, len(G_LT) - 1)]
    return benv - np.nan_to_num(rv, nan=0.0), lo


def lookup_rt3(rt, T, alpha, chi0, f_gev, M, mu_eV, c1, c2, theta_v=0.06, theta_k=0.05, margin=10.0):
    """表からスピンを引く。返り値: (chi, valid, sharp)。
    valid: 表の範囲の中で、重みのある節がどれも tau_q の margin 倍まで解けている。
    sharp: 速い段差の時刻が節のあいだで急に変わる所(印 5)。角の節それぞれの、問いの時刻での値の幅が theta_v を越え、
           しかもその時刻に角の節がいる準位で、段差の時刻が両隣を結ぶ直線から theta_k [dex] 以上ずれている。"""
    T, alpha, chi0, f_gev, M, mu_eV, c1, c2 = np.broadcast_arrays(*[np.asarray(x, float) for x in
                                                                    (T, alpha, chi0, f_gev, M, mu_eV, c1, c2)])
    s = np.asarray(chisr.w_of_chi(chi0, c1, c2), float)
    an, lfn, Mn = rt["A_NODES"], rt["LF_NODES"], rt["M_NODES"]
    lf = np.log10(f_gev)
    tau = mu_eV * T
    inside = ((alpha >= an[0]) & (alpha <= an[-1]) & (s > 0) & (chi0 <= chisr.CHI_TOP) & (M >= Mn[0]) & (M <= Mn[-1])
              & (lf >= lfn[0] - 1e-12) & (lf <= lfn[-1] + 1e-12) & (tau <= rt["TAU_MAX"] / margin))
    chi = chi0.copy()
    valid = np.zeros(len(T), bool)
    sharp = np.zeros(len(T), bool)
    Ln = L_of_M(Mn)
    for piece, snodes in ((0, rt["SA"]), (1, rt["SB"])):
        q = np.flatnonzero(inside & ((s > 1) if piece else (s <= 1)))
        if len(q) == 0:
            continue
        LT, RG, KI = rt["LT"][piece], rt["RG"][piece], rt["KI"][piece]
        ST, TV = rt["status"][piece], rt["tau_valid"][piece]
        sq = s[q]
        sig = np.log(np.maximum(sq - 1, 1e-300)) if piece else np.log(np.maximum(sq, 1e-300))
        c, _, _ = _core2(T[q], alpha[q], sq, f_gev[q], M[q], mu_eV[q], c1[q], c2[q], LT, RG, an, np.log(snodes), lfn, Ln,
                         bool(piece))
        ia, ta = _bracket(np.log(alpha[q]), np.log(an))
        js, ts = _bracket(sig, np.log(snodes))
        kf, tf = _bracket(lf[q], lfn)
        il, tl = _bracket(L_of_M(M[q]), Ln)
        lq = np.log10(np.maximum(tau[q], 1e-300))
        ok = np.ones(len(q), bool)
        vmin = np.full(len(q), np.inf)
        vmax = np.full(len(q), -np.inf)
        kmin = np.full(len(q), 10 ** 9)
        kmax = np.zeros(len(q), int)
        corners = []
        for da in (0, 1):
            for ds in (0, 1):
                for df in (0, 1):
                    wc = (ta if da else 1 - ta) * (ts if ds else 1 - ts) * (tf if df else 1 - tf)
                    a_, s_, f_ = ia + da, js + ds, kf + df
                    pos = wc > 0
                    for il_, wl in ((il, 1 - tl), (il + 1, tl)):
                        st, tv = ST[a_, s_, f_, il_], TV[a_, s_, f_, il_]
                        good = (st == 0) | ((st == 1) & (tv >= margin * tau[q]))
                        ok &= ~((wc * wl > 0) & ~good)
                    bc, kc = _corner_b(LT, RG, a_, s_, f_, il, tl, lq)
                    xc = np.asarray(chisr.chi_of_w(w_of_b(bc, sq), c1[q], c2[q]), float)
                    vmin = np.where(pos, np.minimum(vmin, xc), vmin)
                    vmax = np.where(pos, np.maximum(vmax, xc), vmax)
                    kmin = np.where(pos, np.minimum(kmin, kc), kmin)
                    kmax = np.where(pos, np.maximum(kmax, kc), kmax)
                    corners.append((a_, s_, f_, pos))
        b0 = np.clip(kmin // KI_BLOCK, 0, KI.shape[-1] - 1)
        b1 = np.clip((kmax + 1) // KI_BLOCK, 0, KI.shape[-1] - 1)
        kmaxv = np.zeros(len(q))
        nblk = int(np.max(b1 - b0)) + 1 if len(q) else 0
        for (a_, s_, f_, pos) in corners:
            for il_ in (il, il + 1):
                for d in range(nblk):
                    bb = np.minimum(b0 + d, b1)
                    v = KI[a_, s_, f_, il_, bb] * KI_STEP
                    kmaxv = np.where(pos, np.maximum(kmaxv, v), kmaxv)
        chi[q] = c
        valid[q] = ok & np.isfinite(c)
        sharp[q] = ((vmax - vmin) > theta_v) & (kmaxv > theta_k)
    return chi, valid, sharp
