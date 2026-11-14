"""chi_SR: 自己相互作用のある ALP の超放射で、年齢 T のあとの BH のスピン(単体版、numpy と scipy だけ)。

使い方:
    from chi_sr import ChiSR
    model = ChiSR("chi_sr_tables.npz")
    chi_T, regime = model(M_msun, chi0, T_yr, mu_eV, f_GeV)   # どれも配列でよい(ブロードキャスト)
M_msun は生まれたときの BH の質量(超放射で減る前)。観測された質量を渡すと、f >= 1e17 GeV で最大 0.04 ずれる(README)。

regime(印):
    0 = 超放射しない(chi0 が 211 と 322 の臨界以下で、433 も年齢のうちに飽和しない)。chi0 を返す
    1 = 自己相互作用の固定点(211・322・544 の三準位)。211 から育つ BH も 322 から育つ BH も、322 の臨界まで
        下がりうる。三準位の外の準位(411、433)は年齢のうちに育たない(422 は調べていない。README の限界 2)
    2 = 自由場(211 が飽和して臨界へ、そのあと 322 の段階)、または 322 から育つ BH で自己相互作用が飽和までに効かない場合
    3 = この模型の外。chi0 を返す(予言ではない): 固定点で alpha > 0.30、値の表の下(alpha < 0.008。固定点では
        alpha < 0.01 で f < 1e8 GeV の点も)で雲が年齢のうちに育ちうる点、alpha > 0.5、211・322 は超放射しないが
        433(m = 3)が年齢のうちに飽和しうる、表の範囲の中なのに値を引けない(M が 1--300 Msun の外、mu T > 0.13 eV yr、
        重みのある節がその年齢まで解けていない)、chi0 > 0.995 で超放射する
    4 = 値は返すが、模型の外の準位が年齢のうちに育ちうるので近似は粗い: 固定点で 411、433、211 の臨界より上の
        544 が育ちうる、自由場で止まったスピンから 433 が飽和しうる、または自由場で 544 の超放射まで進んで
        433 の臨界より下にいる。約束の外
    5 = 値は表の補間だが、表の節のあいだで補間できない所。約束の外: 速い段差の時刻が節のあいだで急に変わる所(E18)、
        段差の時刻が alpha の節のあいだで折れていて値がそれに敏感な所、段のあとの平らな所で角の節の値が割れている所(E19)

約束(照合した範囲は M 2--200 Msun、mu 1e-13--1e-11 eV、f 1e6--1e22 GeV、T 1e5--1.4e10 yr、alpha <= 0.5):
  印 0・1・2 の値は、基準の運動論(三準位、BH の質量と角運動量の変化、重力波の過程)と
  |chi_SR(T) - chi_ref(T')| <= 0.03 を満たす T' が |ln(T'/T)| <= 0.02 の中にある。
  これは標本で確かめた統計的な主張である(照合の点の数、越えたセルの数、最大の誤差は README.md)。

模型の要点(詳しくは README.md):
  値: 超放射する点の値は、基準の運動論を alpha 0.008--0.5、f 1e8--10^20.5 GeV、M 1、300 Msun の節で解いた表から引く
  (進み具合の準位ごとの時刻を節のあいだで補間する)。f < 1e8 GeV は以下の半解析の値。alpha < 0.008 は、雲が育つ前の
  点だけ chi0 を返す。
  固定点は Baryakhtar+ 2011.11646 の式 34--36 に 544 を加えた三準位(Table V・VI の率)。211 の臨界より上は (211, 322)、
  下は 211 が空の (322, 544)。
  dchi/dt = mu (f/M_pl)^2 H(alpha, chi) なので、I(alpha, chi) = int_chi^0.995 dchi'/|H| を表にし、
  I(alpha, chi_T) = I(alpha, chi_0') + mu (f/M_pl)^2 (T - t_grow) を chi_T について解く(両隣の alpha の節で厳密に逆に引く)。
  臨界スピンは、相対論的な率が 0 になるスピンの表(alpha 0.003--0.6)をどこでも使う。
  領域の判定: 準平衡の占有数 eps_eq と行き過ぎの占有数 eps_th が、どちらも自由場で飽和する量 chi0 - chi_c より小さければ準平衡。
  準平衡では、1 個の量子から eps_eq まで育つ時間 t_grow を年齢から引き、行き過ぎで失う eps_th を初期スピンから引く。
超放射の率は Kerr 時空の連分数法(Dolan 0705.2880)の表(211: alpha >= 0.03、322: 0.1 <= alpha <= 0.52)、
範囲の外と 433 は解析式(Baumann+ 1908.10370 の式 2.32--2.34)。自己相互作用の率は主要項(相対論的な補正なし)。

別のプロセスで使うとき(multiprocessing など): ChiSR は pickle できる。中身は表の場所だけで、受け取った側が同じ表を
一度だけ読む(メモリは 1 プロセスあたり約 0.4 GB)。
"""
import os
from math import factorial

import numpy as np
from scipy.interpolate import RegularGridInterpolator


# ---------- 値の表(研究ノート E18): 基準の運動論から作った表の引き方(プロジェクトの calc/reftab.py と同じ) ----------
LT_LO, LT_HI = -17.0, 1.0
LT_STEP = (LT_HI - LT_LO) / 65533
KI_BLOCK = 4
KI_STEP = 0.01
LOW_MARGIN = 7.0      # 値の表で引けない固定点の点を chi0 のままにする条件: 雲が育つ時刻の何 e-fold 前か(E19)
# alpha の方向の折れ(E19): 段差の時刻(log10 tau)が alpha の節のあいだで折れている所を、表そのものから測る
KS_STEP = 0.002                    # [dex]、uint8(255 = 0.51 dex 以上)
KS_MIN = 0.02 / np.log(10.0)       # 約束が許す時刻のずれ 2% [dex]。これより小さい折れは見ない
KS_CAP = 0.3                       # 時刻をずらす幅の上限 [dex]
THETA_SW = 0.015                   # 折れの分だけ時刻をずらしたときの値の変化がこれを越えたら印 5
THETA_P = 0.02                     # 補間が止まった所(平らな所)で、角の節の値の幅がこれを越えたら印 5


def _dequantize_lt(q):
    v = LT_LO + (q.astype(float) - 1) * LT_STEP
    v[q == 0] = -np.inf
    v[q == 65535] = np.nan
    return v


class _LTCodes:
    """log10 tau の表を uint16 の符号のまま持ち、引いた要素だけを実数に直す(実数で持つより 4 倍小さい。値は同じ)。E19"""
    __slots__ = ("q",)

    def __init__(self, q):
        self.q = q

    @property
    def shape(self):
        return self.q.shape

    def __getitem__(self, idx):
        return _dequantize_lt(self.q[idx])


def _compute_ks(LT, x_alpha):
    """alpha の方向の、段差の時刻の折れ KS[セル i, js, kf, iL, block](KS_STEP [dex] 刻みの uint8)。
    セル (i, i+1) の準位ごとに、d1 = 下の 2 節 (i-1, i) を結ぶ直線を節 i+1 まで延ばした値と節 i+1 の log10 tau の差、
    d2 = 上の 2 節 (i+2, i+1) を結ぶ直線を節 i まで延ばした値と節 i の差を測る。折れがこのセルの中にあれば両方が大きく、
    隣のセルにあれば片方だけが大きい(延ばした直線が折れをまたぐ)。そこで小さいほうが KS_MIN を越える準位だけを折れと
    みなし、大きいほうを持つ(セルの中の本当の時刻がどちらかの直線の延長の上にあるなら、節 i と i+1 を結ぶ補間の時刻の
    誤差はこの値以下)。片方しか測れないときはその値、どちらも測れない(外の節が届いていない)ときは上限(255)。
    両方の節が届いている準位だけ測り、準位 KI_BLOCK 個ずつの最大を持つ。"""
    nA, K = LT.shape[0], LT.shape[-1]
    nb = -(-K // KI_BLOCK)
    out = np.zeros((nA - 1,) + tuple(LT.shape[1:-1]) + (nb,), np.uint8)
    node = {}                                            # alpha の節ごとの log10 tau(使い終わった節は捨てる)

    def lt(i):
        if i not in node:
            node[i] = np.asarray(LT[i], float)
        return node[i]

    for i in range(nA - 1):
        node.pop(i - 2, None)
        lo, hi = lt(i), lt(i + 1)
        both = np.isfinite(lo) & np.isfinite(hi)
        d1 = np.full(lo.shape, np.nan)
        d2 = np.full(lo.shape, np.nan)
        with np.errstate(invalid="ignore"):
            if i >= 1:
                l0 = lt(i - 1)
                d1 = np.abs(hi - (lo + (lo - l0) * (x_alpha[i + 1] - x_alpha[i]) / (x_alpha[i] - x_alpha[i - 1])))
            if i + 2 < nA:
                l3 = lt(i + 2)
                d2 = np.abs(lo - (hi + (hi - l3) * (x_alpha[i] - x_alpha[i + 1]) / (x_alpha[i + 1] - x_alpha[i + 2])))
            d1 = np.where(np.isfinite(d1), d1, np.nan)        # 外の節が始めの準位(-inf)のときも「測れない」
            d2 = np.where(np.isfinite(d2), d2, np.nan)
            dmin, dmax = np.fmin(d1, d2), np.fmax(d1, d2)      # 片方が NaN ならもう片方
            code = np.where(np.isfinite(dmax), np.where(dmin > KS_MIN, np.minimum(np.ceil(dmax / KS_STEP), 255), 0), 255)
        code = np.where(both, code, 0)
        pad = np.zeros(code.shape[:-1] + (nb * KI_BLOCK - K,))
        out[i] = np.concatenate([code, pad], axis=-1).reshape(code.shape[:-1] + (nb, KI_BLOCK)).max(-1).astype(np.uint8)
    return out


def _b_of_w(w, s):
    w = np.asarray(w, float)
    s = np.asarray(s, float)
    hi = s > 1
    with np.errstate(divide="ignore", invalid="ignore"):
        bh = np.where(w >= 1, (s - w) / (s - 1), 1 + (1 - w))
        bl = np.where(w >= 0, (s - w) / s, 1 - w)
    return np.where(hi, bh, bl)


def _w_of_b(b, s):
    b = np.asarray(b, float)
    s = np.asarray(s, float)
    hi = s > 1
    wh = np.where(b <= 1, s - b * (s - 1), 1 - (b - 1))
    wl = np.where(b <= 1, s - b * s, -(b - 1))
    return np.where(hi, wh, wl)


def _bracket(x, nodes):
    x = np.clip(x, nodes[0], nodes[-1])
    i = np.clip(np.searchsorted(nodes, x, side="right") - 1, 0, len(nodes) - 2)
    t = (x - nodes[i]) / (nodes[i + 1] - nodes[i])
    return i, t


def _lt_corner(LT, a_, s_, f_, il, tl, k):
    v0, v1 = LT[a_, s_, f_, il, k], LT[a_, s_, f_, il + 1, k]
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        return np.where(tl <= 0, v0, np.where(tl >= 1, v1, np.log10((1 - tl) * 10.0 ** v0 + tl * 10.0 ** v1)))


def _core2(ctx, T, alpha, s, f_gev, M, mu_eV, c1, c2, LT, RG, a_nodes, sig_nodes, lf_nodes, L_nodes, piece_hi):
    """第 2 の形の引き方。返り値: (chi, blocked, Rz)。"""
    T, alpha, s, f_gev, M, mu_eV, c1, c2 = (np.asarray(x, float) for x in (T, alpha, s, f_gev, M, mu_eV, c1, c2))
    n = len(T)
    lq = np.log10(np.maximum(mu_eV * T, 1e-300))
    lf_c = np.clip(np.log10(f_gev), lf_nodes[0], lf_nodes[-1])
    sig = np.log(np.maximum(s - 1, 1e-300)) if piece_hi else np.log(np.maximum(s, 1e-300))
    ia, ta = _bracket(np.log(alpha), np.log(a_nodes))
    js, ts = _bracket(sig, sig_nodes)
    kf, tf = _bracket(lf_c, lf_nodes)
    il, tl = _bracket(ctx.L_of_M(M), L_nodes)
    corners = []
    for da in (0, 1):
        for ds in (0, 1):
            for df in (0, 1):
                wc = (ta if da else 1 - ta) * (ts if ds else 1 - ts) * (tf if df else 1 - tf)
                corners.append((ia + da, js + ds, kf + df, wc))

    def lt_corner(a_, s_, f_, k):
        """節の log10 tau(L について tau を線形に混ぜる)。"""
        v0, v1 = LT[a_, s_, f_, il, k], LT[a_, s_, f_, il + 1, k]
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
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
    benv = ctx.P_LEV2[kk] + fr * (ctx.P_LEV2[nxt] - ctx.P_LEV2[kk])
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
            # 重みが 0 の節では tc が NaN になりうる(準位に届いていない節)。NaN を整数に直した結果は CPU で違い
            # (Apple Silicon は 0、x86 は INT64_MIN で添字の範囲の外になり IndexError)、0 にそろえる(E19。値は変わらない)
            j = np.nan_to_num(np.clip((tc - ctx.G_LT[0]) / (ctx.G_LT[1] - ctx.G_LT[0]), 0, len(ctx.G_LT) - 1.000001), nan=0.0)
            j0 = np.floor(j).astype(int)
            u = j - j0
            r0 = RG[a_, s_, f_, il_, j0]
            r1 = RG[a_, s_, f_, il_, np.minimum(j0 + 1, len(ctx.G_LT) - 1)]
            rv = np.nan_to_num((1 - u) * r0 + u * r1, nan=0.0)
            rz += np.where(pos, wc * wl * rv, 0.0)
    b = benv - rz
    w = _w_of_b(b, s)
    # blocked: 次の準位に届かない角の節があって、補間の進み具合がここで止まる(E19。平らな所の値の判定に使う)
    return np.asarray(ctx.chi_of_w(w, c1, c2), float), np.isnan(l1) & (nxt > kk), rz


def _corner_b(ctx, LT, RG, a_, s_, f_, il, tl, lq):
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
    benv = ctx.P_LEV2[lo] + fr * (ctx.P_LEV2[nxt] - ctx.P_LEV2[lo])
    ilr = np.where(tl < 0.5, il, il + 1)
    j = np.clip((lq - ctx.G_LT[0]) / (ctx.G_LT[1] - ctx.G_LT[0]), 0, len(ctx.G_LT) - 1.000001)
    j0 = np.floor(j).astype(int)
    u = j - j0
    rv = (1 - u) * RG[a_, s_, f_, ilr, j0] + u * RG[a_, s_, f_, ilr, np.minimum(j0 + 1, len(ctx.G_LT) - 1)]
    return benv - np.nan_to_num(rv, nan=0.0), lo


def lookup_rt3(ctx, rt, T, alpha, chi0, f_gev, M, mu_eV, c1, c2, theta_v=0.06, theta_k=0.05, margin=10.0):
    """表からスピンを引く。返り値: (chi, valid, sharp, kink)。
    valid: 表の範囲の中で、重みのある節がどれも tau_q の margin 倍まで解けている。
    sharp: 速い段差の時刻が節のあいだで急に変わる所(印 5)。角の節それぞれの、問いの時刻での値の幅が theta_v を越え、
           しかもその時刻に角の節がいる準位で、段差の時刻が両隣を結ぶ直線から theta_k [dex] 以上ずれている。
    kink: 節のあいだで補間できない所(印 5、E19)。(1) 段差の時刻が alpha の節のあいだで折れている所。問いのセルの折れ KS(_compute_ks、角の節と、その時刻に
          角の節がいる準位での最大)が KS_MIN を越え、しかも時刻を折れの分だけ前後にずらすと値が THETA_SW より変わる。
          (2) 補間が止まった所(次の準位に届かない角の節がある)で、角の節の値の幅が THETA_P を越える。"""
    T, alpha, chi0, f_gev, M, mu_eV, c1, c2 = np.broadcast_arrays(*[np.asarray(x, float) for x in
                                                                    (T, alpha, chi0, f_gev, M, mu_eV, c1, c2)])
    s = np.asarray(ctx.w_of_chi(chi0, c1, c2), float)
    an, lfn, Mn = rt["A_NODES"], rt["LF_NODES"], rt["M_NODES"]
    lf = np.log10(f_gev)
    tau = mu_eV * T
    inside = ((alpha >= an[0]) & (alpha <= an[-1]) & (s > 0) & (chi0 <= ctx.CHI_TOP) & (M >= Mn[0]) & (M <= Mn[-1])
              & (lf >= lfn[0] - 1e-12) & (lf <= lfn[-1] + 1e-12) & (tau <= rt["TAU_MAX"] / margin))
    chi = chi0.copy()
    valid = np.zeros(len(T), bool)
    sharp = np.zeros(len(T), bool)
    kink = np.zeros(len(T), bool)
    Ln = ctx.L_of_M(Mn)
    for piece, snodes in ((0, rt["SA"]), (1, rt["SB"])):
        q = np.flatnonzero(inside & ((s > 1) if piece else (s <= 1)))
        if len(q) == 0:
            continue
        LT, RG, KI, KS = rt["LT"][piece], rt["RG"][piece], rt["KI"][piece], rt["KS"][piece]
        ST, TV = rt["status"][piece], rt["tau_valid"][piece]
        sq = s[q]
        sig = np.log(np.maximum(sq - 1, 1e-300)) if piece else np.log(np.maximum(sq, 1e-300))
        c, blocked, _ = _core2(ctx, T[q], alpha[q], sq, f_gev[q], M[q], mu_eV[q], c1[q], c2[q], LT, RG, an, np.log(snodes), lfn, Ln,
                               bool(piece))
        ia, ta = _bracket(np.log(alpha[q]), np.log(an))
        js, ts = _bracket(sig, np.log(snodes))
        kf, tf = _bracket(lf[q], lfn)
        il, tl = _bracket(ctx.L_of_M(M[q]), Ln)
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
                    bc, kc = _corner_b(ctx, LT, RG, a_, s_, f_, il, tl, lq)
                    xc = np.asarray(ctx.chi_of_w(_w_of_b(bc, sq), c1[q], c2[q]), float)
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
        # alpha の方向の折れ(E19): 問いのセルの KS の最大(角の s・f の節、M の両方の節、角の節がいる準位のブロック)
        dq = np.zeros(len(q))
        for ds in (0, 1):
            for df in (0, 1):
                pos = ((ts if ds else 1 - ts) * (tf if df else 1 - tf)) > 0
                for il_ in (il, il + 1):
                    for d in range(nblk):
                        bb = np.minimum(b0 + d, b1)
                        dq = np.where(pos, np.maximum(dq, KS[ia, js + ds, kf + df, il_, bb] * KS_STEP), dq)
        kq = np.zeros(len(q), bool)
        mk = np.flatnonzero(dq > KS_MIN)
        if len(mk):
            sh = 10.0 ** np.minimum(dq[mk], KS_CAP)
            dv = np.zeros(len(mk))
            for Ts in (T[q][mk] * sh, T[q][mk] / sh):
                cs_, _, _ = _core2(ctx, Ts, alpha[q][mk], sq[mk], f_gev[q][mk], M[q][mk], mu_eV[q][mk], c1[q][mk], c2[q][mk],
                                   LT, RG, an, np.log(snodes), lfn, Ln, bool(piece))
                dv = np.fmax(dv, np.abs(cs_ - c[mk]))
            kq[mk] = ~(dv <= THETA_SW)                    # 値を引けない(NaN)ときは印をつける側
        chi[q] = c
        valid[q] = ok & np.isfinite(c)
        sharp[q] = ((vmax - vmin) > theta_v) & (kmaxv > theta_k)
        # 平らな所の値(E19): 補間は、準位ごとの時刻を節のあいだで混ぜるので、重みのある角の節のうち一つでも次の準位に
        # 届かないと、そこで止まる。止まった所での角の節の値の幅が THETA_P を越えていれば、値は進みのいちばん遅い節の値に
        # 引かれて、その幅までずれうる(BH の質量が減ったあとの臨界の値が、節によって違うときなど)。印 5 にする。
        kink[q] = kq | (blocked & ((vmax - vmin) > THETA_P))
    return chi, valid, sharp, kink


_MODELS = {}


def _load_cached(path):
    """pickle から戻すときに使う: 同じプロセスでは同じ表を一度だけ読む。"""
    if path not in _MODELS:
        _MODELS[path] = ChiSR(path)
    return _MODELS[path]


class ChiSR:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        with np.load(self.path) as z:
            self._read(z)

    def __reduce__(self):
        # pickle には表の場所だけを入れる(表そのものはメモリの上で約 0.4 GB)。multiprocessing などで別のプロセスに
        # 渡すと、受け取った側が同じファイルを一度だけ読む。別の計算機に渡すときは、同じ場所に表が要る。
        return (_load_cached, (self.path,))

    def _read(self, z):
        for k in ["HBAR_EVS", "TSUN", "N0", "YR", "MPL_GEV", "KAPPA_BH", "KAPPA_INF", "K_C", "K_D", "CHI_TOP",
                  "ALPHA_MAX", "F_FP_MAX"]:
            setattr(self, k, float(z[k]))
        self.ALPHAS = z["ALPHAS"]
        self.W = z["W"]
        self.lnI = z["lnI"]
        self.J = z["J"]
        self.thr_alpha = z["thr_alpha"]
        self.thr = {1: z["thr1"], 2: z["thr2"], 3: z["thr3"]}
        self._rel = {}
        self._range = {}
        for key, l in [("211", 1), ("322", 2)]:
            self._rel[key] = RegularGridInterpolator((np.log(z[f"rel{key}_alpha"]), z[f"rel{key}_chi"]), z[f"rel{key}_F"],
                                                     method="cubic", bounds_error=False, fill_value=None)
            self._range[key] = tuple(z[f"rel{key}_range"])
        # 値の表(E18)。ファイルでは log10 tau を準位の差の uint16 で持つ(届いた準位の数 rt_nreach より先は届かない)。
        # メモリの上では、log10 tau は uint16 の符号、戻りは float16 のまま持つ(全体で約 0.4 GB。実数に直して持つと 1.4 GB。E19)
        self.rt = None
        if "rt_LTd" in z.files:
            nr = z["rt_nreach"].astype(np.int64)
            codes = np.cumsum(z["rt_LTd"], axis=-1, dtype=np.uint16)       # 符号は 65534 以下なので、あふれない
            codes[np.arange(codes.shape[-1]) >= nr[..., None]] = 65535
            codes[..., 0] = np.where(nr > 0, 0, 65535)
            self.rt = {"LT": tuple(_LTCodes(c) for c in codes), "RG": z["rt_RG"], "KI": z["rt_KI"],
                       "status": z["rt_status"], "tau_valid": z["rt_tau_valid"]}
            for k in ("A_NODES", "SA", "SB", "LF_NODES", "M_NODES"):
                self.rt[k] = z[f"rt_{k}"]
            self.rt["TAU_MAX"] = float(z["rt_TAU_MAX"])
            self.P_LEV2 = z["rt_P_LEV"]
            self.G_LT = z["rt_G_LT"]
            # 段差の時刻の、alpha の方向の折れ(E19)。表から一度だけ測る
            self.rt["KS"] = [_compute_ks(lt, np.log(self.rt["A_NODES"])) for lt in self.rt["LT"]]

    # 表の引き方に渡す文脈(calc/reftab.py の関数と同じ名前)
    def chi_of_w(self, w, c1, c2):
        return self._chi_of_w(w, c1, c2)

    def w_of_chi(self, chi, c1, c2):
        return self._w_of_chi(chi, c1, c2)

    def L_of_M(self, M):
        return np.log(self.N0 * np.asarray(M, float) ** 2)

    # ---------- 基本の量 ----------
    def alpha(self, M, mu):
        return self.TSUN * M * mu / self.HBAR_EVS

    @staticmethod
    def rplus(chi):
        return 1.0 + np.sqrt(1.0 - chi ** 2)

    @staticmethod
    def _omega_over_mu(a, n, l, m, chi):
        f_nl = -6.0 / (2 * l + 1) + 2.0 / n
        h_l = 16.0 / (2 * l * (2 * l + 1) * (2 * l + 2))
        return 1 - a ** 2 / (2 * n ** 2) - a ** 4 / (8 * n ** 4) + f_nl * a ** 4 / n ** 3 + h_l * chi * m * a ** 5 / n ** 3

    def _gamma_analytic(self, a, chi, n, l, m):
        """1908.10370 の式 2.32 の占有数の率 / mu(= 2 x 振幅の率)。"""
        C = (2.0 ** (4 * l + 1) * factorial(n + l) / (n ** (2 * l + 4) * factorial(n - l - 1))
             * (factorial(l) / (factorial(2 * l) * factorial(2 * l + 1))) ** 2)
        rp = self.rplus(chi)
        w = a * self._omega_over_mu(a, n, l, m, chi)
        g = 1.0
        for k in range(1, l + 1):
            g = g * (k ** 2 * (1 - chi ** 2) + (chi * m - 2 * rp * w) ** 2)
        return 2.0 * 2 * rp * C * g * (m * chi / (2 * rp) - w) * a ** (4 * l + 5) / a

    def g(self, a, chi, level):
        """超放射の率 / mu(占有数の率)。level = "211"、"322"、"433"。"""
        n, l, m = {"211": (2, 1, 1), "322": (3, 2, 2), "433": (4, 3, 3)}[level]
        a, chi = np.broadcast_arrays(np.asarray(a, float), np.asarray(chi, float))
        out = np.array(self._gamma_analytic(a, chi, n, l, m), float)
        lo, hi = self._range.get(level, (np.inf, -np.inf))      # 433 は表を使わず解析式
        inr = (a >= lo) & (a <= hi)
        if inr.any():
            F = self._rel[level](np.stack([np.log(a[inr]), chi[inr]], -1))
            out[inr] = 2.0 * F * a[inr] ** (4 * l + 5) / a[inr]
        return out

    def chic(self, a, m):
        """臨界スピン(m = 1, 2, 3)。率の 0 の表(alpha 0.003--0.6)を補間。下は水素原子の主要項、上は 1(超放射しない)。"""
        a = np.asarray(a, float)
        n = m + 1
        ap = a * (1 - a ** 2 / (2 * n ** 2))
        out = np.where(ap / m < 0.5, 4 * m * ap / (m ** 2 + 4 * ap ** 2), 1.0)
        ta = self.thr_alpha
        inside = np.interp(np.log(np.clip(a, ta[0], ta[-1])), np.log(ta), self.thr[m])
        return np.where(a < ta[0], out, np.where(a > ta[-1], 1.0, inside))

    # ---------- 表の引き方 ----------
    def _w_of_chi(self, chi, c1, c2):
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(chi >= c1, 1 + (chi - c1) / (self.CHI_TOP - c1), (chi - c2) / (c1 - c2))

    def _chi_of_w(self, w, c1, c2):
        return np.where(w >= 1, c1 + (w - 1) * (self.CHI_TOP - c1), c2 + w * (c1 - c2))

    def _alpha_nodes(self, a):
        la = np.log(self.ALPHAS)
        x = np.clip(np.log(a), la[0], la[-1])
        i = np.clip(np.searchsorted(la, x, side="right") - 1, 0, len(la) - 2)
        return i, (x - la[i]) / (la[i + 1] - la[i])

    def _fwd(self, arr, rows, w):
        W = self.W
        j = np.clip(np.searchsorted(W, w, side="right") - 1, 0, len(W) - 2)
        s = (w - W[j]) / (W[j + 1] - W[j])
        return arr[rows, j] + s * (arr[rows, j + 1] - arr[rows, j])

    def _inv(self, rows, y):
        """lnI[row](w) = y となる w(区分線形の補間の厳密な逆。lnI は w について減少)。"""
        lnI, W = self.lnI, self.W
        n = lnI.shape[1]
        y = np.clip(y, lnI[rows, -1], lnI[rows, 0])
        j = np.empty(len(rows), int)
        order = np.argsort(rows, kind="stable")
        r_sorted = rows[order]
        starts = np.flatnonzero(np.r_[True, r_sorted[1:] != r_sorted[:-1]])
        ends = np.r_[starts[1:], len(rows)]
        for st, en in zip(starts, ends):
            idx = order[st:en]
            j[idx] = np.searchsorted(-lnI[r_sorted[st]], -y[idx], side="left") - 1
        j = np.clip(j, 0, n - 2)
        v0, v1 = lnI[rows, j], lnI[rows, j + 1]
        with np.errstate(divide="ignore", invalid="ignore"):
            s = np.where(v0 > v1, (v0 - y) / (v0 - v1), 0.0)
        return W[j] + np.clip(s, 0.0, 1.0) * (W[j + 1] - W[j])

    def _chi_after(self, a, chi0, tau):
        c1, c2 = self.chic(a, 1), self.chic(a, 2)
        w0 = np.clip(self._w_of_chi(chi0, c1, c2), self.W[0], self.W[-1])
        i, t = self._alpha_nodes(a)
        lt = np.log(np.maximum(tau, 1e-300))
        wa = self._inv(i, np.logaddexp(self._fwd(self.lnI, i, w0), lt))
        wb = self._inv(i + 1, np.logaddexp(self._fwd(self.lnI, i + 1, w0), lt))
        out = self._chi_of_w((1 - t) * wa + t * wb, c1, c2)
        keep = (chi0 <= c2) | (a < self.ALPHAS[0]) | (a > self.ALPHAS[-1])
        return np.where(keep, chi0, np.minimum(out, chi0))

    def _flag_fp(self, a, chi_start, chi_T, s2, T_s, mu_s, lnGM2):
        """固定点の領域の印(4 か 1): 三準位の外の準位が年齢のうちに 1 個の量子から育ちうるなら 4。"""
        c1, c2 = self.chic(a, 1), self.chic(a, 2)
        w0 = np.clip(self._w_of_chi(chi_start, c1, c2), self.W[0], self.W[-1])
        wT = np.clip(self._w_of_chi(chi_T, c1, c2), self.W[0], self.W[-1])
        i, t = self._alpha_nodes(a)
        J0 = (1 - t) * self._fwd(self.J, i, w0) + t * self._fwd(self.J, i + 1, w0)
        JT = (1 - t) * self._fwd(self.J, i, wT) + t * self._fwd(self.J, i + 1, wT)
        dJ = np.where(JT > J0 * (1 + 1e-9), JT - J0, 0.0)
        n_desc = dJ / s2
        bottom = wT <= self.W[0] * 1.000001
        g433 = np.maximum(self.g(a, np.minimum(chi_T, 0.995), "433"), 0.0)
        n_bot = np.where(bottom, mu_s * g433 * T_s, 0.0)
        return np.where((n_desc > lnGM2) | (n_bot > lnGM2), 4, 1)

    # ---------- 自由場 ----------
    def _sat_state(self, M, chi, mu, m):
        """飽和したあとの (M_f, chi_f)。BH の質量の減少を入れ、臨界は self.chic。Illinois 法で解く。"""
        n = m + 1
        a_sun = self.alpha(1.0, mu)
        J0 = chi * M ** 2

        def F(Mf):
            a = a_sun * Mf
            return Mf - M + a_sun * (1 - a ** 2 / (2 * n ** 2)) * (J0 - self.chic(a, m) * Mf ** 2) / m

        lo, hi = 0.5 * np.asarray(M, float), np.array(M, float)
        Flo, Fhi = F(lo), F(hi)
        side = np.zeros(lo.shape, int)
        x = hi.copy()
        for _ in range(60):
            with np.errstate(divide="ignore", invalid="ignore"):
                x = np.where(Fhi != Flo, hi - Fhi * (hi - lo) / (Fhi - Flo), 0.5 * (lo + hi))
            x = np.clip(x, lo, hi)
            Fx = F(x)
            neg = Fx < 0
            Fhi = np.where(neg & (side == -1), 0.5 * Fhi, Fhi)
            Flo = np.where(~neg & (side == 1), 0.5 * Flo, Flo)
            lo, Flo = np.where(neg, x, lo), np.where(neg, Fx, Flo)
            hi, Fhi = np.where(neg, hi, x), np.where(neg, Fhi, Fx)
            side = np.where(neg, -1, 1)
            if np.all((hi - lo) <= 1e-13 * hi) or np.all(np.abs(Fx) <= 1e-13 * hi):
                break
        return x, self.chic(a_sun * x, m)

    def _free(self, M, chi0, a, mu, mu_s, T, f_gev, g1):
        M1, c1 = self._sat_state(M, chi0, mu, 1)
        c1 = np.minimum(c1, chi0)
        lnN1 = np.log(np.maximum(chi0 * M ** 2 - c1 * M1 ** 2, 1e-300)) + np.log(self.N0)
        tsat1 = lnN1 / g1
        M2, c2 = self._sat_state(M1, c1, mu, 2)
        lnN2 = np.log(np.maximum((c1 * M1 ** 2 - c2 * M2 ** 2) / 2, 1e-300)) + np.log(self.N0)
        g2a = np.maximum(self.g(a, np.minimum(chi0, 0.995), "322"), 0) * mu_s
        a1 = self.alpha(M1, mu)
        g2b = np.maximum(self.g(a1, c1, "322"), 0) * mu_s
        with np.errstate(divide="ignore", invalid="ignore"):
            t2_first = np.where(g2a > 0, lnN2 / g2a, np.inf)
            t2 = np.where(g2b > 0, tsat1 + np.maximum(lnN2 - g2a * tsat1, 0) / g2b, np.inf)
        res = np.where(T >= tsat1, c1, chi0)
        res = np.where((t2_first < tsat1) & (T >= t2_first), c2, res)
        slow = f_gev > self.F_FP_MAX
        res = np.where(slow & (t2_first >= tsat1) & (c1 > c2) & (T >= t2), c2, res)
        fp = (~slow) & (t2_first >= tsat1) & (T > tsat1)
        if fp.any():
            tau1 = mu_s[fp] * (f_gev[fp] / self.MPL_GEV) ** 2 * (T[fp] - tsat1[fp])
            start = c1[fp]
            drop = start - self._chi_after(a1[fp], start + 1e-15, tau1)
            res[fp] = np.maximum(c1[fp] - drop, c2[fp])
        return res

    # ---------- 本体 ----------
    @staticmethod
    def _check(M, chi0, T_yr, mu_eV, f_gev):
        for name, x in (("M", M), ("chi0", chi0), ("T_yr", T_yr), ("mu_eV", mu_eV), ("f_GeV", f_gev)):
            if not np.all(np.isfinite(x)):
                raise ValueError(f"{name} に有限でない値がある")
        if np.any((chi0 < 0) | (chi0 >= 1)):
            raise ValueError("chi0 は 0 <= chi0 < 1")
        if np.any(M <= 0) or np.any(mu_eV <= 0) or np.any(f_gev <= 0) or np.any(T_yr < 0):
            raise ValueError("M, mu_eV, f_GeV は正、T_yr は 0 以上")

    def __call__(self, M, chi0, T_yr, mu_eV, f_gev):
        M, chi0, T_yr, mu_eV, f_gev = np.broadcast_arrays(*[np.asarray(x, float) for x in (M, chi0, T_yr, mu_eV, f_gev)])
        shape = M.shape
        M, chi0, T_yr, mu_eV, f_gev = (x.ravel() for x in (M, chi0, T_yr, mu_eV, f_gev))
        self._check(M, chi0, T_yr, mu_eV, f_gev)
        a = self.alpha(M, mu_eV)
        mu_s = mu_eV / self.HBAR_EVS
        T = T_yr * self.YR
        lnGM2 = np.log(self.N0 * M ** 2)
        c1 = self.chic(a, 1)
        c2 = self.chic(a, 2)
        big = a > self.ALPHA_MAX
        out = chi0.copy()
        reg = np.zeros(len(M), int)
        reg[big] = 3
        # 値の表(E18): 超放射する点の値を、基準の運動論から作った表で引く
        tv = chi0.copy()
        use_t = np.zeros(len(M), bool)
        out_rt = np.zeros(len(M), bool)
        sharp_t = np.zeros(len(M), bool)
        kink_t = np.zeros(len(M), bool)
        if self.rt is not None:
            cand_t = (chi0 > c2) & ~big
            if cand_t.any():
                it = np.flatnonzero(cand_t)
                lfmax = self.rt["LF_NODES"][-1]
                ff = np.minimum(f_gev[it], 10 ** lfmax)
                v, ok, sh, kk = lookup_rt3(self, self.rt, T_yr[it], a[it], chi0[it], ff, M[it], mu_eV[it], c1[it], c2[it])
                tv[it] = v
                use_t[it] = ok
                sharp_t[it] = sh & ok
                kink_t[it] = kk & ok
                below = (np.log10(f_gev[it]) < self.rt["LF_NODES"][0]) | (a[it] < self.rt["A_NODES"][0])
                out_rt[it] = ~ok & ~below
        sr211 = (chi0 > c1) & ~big
        g1 = np.zeros(len(M))
        i1 = np.flatnonzero(sr211)
        g1[i1] = self.g(a[i1], np.minimum(chi0[i1], 0.995), "211") * mu_s[i1]
        ginf = self.KAPPA_INF * a ** 8 * (self.MPL_GEV / f_gev) ** 4 * mu_s
        gbh = self.KAPPA_BH * a ** 11 * (self.MPL_GEV / f_gev) ** 4 * self.rplus(np.minimum(chi0, 0.999)) * mu_s
        emax = chi0 - c1
        with np.errstate(divide="ignore", invalid="ignore"):
            eeq = 2 / np.sqrt(3) * np.sqrt(g1 * ginf) / gbh
            eth = np.sqrt(2 * g1 * lnGM2 / gbh)
        si_reg = sr211 & (g1 > 0) & (eeq < emax) & (eth < emax)
        inval = si_reg & (a > self.ALPHAS[-1])
        reg[inval] = 3
        reg[sr211 & (g1 <= 0) & (chi0 > self.CHI_TOP)] = 3
        q = si_reg & ~inval
        if q.any():
            tgrow = (np.log(np.maximum(eeq[q], 1e-300)) + lnGM2[q]) / g1[q]
            Teff = np.maximum(T[q] - tgrow, 0.0)
            tau = mu_s[q] * (f_gev[q] / self.MPL_GEV) ** 2 * Teff
            c0eff = np.maximum(chi0[q] - eth[q], c1[q])
            out[q] = np.where(Teff > 0, self._chi_after(a[q], c0eff, tau), chi0[q])
            out[q] = np.where(use_t[q], tv[q], out[q])
            reg[q] = self._flag_fp(a[q], np.maximum(c0eff, out[q]), out[q], (f_gev[q] / self.MPL_GEV) ** 2, Teff, mu_s[q],
                                   lnGM2[q])
            # alpha < 0.01 で値の表を引けない点(alpha < 0.008 か f < 1e8 GeV): 雲が固定点の占有数に届く時刻の
            # LOW_MARGIN e-fold 前より若ければ chi0 のまま(印 1)、そうでなければ模型の外(印 3)。E19
            iq = np.flatnonzero(q)
            late = (a[q] < self.ALPHAS[0]) & ~use_t[q] & (g1[q] * (tgrow - T[q]) < LOW_MARGIN)
            out[iq[late]] = chi0[iq[late]]
            reg[iq[late]] = 3
        fr = sr211 & ~si_reg & (g1 > 0)
        if fr.any():
            out[fr] = self._free(M[fr], chi0[fr], a[fr], mu_eV[fr], mu_s[fr], T[fr], f_gev[fr], g1[fr])
            out[fr] = np.where(use_t[fr], tv[fr], out[fr])
            reg[fr] = 2
            if self.rt is not None:
                # 値の表の下(alpha < 0.008)の自由場の点: 雲が飽和する時刻の LOW_MARGIN e-fold 前より若い点だけ chi0(印 2)、
                # ほかは模型の外(印 3)。半解析の飽和の時刻は約束の 2% より粗いことがある(E19)
                late = fr & (a < self.rt["A_NODES"][0]) & (g1 * T > lnGM2 + np.log(np.maximum(emax, 1e-300)) - LOW_MARGIN)
                out[late] = chi0[late]
                reg[late] = 3
        only322 = (~sr211) & (chi0 > c2) & ~big
        if only322.any():
            i3 = np.flatnonzero(only322)
            M3, a3, x0, T3, f3, m3, u3 = M[i3], a[i3], chi0[i3], T[i3], f_gev[i3], mu_s[i3], mu_eV[i3]
            c0 = np.minimum(x0, 0.995)
            g2 = self.g(a3, c0, "322") * m3
            s4 = (self.MPL_GEV / f3) ** 4 * m3
            gc = self.K_C * a3 ** 11 * self.rplus(c0) * s4
            gd = self.K_D * a3 ** 8 * s4
            emax2 = (x0 - c2[i3]) / 2
            with np.errstate(divide="ignore", invalid="ignore"):
                eeq2 = 2 / np.sqrt(3) * np.sqrt(np.maximum(g2, 0) * gd) / gc
                eth2 = np.sqrt(2 * np.maximum(g2, 0) * lnGM2[i3] / gc)
                tgrow2 = (np.log(np.maximum(eeq2, 1e-300)) + lnGM2[i3]) / g2
            fpq = (g2 > 0) & (eeq2 < emax2) & (eth2 < emax2)
            ok = fpq & (a3 <= self.ALPHAS[-1])
            res = x0.copy()
            rg = np.where((g2 > 0) | (x0 > self.CHI_TOP), 3, 0)
            if ok.any():
                Teff = np.maximum(T3[ok] - tgrow2[ok], 0.0)
                tau = m3[ok] * (f3[ok] / self.MPL_GEV) ** 2 * Teff
                c0eff = np.maximum(x0[ok] - 2 * eth2[ok], c2[i3][ok] + 1e-12)
                res[ok] = np.where(Teff > 0, self._chi_after(a3[ok], c0eff, tau), x0[ok])
                res[ok] = np.where(use_t[i3][ok], tv[i3][ok], res[ok])
                rg[ok] = self._flag_fp(a3[ok], np.maximum(c0eff, res[ok]), res[ok], (f3[ok] / self.MPL_GEV) ** 2, Teff,
                                       m3[ok], lnGM2[i3][ok])
                jo = np.flatnonzero(ok)
                late = (a3[ok] < self.ALPHAS[0]) & ~use_t[i3][ok] & (g2[ok] * (tgrow2[ok] - T3[ok]) < LOW_MARGIN)
                res[jo[late]] = x0[jo[late]]
                rg[jo[late]] = 3
            free = (g2 > 0) & ~fpq
            if free.any():
                Mf2, cf2 = self._sat_state(M3[free], x0[free], u3[free], 2)
                cf2 = np.minimum(cf2, x0[free])
                lnN2 = np.log(np.maximum((x0[free] * M3[free] ** 2 - cf2 * Mf2 ** 2) / 2, 1e-300)) + np.log(self.N0)
                res[free] = np.where(T3[free] >= lnN2 / g2[free], cf2, x0[free])
                res[free] = np.where(use_t[i3][free], tv[i3][free], res[free])
                rg[free] = 2
                if self.rt is not None:
                    late = free & (a3 < self.rt["A_NODES"][0]) & (g2 * T3 > lnGM2[i3] + np.log(np.maximum(emax2, 1e-300)) - LOW_MARGIN)
                    res[late] = x0[late]
                    rg[late] = 3
            out[i3] = res
            reg[i3] = rg
        # 表の範囲の中なのに値を引けない点は印 3、初期スピンが 0.995 より上で超放射する点も印 3(E18)
        bad_t = out_rt & np.isin(reg, (1, 2, 4))
        out[bad_t] = chi0[bad_t]
        reg[bad_t] = 3
        top = (chi0 > self.CHI_TOP) & np.isin(reg, (1, 2, 4))
        out[top] = chi0[top]
        reg[top] = 3
        # m = 3(433)の超放射: 返すスピンから 433 が年齢のうちに飽和しうるなら、印 0 -> 3、印 2 -> 4(値は変えない)
        c3 = self.chic(a, 3)
        cand = ((reg == 0) | (reg == 2)) & (out > c3) & ~big
        if cand.any():
            i4 = np.flatnonzero(cand)
            g433 = self.g(a[i4], np.minimum(out[i4], 0.995), "433") * mu_s[i4]
            lnN3 = np.log(np.maximum((out[i4] - c3[i4]) / 3, 1e-300)) + lnGM2[i4]
            grows = (g433 > 0) & (g433 * T[i4] > lnN3)
            reg[i4[grows]] = np.where(reg[i4[grows]] == 0, 3, 4)
        # 自由場で、返すスピンが 433 の臨界より下の点(E19): 三準位の運動論で 544(m = 4)の超放射まで進んだ点。433(m = 3)の
        # 率のほうが大きく、同じスピンから先に飽和するはずなので、値は粗い(印 4)
        reg[(reg == 2) & (out < c3) & ~big] = 4
        # 速い段差の時刻が表の節のあいだで急に変わる所(E18): 値は表の補間のまま、印 5
        reg[np.isin(reg, (1, 2)) & sharp_t] = 5
        # 段差の時刻が alpha の節のあいだで折れていて、値がそれに敏感な所(E19): 値は表の補間のまま、印 5
        reg[np.isin(reg, (1, 2)) & kink_t] = 5
        return out.reshape(shape), reg.reshape(shape)
