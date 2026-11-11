"""chi_SR: 自己相互作用のある ALP の超放射で、年齢 T のあとの BH のスピン(単体版、numpy と scipy だけ)。

使い方:
    from chi_sr import ChiSR
    model = ChiSR("chi_sr_tables.npz")
    chi_T, regime = model(M_msun, chi0, T_yr, mu_eV, f_GeV)   # どれも配列でよい(ブロードキャスト)

regime(印):
    0 = 超放射しない(chi0 が 211 と 322 の臨界以下。m >= 3 は入れていない)
    1 = 自己相互作用の固定点(211・322・544 の三準位)。211 から育つ BH も 322 から育つ BH も、322 の臨界まで
        下がりうる。三準位の外の準位(411、433 など)は年齢のうちに育たない
    2 = 自由場(211 が 1 個の量子から飽和して臨界へ。そのあと f <= 1e17 GeV は固定点の降下、それより上は 322 を
        1 個の量子から)、または 322 から育つ BH で自己相互作用が飽和までに効かない場合(どちらも BH の質量の減少を入れる)
    3 = この模型の外。chi0 を返す: 固定点で alpha < 0.01 または > 0.30、alpha > 0.5、
        211・322 は超放射しないが 433(m = 3)が年齢のうちに飽和しうる
    4 = 値は計算して返すが、模型の外の準位が年齢のうちに育ちうるので近似は粗い: 固定点で 411、433、211 の臨界より上の
        544 が育ちうる(有限 alpha の補正込みで判定)、または自由場で止まったスピンから 433 が飽和しうる

模型の要点(詳しくは README.md):
  固定点は Baryakhtar+ 2011.11646 の式 34--36 に 544 を加えた三準位(Table V・VI の率)。211 の臨界より上は (211, 322)、
  下は 211 が空の (322, 544)。
  dchi/dt = mu (f/M_pl)^2 H(alpha, chi) なので、I(alpha, chi) = int_chi^0.995 dchi'/|H| を表にし、
  I(alpha, chi_T) = I(alpha, chi_0') + mu (f/M_pl)^2 (T - t_grow) を chi_T について解く(両隣の alpha の節で厳密に逆に引く)。
  臨界スピンは、相対論的な率が 0 になるスピンの表(alpha 0.003--0.6)をどこでも使う。
  領域の判定: 準平衡の占有数 eps_eq と行き過ぎの占有数 eps_th が、どちらも自由場で飽和する量 chi0 - chi_c より小さければ準平衡。
  準平衡では、1 個の量子から eps_eq まで育つ時間 t_grow を年齢から引き、行き過ぎで失う eps_th を初期スピンから引く。
超放射の率は Kerr 時空の連分数法(Dolan 0705.2880)の表(211: alpha >= 0.03、322: 0.1 <= alpha <= 0.52)、
範囲の外と 433 は解析式(Baumann+ 1908.10370 の式 2.32--2.34)。自己相互作用の率は主要項(相対論的な補正なし)。
"""
from math import factorial

import numpy as np
from scipy.interpolate import RegularGridInterpolator


class ChiSR:
    def __init__(self, path):
        z = np.load(path)
        self.z = z
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
        inval = si_reg & ((a > self.ALPHAS[-1]) | (a < self.ALPHAS[0]))
        reg[inval] = 3
        reg[sr211 & (g1 <= 0) & (chi0 > self.CHI_TOP)] = 3
        q = si_reg & ~inval
        if q.any():
            tgrow = (np.log(np.maximum(eeq[q], 1e-300)) + lnGM2[q]) / g1[q]
            Teff = np.maximum(T[q] - tgrow, 0.0)
            tau = mu_s[q] * (f_gev[q] / self.MPL_GEV) ** 2 * Teff
            c0eff = np.maximum(chi0[q] - eth[q], c1[q])
            out[q] = np.where(Teff > 0, self._chi_after(a[q], c0eff, tau), chi0[q])
            reg[q] = self._flag_fp(a[q], c0eff, out[q], (f_gev[q] / self.MPL_GEV) ** 2, Teff, mu_s[q], lnGM2[q])
        fr = sr211 & ~si_reg & (g1 > 0)
        if fr.any():
            out[fr] = self._free(M[fr], chi0[fr], a[fr], mu_eV[fr], mu_s[fr], T[fr], f_gev[fr], g1[fr])
            reg[fr] = 2
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
            ok = fpq & (a3 >= self.ALPHAS[0]) & (a3 <= self.ALPHAS[-1])
            res = x0.copy()
            rg = np.where((g2 > 0) | (x0 > self.CHI_TOP), 3, 0)
            if ok.any():
                Teff = np.maximum(T3[ok] - tgrow2[ok], 0.0)
                tau = m3[ok] * (f3[ok] / self.MPL_GEV) ** 2 * Teff
                c0eff = np.maximum(x0[ok] - 2 * eth2[ok], c2[i3][ok] + 1e-12)
                res[ok] = np.where(Teff > 0, self._chi_after(a3[ok], c0eff, tau), x0[ok])
                rg[ok] = self._flag_fp(a3[ok], c0eff, res[ok], (f3[ok] / self.MPL_GEV) ** 2, Teff, m3[ok], lnGM2[i3][ok])
            free = (g2 > 0) & ~fpq
            if free.any():
                Mf2, cf2 = self._sat_state(M3[free], x0[free], u3[free], 2)
                cf2 = np.minimum(cf2, x0[free])
                lnN2 = np.log(np.maximum((x0[free] * M3[free] ** 2 - cf2 * Mf2 ** 2) / 2, 1e-300)) + np.log(self.N0)
                res[free] = np.where(T3[free] >= lnN2 / g2[free], cf2, x0[free])
                rg[free] = 2
            out[i3] = res
            reg[i3] = rg
        # m = 3(433)の超放射: 返すスピンから 433 が年齢のうちに飽和しうるなら、印 0 -> 3、印 2 -> 4(値は変えない)
        c3 = self.chic(a, 3)
        cand = ((reg == 0) | (reg == 2)) & (out > c3) & ~big
        if cand.any():
            i4 = np.flatnonzero(cand)
            g433 = self.g(a[i4], np.minimum(out[i4], 0.995), "433") * mu_s[i4]
            lnN3 = np.log(np.maximum((out[i4] - c3[i4]) / 3, 1e-300)) + lnGM2[i4]
            grows = (g433 > 0) & (g433 * T[i4] > lnN3)
            reg[i4[grows]] = np.where(reg[i4[grows]] == 0, 3, 4)
        return out.reshape(shape), reg.reshape(shape)
