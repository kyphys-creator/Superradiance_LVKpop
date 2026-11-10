"""chi_SR: 自己相互作用のある ALP の超放射で、年齢 T のあとの BH のスピン(単体版、numpy と scipy だけ)。

使い方:
    from chi_sr import ChiSR
    model = ChiSR("chi_sr_tables.npz")
    chi_T, regime = model(M_msun, chi0, T_yr, mu_eV, f_GeV)   # どれも配列でよい(ブロードキャスト)

regime(印):
    0 = 超放射しない(chi0 が 211 と 322 の臨界以下。m >= 3 は入れていない)
    1 = 自己相互作用の準平衡(211・322 の二準位の固定点、alpha <= 0.15)。322 の臨界まで下がりうる
    2 = 自由場(211 が 1 個の量子から飽和して臨界へ。そのあと f <= 1e17 GeV は固定点の降下、それより上は 322 を
        1 個の量子から。BH の質量の減少を入れる)、または 322 から先に育つ BH で自己相互作用が飽和までに効かない場合
    3 = この模型の外。chi0 を返す: 322 から先に育ち自己相互作用が効く BH、準平衡で alpha < 0.01 または > 0.30、alpha > 0.5
    4 = 準平衡だが 0.15 < alpha <= 0.30(433・411 などが育つので二準位の近似は粗い)。値は計算して返す

模型の要点(詳しくは README.md):
  準平衡は Baryakhtar+ 2011.11646 の式 34--36 の右辺を 0 とおいた固定点で、322 の超放射も入れる。
  dchi/dt = mu (f/M_pl)^2 H(alpha, chi) なので、I(alpha, chi) = int_chi^0.995 dchi'/|H| を表にし、
  I(alpha, chi_T) = I(alpha, chi_0') + mu (f/M_pl)^2 (T - t_grow) を chi_T について解く。
  領域の判定: 準平衡の占有数 eps_eq と行き過ぎの占有数 eps_th が、どちらも自由場で飽和する量 chi0 - chi_c より小さければ準平衡。
  準平衡では、1 個の量子から eps_eq まで育つ時間 t_grow を年齢から引き、行き過ぎで失う eps_th を初期スピンから引く。
超放射の率は Kerr 時空の連分数法(Dolan 0705.2880)の表(211: alpha >= 0.03、322: 0.1 <= alpha <= 0.52)、
範囲の外は解析式(Baumann+ 1908.10370 の式 2.32--2.34)。自己相互作用の率は主要項(相対論的な補正なし)。
"""
from math import factorial

import numpy as np
from scipy.interpolate import RegularGridInterpolator


class ChiSR:
    def __init__(self, path):
        z = np.load(path)
        self.z = z
        for k in ["HBAR_EVS", "TSUN", "N0", "YR", "MPL_GEV", "KAPPA_BH", "KAPPA_INF", "KAPPA_BH_544", "CHI_TOP",
                  "ALPHA_2LEVEL", "ALPHA_MAX", "F_FP_MAX"]:
            setattr(self, k, float(z[k]))
        self.ALPHAS = z["ALPHAS"]
        self.W = z["W"]
        la = np.log(self.ALPHAS)
        self._lnI = RegularGridInterpolator((la, self.W), z["lnI"], bounds_error=False, fill_value=None)
        self._winv = RegularGridInterpolator((la, z["Y"]), z["winv"], bounds_error=False, fill_value=None)
        self.Y = z["Y"]
        self.chic_tab = {1: z["chic"], 2: z["chic2"]}
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
        """超放射の率 / mu(占有数の率)。level = "211" か "322"。"""
        n, l, m = (2, 1, 1) if level == "211" else (3, 2, 2)
        a, chi = np.broadcast_arrays(np.asarray(a, float), np.asarray(chi, float))
        out = np.array(self._gamma_analytic(a, chi, n, l, m), float)
        lo, hi = self._range[level]
        inr = (a >= lo) & (a <= hi)
        if inr.any():
            F = self._rel[level](np.stack([np.log(a[inr]), chi[inr]], -1))
            out[inr] = 2.0 * F * a[inr] ** (4 * l + 5) / a[inr]
        return out

    def chic(self, a, m, n):
        """臨界スピン。表の範囲では表(相対論的な率の 0)、外は水素原子の主要項。alpha_E/m >= 1/2 では 1。"""
        a = np.asarray(a, float)
        ap = a * (1 - a ** 2 / (2 * n ** 2))
        out = np.where(ap / m < 0.5, 4 * m * ap / (m ** 2 + 4 * ap ** 2), 1.0)
        if m in (1, 2):
            inr = (a >= self.ALPHAS[0]) & (a <= self.ALPHAS[-1])
            out = np.where(inr, np.interp(np.log(np.clip(a, self.ALPHAS[0], self.ALPHAS[-1])),
                                          np.log(self.ALPHAS), self.chic_tab[m]), out)
        return out

    # ---------- 準平衡(固定点) ----------
    def _chi_of_w(self, w, c1, c2):
        return np.where(w >= 1, c1 + (w - 1) * (self.CHI_TOP - c1), c2 + w * (c1 - c2))

    def _chi_after(self, a, chi0, tau):
        la = np.log(self.ALPHAS)
        lal = np.clip(np.log(a), la[0], la[-1])
        c1 = np.interp(lal, la, self.chic_tab[1])
        c2 = np.interp(lal, la, self.chic_tab[2])
        with np.errstate(divide="ignore", invalid="ignore"):
            w0 = np.where(chi0 >= c1, 1 + (chi0 - c1) / (self.CHI_TOP - c1), (chi0 - c2) / (c1 - c2))
        w0 = np.clip(w0, self.W[0], self.W[-1])
        lnI0 = self._lnI(np.stack([lal, w0], -1))
        yT = np.clip(np.logaddexp(lnI0, np.log(np.maximum(tau, 1e-300))), self.Y[0], self.Y[-1])
        wT = np.clip(self._winv(np.stack([lal, yT], -1)), self.W[0], self.W[-1])
        out = self._chi_of_w(wT, c1, c2)
        keep = (chi0 <= c2) | (np.log(a) < la[0]) | (np.log(a) > la[-1])
        return np.where(keep, chi0, np.minimum(out, chi0))

    # ---------- 自由場 ----------
    def _sat_state(self, M, chi, mu, m, n):
        a_sun = self.alpha(1.0, mu)
        J0 = chi * M ** 2

        def crit(Mf):
            a = a_sun * Mf
            ap = a * (1 - a ** 2 / (2 * n ** 2))
            return a, np.where(ap / m < 0.5, 4 * m * ap / (m ** 2 + 4 * ap ** 2), 1.0)

        def F(Mf):
            a, cf = crit(Mf)
            return Mf - M + a_sun * (1 - a ** 2 / (2 * n ** 2)) * (J0 - cf * Mf ** 2) / m

        lo, hi = 0.5 * M, M.copy()
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            neg = F(mid) < 0
            lo = np.where(neg, mid, lo)
            hi = np.where(neg, hi, mid)
        Mf = 0.5 * (lo + hi)
        return Mf, crit(Mf)[1]

    def _free(self, M, chi0, a, mu, mu_s, T, f_gev, g1):
        M1, c1 = self._sat_state(M, chi0, mu, 1, 2)
        c1 = np.minimum(c1, chi0)
        lnN1 = np.log(np.maximum(chi0 * M ** 2 - c1 * M1 ** 2, 1e-300)) + np.log(self.N0)
        tsat1 = lnN1 / g1
        M2, c2 = self._sat_state(M1, c1, mu, 2, 3)
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
            start = self.chic(a1[fp], 1, 2)
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
        c1 = self.chic(a, 1, 2)
        c2 = self.chic(a, 2, 3)
        big = a > self.ALPHA_MAX
        out = chi0.copy()
        reg = np.zeros(len(M), int)
        reg[big] = 3
        sr211 = (chi0 > c1) & ~big
        g1 = np.where(sr211, self.g(a, np.minimum(chi0, 0.995), "211"), 0.0) * mu_s
        ginf = self.KAPPA_INF * a ** 8 * (self.MPL_GEV / f_gev) ** 4 * mu_s
        gbh = self.KAPPA_BH * a ** 11 * (self.MPL_GEV / f_gev) ** 4 * self.rplus(np.minimum(chi0, 0.999)) * mu_s
        emax = chi0 - c1
        with np.errstate(divide="ignore", invalid="ignore"):
            eeq = 2 / np.sqrt(3) * np.sqrt(g1 * ginf) / gbh
            eth = np.sqrt(2 * g1 * lnGM2 / gbh)
        si_reg = sr211 & (g1 > 0) & (eeq < emax) & (eth < emax)
        inval = si_reg & ((a > self.ALPHAS[-1]) | (a < self.ALPHAS[0]))
        reg[inval] = 3
        q = si_reg & ~inval
        if q.any():
            tgrow = (np.log(np.maximum(eeq[q], 1e-300)) + lnGM2[q]) / g1[q]
            Teff = np.maximum(T[q] - tgrow, 0.0)
            tau = mu_s[q] * (f_gev[q] / self.MPL_GEV) ** 2 * Teff
            c0eff = np.maximum(chi0[q] - eth[q], c1[q])
            out[q] = np.where(Teff > 0, self._chi_after(a[q], c0eff, tau), chi0[q])
            reg[q] = np.where(a[q] > self.ALPHA_2LEVEL, 4, 1)
        fr = sr211 & ~si_reg & (g1 > 0)
        if fr.any():
            out[fr] = self._free(M[fr], chi0[fr], a[fr], mu_eV[fr], mu_s[fr], T[fr], f_gev[fr], g1[fr])
            reg[fr] = 2
        only322 = (~sr211) & (chi0 > c2) & ~big
        if only322.any():
            i3 = np.flatnonzero(only322)
            c0 = np.minimum(chi0[i3], 0.995)
            g2 = self.g(a[i3], c0, "322") * mu_s[i3]
            g1d = np.abs(self.g(a[i3], c0, "211")) * mu_s[i3]
            gbh2 = self.KAPPA_BH_544 * a[i3] ** 11 * (self.MPL_GEV / f_gev[i3]) ** 4 * self.rplus(c0) * mu_s[i3]
            emax2 = (chi0[i3] - c2[i3]) / 2
            with np.errstate(divide="ignore", invalid="ignore"):
                eth2 = np.sqrt(2 * g2 * lnGM2[i3] / gbh2)
                epump = np.sqrt(g1d / ginf[i3])
                tsat2 = (np.log(np.maximum(emax2, 1e-300)) + lnGM2[i3]) / g2
            free2 = (g2 > 0) & (eth2 >= emax2) & (epump >= emax2)
            out[i3] = np.where(free2 & (T[i3] >= tsat2), c2[i3], chi0[i3])
            reg[i3] = np.where(free2, 2, 3)
        return out.reshape(shape), reg.reshape(shape)
