"""BH の年齢(生まれてから合体までの時間)の分布の模型: 1 世代目と 2 世代目(numpy だけで動く)。

使い方:
    from ages import AgeModel
    am = AgeModel()                          # 既定の模型(下の「模型」)
    T = am.nodes(2, z_merge, K=32)           # 2 世代目の年齢 [yr]。形は z_merge の形 + (K,)。等確率の分位点(求積の節)
    chi, regime = chi_sr_model(M[:, None], chi0[:, None], T, mu, f)     # 年齢について平均すれば、年齢の周辺化になる
ほかに am.quantile、am.cdf、am.pdf_lnT、am.sample、am.tmax、am.cosmic_time。時間と年齢の単位はすべて年。

模型(根拠と数値は研究ノート E22、README.md):
  年齢 T は、BH が生まれてから、その BH を含む連星が合体するまでの時間。合体の赤方偏移 z_m での宇宙の年齢を t_m とする。
  1 世代目(gen = 1): T は星形成から合体までの遅れ(BH は星形成の 3--20 Myr あとに生まれるが、その差は無視する)。
      p_1(T | z_m) ∝ T^(-kappa1) psi(t_m - T),    tmin1 <= T <= t_m - t_start
    T^(-kappa1) は遅れの分布(kappa1 = 1、tmin1 = 10 Myr が孤立連星の標準の形)、psi は宇宙の星形成率
    (Madau & Dickinson 2014 の式 15)、t_start は星形成の始まり(z_start = 20)。kappa1 = 0 は古い星団のような年齢に寄る。
  2 世代目(gen = 2): T は、親の合体(この BH ができた時刻)から次の合体までの時間。
      p_2(T | z_m) ∝ p_D(T) R_1(t_m - T),    t2lo <= T <= min(t2hi, t_m - t_start - tmin1)
    p_D は時間の分布の形で、family2 = "loguniform"(1/T。既定)か "lognormal"(ln T が中央値 t2med、幅 s2 dex の正規分布)。
    R_1(t) = int T'^(-kappa1) psi(t - T') dT' は、同じ 1 世代目の模型の合体率(親の合体が起きた時刻の重み)。
  sfr_weight = False で psi = 1、parent_weight = False で R_1 = 1(重みを外した形。先行研究の「年齢は対数一様」を再現するとき)。
  宇宙論は平坦な LCDM(H0 = 67.4 km/s/Mpc、Omega_m = 0.315。Planck 2018 の要旨の値)。

分位点と累積分布は、合体の時刻と ln T の表から補間する(誤差は累積分布で 3e-4 以下、分位点の年齢で 0.5% 以下。年齢の模型そのものの不定性よりずっと小さい)。
"""
import numpy as np

GYR_PER_INV_KMSMPC = 977.7922216807891      # 1 / (1 km/s/Mpc) [Gyr]
N_M, N_X, N_P = 160, 800, 801               # 表の大きさ: 合体の時刻、ln T、確率


def _cumtrapz(y, x):
    """最後の軸に沿った台形の累積(先頭は 0)。"""
    d = 0.5 * (y[..., 1:] + y[..., :-1]) * np.diff(x, axis=-1)
    return np.concatenate([np.zeros(y.shape[:-1] + (1,)), np.cumsum(d, axis=-1)], axis=-1)


class AgeModel:
    def __init__(self, kappa1=1.0, tmin1=1e7, family2="loguniform", t2lo=1e6, t2hi=1e10, t2med=1e8, s2=1.0,
                 sfr_weight=True, parent_weight=True, z_start=20.0, H0=67.4, Om=0.315):
        if family2 not in ("loguniform", "lognormal"):
            raise ValueError("family2 は 'loguniform' か 'lognormal'")
        if not (tmin1 > 0 and t2lo > 0 and t2hi > t2lo and t2med > 0 and s2 > 0):
            raise ValueError("tmin1, t2lo, t2med, s2 は正、t2hi > t2lo")
        self.kappa1, self.tmin1 = float(kappa1), float(tmin1)
        self.family2, self.t2lo, self.t2hi, self.t2med, self.s2 = family2, float(t2lo), float(t2hi), float(t2med), float(s2)
        self.sfr_weight, self.parent_weight = bool(sfr_weight), bool(parent_weight)
        self.z_start, self.H0, self.Om = float(z_start), float(H0), float(Om)
        self._hub = self.H0 / GYR_PER_INV_KMSMPC * 1e-9          # H0 [1/yr]
        self.t_start = float(self.cosmic_time(self.z_start))
        self.t0 = float(self.cosmic_time(0.0))
        self._build()

    def describe(self):
        """模型のパラメータ(記録用)。"""
        return {k: getattr(self, k) for k in ("kappa1", "tmin1", "family2", "t2lo", "t2hi", "t2med", "s2", "sfr_weight",
                                              "parent_weight", "z_start", "H0", "Om")}

    # ---------- 宇宙論 ----------
    def cosmic_time(self, z):
        """赤方偏移 z での宇宙の年齢 [yr](平坦な LCDM、放射は無視)。"""
        ol = 1.0 - self.Om
        return 2.0 / (3.0 * self._hub * np.sqrt(ol)) * np.arcsinh(np.sqrt(ol / self.Om) * (1.0 + np.asarray(z, float)) ** -1.5)

    def z_at(self, t):
        """宇宙の年齢 t [yr] での赤方偏移(cosmic_time の逆)。"""
        ol = 1.0 - self.Om
        return (np.sqrt(self.Om / ol) * np.sinh(1.5 * self._hub * np.sqrt(ol) * np.asarray(t, float))) ** (-2.0 / 3.0) - 1.0

    def sfr(self, t):
        """宇宙の年齢 t [yr] での星形成率(規格化は任意)。Madau & Dickinson 2014(arXiv:1403.0007)の式 15。"""
        if not self.sfr_weight:
            return np.ones_like(np.asarray(t, float))
        z1 = 1.0 + self.z_at(t)
        return z1 ** 2.7 / (1.0 + (z1 / 2.9) ** 5.6)

    # ---------- 分布の形 ----------
    def _w1(self, T, tm):
        """1 世代目の d P / d ln T(規格化なし)。"""
        return T ** (1.0 - self.kappa1) * self.sfr(tm - T)

    def _shape2(self, T):
        """2 世代目の時間の分布の形 T p_D(T)(d / d ln T、規格化なし)。"""
        if self.family2 == "loguniform":
            return np.ones_like(T)
        return np.exp(-0.5 * (np.log(T / self.t2med) / (self.s2 * np.log(10.0))) ** 2)

    def _r1(self, t):
        """1 世代目の合体率 R_1(t)(規格化なし)。表を ln(t - t_start) について補間する。"""
        if not self.parent_weight:
            return np.ones_like(np.asarray(t, float))
        return np.exp(np.interp(np.log(np.maximum(t - self.t_start, self.tmin1)), self._ur, self._lnr1))

    def _w2(self, T, tm):
        return self._shape2(T) * self._r1(tm - T)

    def tmax(self, gen, z):
        """年齢の上限 [yr]。1 世代目は t_m - t_start、2 世代目は min(t2hi, t_m - t_start - tmin1)。"""
        d = self.cosmic_time(z) - self.t_start
        return d if gen == 1 else np.minimum(self.t2hi, d - self.tmin1)

    def _range(self, gen, d):
        """ln T の範囲 (下端, 上端)。d = t_m - t_start。"""
        if gen == 1:
            return np.log(self.tmin1) * np.ones_like(d), np.log(d)
        return np.log(self.t2lo) * np.ones_like(d), np.log(np.minimum(self.t2hi, d - self.tmin1))

    # ---------- 表 ----------
    def _build(self):
        span = self.t0 - self.t_start
        # 1 世代目の合体率 R_1(t): t - t_start = tmin1 から span まで
        self._ur = np.linspace(np.log(self.tmin1), np.log(span), 1200)
        d = np.exp(self._ur)
        x = np.linspace(0.0, 1.0, 600)
        lnT = np.log(self.tmin1) + x[None, :] * (self._ur[:, None] - np.log(self.tmin1))
        T = np.exp(lnT)
        r1 = _cumtrapz(self._w1(T, self.t_start + d[:, None]), lnT)[:, -1]
        r1[0] = r1[1] * 1e-6                       # t - t_start = tmin1 では 0。対数をとるために小さい値に置き換える
        self._lnr1 = np.log(np.maximum(r1, 1e-300))
        # 合体の時刻の格子: d = t_m - t_start。下端は、どちらの世代も年齢の幅が開く所の少し上
        dmin = {1: self.tmin1 * 1.02, 2: (self.t2lo + self.tmin1) * 1.02}
        self._tab = {}
        p = 0.5 * (1.0 - np.cos(np.pi * np.linspace(0.0, 1.0, N_P)))      # 確率の格子は両端に詰める(裾の分位点のため)
        for gen in (1, 2):
            if dmin[gen] >= span:
                raise ValueError("年齢の幅が開く合体の時刻がない(tmin1 か t2lo が大きすぎる)")
            u = np.linspace(np.log(dmin[gen]), np.log(span), N_M)
            dd = np.exp(u)
            lo, hi = self._range(gen, dd)
            xg = np.linspace(0.0, 1.0, N_X)
            lnT = lo[:, None] + xg[None, :] * (hi - lo)[:, None]
            w = (self._w1 if gen == 1 else self._w2)(np.exp(lnT), (self.t_start + dd)[:, None])
            c = _cumtrapz(w, lnT)
            norm = c[:, -1].copy()
            if not np.all(norm > 0):
                raise ValueError("年齢の分布の規格化が 0 になる合体の時刻がある")
            c = c / norm[:, None]
            c[:, -1] = 1.0
            q = np.empty((N_M, N_P))               # 分位点: 確率 p での ln T
            for i in range(N_M):
                q[i] = np.interp(p, c[i], lnT[i])
            self._tab[gen] = dict(u=u, cdf=c, q=q, lnnorm=np.log(norm))

    def _locate(self, gen, z):
        """合体の時刻の格子での位置 (番号, 重み, d)。格子の外(合体が早すぎて年齢の幅がない)は NaN。"""
        d = np.asarray(self.cosmic_time(z), float) - self.t_start
        u = self._tab[gen]["u"]
        with np.errstate(invalid="ignore", divide="ignore"):
            lu = np.log(np.where(d > 0, d, np.nan))
        ok = (lu >= u[0]) & (lu <= u[-1] * (1 + 1e-12))
        s = np.clip((np.where(ok, lu, u[0]) - u[0]) / (u[1] - u[0]), 0.0, N_M - 1.0)
        i = np.minimum(s.astype(int), N_M - 2)
        return i, s - i, d, ok

    # ---------- 呼び出し ----------
    def quantile(self, gen, z, p):
        """確率 p の分位点の年齢 [yr]。z と p はブロードキャストする。"""
        z, p = np.broadcast_arrays(np.asarray(z, float), np.asarray(p, float))
        i, w, d, ok = self._locate(gen, z)
        q = self._tab[gen]["q"]
        pc = np.clip(p, 0.0, 1.0)
        s = np.arccos(1.0 - 2.0 * pc) / np.pi * (N_P - 1)                   # _build の確率の格子(両端に詰めたもの)での位置
        j = np.minimum(s.astype(int), N_P - 2)
        pj = 0.5 * (1.0 - np.cos(np.pi * j / (N_P - 1)))
        pk = 0.5 * (1.0 - np.cos(np.pi * (j + 1) / (N_P - 1)))
        v = (pc - pj) / (pk - pj)                                            # 節のあいだは確率について線形
        out = ((1 - w) * ((1 - v) * q[i, j] + v * q[i, j + 1]) + w * ((1 - v) * q[i + 1, j] + v * q[i + 1, j + 1]))
        return np.where(ok, np.exp(out), np.nan)

    def nodes(self, gen, z, K=32):
        """年齢の求積の節 [yr]: 確率 (k + 1/2)/K の分位点(k = 0..K-1)。形は z の形 + (K,)。重みはどれも 1/K。"""
        return self.quantile(gen, np.asarray(z, float)[..., None], (np.arange(K) + 0.5) / K)

    def sample(self, gen, z, rng, size=None):
        """年齢 [yr] の標本。size を渡さなければ z の形。"""
        z = np.asarray(z, float)
        return self.quantile(gen, z, rng.random(z.shape if size is None else size))

    def cdf(self, gen, z, T):
        """年齢が T [yr] 以下である確率。"""
        z, T = np.broadcast_arrays(np.asarray(z, float), np.asarray(T, float))
        i, w, d, ok = self._locate(gen, z)
        lo, hi = self._range(gen, np.where(ok, d, self.t0 - self.t_start))
        with np.errstate(invalid="ignore", divide="ignore"):
            x = np.clip((np.log(T) - lo) / (hi - lo), 0.0, 1.0)
        s = x * (N_X - 1)
        j = np.minimum(s.astype(int), N_X - 2)
        v = s - j
        c = self._tab[gen]["cdf"]
        out = ((1 - w) * ((1 - v) * c[i, j] + v * c[i, j + 1]) + w * ((1 - v) * c[i + 1, j] + v * c[i + 1, j + 1]))
        return np.where(ok, out, np.nan)

    def pdf_lnT(self, gen, z, T):
        """確率密度 d P / d ln T(年齢 T [yr]。範囲の外は 0)。"""
        z, T = np.broadcast_arrays(np.asarray(z, float), np.asarray(T, float))
        i, w, d, ok = self._locate(gen, z)
        dd = np.where(ok, d, self.t0 - self.t_start)
        lo, hi = self._range(gen, dd)
        ln = self._tab[gen]["lnnorm"]
        norm = np.exp((1 - w) * ln[i] + w * ln[i + 1])
        with np.errstate(invalid="ignore", divide="ignore"):
            inside = (np.log(T) >= lo) & (np.log(T) <= hi)
        Tc = np.where(inside, T, np.exp(lo))
        val = (self._w1 if gen == 1 else self._w2)(Tc, self.t_start + dd) / norm
        return np.where(ok, np.where(inside, val, 0.0), np.nan)
