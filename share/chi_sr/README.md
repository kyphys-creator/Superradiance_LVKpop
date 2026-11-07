# chi_SR: 自己相互作用のある ALP の超放射で、年齢 T のあとの BH のスピン

集団解析(モンテカルロの標本ごとに呼ぶ)のための速い関数。numpy と scipy だけで動く。未発表の研究の一部なので、共有は共同研究者のあいだに限る。

## 使い方

```python
from chi_sr import ChiSR
model = ChiSR("chi_sr_tables.npz")
chi_T, regime = model(M_msun, chi0, T_yr, mu_eV, f_GeV)   # どれも配列でよい(ブロードキャスト)
```

| 入力 | 単位 | 意味 |
|---|---|---|
| `M_msun` | $M_\odot$ | BH の質量 |
| `chi0` | — | 生まれたときのスピン |
| `T_yr` | 年 | 年齢(生まれてから合体まで) |
| `mu_eV` | eV | ALP の質量 $\mu$ |
| `f_GeV` | GeV | ALP の崩壊定数 $f$(自己相互作用 $\lambda=\mu^2/f^2$、引力) |

出力の `chi_T` は年齢 $T$ のあとのスピン、`regime` はどの領域で計算したかの印である。

| `regime` | 意味 |
|---|---|
| 0 | 超放射しない($\chi_0$ が 322 の臨界以下)。$\chi_T=\chi_0$ |
| 1 | 自己相互作用の準平衡(211・322 の二準位) |
| 2 | 自由場(211 が飽和して臨界へ、時間があれば 322 も) |
| 3 | この模型の外。$\chi_T=\chi_0$ を返す(下の「限界」) |

速さ: $10^6$ 個の標本で約 0.4 秒(Mac Studio、numpy)。

## 模型

- **重力の微細構造定数**: $\alpha=GM\mu/\hbar c=0.0748\,(M/10M_\odot)(\mu/10^{-12}\,\mathrm{eV})$。
- **超放射の率**: Kerr 時空の Klein–Gordon 方程式を連分数法(Dolan, arXiv:0705.2880)で解いた表を使う。範囲は 211 で $\alpha\ge0.03$、322 で $0.1\le\alpha\le0.52$。範囲の外は解析式(Baumann+ arXiv:1908.10370 の式 2.32–2.34)。占有数の率で、振幅の率の 2 倍。
- **自己相互作用の過程**: $211\times211\to322\times\mathrm{BH}$ と $322\times322\to211\times\infty$(Baryakhtar+ arXiv:2011.11646、以下 BGLS)。占有数を $\varepsilon=N/(GM^2/\hbar c)$ で数えると、率は次のとおり($M_\mathrm{Pl}=1.22\times10^{19}$ GeV)。
  - $\gamma_\mathrm{BH}=4.3\times10^{-7}\alpha^{11}(M_\mathrm{Pl}/f)^4\tilde r_+\mu$
  - $\gamma_\infty=1.1\times10^{-8}\alpha^8(M_\mathrm{Pl}/f)^4\mu$
- **準平衡**(BGLS の式 55、36): $\varepsilon_{211}^\mathrm{eq}=(2/\sqrt3)\sqrt{\gamma_\mathrm{SR}\gamma_\infty}/\gamma_\mathrm{BH}$ で、$\dot\chi=-\gamma_\mathrm{SR}\varepsilon_{211}^\mathrm{eq}=\mu(f/M_\mathrm{Pl})^2H(\alpha,\chi)$。
  - $f$ と $\mu$ は前の因子にしか入らない。そこで $I(\alpha,\chi)=\int_\chi^{0.995}d\chi'/|H|$ を表にし、$I(\alpha,\chi_T)=I(\alpha,\chi_0)+\tau$($\tau=\mu(f/M_\mathrm{Pl})^2T$)を解く。
  - スピンは $(\alpha,\chi_0,\tau)$ の 3 つの量だけで決まる。
- **領域の判定**(BGLS の式 53・56): 次の二つが、どちらも自由場で飽和する量 $\chi_0-\chi_c$ より小さければ準平衡、そうでなければ自由場とする。
  - 準平衡の占有数 $\varepsilon_\mathrm{eq}$
  - 行き過ぎの占有数 $\varepsilon_\mathrm{th}=\sqrt{2\gamma_\mathrm{SR}\ln(GM^2)/\gamma_\mathrm{BH}}$
- **準平衡の側**: 1 個の量子から $\varepsilon_\mathrm{eq}$ まで育つ時間を年齢から引き、行き過ぎで失う $\varepsilon_\mathrm{th}$ を初期スピンから引く。
- **自由場の側**: 1 個の量子から $\ln N_\mathrm{sat}\approx180$ 回の $e$-fold で飽和し、臨界スピン($m\Omega_H=E_n$)まで落ちる。BH の質量の減少を入れる。322 は $t=0$ から $e$-fold を積む。

## 照合

研究ノート E11 にまとめたもの。
- **運動論の ODE**(BGLS の式 34–36 を 1 個の量子から解いたもの)との比較: 対象は GW241011 の主星($19.6M_\odot$、$\chi_0=0.78$、$\mu=0.71\times10^{-12}$ eV)で、$f=10^{11}$–$10^{15}$ GeV、$T=10^6$–$10^8$ 年の範囲。差は最大 $1.6\times10^{-4}$。$15$–$25M_\odot$、$f=10^{12}$–$10^{13}$ GeV では、$10^{10}$ 年まで小数 3 桁で一致した。
- **準平衡の表と直接の積分**の比較: 差は最大 $1.8\times10^{-4}$。
- **自由場の半解析**との比較($f=10^{19}$ GeV、無作為な 400 点): 差の中央値は $3\times10^{-7}$、93.5% が 0.01 以内、最大 0.034。残りの差は、臨界スピンを解くときのエネルギーの近似の違いによる。

## 限界(使う前に読む)

1. **重い BH(印 3)**: 211 が超放射せず 322 から育つ BH の自己相互作用は、模型に入っていない。$\mu=0.71\times10^{-12}$ eV では $\gtrsim35M_\odot$($\alpha\gtrsim0.19$)の BH にあたり、GWTC-5 の $45M_\odot$ 以上の速く回る群が含まれる。準平衡で $\alpha>0.30$ の BH も印 3 にした。
2. **領域の境目**: $f\approx10^{16}$ GeV 付近(BGLS の領域 B)では、自由場の側の 322 の段階と準平衡の側の扱いが合わず、スピンが最大 0.2 ほど飛ぶ。窓($f\sim10^{12}$ GeV)には効かない。
3. **自己相互作用の率は主要項**: Kerr 時空での補正は、$211\times211\to322\times\mathrm{BH}$ で $\alpha=0.1$ のとき約 2.5 倍、$0.2$ のとき約 5 倍(Witte & Mummery arXiv:2412.03655)で、まだ入れていない。入れると回転を失う時間はその分長くなる。
4. **その他の近似**: 二準位の近似は $\alpha\lesssim0.15$ で信頼できる。準平衡の側では BH の質量を固定した(補正は回転を失う時間で 2–5%)。ガスの降着は入れていない。
5. **入力は呼ぶ側が決める**: 生まれたときのスピンと年齢の分布(世代ごと)は、集団解析の側で決める。

## 作り直し

表はプロジェクトの `research/calc/export_chisr.py` が作る。その際、単体版とプロジェクト版が一致することも確かめる。
