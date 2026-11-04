# Astrophysical inputs: 2G natal spin and mass, and BH ages by generation (full-text check, 2026-11-03)

Every number below comes from an arXiv **full text** (LaTeX source and/or the compiled arXiv PDF) downloaded on
2026-11-03 to `Refs/arxiv/<ID>/` (`src/` = extracted e-print, `<ID>.pdf`, `<ID>.txt` = `pdftotext -layout`, `raw.txt` = plain
`pdftotext`). Equation, figure, table and section numbers follow the arXiv PDF unless marked "(source order)".
`2510.26931/` and `2107.11730/` already existed; they were read as is (no compiled PDF is present in `2510.26931/`, so its
numbers are "(source order)").

Legend:
- no tag = stated in the text of the paper at the location given.
- **[fig]** = read by eye from a figure (rendered from the source PDF/EPS); approximate (factor ~1.5 in log-axis readings).
- **[inference]** = my reading/derivation from the paper's own statements; not stated by the paper.
- **unconfirmed** = not found in the full text.
- Quotes are short verbatim fragments in "…".

Definitions used here:
- $\Delta t_{\rm 2G}$ ("2G age at merger") $\equiv t_{\rm merge}(\text{binary containing the 2G BH}) - t_{\rm merge}(\text{its 1G+1G parent binary})$.
- $\tau_{\rm 1G}$ ("1G age at merger") $\equiv t_{\rm merge} - t_{\rm BH\ formation}$. Most papers give the delay from star/cluster formation;
  BH formation lags star formation by only $\sim3$–20 Myr (Mapelli+ 2103.05016 Sec. 6.3: original binaries form the second BH at
  "$t_{\rm form}\approx3$ and $\approx20$ Myr"; Baibhav+ 2004.00650 Sec. II A (source order): stellar lifespans "of order $\mathcal O(1)$ Myr"
  neglected).

---

## 要点 (Japanese summary)

- **終スピンの式**: HBR 2016 (1605.01938) の一般歳差スピン版 Eqs. (12)–(18) と Table 1 上段 ($n_M=1,n_J=2$、6 係数、$\xi=0.41616$)。等質量・無スピンで厳密に 0.68646 (Eq. 11)。終質量は Barausse+ 2012 (1206.3803) Eq. (18)(同じ $\tilde a=a_{\rm tot}(\beta,\gamma,q)$ を使うので HBR と整合)。Baibhav+ 2020 も同じ組み合わせ (質量 Barausse 2012、スピン HBR) を使っている。
- **2G の生まれたときのスピン**: 全残骸では $\chi\simeq0.69$ に鋭いピーク、裾は $0.4\lesssim\chi\lesssim0.9$ (Fishbach+ 2017, $q\ge0.7$)。CMC で 1G 無スピンなら 90% が 0.67–0.69 (Rodriguez+ 2019 Fig. 8)。1G スピン 0.5 なら中央値 0.62、90% が 0.53–0.72。ただし**キックで残る 2G に限ると分布が歪む** (Borchers+ 2025): 1G スピン $a_i\gtrsim0.3$ でピークが 0.69 より上 ($a_i=0.6$ で $\sim0.8$)、脱出速度が低いと $\sim0.6$ と $\sim0.8$ の二峰。$q\sim1/3$ の 1G+1G は $\chi\approx0.55$。
- **2G の年齢 $\Delta t_{\rm 2G}$ の分布を直接示した論文は見つからなかった。** 最も近いもの: (i) Baibhav+ 2020 は「2g の赤方偏移分布は 1g+1g とほぼ同じ」「重いのでごく短時間で合体」とし、時間を Eq. (26) で $t_{\rm MS}\to t_{\rm MS}(v_{\rm kick}/v_{\rm esc})^3$ として計算 (図 2 から $(50+50)M_\odot$ で最短 $\sim5\times10^7$ yr [fig])。(ii) Mapelli+ 2021 Fig. 1–2 の $t_{\rm dyn,ng}$ は 0.3–0.6 Myr から 7 Gyr、ピーク 10–100 Myr [fig; 増分時間と解釈、inference]。(iii) N 体の具体例: Dragon-II で 2g IMBH が「わずか 25.5 Myr」で再合体。(iv) AGN: 円盤寿命 (中央値 $\sim$1–2 Myr、0.1–10 Myr) 以内で全世代が合体、2g は $<0.1$ Myr で出現 (McFACTS)、Ng の移動時間ピーク $10^{-2}$ Myr (Vaccaro+ 2024)。
- **1G の年齢**: 古い GC 起源の低赤方偏移合体は $\gtrsim6$ Gyr (Kremer+ 2020 Fig. 15)。放出連星は中央値 $\sim10$ Gyr、クラスター内合体は即座 (Fishbach & Kalogera 2021 の総説部分)。孤立連星は $p(\tau)\propto\tau^{-1}$、$\tau_{\min}\sim10$ Myr が典型予言、GWTC-2 からは $\tau_{50\%}=2.8^{+3.3}_{-2.6}$ Gyr。
- **GW241011 (2510.26931)**: 環境は「球状・核・若い大質量星団」、金属量 $Z\gtrsim0.1Z_\odot$ の星団と整合。2G の年齢や合体時刻の推定は**書かれていない**。超放射の制限は「控えめに」BH 年齢 $10^5$ yr を仮定。AGN は離心率の文脈でしか触れていない。

---

## 0. Sources read (all arXiv metadata checked via export.arxiv.org on 2026-11-03)

| arXiv (version) | Short ref | Topic | Local dir |
|---|---|---|---|
| 1605.01938v2 | Hofmann, Barausse, Rezzolla 2016 (HBR) | final spin, generic spins | `Refs/arxiv/1605.01938/` |
| 1206.3803v3 | Barausse, Morozova, Rezzolla 2012 | radiated energy / final mass, generic spins | `Refs/arxiv/1206.3803/` |
| 1611.00332v2 | Jiménez-Forteza+ 2017 (UIB) | final spin & $E_{\rm rad}$, aligned spins only | `Refs/arxiv/1611.00332/` |
| 1703.06223v3 | Gerosa & Berti 2017 | 1g/2g populations | `Refs/arxiv/1703.06223/` |
| 1703.06869v3 | Fishbach, Holz & Farr 2017 | hierarchical spin distribution | `Refs/arxiv/1703.06869/` |
| 1906.10260v2 | Rodriguez+ 2019 | CMC 2G mergers | `Refs/arxiv/1906.10260/` |
| 2005.04239v2 | Rodriguez+ 2020 | GW190412 as 3G (CMC) | `Refs/arxiv/2005.04239/` |
| 2503.21278v1 | Borchers, Ye, Fishbach 2025 | kick-selected 2G spins | `Refs/arxiv/2503.21278/` |
| 1911.00018v2 | Kremer+ 2020 | CMC Cluster Catalog | `Refs/arxiv/1911.00018/` |
| 2507.07183v2 | Ye, Fishbach, Kremer, Reina-Campos 2025 | CMC + formation histories | `Refs/arxiv/2507.07183/` |
| 1811.03640v2 | Antonini, Gieles, Gualandris 2019 | semi-analytic hierarchical growth | `Refs/arxiv/1811.03640/` |
| 2208.01081v2 | Antonini, Gieles, Dosopoulou, Chattopadhyay 2023 (cBHBd) | GC hierarchical mergers | `Refs/arxiv/2208.01081/` |
| 2103.05016v2 | Mapelli+ 2021 (FASTCLUSTER) | hierarchical in YSC/GC/NSC | `Refs/arxiv/2103.05016/` |
| 2109.06222v2 | Mapelli+ 2022 | multi-channel incl. isolated | `Refs/arxiv/2109.06222/` |
| 2401.14837v2 | Torniamenti+ 2024 | FASTCLUSTER + cluster evolution | `Refs/arxiv/2401.14837/` |
| 2210.10055v3 | Kritos+ 2022/24 (Rapster) | rapid cluster code | `Refs/arxiv/2210.10055/` |
| 2004.00650v2 | Baibhav, Gerosa, Berti+ 2020 | mass gap / spin gap; 2g time delay | `Refs/arxiv/2004.00650/` |
| 2006.01867v2 | Fragione & Silk 2020 | NSC repeated mergers | `Refs/arxiv/2006.01867/` |
| 2211.09670v4 | Atallah+ 2023 | NSC successive mergers (few-body) | `Refs/arxiv/2211.09670/` |
| 2308.10884v2 | Chattopadhyay, Stegmann, Antonini+ 2023 | NSC (cBHBd) | `Refs/arxiv/2308.10884/` |
| 2307.04806v1 / 2307.04807v1 | Arca Sedda+ 2023 (Dragon-II II / III) | direct N-body | `Refs/arxiv/2307.0480{6,7}/` |
| 2012.00011v2 | Tagawa+ 2021 | AGN hierarchical | `Refs/arxiv/2012.00011/` |
| 2311.18548v1 | Vaccaro+ 2024 | AGN hierarchical, timescales | `Refs/arxiv/2311.18548/` |
| 2410.16515v2 | McKernan+ 2025 (McFACTS I) | AGN Monte Carlo | `Refs/arxiv/2410.16515/` |
| 2105.03439v2 | Gerosa & Fishbach 2021 | review | `Refs/arxiv/2105.03439/` |
| 2105.06491v2 | Fishbach & Kalogera 2021 | BBH delay-time distribution | `Refs/arxiv/2105.06491/` |
| 2107.11730 (existing) | Payne+ 2021 | superradiance + hierarchical | `Refs/arxiv/2107.11730/` |
| 2510.26931 (existing) | LVK 2025 | GW241011 / GW241110 | `Refs/arxiv/2510.26931/` |

---

## 1. NR remnant fits (2G natal spin and mass)

### 1a. Final spin: Hofmann, Barausse & Rezzolla 2016 (1605.01938) — recommended

Conventions [Sec. 2]: $M_1\ge M_2$, $q\equiv M_2/M_1\le1$, $\nu\equiv q/(1+q)^2$, $G=c=1$; $\boldsymbol a_{1,2}$ dimensionless spin vectors;
$\alpha$ = angle between the spins, $\beta$ ($\gamma$) = angle between $\boldsymbol a_1$ ($\boldsymbol a_2$) and $\boldsymbol L$, all defined at the
initial separation $r_{\rm in}$ [Eq. (17)]:
$$\cos\alpha\equiv\hat{\boldsymbol a}_1\cdot\hat{\boldsymbol a}_2|_{r_{\rm in}},\quad \cos\beta\equiv\hat{\boldsymbol L}\cdot\hat{\boldsymbol a}_1|_{r_{\rm in}},\quad \cos\gamma\equiv\hat{\boldsymbol L}\cdot\hat{\boldsymbol a}_2|_{r_{\rm in}}.$$

Kerr ISCO functions [Eqs. (2)–(6)]:
$$E_{\rm ISCO}(a)=\sqrt{1-\frac{2}{3r_{\rm ISCO}(a)}},\qquad L_{\rm ISCO}(a)=\frac{2}{3\sqrt3}\left[1+2\sqrt{3r_{\rm ISCO}(a)-2}\right],$$
$$r_{\rm ISCO}(a)=3+Z_2-\frac{a}{|a|}\sqrt{(3-Z_1)(3+Z_1+2Z_2)},\quad Z_1=1+(1-a^2)^{1/3}\left[(1+a)^{1/3}+(1-a)^{1/3}\right],\quad Z_2=\sqrt{3a^2+Z_1^2}.$$

Aligned-spin ansatz [Eq. (7)], with $a_{\rm tot}=(a_1+a_2q^2)/(1+q)^2$ and $a_{\rm eff}=a_{\rm tot}+\xi\nu(a_1+a_2)$ [text below Eq. (7)]:
$$a_{\rm fin}=a_{\rm tot}+\nu\left[L_{\rm ISCO}(a_{\rm eff})-2a_{\rm tot}\big(E_{\rm ISCO}(a_{\rm eff})-1\big)\right]+\sum_{i=0}^{n_M}\sum_{j=0}^{n_J}k_{ij}\,\nu^{2+i}a_{\rm eff}^{\,j}.$$

**Generic (precessing) spins** — the formula to use:
- Projections [Eqs. (14)–(15)]:
  $$a_{\rm tot}(\beta,\gamma,q)\equiv\frac{|\boldsymbol a_1|\cos\beta+|\boldsymbol a_2|\cos\gamma\,q^2}{(1+q)^2},\qquad a_{\rm eff}(\beta,\gamma,q)\equiv a_{\rm tot}(\beta,\gamma,q)+\xi\nu\left(a_1\cos\beta+a_2\cos\gamma\right).$$
- Orbital term [Eq. (13)] (the PDF prints the inner upper limit as $n_S$; it is $n_J$ of Eq. (7)):
  $$|\boldsymbol\ell|=\left|L_{\rm ISCO}(a_{\rm eff})-2a_{\rm tot}\big(E_{\rm ISCO}(a_{\rm eff})-1\big)+\sum_{i=0}^{n_M}\sum_{j=0}^{n_J}k_{ij}\,\nu^{1+i}a_{\rm eff}^{\,j}\right|.$$
- Final-spin magnitude [Eq. (16)]:
  $$|\boldsymbol a_{\rm fin}|=\frac{1}{(1+q)^2}\Big[|\boldsymbol a_1|^2+|\boldsymbol a_2|^2q^4+2|\boldsymbol a_1||\boldsymbol a_2|q^2\cos\alpha+2\big(|\boldsymbol a_1|\cos\beta+|\boldsymbol a_2|q^2\cos\gamma\big)|\boldsymbol\ell|q+|\boldsymbol\ell|^2q^2\Big]^{1/2}.$$
- Vector form [Eq. (12)]: $\boldsymbol a_{\rm fin}=\boldsymbol a_{\rm tot}+\boldsymbol\ell\nu$, $\boldsymbol a_{\rm tot}=(\boldsymbol a_1+\boldsymbol a_2q^2)/(1+q)^2$.
- Optional bias correction [Eq. (18)]: replace $\Theta=\alpha,\beta,\gamma$ by $\Theta^*\equiv2\arctan[(1+\epsilon_\Theta)\tan(\Theta/2)]\approx\Theta+\epsilon_\Theta\sin\Theta$,
  with $\epsilon_\alpha\approx0$, $\epsilon_\beta=\epsilon_\gamma\approx0.024$ [text below Eq. (18)].
- Final-spin direction [Eq. (19)]: $\cos\theta_{\rm fin}=\hat{\boldsymbol J}(r_{\rm in})\cdot\hat{\boldsymbol L}(r_{\rm in})$; checked against 157 NR runs: within
  $\sim18^\circ$ worst case, within $4^\circ$ ($6^\circ$) in 64% (78%) of cases [text after Eq. (19), Fig. 3].

**Coefficients** [Table 1, top block, $n_M=1$, $n_J=2$ — the paper's recommended set, "the optimal choice to avoid overfitting appears to be $n_M=1,\ n_J=2$" (text after Table 2)]:

| $k_{01}$ | $k_{02}$ | $k_{10}$ | $k_{11}$ | $k_{12}$ | $\xi$ |
|---|---|---|---|---|---|
| −1.2019 | −1.20764 | 3.79245 | 1.18385 | 4.90494 | 0.41616 |

$k_{00}$ is fixed by the equal-mass non-spinning constraint [Eq. (11)], using $a_{\rm fin}=0.68646\pm0.00004$ (Scheel+ 2009):
$$\frac{\sqrt3}{2}+\sum_{i=0}^{n_M}\frac{k_{i0}}{4^{2+i}}=0.68646\pm0.00004\;\Rightarrow\;k_{00}=16\left(0.68646-\tfrac{\sqrt3}{2}\right)-\tfrac{k_{10}}{4}\quad(n_M=1),$$
and the text states "the fitted value is $k_{00}\approx-3.82$" [text before Eq. (11)]. (The perturbative expectation is $k_{00}\approx-3.87$, Eq. (10).)
Table 1 middle/bottom blocks give the 16- and 20-coefficient sets ($\xi=0.463926$, $0.474046$); the text notes these "are probably overfitting".

**Accuracy / validity** [abstract; Table 2; text after Eq. (18)]:
- Calibration: 619 NR simulations (246 parallel-spin used for the aligned fit); blind test on 248 GaTech runs (83 aligned + 165 precessing).
- Aligned spins, $n_M=1,n_J=2$: residual mean $-0.000215$, rms $0.00198$, $\chi^2_{\rm red}=0.985$ [Table 2]. Abstract: "$\sigma\approx0.002$ for aligned/anti-aligned binaries and $\sigma\approx0.006$ for generic binaries".
- Generic spins: residual $\mu\approx-0.005$, rms $\approx0.007$ (formula overpredicts); with Eq. (18) angles $\mu\approx-0.001$, $\sigma\approx0.006$ [Fig. 2].
- Blind GaTech test: aligned $\mu\approx-5\times10^{-5}$, $\sigma\approx1.4\times10^{-4}$; precessing $\mu\approx-0.004$ ($-0.0005$ adjusted), $\sigma\approx3.3\times10^{-4}$ ($3.5\times10^{-4}$) [text before Eq. (19)].
- Domain: the EMRL is enforced exactly (Eq. (1)), so the formula is algebraic for **all** $0<q\le1$, $|a_i|\le1$ and any separation $r_{\rm in}$ ("valid also when $r_{\rm in}\gg M_1+M_2$"). Angles are assumed constant during the inspiral, exact only for $q=1$, one spin, or $\nu\approx0$ [text after Eq. (16)]; transitional precession and up–down instability configurations are excluded caveats.
- The explicit $q$ and spin coverage of the 619-run calibration set is **unconfirmed** (not stated in the text).

### 1b. Final mass usable with HBR: Barausse, Morozova & Rezzolla 2012 (1206.3803v3) — recommended

Same projection as HBR: [Eq. (17)]
$$\tilde a\equiv\frac{\hat{\boldsymbol L}\cdot(\boldsymbol S_1+\boldsymbol S_2)}{M^2}=\frac{|\boldsymbol a_1|\cos\beta+q^2|\boldsymbol a_2|\cos\gamma}{(1+q)^2}\;\;(=a_{\rm tot}(\beta,\gamma,q)\text{ of HBR Eq. 14}),$$
with $M=m_1+m_2$, $q=m_2/m_1$, $\nu=m_1m_2/M^2$ [Sec. 2]. Then [Eq. (18)]
$$\frac{E_{\rm rad}}{M}=\left[1-\tilde E_{\rm ISCO}(\tilde a)\right]\nu+4\nu^2\left[4p_0+16p_1\tilde a(\tilde a+1)+\tilde E_{\rm ISCO}(\tilde a)-1\right],\qquad \tilde E_{\rm ISCO}(\tilde a)=\sqrt{1-\frac{2}{3\tilde r^{\rm eq}_{\rm ISCO}(\tilde a)}}\;[\text{Eq. (16)}],$$
$\tilde r^{\rm eq}_{\rm ISCO}$ from Eqs. (3)–(5) (identical to HBR Eqs. (4)–(6)), and [Eq. (12)] $p_0=0.04827\pm0.00039$, $p_1=0.01707\pm0.00032$.
Final mass: $M_{\rm f}/M\equiv1-E_{\rm rad}/M$, with $M$ the mass "at infinite separation" [Sec. 3].
- Check value: for $q=1$, $\boldsymbol a_{1,2}=0$, Eq. (18) reduces to $E_{\rm rad}/M=p_0=0.04827$, i.e. $M_{\rm f}=0.952\,M$ [inference: algebraic identity; equivalently Eq. (11) at $a_1=a_2=0$].
- Accuracy: Eq. (11)–(12) reproduce equal-mass aligned NR data "to within $\sim0.005M$ (except for almost maximal spins)"; fit gives 9.95% at $a_1=a_2=1$ vs 10.95% found by Lovelace+ for $a\approx0.97$ [Sec. 2]. Against all 186 NR runs, residuals "$\lesssim1\%$" and errors "$\lesssim10\%$ relative to the radiated mass" [Sec. 2, Fig. 1], at both small and large ($2\times10^4M$) separations [Sec. 3, Fig. 2]. Coefficients come only from equal-mass aligned runs; in 2012 "spinning binaries with unequal mass ratios are essentially absent" [Sec. 3]. v3 includes an erratum figure (`fig4_erratum.pdf`).
- Precedent for the HBR + Barausse 2012 pairing: Baibhav+ 2020 (2004.00650) compute remnants with "final mass [Barausse 2012], spin [Hofmann 2016]" [Sec. II C 1 "Retention in the cluster" (source order)]; Gerosa & Berti 2017 use Barausse 2012 + Barausse & Rezzolla 2009 [Sec. II B].

### 1c. Alternative (aligned spins only): Jiménez-Forteza+ 2017 (1611.00332)

- Non-precessing only ("quasicircular nonprecessing merging BBHs", abstract); 427 NR runs; uses $\hat S=(m_1^2\chi_1+m_2^2\chi_2)/(m_1^2+m_2^2)$ [Eq. (3)] and $\Delta\chi=\chi_1-\chi_2$.
- $E_{\rm rad}$: Sec. IV, Eqs. (21)–(28): 1D [Eq. (21)] $E_{\rm rad}(\eta,0)=a_4\eta^4+a_3\eta^3+a_2\eta^2+\left(1-\frac{2\sqrt2}{3}\right)\eta$, $a_2=0.5610$, $a_3=-0.847$, $a_4=3.145$;
  [Eq. (22)] $E_{\rm rad}(0.25,\hat S)=\frac{0.0484161\,(0.128b_3\hat S^3+0.211b_2\hat S^2+0.346b_1\hat S+1)}{1-0.212b_5\hat S}$, $b_1=-0.209$, $b_2=-0.197$, $b_3=-0.159$, $b_5=2.985$ [Table VIII];
  2D product ansatz [Eq. (23)], constraints $f_{i2}=16-16f_{i0}-4f_{i1}$ [Eq. (24)], EMRL [Eq. (25)–(26)], spin-difference terms [Eq. (27a–c)], final 3D [Eq. (28)], coefficients $d_{10}=-0.098$, $d_{11}=-3.23$, $d_{20}=0.0112$, $d_{30}=-0.0198$, $d_{31}=-4.92$, $f_{11}=15.7$, $f_{31}=-243.6$, $f_{51}=-0.58$ [Table X]; RMSE $2.62\times10^{-4}$ over 298 cases [Table XI].
- The full assembly of the EMRL coefficients $f_{i0}$ (Table IX) was **not transcribed here**; the $\eta$-fit table number (VII) is **unconfirmed** by the PDF text. Not recommended for precessing 2G parents; Mapelli+ (2103.05016 Sec. 3; 2109.06222 Sec. 2.4) and Torniamenti+ (2401.14837 Sec. 2.2.3) use it for remnants.

---

## 2. Literature on the 2G natal-spin (and mass) distribution

| Source | Setup | Statement / number | Location |
|---|---|---|---|
| Fishbach, Holz & Farr 2017 (1703.06869) | HBR fit, $6.25\times10^6$ BHs, isotropic spins, major mergers $0.7\le q\le1$ | Equal-mass non-spinning → $a=0.6864$; distribution "peaked at $a=0.69$, and with little support below $a\approx0.5$" | Sec. 1 |
| ″ | ″ | "sharply peaked at $a\sim0.7$ with nonzero support over $0.4\lesssim a\lesssim0.9$"; insensitive to 1G spin, to $q\in[0.7,1]$, and to generation | Sec. 2.1, Figs. 1–2 |
| ″ | 1G with isotropic $a=0.99$ | 2G "similarly peaked at $a\sim0.7$ with slightly wider support" | Sec. 2.1, Fig. 1 |
| ″ | always-aligned spins | narrower, "converging to $a\sim0.9$ by the third generation" | Sec. 3 |
| ″ | always anti-aligned | 1G $a\sim0.7$ → $0.43<a<0.54$; maximal anti-aligned 1G → 2G $0.32<a<0.36$ ($q\ge0.7$) | Sec. 3 |
| Gerosa & Berti 2017 (1703.06223) | 1g spins uniform [0,1], isotropic; Barausse 2012 + B&R 2009 fits | "spin magnitudes of 2g BHs should cluster (on average) around … $\sim0.7$"; Fig. 3: remnant spins "clustered around $\sim0.7$, quite independently of the progenitor parameters". No width quoted | Sec. I item (ii); Sec. II E; Fig. 3 |
| Rodriguez+ 2019 (1906.10260), CMC, 96 GC models | $\chi_{\rm birth}=0$ | 2G spins "strongly concentrated at $\chi_{\rm final}\approx0.69$, with 90% of sources lying between 0.67 and 0.69" | Sec. IV D, Fig. 8 |
| ″ | $\chi_{\rm birth}=0.5$ | median 0.62, 90% between 0.53 and 0.72 | Sec. IV D, Fig. 8 |
| ″ | $\chi_{\rm birth}=0$, $q$ dependence | kick max "$\sim175\pm11$ km/s when the mass ratio is $q\sim1/3$", remnant spin "decreases to $\chi_{\rm final}\approx0.55$" | Sec. IV D |
| ″ | retention | ~60% of merger products retained for $\chi_{\rm birth}=0$; <3% for 0.5; 1G+2G and 2G+2G products "virtually never retained" (1 case in ~$10^4$ mergers) | Sec. III A, Fig. 1 |
| ″ | 2G fractions ($\chi_{\rm birth}=0$) | 1G+2G ≈13% of GC mergers (18% of detected); 2G+2G 1% (3% detected) | Sec. III B, Fig. 3 |
| ″ | 2G mass / $q$ | 2G "nearly twice the mass of the most massive 1G BHs"; 1G+2G mass ratios peak at $\sim0.5$ | Sec. IV B, Fig. 6 |
| Rodriguez+ 2020 (2005.04239) | 3G from 1G+2G vs 2G+2G | median 3G spin 0.62 / 0.68; only kicks <100 km/s → 0.43 / 0.65; 1G+2G minimum kick $\sim90$ km/s | Sec. 2, Fig. 2 |
| ″ | CMC merger trees | 2G spins in the trees 0.68–0.69 (no merger times given) | Fig. 4 |
| Borchers, Ye & Fishbach 2025 (2503.21278) | NRSur7dq4Remnant; CMC $q$ distribution and $v_{\rm esc}$ profile (peak ~50 km/s, range ~5–145 km/s, Fig. 2) | all remnants centred at $a=0.69$; **retained** remnants peak above 0.69 once $a_i\gtrsim0.3$; "for initial spins of $a_i=0.6$ the distribution peaks at $a_f\sim0.8$"; $a_i=0.1$ peaks at $\sim0.7$ | Sec. 3.1, Fig. 8; Sec. 5 |
| ″ | fixed $q$ | $q=1$: retained distribution centred at $a_f\sim0.7$ independent of $a_i$; $q=2$: above 0.7 ($a_i=0.3$) and above 0.8 ($a_i=0.8$) | Sec. 3.1, Fig. 10 |
| ″ | low $v_{\rm esc}$ | bimodal; at low escape velocity "two peaks around $a_f\sim0.6$ and $a_f\sim0.8$"; e.g. $v_{\rm esc}=20$ km/s, $a_i=0.3$ | Sec. 3.2, Figs. 11–12; Sec. 5 |
| ″ | higher generations | no convergence; "$a_f\in(0.4,1)$" | Sec. 3.3, Fig. 13; abstract |
| ″ | 25% near-aligned 1G+1G, $v_{\rm esc}=50$ km/s | peaks at 0.69 and 0.72 ($a_i=0.1$); second peak 0.77 ($a_i=0.3$) | Sec. 4, Fig. 14 |
| ″ | kick-model systematics | PRECESSION (used in CMC) kicks ~10× larger on average than NRSur7dq4Remnant | Sec. 2.3, Fig. 7 |
| Ye+ 2025 (2507.07183) | CMC + SFR-like formation history | ~10% of cluster mergers contain higher-generation BHs; their spins span "about 0.4 to 1", "major peak above 0.7" (citing Borchers+) | Sec. 3.2 (source order) |
| Mapelli+ 2021 (2103.05016) | JF2017 remnant fit | "final mass is $\approx0.9$ the total mass", "final spin magnitude clusters around $\chi_{\rm f}\approx0.7$–0.9" | Sec. 3 |
| Baibhav+ 2020 (2004.00650) | HBR + Barausse 2012 | "spin distribution peaked at $\chi\sim0.7$"; 74% (50%) of 2g mergers in the spin gap for $\chi_{\max}=0.1$ (0.2) | Sec. III B (source order) |
| Gerosa & Fishbach 2021 (2105.03439) | review | equal-mass non-spinning $\chi\sim0.69$; generic mergers "strongly peaked at the characteristic value $\chi\simeq0.7$" | Sec. 2 (source order) |
| Tagawa+ 2021 (2012.00011), AGN | gas accretion on | mergers "tend to have high spin around $\sim0.7$"; accretion "significantly enhances the BH spin" and partly aligns | Sec. 3.1, Fig. 1b |
| Atallah+ 2023 (2211.09670), NSC seeds | — | "First and second generation mergers fall between $\chi_{\rm eff}=0.5-0.75$" (wording as printed) | Sec. 5 (source order) |
| LVK 2510.26931 | — | remnant spins "generically and robustly concentrated about $\chi\approx0.7$" | Sec. 4.2 (source order) |

**Bottom line for 2G natal spin** [inference, combining the rows above]: for low-spin 1G parents ($\chi\lesssim0.1$) and near-equal masses, the
2G natal spin is $0.69$ with a 90% width of only $\sim0.02$ (CMC: 0.67–0.69). The width grows with 1G spin (CMC, $\chi_{\rm birth}=0.5$:
0.53–0.72). Unequal parents ($q\sim1/3$) give $\approx0.55$. The retention (kick) selection shifts the retained 2G distribution upward
(peak $\sim0.8$ for $a_i=0.6$) and can make it bimodal ($\sim0.6$/$\sim0.8$) at low $v_{\rm esc}$. In AGN disks gas accretion modifies
spins after birth.

**2G natal mass**: $M_{\rm 2G}=(1-E_{\rm rad}/M)(m_1+m_2)$ with Barausse Eq. (18); $\approx0.952(m_1+m_2)$ for equal-mass non-spinning parents (§1b).

---

## 3. Ages (delay times) by generation

### 3a. Is there a published distribution of $\Delta t_{\rm 2G}$?

**No paper read here gives the distribution of $\Delta t_{\rm 2G}$ (time from the parent 1G+1G merger to the 2G's own merger) as an explicit
histogram, median, or table.** Checked: Rodriguez+ 2019, Rodriguez+ 2020, Kremer+ 2020, Antonini+ 2019/2023, Mapelli+ 2021/2022,
Torniamenti+ 2024, Kritos+ (Rapster), Fragione & Silk 2020, Atallah+ 2023, Chattopadhyay+ 2023, Arca Sedda+ 2023 (Dragon-II II/III),
Baibhav+ 2020, Tagawa+ 2021, Vaccaro+ 2024, McFACTS I, Gerosa & Fishbach 2021, Ye+ 2025, Borchers+ 2025. The closest items:

1. **Baibhav+ 2020 (2004.00650)** — the only paper that states how the 1g→2g delay is computed and what it does to the redshift distribution:
   - "We calculate the time delay between a 1g+1g merger and the next merger involving its 2g remnant using the same $t_{\rm tot}$ … but replacing the mass segregation timescale $t_{\rm MS}$ by $t_{\rm MS}\times(v_{\rm kick}/v_{\rm esc})^3$" [Sec. II C 2 "2g+1g or 2g+2g?" (source order)], with [Eq. (26)]
     $$t_{\rm tot}=t_{\rm MS}+\min(t_{1\text-2},t_{\rm 3bb})+t_{\rm hard}+t_{\rm GW}\;\Rightarrow\;\Delta t_{\rm 2G}=t_{\rm MS}\Big(\frac{v_{\rm kick}}{v_{\rm esc}}\Big)^3+\min(t_{1\text-2},t_{\rm 3bb})+t_{\rm hard}+t_{\rm GW}.$$
   - Result: "2g BBHs closely follow the same distribution [as 1g+1g, which follows the Madau SFR], contrary to the expectation that repeated mergers may suffer further time delays", because merger products "merge on very short timescales in our model" [Sec. II C 2, Fig. 7].
   - Component timescales [Eqs. (10)–(25)]: $t_{\rm MS}=(m_\star/m_{\rm BH})\,t_r(r_h)$ [Eq. (10)], $t_r(r_h)=4.2\times10^9\,(15/\ln\Lambda)(r_h/4\,{\rm pc})^{3/2}(M_{\rm cl}/10^7M_\odot)^{1/2}$ yr [Eq. (11)];
     $t_{\rm 3bb}=6.45\times10^9\,(n_{\rm BH}/10^6{\rm pc^{-3}})^{-2}(\sigma_{\rm BH}/10\,{\rm km\,s^{-1}})^9(m_1/10M_\odot)^{-5}$ yr [Eq. (16)];
     $t_{\rm hard}=\frac{3.16\times10^8}{q_3}(\sigma_{\rm BH}/30\,{\rm km\,s^{-1}})(a_{\rm crit}/0.05\,{\rm AU})^{-1}((m_1+m_2)/20M_\odot)^{-1}(n_{\rm BH}/10^6{\rm pc^{-3}})^{-1}$ yr [Eq. (24)];
     $t_{\rm GW}=0.56\times10^8\,[m_1m_2(m_1+m_2)/2\times10^3M_\odot^3]^{-1}(a_{\rm crit}/0.1\,{\rm AU})^4((1-e^2)/0.5)^{7/2}$ yr [Eq. (25)], $a_{\rm crit}=\max(a_{\rm GW},a_{\rm ej})$ with $a_{\rm GW}$ Eq. (23), $a_{\rm ej}$ Eq. (20).
   - [fig] Fig. 2: total delay for $(50+50)M_\odot$ is $\approx1$–$2\times10^8$ yr at $M_{\rm cl}=10^5M_\odot$, minimum $\approx5\times10^7$ yr near $M_{\rm cl}\approx3\times10^6M_\odot$ (hardening-dominated), rising to $\sim10^{10}$ yr at $10^9M_\odot$; for $(10+10)M_\odot$ the minimum is $\approx3\times10^8$ yr near $10^6M_\odot$. Because $t_{\rm MS}$ is suppressed by $(v_{\rm kick}/v_{\rm esc})^3$ for 2G, these totals are upper-side proxies for $\Delta t_{\rm 2G}$ [inference].
2. **Mapelli+ 2021 (2103.05016), FASTCLUSTER** — Nth-generation assembly [Eq. (21)]:
   $$t_{\rm dyn,ng}=t_{\rm merg}+t_{\rm DF}+\min(t_{\rm 3bb},t_{12}),$$
   so $\Delta t_{\rm 2G}=t_{\rm DF}+\min(t_{\rm 3bb},t_{12})+t_{\rm del}$ [inference from the definitions in Sec. 4.2]. Ingredients: $t_{\rm DF}=\frac{3}{4(2\pi)^{1/2}G^2\ln\Lambda}\frac{\sigma^3}{m_{\rm BH}\rho}$ [Eq. (3)];
   $t_{\rm 3bb}=125\,{\rm Myr}\,(10^6{\rm pc^{-3}}/n_{\rm c})^2(\zeta^{-1}\sigma_{\rm 1D}/30\,{\rm km\,s^{-1}})^9(20M_\odot/m_{\rm BH})^5$ [Eq. (4)];
   $t_{12}=3\,{\rm Gyr}\,(0.01/f_{\rm bin})(10^6{\rm pc^{-3}}/n_{\rm c})(\sigma/50\,{\rm km\,s^{-1}})(12M_\odot/(m_{\rm BH}+2m_*))(1\,{\rm AU}/a_{\rm hard})$ [Eq. (6)];
   hardening + GW integration [Eq. (16)]. Cluster priors [Sec. 5]: $t_{\rm SC}=13.6$, 12, 1 Gyr; $\langle\log_{10}M_{\rm tot}\rangle=6.18$, 5.3, 4.3; $\langle\log_{10}\rho\rangle=5$, 3.3, 3.3 (NSC, GC, YSC); $\rho_c=20\rho$.
   - Text [Sec. 6.2]: $t_{\rm dyn}$ (1g assembly) "$\sim3$–20 Myr" in YSCs, "$\sim20$–200 Myr" in GCs, NSCs "peaking at $t_{\rm dyn}\sim20$–200 Myr, but with a considerable fraction of systems with $t_{\rm dyn}>1$ Gyr". Sec. 6.3: 1g dynamical assembly "could be as short as $\sim3$ Myr or as long as $\sim t_{\rm Hubble}$".
   - [fig] Figs. 1 (NSC) and 2 (GC), lower-left panels: the plotted $t_{\rm dyn,ng}$ spans $\approx0.3$ Myr (NSC) / $\approx0.6$ Myr (GC) up to $\approx7$ Gyr, peaking at $\sim10$–100 Myr; the plotted $t_{\rm merg,ng}$ spans $\approx4$ Myr–14 Gyr. Since the plotted $t_{\rm dyn,ng}$ extends **below** the minimum $t_{\rm merg}$ (~4 Myr), it appears to be the increment $t_{\rm DF}+\min(t_{\rm 3bb},t_{12})$ (time from the previous merger to the new pairing), not the cumulative Eq. (21) [inference; the caption calls it the "time for the dynamical formation of nth-generation BBHs" — **interpretation unconfirmed**].
   - Mapelli+ 2022 (2109.06222) uses the same chain [Eq. (8)]: $t_{\rm Ng}=t_{\rm merg}+t_{\rm DF}+\min(t_{\rm 3bb},t_{12})$ [Sec. 2.4]; no $\Delta t$ distribution shown.
3. **Direct N-body example — Arca Sedda+ 2023, Dragon-II II (2307.04806)**: model IBH_Rh0.47f20N300k (re-run with an imposed $v_{\rm GW}=20$ km/s kick): the 2g IMBH "undergoes a series of resonant interactions with two BHs, which drives the IMBH to merge after just 25.5 Myr" [Sec. 4.1.3]. Earlier in the same model the 1g IMBH merged "After only 17 Myr" [Sec. 4.1.3]. This is an IMBH (~270 $M_\odot$), not a stellar-mass 2G. Whether 25.5 Myr is counted from the 2g formation is the natural reading but is not stated explicitly [**partly unconfirmed**].
4. **AGN disks** (all generations within the disk lifetime): see §3c.
5. **Modelling assumptions** (not predictions): Gerosa & Berti 2017 draw $t_D$ "from a flat distribution in $\log(t_D)$ in the range $t_D\in[10^{-4}\,{\rm Gyr},t_L(\tilde z)]$", the lower limit being "very conservative" [Sec. II D]. Payne+ 2021 (2107.11730) state that their 2G population is present with or without the boson because "the majority of black holes in this population are not old enough" [Sec. 4.3 (source order), Fig. 6].

### 3b. Proxy timescales (formulas with the papers' normalizations)

- **Antonini, Gieles & Gualandris 2019 (1811.03640)** (balanced Hénon evolution):
  heating rate [Eq. (3)] $\dot E_0\simeq2.3\times10^5\,M_\odot({\rm km\,s^{-1}})^2{\rm Myr^{-1}}\,M_5^{2/3}\rho_{5,0}^{5/6}$ ($\zeta=0.1$, $\psi=5$, $\bar m=0.6M_\odot$), $\dot E_{\rm bin}=-\dot E$ [Eq. (9)];
  binary lifetime [Eq. (11)] $\tau_{\rm bin}\simeq-\frac{Gm_1m_2}{2a_{\rm m}}\dot E_{\rm bin}^{-1}$, $a_{\rm m}=\max(a_{\rm GW},a_{\rm ej})$ [Eq. (12)], $a_{\rm GW}$ [Eq. (16)];
  $t_{\rm rh,0}\simeq7.5\,{\rm Myr}\,M_5\rho_{5,0}^{-1/2}$, $v_{\rm esc,0}\simeq50\,{\rm km\,s^{-1}}M_5^{1/3}\rho_{5,0}^{1/6}$ [Eqs. (6)–(7)];
  dynamical friction [Eq. (22)] $\tau_{\rm df}\simeq1.65\,r_h^2\frac{\sigma}{\ln\Lambda\,Gm}=4.5\,{\rm Myr}\,(m/20M_\odot)^{-1}M_5\rho_{5,0}^{-1/2}$; a kicked remnant sinks back on $\tau_{\rm df}(r_{\rm in})$ with $r_{\rm in}\simeq r_h\sqrt{v_{\rm esc}^4/(v_{\rm esc}^2-v_{\rm GW}^2)^2-1}$ and "remains 'inactive' until this time has passed"; the time from binary formation to merger is $\tau_{\rm merge}=\tau_{\rm bin}+\tau_{\rm GW}$ [Sec. 3].
  Example numbers [Sec. 3.1 (source order), Fig. 2]: $M_{\rm cl}=10^7M_\odot$; first merger at "100 Myr" ($\rho_{h,0}=10^6M_\odot{\rm pc^{-3}}$) vs "1 Gyr" ($10^4$); for $10^6$: 1873 (1719) mergers in the first 1 Gyr, 617 (1400) remnants retained, for spins $U[0,1)$ ($\chi=0$); for $10^4$: 959 (915) mergers over ~10 Gyr, 51 (459) retained.
  [inference] With ~1800 mergers per Gyr in one cluster, consecutive in-cluster mergers are $\sim0.5$ Myr apart on average; with pairing probability $\propto(m_1+m_2)^4$, a retained (heaviest) remnant should be re-paired soon after $\tau_{\rm df}(r_{\rm in})$, so $\Delta t_{\rm 2G}\sim\tau_{\rm df}+\tau_{\rm bin}+\tau_{\rm GW}$ (Myr–100 Myr scale). Not stated by the paper.
- **Antonini+ 2023 (2208.01081), cBHBd**: first binary at $\tau_{\rm cc}=3.21\,t_{\rm rh,0}$ [Eq. (1)]; retained remnants get $\tau_{\rm df}\simeq1.65\,r_{\rm in}^2\sigma/(\ln\Lambda\,Gm)$ [Eq. (6)] with $r_{\rm in}$ [Eq. (7)] "and only allow the BH to form a new binary after this time"; binary lifetime $\tau_{\rm bin}=m_{\rm ej}/\dot M_{\rm BH}$ [Eq. (8)]; integration to 13 Gyr [Sec. 2]. No $\Delta t$ numbers given.
- **Kritos+ (Rapster, 2210.10055)**: a recoiling binary/remnant that stays bound but leaves the BH subsystem is delayed by the convection time $(\bar m/m_{12})\tau_{\rm rh}$ [Sec. II E]. No $\Delta t$ numbers.
- **Fragione & Silk 2020 (2006.01867), NSC seed growth**: $t_{\rm 3bb}=125\,{\rm Myr}(10^6{\rm pc^{-3}}/n)^2(\xi^{-1}\sigma/30\,{\rm km\,s^{-1}})^9(20M_\odot/M_{\rm MBH})^5$ [Eq. (2)]; $t_1=300\,{\rm Myr}\,\xi^{-1}(10^6{\rm pc^{-3}}/n)(\sigma/30)(30M_\odot/(10M_\odot+M_{\rm MBH}+m_2))(1\,{\rm AU}/a_h)$ [Eq. (4)]; $t_2=200\,{\rm Myr}\,\xi^{-1}(\ldots)(0.05\,{\rm AU}/\max(a_{\rm ej},a_{\rm GW}))(10M_\odot/(M_{\rm MBH}+m_2))^2$ [Eq. (9)]; $T_{\rm GW}=2000\,{\rm Myr}(10^3M_\odot^3/M_{\rm MBH}m_2(M_{\rm MBH}+m_2))(\max(a_{\rm ej},a_{\rm GW})/0.05\,{\rm AU})^4(1-e^2)^{3/2}$ [Eq. (10)]; $\tau_{\rm df}\approx1\,{\rm Myr}(100M_\odot/M_{\rm MBH})(M_{\rm NSC}/10^5M_\odot)(10^5M_\odot{\rm pc^{-3}}/\rho_{\rm NSC})^{1/2}$ [Eq. (18)]; all sampled from $\exp(-t/\tau)$ [Sec. 2].
  [fig] Fig. 2 ($M_{\rm NSC}=10^8M_\odot$, $\rho=10^6M_\odot{\rm pc^{-3}}$): the 50 $M_\odot$ seed's first merger is at $\approx0.4$ Gyr and the next at $\approx1.5$ Gyr; successive mergers are $\sim0.1$–1 Gyr apart, becoming more frequent as the seed grows.
- **Atallah+ 2023 (2211.09670), NSC seeds (50–300 $M_\odot$)**: [fig] Fig. 3 (median merger rate per seed over the first ~100 interactions, $t_{\rm df}$ excluded): $\Gamma_{\rm merge}\approx5$, 1.5, 0.08 Myr$^{-1}$ for $M_{\rm cl}=10^6M_\odot$ and $\epsilon_b=0.005$, 0.01, 0.05; $\approx1.5$, 0.4, 0.02 Myr$^{-1}$ at $10^7$; $\approx0.3$, 0.1, 0.007 at $10^8$; $\approx0.02$, 0.007, 0.0013 at $10^9$. Text [Sec. 5 (source order)]: seeds accumulate mergers "On timescales of ${\sim}0.01-1$ Gyr following their formation"; low-mass NSCs ($10^6$–$10^7M_\odot$) reach escape or runaway in "$<100$ Myr"; escapes mostly "by 1 Gyr … before 0.5 Gyr" [Sec. 3.5].
- **Chattopadhyay+ 2023 (2308.10884), NSC (cBHBd)**: IMBHs grow by hierarchical mergers "in timescales of a few $100$ Myrs to a few Gyrs" (abstract). No per-generation $\Delta t$.

### 3c. Summary by environment

| Environment | 2G age at merger $\Delta t_{\rm 2G}$ | 1G age at merger $\tau_{\rm 1G}$ | Shortest plausible 2G age |
|---|---|---|---|
| **GCs** | No explicit distribution. Qualitative: 2G "most likely to rapidly merge again before the cluster can form another 2G BH" (Rodriguez+ 2019 Sec. IV B); of 31 retained 2G in one model, "16 are later ejected … as binaries while 11 merge again inside the cluster" (Sec. IV A). Baibhav+ 2020: 2g follow 1g in redshift, "very short timescales" (Sec. II C 2, Fig. 7); total delay for $(50+50)$ min $\approx5\times10^7$ yr [fig, Fig. 2]. FASTCLUSTER GC pairing increment peaks $\sim30$–100 Myr, range $\approx0.6$ Myr–7 Gyr [fig, Fig. 2; inference]. | 1G formation follows cluster formation: z peak $\sim4$ (Antonini+ 2023 Sec. 2.1 (source order)), z\_peak = 2.2 for SFR-like histories (Ye+ 2025 Eq. (1)). Low-z mergers from old GCs: "$z<0.5$ ($t_{\rm mergers}\gtrsim6$ Gyr, assuming cluster ages of 12 Gyr" (Kremer+ 2020 Fig. 15). Ejected (~50%) "median of $\sim10$ Gyr", in-cluster mergers "extremely promptly" (Fishbach & Kalogera 2021 Sec. 1). Rapster: "the majority of mergers occur early"; late-time tail $\propto1/t_{\rm me}$ for $r_v>0.5$ pc (Sec. III B 3, Fig. 9). FASTCLUSTER GC 1g $t_{\rm merg}$: $\approx4$ Myr–14 Gyr, broad, highest at $\sim1$–5 Gyr [fig, Fig. 2]. | ~1 Myr–few×10 Myr [inference]: Antonini 2019 $\tau_{\rm df}=4.5$ Myr$(m/20)^{-1}M_5\rho_5^{-1/2}$ (Eq. 22) + hardening; FASTCLUSTER increment down to $\approx0.6$ Myr [fig]. No paper quotes a minimum. |
| **NSCs** | No explicit distribution. FASTCLUSTER NSC increment peaks $\sim10$–100 Myr, range $\approx0.3$ Myr–7 Gyr [fig, Fig. 1; inference]. Fragione & Silk: successive mergers $\sim0.1$–1 Gyr apart in $10^8M_\odot$ NSC [fig, Fig. 2]. Atallah: per-seed merger intervals $\sim0.2$ Myr ($10^6M_\odot$, compact) to $\sim1$ Gyr ($10^9M_\odot$) [fig, Fig. 3]; "${\sim}0.01-1$ Gyr" (Sec. 5). Baibhav Fig. 2: delays grow for $M_{\rm cl}\gtrsim10^7M_\odot$ ($t_{\rm 3bb}\propto M_{\rm cl}^{5/2}$, Sec. II C 1). | NSC 1g $t_{\rm dyn}$ peaks 20–200 Myr with a ">1 Gyr" fraction (Mapelli+ 2021 Sec. 6.2); NSC $t_{\rm SC}=13.6$ Gyr assumed (Sec. 5). FASTCLUSTER NSC 1g $t_{\rm merg}$ roughly flat from ~100 Myr to 14 Gyr [fig, Fig. 1]. Absolute stellar ages of NSC 1G BHs at LVK redshifts: **unconfirmed** (not given). | ~0.3–1 Myr [fig/inference] in compact, low-mass NSCs; Myr–100 Myr typical. |
| **YSCs** (extra) | Dragon-II: 2g IMBH re-merges "after just 25.5 Myr" (2307.04806 Sec. 4.1.3). Hierarchical mergers rare in YSCs ("YSCs have the lowest efficiency of hierarchical mergers", Mapelli+ 2021 Sec. 6.7 (source order); "none of the mergers with component masses $<100\,M_\odot$ undergo multiple mergers", 2307.04807 Sec. 2.2.1). | YSC delay times peak "$\sim100$ Myr", $\tau^{-1}$ above 400 Myr (Fishbach & Kalogera Sec. 1, citing Di Carlo+ 2020); Dragon-II mergers peak at $t_{\rm GW}\simeq0.5$, 1.5, 10 Gyr, i.e. $t_{\rm GW}/t_{\rm rlx}\simeq8$–30 (2307.04807 Sec. 3.1.2, Fig. 5). YSC 1g $t_{\rm dyn}\sim3$–20 Myr (Mapelli+ 2021 Sec. 6.2). | ~10 Myr (N-body example) |
| **AGN disks** | All mergers occur within the disk lifetime: AGN lifetime log-normal with "$\langle\log\tau/{\rm Myr}\rangle=0.22\pm0.80$" and Ng requires $t^{(N)}_{\rm merg}<\tau$ (Vaccaro+ 2024 Secs. 2.1, 2.5); Ng timescales (gas hardening): damping peak $10^{-5}$ Myr, migration peak $10^{-2}$ Myr, $t_{\rm del}$ peak "$\sim1$ yr" (Sec. 3.4, Fig. 12). McFACTS: $\tau_{\rm AGN}$ "plausibly $\sim0.1-10$ Myr" (Sec. 2.1); "high generation mergers (2g,3g) quickly appear in $<0.1$ Myr" (Sec. 4.2, Fig. 4); with dynamics 3g-ng "first appear after $\sim0.4$ Myr vs $\sim0.2$ Myr" (Sec. 4.3, Fig. 6). Tagawa+ 2021: max merging mass 130, 185, 185 $M_\odot$ and 260, 1600, 7300 mergers at 1, 10, 100 Myr, from 1g BHs $\le15M_\odot$ (App. B; Sec. 2.3) → several generations within ~1 Myr [inference]. Fishbach & Kalogera Sec. 1: AGN delays "$\tau<100$ Myr". | 1g BHs are pre-existing NSC BHs captured by the disk (Tagawa Sec. 2.1; Vaccaro Sec. 2.3; McFACTS Sec. 2.2). Their stellar age is **unconfirmed** (none of these papers states it). | $\lesssim10^{-2}$ Myr (Vaccaro Ng migration peak) to $<0.1$ Myr (McFACTS) |
| **Isolated binaries** (1G only) | — | "$p(\tau)\propto\tau^\alpha$", typical "$\alpha=-1$" (Fishbach & Kalogera Sec. 1); $\tau_{\rm min}=10$ Myr "as predicted for the isolated binary evolution channel" (Sec. 4). GWTC-2 inference: "43–100% of mergers have delay times $<4.5$ Gyr", $\tau_{50\%}=2.8^{+3.3}_{-2.6}$ Gyr, $\tau_{\min}<2.2$ Gyr (90%) for SFR-following formation (abstract; Sec. 3.3, Sec. 5); posterior on $\tau_{\min}$ peaks at 10 Myr, the prior edge (Sec. 3.3). Payne+ 2021: DTD "has little support for systems merging within 10 Myr of formation" (Fig. 1 caption, source order). Mapelli+ 2022: isolated delay times from MOBSE, no numbers quoted (Sec. 2.1). | — |

### 3d. Takeaway for the 1G–2G age contrast [inference]

At an LVK merger at $z\lesssim0.3$: a 1G BH from an old GC is typically Gyr–10 Gyr old (Kremer Fig. 15; ejected median ~10 Gyr), and an
isolated-binary 1G BH is $\sim10$ Myr–Gyr old ($\tau^{-1}$). A cluster 2G BH is typically $\sim10$–100 Myr old (GC/NSC proxies above),
with a tail to Gyr in massive NSCs and a floor of order Myr. An AGN 2G BH is $\lesssim0.1$–10 Myr old. So the 1G/2G age ratio is
$\gtrsim10$–$10^3$ in clusters and $\gtrsim10^3$ in AGN. **No paper gives these as distributions**; all are proxies.

---

## 4. LVK GW241011/GW241110 (2510.26931): environment and timing of GW241011's primary

(Section, figure and table numbers in source order; macro values read from `pe_macros_241011.tex`, `progenitor_macros_241011.tex`,
`eccentricity_macros.tex`, `superrad_macros.tex`.)
- Source [Sec. 3.1, Table 1]: $m_1=19.6^{+3.6}_{-2.5}M_\odot$, $m_2=5.9^{+0.8}_{-0.8}M_\odot$, $q=0.30^{+0.09}_{-0.08}$, $\chi_1=0.78^{+0.09}_{-0.09}$ (medians, central 90%), "$\chi_1>0.69$ at 95% credibility", $\theta_1=31^{+11}_{-14}$ deg (HPD), $\chi_{\rm eff}=0.49$ (median), $D_L=214^{+44}_{-48}$ Mpc, $D_L<248$ Mpc (90%); "probably the closest binary black hole merger observed to date".
- Environment [Sec. 4.2]: "A more natural interpretation … is that they involve the mergers of second-generation black holes in dense stellar environments, such as globular, nuclear, and young massive star clusters". Spin magnitude and mass ratio of GW241011 "are inconsistent with the properties predicted of first-generation mergers, but lie within ranges predicted by both models [CMC, cBHBd] for higher-generation hierarchical mergers" [Fig. 6]. Both models "predict approximately one higher-generation merger per five first-generation mergers". "clusters with stellar metallicities $Z\approx0.1\,Z_\odot$ and above are predicted to readily yield hierarchical mergers with $\sim20\,M_\odot$ primaries" (citing Ye+ 2025). Conclusions: "strongly suggestive of hierarchical binary black hole mergers in dense stellar environments, such as globular, nuclear, or young stellar clusters".
- Model data [App. "Cluster models"]: CMC catalog (148 models) and cBHBd models with $\rho_{h,0}=10^5M_\odot{\rm pc^{-3}}$ (129 models up to $2\times10^6M_\odot$), each "evolved over one Hubble time".
- Isolated binaries [Sec. 4.1]: primary spins "difficult to explain via isolated binary evolution" (unless large/misaligned natal spins).
- Ancestors (1G parents of GW241011's primary) [Sec. 4.3, Fig. 7; App. progenitor]: Backward approach $m=13.3^{+4.8}_{-3.2}$ and $7.5^{+3.2}_{-3.9}M_\odot$, $\chi_{\rm eff}=0.23^{+0.29}_{-0.28}$, recoils "approximately 100–2000 km/s" (median $750^{+1400}_{-630}$ km/s); Forward (astrophysical prior, near-zero 1G spins) $11.0^{+1.9}_{-1.6}$ and $9.6^{+1.7}_{-1.8}M_\odot$, $\chi_{\rm eff}=0.02^{+0.08}_{-0.06}$, recoils "approximately 10–300 km/s" (median $100^{+190}_{-70}$ km/s). Kick constraints "almost entirely prior dominated". Remnant mapping via NRSur7dq4Remnant (Varma+ 2019).
- Eccentricity [Sec. 4.4, Fig. 8]: $e<0.05$ (90%) at 13.33 Hz for GW241011 under both SEOBNRv5EHM and TEOBResumS-Dalí; "not inconsistent with a dynamical origin".
- AGN: mentioned only in the eccentricity discussion (AGN eccentric fractions $\gtrsim10\%$ to $\sim70\%$), **not** as a proposed environment for GW241011 [Sec. 4.4].
- **Timing / age of the 2G primary: not estimated anywhere in the paper.** Superradiance section [Sec. 5.3, Fig. 11]: exclusions are shown "as a function of the presumed age"; "Conservatively assuming an age of $10^5$ years for the primary black hole", scalars excluded in $[0.3,2.6]\times10^{-12}$ eV and vectors in $[0.1,5.3]\times10^{-12}$ eV; overall "between approximately $10^{-13}$ and $3\times10^{-12}$ eV".

---

## 5. Not confirmed in the full texts

- Any explicit histogram / median / minimum of $\Delta t_{\rm 2G}$ for GCs, NSCs or AGN (§3a). All GC/NSC values above are proxies or figure readings.
- Interpretation of the plotted $t_{\rm dyn,ng}$ in Mapelli+ 2021 Figs. 1–2 as an increment (vs the cumulative Eq. 21).
- Whether Dragon-II's "25.5 Myr" is measured from the 2g IMBH's formation.
- Stellar ages of the 1G BHs that enter AGN disks; absolute NSC 1G ages at LVK redshifts.
- HBR 2016 calibration-set coverage in $q$ and spin (only the EMRL-enforced formal domain is stated).
- UIB 2017: η-fit table number (VII) and the full $f_{i0}$ assembly (not transcribed).
- Isolated-binary DTD from a population-synthesis paper read in full (only Fishbach & Kalogera's summary of Dominik+ 2012 / O'Shaughnessy+ 2010 and their own GWTC-2 inference were read).
