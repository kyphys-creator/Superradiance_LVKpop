# Population-level superradiance with LVK BBH spins: full-text check (2026-11-03)

All facts below come from the arXiv **full texts** (LaTeX source plus compiled PDF), downloaded on 2026-11-03 to
`Refs/arxiv/<ID>/`. Section, equation and figure numbers follow the arXiv PDF. If a number could only be read
from the source order, it is marked "(source order)". "**unconfirmed**" means the full text does not state it.
Quotes are short verbatim fragments in "…". Versions: 2609.02678v1, 2607.01317v1, 2610.11091v1.

## 要点 (Japanese summary)

- **Kou+ (2609.02678)**: GWTC-5 全体で $m_b = 7.1^{+1.8}_{-1.1}\times10^{-13}$ eV(90% CI、Table I)、帰無模型(Gaussian component spins)に対し $\ln\mathcal B\approx7.8$。模型は「$\chi>\chi_{\rm SR}$ の質量を幅 $\sigma_{\rm SR}$ の切断正規分布に寄せる」もの(Eq. 1)。$\chi_{\rm SR}$ は SuperRad で計算する(相対論的なのは $m\le2$、半解析は $m\le15$)。$T_{\rm age}$ は事象ごとに log-uniform [1 Myr, 2.2 Gyr] で周辺化する。選択効果は LVK O1–O4 の注入で評価する。**自己相互作用なし**(minimal coupling)、**階層合体の模型化なし**。
- **Ning+ (2607.01317)**: 単一の、質量に依らない Beta 分布の自然スピン(確認済み、Eq. 2)。速度は解析的な水素原子近似(+$\alpha^2$ 補正)で $m=1$–6 を使い、Witte–Mummery のコードは系統誤差の検査だけに使う。**選択関数なし**、質量分布は固定(Salpeter、一様な $q$)。自己相互作用は QCD の $f_a$ に限って Witte のコードで検査し、「最終スピンは場合により最大 ~0.05 変わる、無視できる」(SM「Superradiance for higher α」)。ほかの $f_a$ は future work と書かれている。
- **新規性**: 2025-01 から今日までに、**LVK のスピン集団で $(\mu,f)$ を走査した解析は見つからない**。レビュー 2610.11091 の Sec. V は、集団解析は自由場のままであると明記している。近い研究は (i) 2609.12052(GWTC-5 で $(\mu,1/f)$ 面を扱うが、スピンではなく $q<0.3$ の合体率を使う)、(ii) 個別事象の $(\mu,f)$(2507.20979、2507.21788、2503.10347)。2G のスピン ≈0.69–0.7 を超放射の較正に使った論文は見つからない。2021 年の Payne+ (2107.11730) は「2G は若すぎて超放射の影響を受けない」としている(Sec. 4.3)。

---

## 1. Kou, Mandic, Ding, Tian — arXiv:2609.02678 (v1, 2026-09-02)

"Ultralight Bosons Explain the Mass–Spin Correlations in the Merging Binary Black Hole Population".
The arXiv metadata check is OK. Local files: `Refs/arxiv/2609.02678/{main_v2.tex, 2609.02678.pdf, .txt}`.
Layout: an untitled PRL-style main text; End Matter Secs. I–III; Supplemental Material Secs. S1–S4. Main-text
equations are (1)–(3) and SM equations (S1)–(S8). Figures 1–6, S1–S8. Table I.

### 1a. Boson-mass result (exact)
- GWTC-5.0, full catalog: $m_b = 0.71^{+0.18}_{-0.11}\times10^{-12}$ eV, $\sigma_{\rm SR}=0.18^{+0.05}_{-0.03}$,
  $\ln\mathcal B = 7.8$ [Table I]. The **90%** credibility comes from the Table I caption: "Inferred 90% credible intervals of boson mass $m_b$ …".
  The Introduction (last paragraph) gives the same number: "$m_b=7.1_{-1.1}^{+1.8}\times10^{-13}$ eV yields a Bayes factor of $\ln\mathcal B\approx7.8$".
  The abstract states only "$m_b\sim10^{-12}$ eV" and "$\ln\mathcal B\approx7.8$ for GWTC-5.0".
- The other rows of Table I ($m_b$ in $10^{-12}$ eV; $\sigma_{\rm SR}$; $\ln\mathcal B$):
  - GWTC-3.0: $0.74^{+1.75}_{-0.43}$; $0.18^{+0.12}_{-0.07}$; 2.8
  - GWTC-4.0: $0.74^{+0.28}_{-0.22}$; $0.16^{+0.04}_{-0.03}$; **9.1**
  - GWTC-5.0 without GW241011: $0.69^{+0.22}_{-0.11}$; 8.1
  - GWTC-5.0 without GW250114: $0.72^{+0.16}_{-0.13}$; 8.4
  - GWTC-5.0 without both: $0.71^{+0.15}_{-0.11}$; $0.16^{+0.04}_{-0.02}$; **9.5**
- Fig. 1 shows the $(m_b,\sigma_{\rm SR})$ posteriors for GWTC-3/4/5 (68/95% contours). Fig. 3 (End Matter I) plots the
  illustrative trajectory at $m_b=7.1\times10^{-13}$ eV and $T_{\rm age}=10^7$ yr.
- A commented-out draft line in the source (not in the PDF) reads $\log_{10}(m_b/{\rm eV})=-12.10\pm0.14$, $\ln\mathcal B\approx8.7$. **Do not cite it.**

### 1b. Spin model ("environment-induced spin variation")
- Natal (baseline) spin model: the LVK *Gaussian component spin model*, a truncated normal $\mathcal N_{[0,1]}(\mu_\chi,\sigma_\chi)$
  "shared by both BHs" [Baseline model paragraph]. Mass model: *Broken Power Law + 2 Peaks*. Redshift model: $\mathcal R(z)\propto(1+z)^\kappa$.
  All sectors are inferred **jointly** [Baseline model; S4 calls it "the first fully joint hierarchical inference"].
- Superradiance spin model [Eq. (1)]:
  $$p(\chi\,|\,M_f,T_{\rm age},m_b,\sigma_{\rm SR},\Lambda_s)=\pi_{\rm base}(\chi|\Lambda_s)\,\Theta(\chi_{\rm SR}-\chi)+C\,\mathcal N_{[0,1]}(\chi_{\rm SR},\sigma_{\rm SR})$$
  $C$ "normalizes the total distribution to unity". By normalization this implies $C=\int_{\chi_{\rm SR}}^1\pi_{\rm base}\,d\chi$
  (the paper does not write this form). In words, all baseline probability above $\chi_{\rm SR}$ is moved into a
  truncated Gaussian centred on $\chi_{\rm SR}$ with free width $\sigma_{\rm SR}$. The width stands for accretion and environmental
  spin-up or scatter [paragraphs before Eq. (1); footnote 1: mass accretion is "expected" to be absorbed by $\sigma_{\rm SR}$].
- "For a given final mass $M_f$, any BH with a spin exceeding the maximum value $\chi_{\rm SR}(M_f,T_{\rm age},m_b)$ is reassigned to this maximum spin."
  $M_f$ is the post-superradiance (observed) mass. This matches SuperRad `max_spin`, which bins by final BH mass (read from the code, see §4).
- The model is applied to each component spin separately (S2: "treats the primary and secondary spins separately").
- Tilts: the Baseline-model paragraph says "we inherit the isotropic tilt-angle distribution". However, Fig. S8 (the validation of
  the baseline) uses the five-parameter Gaussian Component Spins model, which includes the tilt mixing fraction $\xi$.
  Whether the superradiance runs use isotropic tilts or the LVK tilt mixture is **ambiguous (unconfirmed)**.
- Priors: $\log_{10}(m_b/{\rm eV})\sim\mathcal U(-14.5,-10)$ and $\sigma_{\rm SR}\sim\mathcal U(0,0.35)$. The rest follow the LVK GWTC-4/5 population papers.
- The delta-function variant [End Matter III, Eq. (2)] replaces the Gaussian with $C\,\delta(\chi-\chi_{\rm SR})$, which is the Ng+/Ning+ form.

### 1c. Null model(s) and Bayes factors
- The only null model is the **Gaussian component spin model** (the LVK baseline). It was re-run in their own pipeline with dynesty
  $n_{\rm live}=1000$, so that the evidence comparison is like for like [S4].
- Bayes factors: Table I (see 1a). Fixing $\sigma_{\rm SR}=0.1$ makes "$\ln\mathcal B$ decreases by an additional $\sim2$" [Discussions, para. 3].
- No Bayes factor is given for the delta-function model. Its posterior has "a region of vanishing support, approximately between
  $[3.5\times10^{-14},1.6\times10^{-12}]$ eV" (consistent with Ning+). The peak returns when GW241011 is removed:
  "the excluded region arises solely from this event" [End Matter III, Fig. 6].
- **No hierarchical-merger or astrophysical mass–spin-correlation null model is compared.** The authors warn:
  "caution is needed when interpreting our results from the perspective of model selection" [Discussions, last paragraph].

### 1d. Events and cuts
- "69 BBHs from GWTC-3.0, 153 BBHs from GWTC-4.0 and 259 events from GWTC-5.0 with false alarm rates $\le1\,{\rm yr}^{-1}$"
  [paragraph before Results]. Ning+ have 257, because they also cut on median component mass $\ge3\,M_\odot$.
- PE samples [S1, after Eq. (S4)]: O4b uses `IMRPhenomXPHM_SpinTaylor`. O4a uses `NRSur7dq4` where available, else `Mixed`
  (XPHM_SpinTaylor + SEOBNRv5PHM). Pre-O4a uses the GWTC-3.0/2.1 `Mixed` samples.

### 1e. Selection effects
- Yes. The hierarchical likelihood (S1) includes $\xi(\Lambda)$, estimated from "the publicly available set of BBH injections
  recovered by the search pipelines across all four observing runs" (Essick+ 2025, arXiv:2508.10638).
- Each injection is given an independent $T_{\rm age}$ draw [Eq. (S6)].
- Monte Carlo control [Eqs. (S7)–(S8)]: $N_{\rm eff}^{\rm inj}\ge4N_{\rm det}$ and $\min_i N^{\rm event}_{{\rm eff},i}\ge20$, applied as
  smooth penalties $\ln S(x)$ with $S(x)=1/(1+x^{-30})$. They use $N_{\rm eff}$ and not the LVK $\sigma^2_{\ln\hat{\mathcal L}}$
  threshold; S4 discusses the difference.
- Mock tests [S2] use GWMockCat injections ($7.5\times10^8$, $\rho_{\rm obs}\ge12$).

### 1f. Superradiance spin-down computation
- Code: "the public package SuperRad" [refs: Siemonsen+ 2211.03845; May+ 2410.21442]. In the paper's words, it uses
  "fully relativistic templates" for $m\le2$ and "semi-analytical formulae" up to $m\le15$ [paragraph before Results].
  $\chi_{\rm SR}$ "is the single spin value remaining after $T_{\rm age}$". **Scalar only.**
- BH age: $T_{\rm age}$ is a **per-event** source parameter with a log-uniform prior over [1 Myr, 2.2 Gyr]. It is marginalized by appending
  one independent draw to every PE sample and every injection [Eqs. (S2)–(S6)]. The lower bound is the minimum delay-time
  floor. The upper bound is age-of-Universe consistency for the highest-z event [paragraph before "Our dataset…"].
- Fixing $T_{\rm age}=10^6,10^7,10^8$ yr moves $m_b$ monotonically, by ~0.15 dex over two decades [S1, Fig. S1 right].
- Cloud depletion (GW emission or annihilation after saturation) and level interactions: **not discussed in the paper (unconfirmed).**
  The current SuperRad `spin_down.spindown()` (Bitbucket HEAD, read 2026-11-03) grows modes one after another. Each mode saturates after
  $t_c=\tau_I\ln(M_c/\mu)$, starting from a single boson, and the code then moves to $m+1$. It removes no GW dissipation and has no
  self-interaction (my reading of the code, not a statement by Kou+).
- Self-interactions: **not included**. "We also assume that the ultralight boson is only minimally coupled to gravity." Exploring
  interaction models "lies beyond the scope of this study" [Discussions, para. on vectors/interactions; cites Baryakhtar+ 2011.11646, Witte & Mummery 2412.03655].

### 1g. Hierarchical mergers / 2G BHs
- **Not modeled.** They appear only in the Introduction, as the standard explanation of the high-spin subpopulation at large masses
  (hierarchical mergers in clusters and AGN disks). There is no 2G component, no 2G natal spin ≈0.69, and no 2G-specific $T_{\rm age}$.

### 1h. Robustness and future work (as stated)
- Mock catalogs [S2, Figs. S4–S5]: 10 catalogs of 200 events each, with $q=1$, $\chi_2=0$ and aligned $\chi_1$. The injected
  $\log_{10}m_b=-12.1$ is recovered in all 10. Shifting the mass-function peak to 20 $M_\odot$ gives the same recovery.
  $\mu_\chi,\sigma_\chi$ are poorly recovered and depend on the realization [Results; S2].
- Leave-out tests: GW241011, GW250114, GW241113 [End Matter I]; GW190728 [End Matter II, Fig. 5]. $m_b$ is stable under all of them.
- Future work, as stated: a vector analysis [Discussions]; non-gravitational interactions ("beyond the scope"); a "modified BBH catalog
  that incorporates the effects of boson clouds"; extending the truncated-Gaussian-mixture likelihood (Hussain+) [End Matter III];
  mocks with unequal masses and two spins [S2 end]; O5.
- Data: the Zenodo DOI 10.5281/zenodo.22242056 (CC-BY-4.0, 146 MB) holds bilby HDF5 results for the extended, delta and baseline runs.
  The README says $\chi_{\rm SR}$ comes from SuperRad. **The pipeline code is not released (unconfirmed / not stated).**

---

## 2. Ning, Safdi, Welch — arXiv:2607.01317 (v1, 2026-07-01)

"No Evidence for Superradiant Axions in LIGO-Virgo-KAGRA GWTC-5 Binary Black Hole Spins".
Local files: `Refs/arxiv/2607.01317/{main.tex, 2607.01317.pdf, .txt}`.
Layout: an unnumbered PRL main text; End Matter (EM) "SR Spin-down Rates"; an SM with unnumbered sections. Equations (1)–(6) and
(S1)–(S5). Figures 1–5 and S1–S10.

### 2a. Result
- 95% HPD exclusion for GWTC-5: "$\log_{10}(m_a/{\rm eV})\in[-13.8,-11.5]$", i.e. $m_a\in[1.7\times10^{-14},3.3\times10^{-12}]$ eV
  [Results, Fig. 3].
- Bayes factor: "$\log B^B_A=0.22\pm0.06$" (the log base is not stated). Jeffreys prior: $[2.5\times10^{-14},2.2\times10^{-12}]$ eV [Results; SM Jeffreys, Fig. S3].
- There is a local posterior maximum at $\log_{10}m_a\approx-11.4$, near $m=3,4$ saturation curves that pass through several high-spin events [Results].
- Natal-spin posterior medians: $\alpha=0.69^{+0.21}_{-0.12}$, $\beta=2.28^{+0.62}_{-0.43}$ (the interval credibility is **unconfirmed**) [Results, last paragraph].

### 2b. Spin-down mapping
- **Fiducial**: analytic hydrogenic rate [EM Eq. (4)]
  $$\Gamma^{[n\ell m]}_{\rm SR}=2\alpha^{4\ell+5}M^{-1}(m\chi-2m_ar_+)\,C_{n\ell m}\,\mathcal K(\chi),\quad C_{011}=1/24,$$
  with the leading $\mathcal O(\alpha^2)$ correction in $\omega_R$ kept [main text "Axion SR…" para. 4; SM "Superradiance for higher α"].
- Modes: $\ell=m$, $n=0$, for $m=1$–6: "for our fiducial analysis we include modes up to $m=6$" [main text; EM].
- Saturation spin [EM Eq. (6)]: $\Gamma^{[n\ell m]}_{\rm SR}(m_a,M,\chi_{\rm sat})=N_{\rm efolds}/\tau_m$, with $N_{\rm efolds}=180$
  (also tested at 90 and 360, Fig. S10). It applies if $\tau_m>\tau_{\rm grow}$. The BH spins down successively to the saturated
  spin of the highest mode that can grow within $\tau_m$. The critical spin is EM Eq. (5).
- Forward map [SM "Hierarchical Likelihood", Eq. (S1)]:
  $p(\chi_M|\Lambda,\mathcal H_B)=p(\chi_M|\alpha,\beta)\Theta(\chi_{\rm sat}-\chi_M)+f_{\rm SR}\,\delta(\chi_M-\chi_{\rm sat})$.
  The likelihood terms are $\mathcal L^{\rm OFF}$ and $\mathcal L^{\rm SAT}$ [Eqs. (S2)–(S3)]. $\mathcal L^{\rm SAT}$ uses the single-event
  posterior **density** at $\chi_{\rm sat}$.
- **Witte & Mummery code (2412.03655)**: used only as a **systematic test**. Ning+ say they "use the code package released by
  Ref. [Witte:2024drg] to numerically solve the SR equations" without the small-α approximation (Leaver), up to $m=6$.
  The posteriors are "similar" [main text; SM "Superradiance for higher α", Fig. S6]. SuperRad is **not** used.
- BH age: $\tau_m$ is per-event, with a log-uniform DTD $p(\tau_m)\propto\tau_m^{-1}$ over $[10^7,\tau_H]$ yr and $\tau_H\approx1.4\times10^{10}$ yr,
  marginalized as a nuisance [main text]. The hyperparameter list has $\tau_m^1\ldots\tau_m^{2N}$, which is per **BH**; the text says
  "each BBH event". **This is inconsistent within the paper.**
- Alternatives [SM "Astrophysical Systematic Tests", Figs. S7–S8]: a shared $\tau_m$ with a log-uniform prior on [$10^5,10^9$] yr, and
  mass-dependent DTDs split at $M_s=45\,M_\odot$.

### 2c. Natal spin model
- **Confirmed**: the fiducial model is a single, mass-independent Beta distribution [Eq. (2)]: $p(\chi|\alpha,\beta)\propto\chi^{\alpha-1}(1-\chi)^{\beta-1}$,
  with priors $\log_{10}\alpha,\log_{10}\beta\in[-1,1]$. Both components use the same Beta (the null mocks draw "both components
  independently from the population Beta").
- Alternatives (all gave robust $m_a$) [SM "Astrophysical Systematic Tests", Figs. S7–S8, Eq. (S5)]:
  - a truncated Gaussian;
  - a two-component Beta mixture;
  - a mass-split Beta at $M_s=45\,M_\odot$;
  - a logistic mass-dependent Beta mixture, following Berti+ 2512.03152.

### 2d. Hierarchical mergers
- **Not explicitly modeled.** They appear only as the motivation for the mass-dependent tests ("plausibly hierarchical high-mass primary
  population" near 45 $M_\odot$ [SM "Mass-Dependent Distributions"]). There is no 2G spin ≈0.69 component and no 2G-specific age.

### 2e. Selection effects
- **None.** The hyperposterior, Eq. (3), has no $\xi(\Lambda)$ term. The full text never mentions selection functions or sensitivity injections
  (searched for "selection", "injection" and "detectab"). Its "Signal Injection Test" means synthetic SR catalogs.
- The mass prior is fixed in importance sampling: Salpeter $\pi(M_1)\propto M_1^{-2.35}$ and uniform $q$, following Ng+ [SM "Hierarchical Likelihood"].
- Data: GWTC-3/4/5 `Mixed` PE, with $K=2000$ samples per event. Cuts: FAR $<1\,{\rm yr}^{-1}$, BBH, $M_{1,2}^{\rm med}\ge3\,M_\odot$, no upper mass cut.
  This gives N = 69 / 153 / 257 [SM "GWTC Catalog Details"]. Sampler: dynesty with $n_{\rm live}=2000$.

### 2f. Self-interaction check (exact)
- Location: SM section "**Superradiance for higher α**" (unnumbered), last part of the paragraph discussing Fig. S6. **No figure shows the self-interaction run.**
- Method: the same Witte code base, with self-interactions switched on. $f_a$ is "the QCD axion decay constant", i.e. the QCD relation
  $m_a\approx0.57\,{\rm neV}(10^{16}\,{\rm GeV}/f_a)$ [Introduction]. Which levels ($n_{\max}$) were used is **unconfirmed**.
- Result: "generally negligible effect on the final spins when assuming the QCD axion decay constant". The number **0.05 is confirmed**:
  "the BH's final spin $\chi$ may be altered by up to $\sim$0.05 in some cases". Self-interactions there do not prevent modes that grow
  in the linear theory from doing so.
- Other $f$: "a more intensive analysis for other $f_a$ values could be an avenue for future work" [same paragraph].
  Discussion: for 3G detectors above $10^{-11}$ eV, "further work understanding the effects of axion self-interactions is needed".
- Also future work: tidal mixing beyond $m=2$ [SM "Tidal Perturbations", Fig. S9]; O4c and O5 [Discussion].

### 2g. Code and data availability
- **No code or data availability statement** in the full text (no GitHub or Zenodo link). The Witte code they used is public (see §4).

---

## 3. Novelty re-check (2025-01 → 2026-11-03)

Searches:
- INSPIRE `refersto:recid` for 2607.01317 (11 citing papers), 2609.02678 (3), 2011.11646 (322, of which 2025+ ≈ 95), 2412.03655 (69),
  2507.21788 (37) and 2610.11091 (0). The `refersto:arxiv:` form returned the whole database, so record IDs were used instead.
- INSPIRE `refersto` for 2011.06010, 1908.02312, 2201.11338, 1911.07862, 2608.11620, 2507.20979, 2503.10347 and 2107.11730.
- INSPIRE abstract and full-text keyword queries: superradiance × {population, self-interaction, hierarchical, GWTC, decay constant, 0.69}.
- arXiv API abstract searches.

All IDs below were confirmed with the arXiv API. "FT" means the full text was checked.

**Bottom line: no paper does the planned analysis.** That means no paper combines an LVK BBH spin population, hierarchical inference,
superradiance with quartic self-interaction, and a scan or exclusion in $(\mu,f)$. No paper uses the 2G natal spin as a calibrator either.
The newest review states the gap directly.

### 3a. Population-level, LVK spins (free field only)
| arXiv | Title | What it does / overlap |
|---|---|---|
| 2607.01317 | No Evidence for Superradiant Axions in LVK GWTC-5 BBH Spins (Ning, Safdi, Welch) | Free-field population exclusion. Self-interactions checked only at the QCD $f_a$ (Δχ ≲ 0.05). **Overlaps on data and method, but has no $f$ scan.** FT |
| 2609.02678 | Ultralight Bosons Explain the Mass–Spin Correlations… (Kou+) | Free-field population "evidence" model with joint mass, z, spin and selection. Self-interactions "beyond the scope". **Overlaps on data, no $f$.** FT |
| 2610.11091 | Systematic Uncertainties in Ultralight Boson Constraints from Superradiance (Aswathi, Sun) | **Review**, with no new analysis. Sec. V says the GWTC-2 and GWTC-5 population analyses "model superradiant spin-down in the free-field limit, without explicitly including self-interactions". It adds that extending them to interacting models would help. **This confirms the gap.** Sec. III: Kou+'s preferred mass lies inside the regions excluded by Ning+ and GW241011. FT |

### 3b. Population-level, a different observable, but with $f$ (partial overlap)
| arXiv | Title | Overlap |
|---|---|---|
| 2609.12052 | Constraining Ultralight Scalars with BBH Mergers in Galactic Nuclei (Khalaf, Kuflik, Lenoci, Stone, Telem) | Uses GWTC-5 and the **$(\mu,1/f)$ plane** (Fig. 1, source order). The observable is the cloud-enhanced GW-capture **merger rate** of $q<0.3$ events, **not spins**. It excludes $[4\times10^{-14},10^{-13}]$ eV for QCD-like $f$. "A full Bayesian analysis of the events is left for future work." **Same plane, different physics. Cite it as complementary.** FT |
| 2609.33596 | SGWB from self-interacting superradiant clouds (Guo, Brito, Yuan) | Self-interacting cloud evolution (Baryakhtar+ formalism) applied to BH populations, including BBH merger remnants, to compute an SGWB. Finds LIGO design is insensitive for $f_s\lesssim3\times10^{17}$ GeV. **No spin inference. Useful as a physics reference.** FT (grep) |
| 2602.23415 | Black hole scalar sirens in the Milky Way (Gavilan-Martin+) | Milky Way isolated-BH population with self-interaction ($f$) driven scalar emission. Not LVK. No overlap. (abstract) |
| 2603.15734 | Crowdsourcing GWs from Superradiant Axions (Ellis, Ning, Rodd, Schütte-Engel) | Population GW (Milky Way + SGWB) projections. Not a spin inference. No overlap. (abstract) |

### 3c. Individual-event $(\mu,f)$ constraints from LVK spins (not population)
| arXiv | Title | Overlap |
|---|---|---|
| 2507.20979 | Ultralight boson constraints from GW observations of spinning BBHs (Aswathi, East, Siemonsen, Sun, Jones) | GW231123 + GW190517, using SuperRad. The free-field exclusion holds for "$f\gtrsim10^{15}$ GeV" at $T_{\rm age}=10^5$ yr [Fig. 2 left; App. D "Scalar self-interactions", Eq. (D8)]: $\tau_{\rm sd}$ from Baryakhtar+ is only valid for $\alpha\le0.2$, and $f\propto T_{\rm age}^{-1/2}$. FT |
| 2507.21788 | Superradiance Constraints from GW231123 (Caputo, Franciolini, Witte) | GW231123 only. Exclusion in $(\mu,f_\Phi)$ with the Witte code at $N_{\max}=5$: $\mu\sim[0.6,5]\times10^{-13}$ eV and $f_\Phi\gtrsim10^{14}$ GeV. The age comes from the formation channel, with a hierarchical origin discussed. FT |
| 2503.10347 | Self-interaction effects on Kerr BH superradiance… (Xie, Huang) | Single-mode, relativistic, self-interacting cloud. Applied to GW190412 and GW190517 [Sec. "Observational consequences…", Table and Fig. for GW190517]. FT (grep) |
| 2510.26931 | GW241011 and GW241110 (LVK) | Free-field SuperRad constraint from one event. Its astrophysics section notes that 2G spins concentrate "about χ≈0.7". No self-interaction, no population. FT (grep) |
| 2608.11620 | LVK O4b/O4c ultralight-boson constraints | Free-field. HMM remnant searches plus high-spin events GW240515, GW241113, GW241225_08. No self-interaction. FT (grep) |

### 3d. 2G (hierarchical) natal spin ≈0.69 as a superradiance calibrator
- **None found in 2025–2026.** Related pre-2025 work: **2107.11730** (Payne, Sun, Kremer, Lasky, Thrane 2021, "The imprint of superradiance on
  hierarchical BH mergers"). There, superradiance lowers kicks and so raises hierarchical retention. Its Sec. 4.3 (Fig. 6, source order) says
  "a distinct 2G population with spins ∼0.7 is observed", present with or without the boson, because those BHs are "not old enough" for
  significant cloud growth. **This is a key caveat for a 2G-calibrator design**, since a 2G BH's age is the time since its previous merger. FT
- Population spin-model papers that motivate a 2G component but contain no superradiance: 2512.03152 (Berti+, cites superradiance only as
  motivation), 2605.24281 (Hussain+, finding mass-dependent spin subpopulations) and 2607.15431 (Zeeshan & O'Shaughnessy). Ning+ use the
  Berti+ logistic mixture only as a robustness test.

### 3e. Pre-2025 context (outside the window, for citation)
- 2011.06010 (Ng+ GWTC-2) says its bound is "valid for bosons with negligible self-interaction", with $f_a\gtrsim10^{14}$ GeV
  [Sec. "Constraints from GWTC-2"; the estimate uses the Yoshino–Kodama nonlinear condition]. It is free-field only. FT
- 1911.07862 (Fernandez+ 2019): LIGO data treated with a non-self-interacting boson; X-ray data with large self-interactions (abstract).
- 1908.02312 (Ng+ 2019, forecast) and 2201.11338 (Cheng+ 2022, $T_{\rm M}$-marginalized, μ only): abstracts only.

---

## 4. Code availability

- **Witte & Mummery (2412.03655)**: public. The paper's footnote in Sec. VI says "All code used to produce the results shown here can be found
  at [WitteGit]", which is https://github.com/SamWitte/Axion_SR_TDE. GitHub now redirects this to **https://github.com/SamWitte/Axion_SR**.
  The repo is public, **has no license file** (GitHub API `license: null`) and is written in Julia (plus Python and notebooks); last push 2026-11-03.
  Per the README, it computes Leaver/Heun eigenvalues, self-interaction transition rates, and coupled ODEs for $M,a$ and the level
  occupations. It "supports" $N_{\max}=8$ and has GW-event limit scripts (`Ligo_limit.jl`).
  Paper scope: coupled evolution of all $n\le5$ states, relativistic, with $\lambda=\mu^2/f_a^2$ [Sec. II.C–D].
  The fiducial X-ray limit uses only the $m\le2$ region [Fig. 1 caption]. Level mixing is negligible for $m=1$ at $\alpha\lesssim0.15$ [Sec. VI].
- **SuperRad (2211.03845)**:
  - Install: `pip install superrad`. The latest is v0.3.0 (uploaded 2026-01-07), depending on numpy and scipy. Home: https://bitbucket.org/weast/superrad
  - License: **BSD 3-Clause** (LICENSE file and PyPI classifier).
  - Coverage per the paper abstract: scalar and vector, the "lowest two azimuthal numbers ($m=1$ and 2)", relativistic plus fits to analytic
    results. Cloud frequency including the self-gravity shift [Sec. IV]. GW power and strain from linearized metric perturbations [Sec. V.B].
    Quasi-adiabatic cloud growth and decay [Sec. VI, Eq. (36)].
  - Self-interactions: "we neglect field self-interactions and non-minimal couplings" [Sec. III]. **No self-interactions.**
  - Current code: the scalar `NonrelScalar` model sets `_max_m = 15`, and `HybridCloudModel` combines it with the relativistic model.
    This matches Kou+'s "$m\le2$ relativistic, up to $m\le15$ semi-analytic".
  - Entry points: `spin_down.spindown(mu, Mbh0, abh0, tage, cloud_model)` and `max_spin(...)`.
- **Code used by 2609.02678**: SuperRad and GWPopulation/Bilby/dynesty, all public. The authors' own pipeline is **not released**
  (only posterior samples, on Zenodo 10.5281/zenodo.22242056).
- **Code used by 2607.01317**: their own hierarchical code is **not released, with no statement**. For the systematic and self-interaction
  checks they used the public Witte code (above).

## Downloaded sources (`Refs/arxiv/`)
2609.02678, 2607.01317, 2412.03655, 2610.11091, 2211.03845, 2507.20979, 2507.21788, 2503.10347, 2608.11620, 2510.26931,
2512.03152, 2607.15431, 2609.33596, 2609.12052, 2107.11730, 2011.06010. Each directory holds `src.tar.gz` and the extracted files, plus a PDF
and text export for the main ones.
