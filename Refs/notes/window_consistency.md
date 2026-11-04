# Is the window μ ≈ (0.6–0.9)×10⁻¹² eV, f ≈ 10¹²–10¹³ GeV already excluded? (full-text check, 2026-11-03)

Question: a first estimate found a window where an attractive quartic self-interaction (λ = μ²/f², axion-like)
keeps a young 2G BH (GW241011 primary, 19.6 M☉, χ = 0.78) spinning for ~10⁷ yr. This note checks every probe
against the **full texts** (LaTeX source + arXiv PDF) in `Refs/arxiv/<ID>/`.

Legend
- **[confirmed]**: stated in the full text (equation, figure, table or section given).
- **[figure read-off]**: number read from the figure image by pixel calibration of the log axes against the
  detected major ticks (temporary scripts, not saved; reading error about ±0.05 dex in f). Read-offs are not statements
  by the authors.
- **[my estimate]**: my own arithmetic with formulas from the papers; not a statement of any paper.
- Fine-structure constant: α = GMμ = 7.48×10⁹ (M/M☉)(μ/eV) **[my estimate]**.
- f always means the decay constant with λ = μ²/f²; "excludes f ≳ X" means the region from X up to the free-field
  limit (f → ∞) is excluded, while f < X (stronger self-interaction) is allowed.

Target box used below: **μ ∈ [0.6, 0.9]×10⁻¹² eV, f ∈ [10¹¹, 10¹⁴] GeV**; window claimed: f ∈ [10¹², 10¹³] GeV.

---

## 0. Bottom line

- **One probe cuts deep into the window: the X-ray-binary spin limit of Witte & Mummery (WM, 2412.03655), driven by
  GRS 1915+105.** Its fiducial limit (Fig. 1) excludes f ≳ 2.7, 2.3, 1.9 and 1.6 ×10¹² GeV at μ = 6, 7, 8 and 9 ×10⁻¹³ eV.
  Only **f ≲ (1.6–2.7)×10¹² GeV survives**, so roughly the lower quarter-decade of the claimed window survives. If the
  Cyg X-1 "fiducial" spin (ã > 0.9985) is trusted, Cyg X-1 alone excludes f ≳ (0.5–1.0)×10¹² GeV. Then essentially
  the whole claimed window is gone, and only f ≲ 10¹² GeV survives (WM Figs. 14 and 16).
- Older or weaker spin-based limits also cut into the box, but only above 10¹³ GeV. Baryakhtar+ (BGLS) Fig. 10 gives
  f ≳ (1.4–3.8)×10¹³ GeV. Hoof+ (M33 X-7, quasi-equilibrium) gives f ≳ (1.0–2.6)×10¹³ GeV.
- None of the following cuts into f ≤ 10¹⁴ GeV at this μ:
  - the GW231123 analyses;
  - the GW241011 LVK constraint (free field);
  - the O4b/O4c spin and continuous-wave (CW) results (free field or vector);
  - all continuous-wave and stochastic-background (SGWB) searches;
  - the merger-rate bound (Khalaf+);
  - any non-superradiance probe (none found).
- My scaling check gives the same picture **[my estimate]**: WM's GRS limit and a "Δχ ≲ 0.1 in 10⁷ yr" requirement for
  GW241011 sit at almost the same f (≈ 2–3×10¹² GeV at 7×10⁻¹³ eV). The surviving window is therefore
  **f ≈ (lower edge of the window) – 2×10¹² GeV**, not 10¹²–10¹³ GeV. This must be redone with the project's own kinetic
  model before it goes into the note (project rule).

---

## 1. X-ray binary spins with self-interactions

### 1a. Witte & Mummery, arXiv:2412.03655v3 (7 Apr 2025), "Stepping Up Superradiance Constraints on Axions"
Local: `arxiv/2412.03655/` (existing; read as is). The figure images are the JPEGs in that directory.

**Method [confirmed].** The paper uses a Bayesian MCMC over (f_a, m_a, M_i, ã_i) [Sec. III].
- Priors: log-flat in m_a ∈ [10⁻¹⁴, 10⁻¹⁰] eV and f_a ∈ [10⁹, 10¹⁸] GeV. Initial-spin and initial-mass priors are
  Eqs. (76)–(77); for lower-limit spins the initial-spin prior is flat up to ã_max = 0.998.
- Each system is evolved from t₀ − τ_max to t₀, with τ_max = Min[τ_bh, τ_Sal], τ_Sal ≈ 4.5×10⁷ yr, and the *lower*
  end of τ_bh used. The likelihood is Eq. (78). For lower-limit spins it is an asymmetric Gaussian centred on the
  2σ limit with a lower width of 10⁻³, so any present-day spin below the limit is excluded outright.
- The limit is the 2σ (95%) marginal posterior per m_a bin, with each BH analysed independently.
- Relativistic (Green's-function) scattering rates are used for all processes with n ≤ 5. Growth is capped at the
  bosenova threshold.
- Companion tidal mixing is checked with Eq. (79) and found irrelevant.

**BH inputs, Table III [confirmed]:**

| BH | Mass [M☉] (1σ) | ã (2σ, or 3σ threshold for limits) | τ_bh [yr] (τ used) | P_orb [d] | M_comp [M☉] |
|---|---|---|---|---|---|
| **Cyg X-1 [Fid.]** | 21.2 ± 2.2 | > 0.9985 | 4.8–7.6×10⁶ (**4.8×10⁶**) | 5.599829 | 40 |
| **Cyg X-1 [Cnsrv.]** | 21.2 ± 2.2 | 0.92⁺⁰·⁰⁵₋₀.₀₇ (warm-corona spectral model, Zdziarski+24) | 4.8–7.6×10⁶ (**4.8×10⁶**) | 5.599829 | 40 |
| **GRS 1915+105** | 12.4⁺²·⁰₋₁.₀ | > 0.95 | 3–5×10⁹ (**τ_Sal ≈ 4.5×10⁷**) | 33.85 | 0.47 ± 0.27 |
| M33 X-7 | 15.65 ± 1.45 | 0.84 ± 0.1 | 2–3×10⁶ | 3.4530 | ≳ 20 |
| LMC X-1 | 10.91 ± 1.4 | 0.92⁺⁰·⁰⁶₋₀.₁₈ | 5–6×10⁶ | 3.9093 | 31.79 ± 3.48 |
| GRO J1655-40 | 6.3 ± 0.5 | 0.72⁺⁰·¹⁶₋₀.₂₄ | 3.4–10×10⁸ | 2.622 | 2.3–4 |

The Table III footnote rejects the "conservative" Cyg X-1 age of ~0.05 Myr from the jet-shock age (a lower limit only).

**Fig. 1 (fiducial; "2σ"; Cyg X-1 Cnsrv. + GRS 1915+105; n ≤ 5; m ≤ 2 spin-down region in dark blue).**
The y-axis is 1/f_a. The excluded region runs from the curve down to 1/f = 10⁻¹⁸ GeV⁻¹, i.e. from the f on the curve
up to the top of the prior and the free-field limit. **[figure read-off]**

| m_a [eV] | 6×10⁻¹³ | 7×10⁻¹³ | 8×10⁻¹³ | 9×10⁻¹³ | 10⁻¹² |
|---|---|---|---|---|---|
| lowest f excluded (Fig. 1) [GeV] | 2.7×10¹² | 2.3×10¹² | 1.9×10¹² | 1.6×10¹² | 1.3×10¹² |

The text gives no number at this mass. It says only that the Fig. 1 limit is the main result [Sec. VI] and that the
procedure is ~10× more sensitive in f than BGLS for Cyg X-1 and ~2× for M33 X-7 [Sec. III, last paragraph; Fig. 18].

**Per-BH curves [figure read-off]**, lowest f excluded [GeV]:

| m_a [eV] | 6×10⁻¹³ | 7×10⁻¹³ | 8×10⁻¹³ | 9×10⁻¹³ | source |
|---|---|---|---|---|---|
| GRS 1915+105, n ≤ 3 | 3.0×10¹² | 2.2×10¹² | 1.7×10¹² | 1.5×10¹² | Fig. 14 (dark red) |
| GRS 1915+105, n ≤ 5 | 1.9×10¹² | 1.5×10¹² | 1.3×10¹² | 0.96×10¹² | Fig. 16 right (green) |
| Cyg X-1 Fid., n ≤ 3 | 1.0×10¹² | 0.78×10¹² | 0.67×10¹² | 0.56×10¹² | Fig. 14 (gold solid) |
| Cyg X-1 Fid., n ≤ 5 | 0.92×10¹² | 0.75×10¹² | 0.63×10¹² | 0.51×10¹² | Fig. 16 left (green) |
| Cyg X-1 Cnsrv., n ≤ 3 | 1.1×10¹³ | 8.3×10¹² | 6.7×10¹² | 6.4×10¹² | Fig. 14 (gold dotted) |
| M33 X-7, n ≤ 3 | 2.8×10¹³ | 2.0×10¹³ | 1.8×10¹³ | 1.4×10¹³ | Fig. 14 (navy) |
| LMC X-1, n ≤ 3 | 3.8×10¹³ | 2.7×10¹³ | 1.9×10¹³ | 1.5×10¹³ | Fig. 14 (green) |

- The Fig. 1 boundary at this μ is set by **GRS 1915+105**; the Cnsrv. Cyg X-1 curve is ~4× weaker in f. The Fig. 1 values
  match the n ≤ 3 GRS curve of Fig. 14 and are 0.15–0.2 dex above the n ≤ 5 GRS curve of Fig. 16. These are different
  runs; WM state the n ≤ 3, 4 and 5 results "converge in the m = 1 spin down region" [Sec. III], so the spread is
  MCMC/binning noise. Quote 1.5–2.7×10¹² GeV as the GRS range at 6–9×10⁻¹³ eV.
- The Fig. 14 caption says these per-BH curves are 95% limits.
- Fig. 18 left reruns Cyg X-1 with n ≤ 3 and the old mass M = 14.8 M☉ ("This work"). It excludes f ≳ 1.6, 1.1, 0.93 and
  0.72 ×10¹² GeV at 6, 7, 8 and 9 ×10⁻¹³ eV **[figure read-off]**.
- WM's bosenova-only treatment ("BN Cut", Fig. 15) extends the limits to "as much as an order of magnitude smaller f_a"
  [Sec. III, text on Fig. 15]. Level mixing is therefore the *conservative* choice.

**m = 1 or m = 2, and α.** At 6–9×10⁻¹³ eV the limit is the **m = 1 (|211⟩) spin-down region**, with |322⟩ as the
self-interaction sink.
- α values **[my estimate]**: GRS 1915+105 0.056–0.083; Cyg X-1 (21.2 M☉) 0.095–0.143. These are below the m = 1
  superradiance threshold for both BHs (≈ 0.33 for ã = 0.92; 0.36–0.47 for ã = 0.95–0.998).
- WM: higher-order mixing is "negligible ... in the m = 1 spin down at α ≲ 0.15" [Sec. VI]. Fig. 16 shows n ≤ 3, 4, 5
  coincide below ~2×10⁻¹² eV.
- The m = 2 region (and the dip between the m = 1 and m = 2 regions) lies at m_a ≳ 2.5×10⁻¹² eV in Figs. 1 and 16.
  Fig. 12 (Cyg X-1 parameters, τ = 4×10⁶ yr) marks the approximate m = 1 spin-down in green at α ≲ 0.45.

**Caveats that matter for the window:**
1. The X-ray spin systematics are discussed at length in Sec. IV. For Cyg X-1 the spin choice moves the limit by
   ~10× in f (Fid. vs Cnsrv.). No "conservative" alternative is given for GRS 1915+105 (ã > 0.95, Reid+ 2014).
2. τ_max = 4.5×10⁷ yr for GRS, with no accretion torque during that time. WM note that GRS 1915+105 "has reached rates
   close to its Eddington limit during its current 30 year outburst" [Sec. III], but treat this only as a caveat.
   **[my estimate]** In the regime that matters here the spin-down is linear in time. Excluding f only needs
   Δã ~ (f/f_CD)² ~ 0.04–0.09 over 4.5×10⁷ yr, so the bound is sensitive to any spin-up by accretion over that time.
   This is not quantified in WM.
3. It is a Bayesian 95% limit with log-flat priors in (m_a, f_a); Hoof+ show such limits can depend on the prior.
4. WM give no compact f(μ) formula; the numbers above are figure read-offs only.

### 1b. Baryakhtar, Galanis, Lasenby, Simon (BGLS), arXiv:2011.11646v2 (16 Feb 2021), Fig. 10
Local: `arxiv/2011.11646/` (existing; `src/fabounds.pdf` rendered for reading).

**Method [confirmed].** Fig. 10 caption: for each BH the solid line encloses "the intersection of the m = 1 spin-down
regions for different BH parameters (mass, spin, lifetime, binary period, and mass of the binary companion) within the
observational error intervals". The light shading is the union and is *not* claimed as a constraint. For m = 2 only
regimes A and B are included.
- Spin-down follows the two-level (211–322) quasi-equilibrium, with τ_sd Eq. (60) and the f_CD criterion Eq. (62).
- WM, citing a private communication, say BGLS equated the BH age with the quasi-equilibrium spin-depletion time and
  truncated m = 2 at f ≲ 10¹⁷ GeV.
- The low-μ edge is set by companion tidal mixing [App. I]. α ≲ 0.2 is required [Sec. VI].

**Numbers [figure read-off]**: Cyg X-1 solid line (the strongest BH at this μ), lowest f excluded:
3.8×10¹³ (6×10⁻¹³ eV), 2.6×10¹³ (7×10⁻¹³ eV), 1.9×10¹³ (8×10⁻¹³ eV), 1.4×10¹³ (9×10⁻¹³ eV) GeV. This agrees with WM's
replot of BGLS (WM Fig. 1, white dashed; WM Fig. 18 left, green): 4.1, 2.8, 2.0 and 1.5 ×10¹³ GeV.

**Cyg X-1 parameters in BGLS.** **BGLS do not list them.** The text and Fig. 10 caption only cite Miller & Miller 2014
and McClintock+ 2014. WM say BGLS used the older mass "M ∼ 14.8 ± 1.0 M☉, and a spin roughly consistent with the
fiducial model" [WM Sec. III].
- BGLS cite Baryakhtar, Lasenby & Teo 2017 (arXiv:1704.05081, downloaded) for the same five BHs. Its Table I gives
  Cyg X-1: **14.8 ± 1.0 M☉, a* > 0.99, age 4.8–7.6×10⁶ yr, P = 5.599829 d, M_comp = 17.8 M☉**.
  For GRS 1915+105 it uses τ_bh = τ_Salpeter/10.
- That BGLS used this table is my inference (same group, same BH list, same parameter columns as the Fig. 10
  caption), not a BGLS statement.
- α(14.8 M☉) = 0.066–0.100 at 6–9×10⁻¹³ eV **[my estimate]**, so the limit is m = 1 (211) in the two-level regime.

### 1c. Hoof, Marsh, Sisk-Reynés, Matthews, Reynolds, arXiv:2406.10337v2 (9 Oct 2025) — M33 X-7
Downloaded to `arxiv/2406.10337/`. This is pre-2025, but the review 2610.11091 treats it as a (μ, f) X-ray limit.
- Method [confirmed]: Bayesian timescale analysis with |211⟩ only (n = 2). M33 X-7 (M, a*) posterior samples come from
  digitised Liu+ 2008, with τ_BH = τ_age = 3×10⁶ yr [Sec. 3]. BGLS rates are used for the "equilibrium" scenario and
  the ABH15 bosenova condition for the "bosenova" scenario [Sec. 2]. Limits are 95% HPD.
- Fig. 4 right **[figure read-off]**, lowest f excluded:
  - equilibrium: 2.6×10¹³, 1.9×10¹³, 1.4×10¹³ and 1.2×10¹³ GeV at 6, 7, 8 and 9 ×10⁻¹³ eV;
  - bosenova: 3.6×10¹², 2.2×10¹², 1.9×10¹² and 1.7×10¹² GeV.
- The equilibrium values agree with WM Fig. 14 (M33 X-7: 2.0×10¹³ GeV at 7×10⁻¹³ eV).

### 1d. Summary for probe 1
- At μ = 6–9×10⁻¹³ eV, WM (fiducial, Fig. 1) excludes f ∈ [~2×10¹², 10¹⁸] GeV. Of the target box, **f ∈ [10¹¹,
  (1.6–2.7)×10¹²] GeV survives** and f ∈ [(1.6–2.7)×10¹², 10¹⁴] GeV is excluded.
- With Cyg X-1 Fid., only f ≲ (0.5–1.0)×10¹² GeV survives.
- BGLS alone would leave f ≲ 1.4–3.8×10¹³ GeV.

---

## 2. GW231123

### 2a. Caputo, Franciolini, Witte (CFW), arXiv:2507.21788v2 (27 Mar 2026)
Local: `arxiv/2507.21788/` (existing).
- Assumptions [confirmed]:
  - NRSur PE samples, M₁ = 137⁺²²₋₁₇ M☉, ã₁ = 0.90⁺⁰·¹⁰₋₀.₁₉, M₂ = 103⁺²⁰₋₅₂ M☉, ã₂ = 0.80⁺⁰·²⁰₋₀.₅₁ [Introduction].
  - WM machinery (Eqs. (1)–(3), WM rate tables), fiducial N_max = 5.
  - τ_max = 10⁵ (shaded), 10⁶ (dashed) and 10⁷ yr (dot-dashed), motivated by hierarchical merger times
    [Sec. "Determination of merger timescales"].
- Abstract [confirmed]: "conservative constraints on axion masses μ ∼ [0.6–5]×10⁻¹³ eV and decay constants
  f_Φ ≳ 10¹⁴ GeV". The 10¹⁴ GeV is the lowest f reached anywhere, around μ ~ 1.3×10⁻¹³ eV, not at our μ.
- **At our μ (Fig. 1, N_max = 5) [figure read-off]:** the N_max = 5 region extends past the abstract's 5×10⁻¹³ eV, in
  the m ≥ 2 regime (α(137 M☉) ≈ 0.6–0.7).
  - At 6×10⁻¹³ eV the curves sit at f ≈ 0.8×10¹⁶ (τ = 10⁷ yr) to ≈ 3×10¹⁶ GeV (τ = 10⁵ yr).
  - The shaded (10⁵ yr) region ends at μ ≈ 6.95×10⁻¹³ eV; the 10⁶ and 10⁷ yr curves end at ≈ 7.1–7.2×10⁻¹³ eV.
  - **At 7×10⁻¹³ eV CFW exclude only f ≳ 10¹⁶ GeV, and nothing above ~7.2×10⁻¹³ eV.**
  - With N_max = 3 (Fig. 3 left) the region stops at ≈ 5×10⁻¹³ eV.
- CFW footnote (Introduction) on GW241011: its mass lies in the X-ray-binary range, so "the derived constraints are
  sub-dominant to pre-existing limits".
- **Verdict: no cut into f ≤ 10¹⁴ GeV at 0.6–0.9×10⁻¹² eV.**

### 2b. Aswathi, East, Siemonsen, Sun, Jones, arXiv:2507.20979v2 (6 Jan 2026), App. D
Local: `arxiv/2507.20979/` (existing).
- **Free-field result [confirmed]**: SuperRad, relativistic for scalar m ≤ 2, primary BH only, 90% (P < 0.1P′).
  At T_age = 10⁵ yr it excludes [0.55, 3.3]×10⁻¹³ eV (GW231123) and **[2.2, 11]×10⁻¹³ eV (GW190517)**
  [Sec. "Constraints", Fig. 1]. **GW190517 therefore covers 6–9×10⁻¹³ eV, but only for a free field.**
- **Self-interacting extension, App. D "Scalar self-interactions" [confirmed]**:
  - m = 1 only, with τ_sd from BGLS, Eq. (D8):
    τ_sd = (√6/m_b) Γ^{322×BH}_{211×211} / [(ω_I/m_b)^{3/2} (Γ^{211×∞}_{322×322})^{1/2}].
  - It is used only for α ≤ 0.2, with τ_sd set to 10⁵ yr.
  - Main text: the bounds "apply to axion models with energy scale f ≳ 10¹⁵ GeV", and f ∝ T_age^{−1/2}.
- Fig. 2 left **[figure read-off]**:
  - The GW231123 band covers ≈ 0.55–1.6×10⁻¹³ eV.
  - The GW190517 band covers ≈ 2.1–5.0×10⁻¹³ eV; its lowest f is ≈ 2.5×10¹⁴ GeV near 4.5×10⁻¹³ eV, and ≈ 4×10¹⁴ GeV at
    its upper edge.
  - The bands stop because of the α ≤ 0.2 cut. The caption says the dotted lines mark the free-field mass range that a
    better treatment of interactions could eventually cover.
- **Verdict: nothing at 6–9×10⁻¹³ eV in the self-interacting case. In the free-field limit GW190517 excludes this μ.**

---

## 3. GW-signal searches (CW, directed, SGWB)

| Search | Full text | What it says near 7×10⁻¹³ eV | Self-interaction treatment |
|---|---|---|---|
| **Sun, Brito, Isi**, Cyg X-1 directed CW, O2 (1909.11267v3) | downloaded | Free field: disfavors 6.4–8.0×10⁻¹³ eV (t_age = 5×10⁶ yr) and 6.3–13.2×10⁻¹³ eV (t_age = 10⁵ yr) [Fig. 3; Sec. "Disfavored boson mass"]. Assumes M = 14.8 M☉, χ_i = 0.99 and *unknown* post-superradiance spin (χ_f ≲ 0.6, i.e. it does not use the measured spin). Bosenova "string axiverse" model: f_a ~ 10¹⁵ GeV excluded at 9.6–15.5×10⁻¹³ eV [Fig. 4] | Bosenova picture (Yoshino–Kodama), later superseded (next row) |
| **Collaviti, Sun, Galanis, Baryakhtar** (2407.04304v2) | downloaded | Reanalysis of the same O2 Cyg X-1 limits with BGLS two-level kinetics (Table 1: M = 14.8 M☉, χ_i = 0.99, ages 10⁵ and 10⁶ yr, α ≤ 0.2). Disfavored range ≈ [0.6, 1.4]×10⁻¹² eV "mostly in the gravitational regime" [Sec. 4.2]. Fig. 7 (blue) **[figure read-off]**: lowest f ≈ 3.8×10¹⁷, 6.7×10¹⁷, 9×10¹⁷ and 1.2×10¹⁸ GeV at 6, 7, 8 and 9 ×10⁻¹³ eV. The f ~ 10¹⁵ GeV bosenova exclusion is withdrawn: "we do not find that ultralight scalars with ... f ∼ 10¹⁵ GeV ... could be detected" [Sec. 4.2]. No observable GW in the strong or no-spin-down regimes even for ET/CE [Sec. 5] | Full BGLS two-level treatment |
| **LVK O3 all-sky boson-cloud search** (2111.15507v2) | downloaded | 20–610 Hz, so μ up to ~1.26×10⁻¹² eV. Exclusions are given only in the (m_b, M_BH) plane for fixed D, χ_i, t_age [Figs. 6–7], plus D_max for a Galactic population [Fig. 8]: 2–8×10⁻¹³ eV excluded within 1 kpc for M_max = 50 M☉. Probes α ≈ 0.02–0.13 [Sec. V] | Covers the "small self-interaction" regime F_b ≳ 2×10¹⁸ GeV [Sec. II, Eq. (6)] and partly the intermediate one [Eq. (7)]; the strong regime is not covered [Sec. V, last paragraph] |
| O4 all-sky scalar boson-cloud search | — | **None found** (arXiv/INSPIRE searches, 2026-11-03). The O4a all-sky isolated-NS CW search (2603.14168, abstract only) has no boson interpretation in its abstract | — |
| LVK O4a directed searches (2509.07352; abstract only) | abstract | **Vector bosons only** (Cyg X-1: [0.85, 1.59]×10⁻¹³ eV, 95%) | — |
| **LVK O4b/O4c**, 2608.11620 | existing | Free-field spin limits (GW240515, GW241113, GW241225_08; SuperRad m ≤ 2 relativistic; T_age 10⁵–10⁷ yr): scalar [1.39, 6.94]×10⁻¹³ eV at 10⁵ yr [abstract; Fig. 4]. The upper edge (GW241113) stays at ≈ 6.9–7.1×10⁻¹³ eV for all ages; an extra island at ~1.3–2.3×10⁻¹² eV appears for t ≳ 2×10⁶ yr [figure read-off]. The CW HMM remnant searches are vector-only, [2.80, 3.95]×10⁻¹³ eV | None ("if gravitational interactions dominate", Introduction; Sec. IV uses SuperRad) |
| **Yuan, Jiang, Huang**, SGWB O1–O3 (2204.03482v2) | downloaded | "minimally coupled" scalars: excludes [1.5, 15], [1.8, 8.1] and [1.3, 17]×10⁻¹³ eV for natal spins uniform in [0,1], [0,0.5] and [0.5,1] [abstract; Table I] | None |
| Tsukada+ 2019 SGWB O1 (1812.09622) | not read (via review) | [2.0, 3.8]×10⁻¹³ eV, optimistic spins; does not reach 6×10⁻¹³ eV | None |
| **Zhu+**, Galactic ensemble vs O2 all-sky (2003.03359v2) | downloaded | Non-detection "disfavors" 2×10⁻¹³–2.5×10⁻¹² eV for a 10⁸-BH population (assumption-dependent) [abstract; Sec. V] | "We assume the bosons only interact via gravity" [Sec. II]; the footnote says this holds for the QCD axion (f ≳ 10¹⁸ GeV) |
| **Guo, Brito, Yuan**, SGWB with self-interaction (2609.33596) | existing | Two-mode BGLS kinetics, α ≤ 0.22 (0.15 tested). **"LIGO detectors at design sensitivity are insensitive to a SGWB produced by scalar fields with a decay constant f_s ≲ 3×10¹⁷ GeV, independently of the scalar field mass"** [abstract; Sec. IV C; Fig. 3 top]. ET and CE need f_s ≳ 6×10¹⁶ and 3×10¹⁶ GeV | Yes. Same sign convention, V ⊃ −λΦ⁴/4!, λ = μ²/f_s² [Sec. II] |
| **BGLS, Sec. VII A ("Annihilations")** (2011.11646) | existing | In the strong-coupling regime the peak annihilation strain drops ∝ (f/f_BC)². "self-interactions strong enough to avoid BH spin-down constraints ... also render GW annihilation signals undetectable". For f ≲ f_BC ≈ 10¹⁶ GeV the expected number of blind-search events is ≲ 10⁻³. (The source writes "f > f_BC" for the strong regime, a typo for f < f_BC; f_BC > f_CD is stated.) | Theory statement |

**Verdict for probe 3:** no GW-signal search constrains f ≲ 10¹⁴ GeV at 0.6–0.9×10⁻¹² eV.
- The lowest f reached by any of them at this μ is ≈ 4×10¹⁷ GeV (Collaviti+ Fig. 7).
- The self-interacting SGWB paper puts current LIGO's reach at f ≳ 3×10¹⁷ GeV.
- This is consistent with BGLS's argument that annihilation GWs are suppressed in the strong-coupling regime.

---

## 4. GW241011's own boson constraint — LVK, arXiv:2510.26931
Local: `arxiv/2510.26931/` (existing; LaTeX only, no PDF, so figure numbers follow source order).
- Event [confirmed]: m₁ = 19.6⁺³·⁶₋₂.₅ M☉, χ₁ = 0.78 ± 0.09 (90%), 90% lower limit 0.69 (`pe_macros_241011.tex`).
  The paper discusses a hierarchical (2G) origin [Sec. "GW241011 and GW241110 as hierarchical mergers"].
- Constraint [confirmed]: SuperRad, all relevant azimuthal modes, linear (free-field) instability. The excluded
  region is P < 0.1P′, i.e. ≥ 90% and likelihood-driven [App. "Constraining ultralight bosons through superradiance"].
  At T_age = 10⁵ yr: **scalar [0.3, 2.6]×10⁻¹² eV**, vector [0.1, 5.3]×10⁻¹² eV [Sec. "Ultralight bosons";
  `superrad_macros.tex`; Fig. 11, source order].
- Fig. 11 top **[figure read-off]**: at T_age = 10⁷ yr the scalar range is ≈ [0.19, 4.0]×10⁻¹² eV. A notch at ≈ 2.6–3.8×10⁻¹²
  eV is not excluded for ages below ~4×10⁵ yr.
- Assumptions [confirmed]:
  - "This analysis assumes only a minimally coupled boson with gravitational interactions, but will also apply to
    scalar or vector bosons with sufficiently weak interactions so as to not disrupt the black hole spindown."
  - The companion's effect is neglected.
  - The age is the time since formation or since the BH gained its angular momentum.
- **Verdict:** this is the free-field exclusion the window is built to evade. It says nothing about f ≲ 10¹⁴ GeV.
  CFW (with Witte) call GW241011's constraints "sub-dominant to pre-existing limits" (see 2a).

---

## 5. Merger rates in galactic nuclei — Khalaf, Kuflik, Lenoci, Stone, Telem, arXiv:2609.12052
Local: `arxiv/2609.12052/` (existing).
- Observable [confirmed]: GW-capture binary formation in galactic nuclei is enhanced by the cloud's energy
  dissipation. The predicted LLO rate of q < 0.3 detections is compared with GWTC-5 counts at 90% [Fig. 1 caption;
  Fig. 2]. The BH mass function in nuclei ∝ M^−β with β = 0.5, 1, 1.5, and M_max = 100 M☉ in Fig. 2. "A full Bayesian
  analysis of the events is left for future work."
- Excluded region, Fig. 1 in the (μ, 1/f = √λ/μ) plane **[figure read-off]**: μ ≈ 3.5×10⁻¹⁴–1.3×10⁻¹³ eV (β = 0.5) or
  ≈ 5×10⁻¹⁴–1.1×10⁻¹³ eV (β = 1.5), for 1/f ≲ 1.3×10⁻¹⁷ GeV⁻¹, i.e. **f ≳ 8×10¹⁶ GeV** (≳ 3×10¹⁷ for β = 1.5).
  Abstract: μ ∈ [4×10⁻¹⁴, 10⁻¹³] eV. Smaller f weakens the bound because the cloud is lighter [text after Fig. 2].
- **Verdict: does not reach 7×10⁻¹³ eV (it stops a factor ~5 lower in μ) and does not reach f ≤ 10¹⁴ GeV.**

---

## 6. Anything else

### 6a. INSPIRE citations of 2412.03655 (recid 2856292) and 2011.11646 (recid 1832641), 2025-01 → 2026-11-03
Queries: `refersto:recid:2856292 and de > 2024` and the same for 1832641; 69 and 322 total citations, of which the 2025–2026
lists were screened by title and abstract, and by full text where relevant. Papers that set or discuss (μ, f) limits from
stellar-mass BHs near 10⁻¹² eV:
- **2507.21788 (CFW), 2507.20979 (Aswathi+), 2609.12052 (Khalaf+), 2609.33596 (Guo+)**: covered above; none cuts into
  f ≤ 10¹⁴ GeV at 0.6–0.9×10⁻¹² eV.
- **2503.10347 (Xie & Huang)**, local full text: single-mode, relativistic, self-interacting cloud; no exclusion.
  Self-interaction stops |211⟩ growth at Ñ_c ≪ N_c. For GW190517 the primary can keep its spin within the 95% region
  if f_a ≲ 1.91×10²⁰ GeV [Sec. V A; Fig. 7, source order]. This only *relaxes* free-field limits.
- **2607.01317 (Ning, Safdi, Welch)**: free-field GWTC-5 population, excluding [1.7×10⁻¹⁴, 3.3×10⁻¹²] eV (95% HPD).
  Self-interactions are tested only at the QCD f_a (≈ 8×10¹⁸ GeV at 7×10⁻¹³ eV **[my estimate]** from their
  m_a ≈ 0.57 neV (10¹⁶ GeV/f_a)), with Δχ ≲ 0.05; other f_a are "future work" [SM "Superradiance for higher α"; see
  `population_papers.md` §2f].
- **2609.02678 (Kou+)**: free-field population *evidence* at m_b = 7.1⁺¹·⁸₋₁.₁×10⁻¹³ eV; self-interactions "beyond the
  scope" (`population_papers.md` §1).
- **2610.11091 (Aswathi & Sun review)**, local full text: Sec. V confirms that population and SGWB/CW interpretations are
  free-field. Its 2-D (mass–coupling) entries at this μ are WM, CFW, Aswathi+, Hoof+, Collaviti+ and Sun–Brito–Isi,
  all covered here [Sec. V; Appendix table].
- Abstract-only, no limits:
  - 2602.23415 ("BH scalar sirens"): for f ≲ 10¹⁴–10⁹ GeV, Milky Way BHs emit scalars persistently; a detection target
    needing SM couplings, not a limit. It is relevant to the window physics (the BGLS Sec. VIII axion waves).
  - 2603.15734 (Ellis+): LIGO population projections.
  - 2608.18200: needs g_ae.
  - 2510.17967: scalar density around binaries; no f.
  - 2604.21239, 2606.10776, 2512.17887, 2512.23279, 2608.02051, 2512.06790, 2509.10892: theory only.
- **No 2025–2026 paper other than WM's own limit (re-plotted as "X rays" in CFW Fig. 1 and Khalaf+ Fig. 1) puts a
  stellar-mass (μ, f) limit inside the target box.**

### 6b. Non-superradiance bounds on a pure-gravity ALP with μ ~ 10⁻¹² eV, f ~ 10¹² GeV
**None found.** What exists:
- **SM-coupling bounds** need extra couplings, and even then reach only small f. Collaviti+ Fig. 7 caption: "Axion
  interactions with Standard Model photons, electrons, and nuclei result in exclusions in typical models of
  GeV/f ≳ 10⁻⁸–10⁻⁹", i.e. f ≲ 10⁸–10⁹ GeV, far below the window.
- **Dark-matter abundance is a consistency condition, not an exclusion.** BGLS Sec. VI A and the Fig. 10 "ALP DM" band:
  ρ_DM ∝ μ^{1/2} θ_init² f². Read off Fig. 10 at 7×10⁻¹³ eV **[figure read-off]**, misalignment gives all of the DM
  only for:
  - f ≈ (1.5–4)×10¹⁴ GeV for θ_init ∈ (1, π−1);
  - f ≈ 4×10¹⁵ GeV for θ_init = 0.1;
  - f ≈ 2×10¹³ GeV for θ_init = π − 10⁻⁶ (top edge, "large misalignment").

  For f ~ 10¹²–2×10¹² GeV the scalar would therefore be a **subdominant DM component** (or not DM) unless it is produced
  otherwise. WM Fig. 1 shows the same "ALP Misalign." band at 1/f ~ 10⁻¹⁵·⁵–10⁻¹⁴·⁵ GeV⁻¹ near 10⁻¹² eV. This matters
  only if the paper wants the scalar to be all of the DM.
- **Late-time DM self-scattering** is irrelevant. BGLS Eq. (96): Γ ~ 3×10⁻²⁶ yr⁻¹ (ρ/GeV cm⁻³)² (10¹¹ GeV/f)⁴
  (10⁻¹² eV/μ)³ (10⁻³/v)².
- Isocurvature, strings and domain walls depend on the cosmological history (whether the symmetry breaks before or
  after inflation, and on H_inf). They constrain the inflation scale or the DM fraction, not (μ, f) by themselves.
  No paper found sets such a bound at this (μ, f).

---

## 7. Consistency of the window with the X-ray limit — scaling check **[my estimate]**
Leading-order BGLS κ form of Eq. (60): τ_sd = (√3/2) κ^BH r̃₊ (M_pl/f)² / [α⁵ μ √κ^∞ (κ^SR(a − 2α r̃₊))^{3/2}],
with κ^SR = 1/24, κ^BH = 4.3×10⁻⁷, κ^∞ = 1.1×10⁻⁸ and M_pl = 1.22×10¹⁹ GeV. Setting τ_sd(a₀) = T gives f_CD.
For f ≪ f_CD, |Δa| ≃ T/τ_sd = (f/f_CD)² (BGLS Eq. (63)).

| object (M, a₀, T) | f_CD at 6, 7, 8, 9 ×10⁻¹³ eV [GeV] |
|---|---|
| GRS 1915+105 (12.4, 0.97, 4.5×10⁷ yr) | 1.2×10¹³, 7.8×10¹², 5.3×10¹², 3.8×10¹² |
| Cyg X-1 Cnsrv. (21.2, 0.92, 4.8×10⁶ yr) | 1.2×10¹³, 8.2×10¹², 5.8×10¹², 4.3×10¹² |
| GW241011 primary (19.6, 0.78, 10⁷ yr) | 1.4×10¹³, 9.4×10¹², 6.9×10¹², 5.3×10¹² |
| GW241011 primary (19.6, 0.78, 10⁶ yr) | 4.4×10¹³, 3.0×10¹³, 2.2×10¹³, 1.7×10¹³ |

- WM's GRS limit (1.5–2.3×10¹² GeV at 7×10⁻¹³ eV) is ≈ 0.2–0.3 f_CD(GRS). This matches their likelihood, under which a
  present-day spin below 0.95 is excluded, i.e. |Δã| of only a few ×10⁻².
- "Keeps its spin for 10⁷ yr" with f ≈ f_CD ≈ 10¹³ GeV means no *full* spin-down. But GW241011 is measured at
  χ = 0.78 ± 0.09, and a 2G BH is born at ~0.7. Requiring |Δχ| ≲ 0.1 gives f ≲ 0.3 f_CD ≈ 3×10¹² GeV at 7×10⁻¹³ eV
  (10⁷ yr). That is about the same f as WM's GRS exclusion.
- So the part of the window that survives WM is f ≲ 2×10¹² GeV. Whether a non-empty window remains depends on the
  window's lower edge, which this note does not examine.
- This is leading-order and two-level, with BGLS printed-prefactor caveats (`theory_rates.md` §4.5). The project's own
  kinetic model must redo this comparison before anything is written in the note.

---

## 8. Summary table

| Probe | Excluded region near μ ≈ 7×10⁻¹³ eV | Assumptions | Does the window survive? |
|---|---|---|---|
| **WM 2412.03655 Fig. 1** (Cyg X-1 Cnsrv. + GRS 1915+105; driven by GRS) | f ≳ 2.7, 2.3, 1.9, 1.6 ×10¹² GeV at 6, 7, 8, 9 ×10⁻¹³ eV, up to the free field [figure read-off]. GRS alone: 1.5–3×10¹² (Figs. 14, 16) | n ≤ 5, relativistic rates; m = 1 (211, α_GRS = 0.056–0.083); Bayesian 2σ; GRS 12.4 M☉, ã > 0.95, τ = τ_Sal ≈ 4.5×10⁷ yr; no accretion torque | **Partly.** f ∈ [10¹¹, ~2×10¹²] survives; **[~2×10¹², 10¹⁴] excluded** (most of the claimed 10¹²–10¹³) |
| WM, Cyg X-1 Fid. (Figs. 14, 16 left) | f ≳ (0.5–1.0)×10¹² GeV | 21.2 ± 2.2 M☉, ã > 0.9985, τ = 4.8×10⁶ yr; m = 1 (α = 0.095–0.143) | **No** above ~10¹²; only f ≲ 10¹² survives, if this spin is trusted |
| WM, Cyg X-1 Cnsrv. (Fig. 14 dotted) | f ≳ 1.1×10¹³ – 6.4×10¹² GeV | ã = 0.92⁺⁰·⁰⁵₋₀.₀₇, otherwise as above; n ≤ 3 | Yes for f ≲ 6–11×10¹² |
| **BGLS 2011.11646 Fig. 10** (Cyg X-1 solid) | f ≳ 3.8, 2.6, 1.9, 1.4 ×10¹³ GeV | Two-level quasi-equilibrium, "age = spin-depletion time"; intersection over error ranges; Cyg 14.8 ± 1.0 M☉, a* > 0.99, 4.8–7.6 Myr (BLT17 Table I, inferred); m = 1, α ≈ 0.08 | **Yes** (cuts only f ≳ 1.4×10¹³) |
| Hoof+ 2406.10337 Fig. 4 (M33 X-7) | Equilibrium: f ≳ 2.6–1.2×10¹³; bosenova: f ≳ 3.6–1.7×10¹² GeV | 211 only; τ = 3 Myr; 95% HPD | Yes (equilibrium); bosenova variant similar to WM |
| **CFW 2507.21788 Fig. 1** (GW231123) | Only f ≳ (0.8–3)×10¹⁶ GeV; region ends at ≈ 7.0–7.2×10⁻¹³ eV | N_max = 5, WM rates; τ_max = 10⁵–10⁷ yr (hierarchical); NRSur PE | **Yes** |
| Aswathi+ 2507.20979 App. D (GW231123, GW190517) | Interacting: nothing above ~5×10⁻¹³ eV (α ≤ 0.2 cut; lowest f ≈ 2.5×10¹⁴ GeV at 4.5×10⁻¹³). Free field: GW190517 excludes 2.2–11×10⁻¹³ eV | SuperRad; BGLS τ_sd (D8); T_age = 10⁵ yr; m = 1 only for interacting | **Yes** |
| Sun–Brito–Isi 1909.11267 (Cyg X-1 CW, O2) | Free field 6.3–13.2×10⁻¹³ eV (10⁵ yr) or 6.4–8.0×10⁻¹³ eV (5×10⁶ yr); bosenova f ~ 10¹⁵ GeV at 9.6–15.5×10⁻¹³ eV | 14.8 M☉, χ_i = 0.99, unknown χ_f; Yoshino–Kodama bosenova | **Yes** |
| Collaviti+ 2407.04304 Fig. 7 (Cyg X-1 CW reanalysis) | f ≳ 3.8×10¹⁷ – 1.2×10¹⁸ GeV at 0.6–0.9×10⁻¹² eV; f ~ 10¹⁵ withdrawn | BGLS two-level; 14.8 M☉, χ_i = 0.99, 10⁵ yr | **Yes** |
| LVK O3 all-sky 2111.15507 | Only (m_b, M_BH) exclusions for fixed D, χ_i, t_age; no population limit on μ | Small self-interaction regime, F_b ≳ 2×10¹⁸ GeV [Eq. (6)] | **Yes** |
| LVK O4b/O4c 2608.11620 | Free-field scalar 1.39–6.94×10⁻¹³ eV (10⁵ yr), edge at ≈ 7×10⁻¹³; CW searches vector-only | SuperRad, minimally coupled | **Yes** |
| SGWB: Yuan+ 2204.03482; Guo+ 2609.33596 | Yuan+: free field 1.3–17×10⁻¹³ eV. Guo+: LIGO design blind for f ≲ 3×10¹⁷ GeV at any μ | Population models; BGLS two-mode | **Yes** |
| BGLS Sec. VII A (annihilation GWs) | Blind-search events ≲ 10⁻³ for f ≲ 10¹⁶ GeV | Two-level kinetics | **Yes** |
| **LVK 2510.26931** (GW241011) | Free field: scalar 0.3–2.6×10⁻¹² eV (10⁵ yr); ≈ 0.19–4.0×10⁻¹² eV (10⁷ yr) | SuperRad, minimally coupled, all m; m₁ = 19.6 M☉, χ₁ = 0.78 | **Yes** (the window is designed to evade it) |
| **Khalaf+ 2609.12052 Fig. 1** | μ ≈ 0.35–1.3×10⁻¹³ eV, f ≳ 8×10¹⁶ GeV | GW-capture in nuclei; q < 0.3 LLO rate vs GWTC-5; β = 0.5–1.5 | **Yes** (does not reach 7×10⁻¹³ eV) |
| Ning+ 2607.01317 (GWTC-5 population) | Free field 1.7×10⁻¹⁴–3.3×10⁻¹² eV | Self-interaction only at the QCD f_a (~10¹⁹ GeV) | **Yes** (no statement for f ≤ 10¹⁴) |
| Non-superradiance (SM couplings, DM, halos) | None without extra couplings. With SM couplings: f ≲ 10⁸–10⁹ GeV. DM-by-misalignment needs f ≳ 2×10¹³ GeV here | — | **Yes** (if DM: only a subdominant fraction) |

**Overall:** inside f ∈ [10¹¹, 10¹⁴] GeV at μ ∈ [0.6, 0.9]×10⁻¹² eV, the only existing exclusions are spin-based X-ray limits.
- WM (GRS 1915+105) excludes f ≳ 1.6–2.7×10¹² GeV, and f ≳ 0.5–1×10¹² with Cyg X-1 Fid.
- BGLS (Cyg X-1) and Hoof+ (M33 X-7) exclude f ≳ 1–4×10¹³ GeV.
- The claimed window f ≈ 10¹²–10¹³ GeV is therefore **mostly excluded**. A strip f ≈ 10¹²–2×10¹² GeV survives the
  fiducial WM limit but not the Cyg X-1 fiducial-spin limit.

---

## 9. Files read

Existing directories, read as is: 2412.03655, 2011.11646, 2507.21788, 2507.20979, 2608.11620, 2510.26931, 2609.12052,
2609.33596, 2610.11091, 2503.10347, 2607.01317 and 2609.02678 (via `population_papers.md`).

Downloaded today (PDF + e-print + pdftotext; arXiv IDs checked with the export API):
- 1704.05081 (Baryakhtar–Lasenby–Teo, BH table)
- 1909.11267 (Sun–Brito–Isi)
- 2407.04304 (Collaviti+)
- 2111.15507 (LVK O3 all-sky)
- 2204.03482 (Yuan+)
- 2406.10337 (Hoof+)
- 2003.03359 (Zhu+)

Abstract only, not downloaded: 2509.07352, 2204.04523, 2603.14168, 2603.15734, 2602.23415, 2608.18200, 2608.02051,
2512.06790. Not read: 1812.09622 and 1801.01420 (cited via the review table).
