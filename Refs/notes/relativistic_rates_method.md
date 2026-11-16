<!-- 2026-11-16: 相対論的な自己相互作用の率の計算方法の下調べ(値の表の作り直し D18 の手順 3)。[inferred] は原文から読んだのではなく推論した所。実装は research/calc/relrates.py。BH が吸収する過程の式は、n00 の減衰の恒等式と alpha -> 0 の極限(hyd.py)で確かめた(SESSION.md の進み具合)。 -->

# Relativistic (Kerr) self-interaction rates: method memo (scoping, 2026-11-16)

Read: WM 2412.03655 `main.tex` and PDF text (equation numbers are the PDF's); BGLS 2011.11646 App. A, B2–B4, E and Fig. 21; `hyd.py`, `leaver.py`, `fp3.py`; E7, E15, E21. I also ran quick scratch checks (about 1 s each).
Units: G=c=ħ=M=1, μ=α, a=χ. Notation: r± = 1±b with b=√(1−χ²); Ω_H = χ/(2r+); Δ = r²−2r+a²; ρ² = r²+a²cos²θ.

## 1. BH absorption, i×j → k×BH

**Fields.** Φ⁽⁰⁾ = Σ_n √N_n (u_n e^{−iω_n t} + c.c.), with u_n = R_n S_n e^{im_nφ} (WM (35); u_n → ψ_n/√(2μ) in the hydrogenic limit). Normalize each u_n to unit Klein–Gordon charge, with ∮|S e^{imφ}|² = 1:
  2∫dr∮dΩ {[(r²+a²)²ω − 2a r m]/Δ − a²ω sin²θ} |R|²|S|² = 1, so that N = εGM² [inferred].
- This differs from the flat norm ∫r²|R|² by O(α²/n²) per state.
- On a Boyer–Lindquist slice the integrand diverges logarithmically at r+, ∝ (ω−mΩ_H)|R(r+)|²/Δ. This is a slicing artefact: use a horizon-penetrating slice, or cut at r+ + δ and check the result does not depend on δ.

**Source.** WM (34): (□−μ²)Φ = −(λ/6)Φ³. The part at frequency ω_X = Re(ω_i+ω_j−ω_k), with m_X = m_i+m_j−m_k, is
  T_X = −(λ/2) s_ij √(N_iN_jN_k) u_i u_j u_k*,  where s_ij = 1 if i=j and 2 otherwise (WM (36)).
Expand Φ⁽¹⁾ = Σ_ℓ R_ℓ S_{ℓm_X}(θ; c_X²) e^{im_Xφ}, with c_X² = a²(ω_X²−μ²). The angular projection (WM (64)) separates exactly:
  𝒯_ℓ(r) = −(λ/2) s_ij √(N_iN_jN_k) R_iR_jR_k* [r²𝒜_ℓ + a²ℬ_ℓ],
  𝒜_ℓ = ∮S_iS_jS_k*S_ℓ*,  ℬ_ℓ = the same integral with an extra cos²θ.

**Green function.** WM (61): (ΔR′)′ + [K²/Δ + 2amω − (a²ω²+μ²r²+Λ)]R = 𝒯, with K = (r²+a²)ω − am. WM (62)–(63) give R_ℓ → Z_ℓ R^in at the horizon, where
  Z_ℓ = W⁻¹∫_{r+}^∞ R^up 𝒯_ℓ dr,  W = Δ(R^upR^in′ − R^inR^up′),
  R^in → e^{−i(ω_X−m_XΩ_H)r*}  (no 1/Δ for s = 0),
  R^up ~ r^{ν_X−1}e^{−κ_Xr},  κ_X = √(μ²−ω_X²),  ν_X = (2ω_X²−μ²)/κ_X → 1/√e (hyd's ν).

**Flux and rate.**
- Klein–Gordon number flux into the horizon: Ṅ_X = 2(ω_X−m_XΩ_H)(r+²+a²) Σ_ℓ|Z_ℓ|² [inferred]. It reduces to E21's A_H|χ_E(0)|² because r+²+a² = 2r+. The BH gains energy ω_XṄ_X and angular momentum m_XṄ_X.
- Each process occurrence moves one quantum (Ṅ_k = Ṅ_X, as in the E21 derivation), so Γ = Ṅ_X/(N_iN_jN_k) and γ = Γ(GM²)², using λ²(GM²)² = α⁴(M_pl/f)⁴.
- In code: Γ̂ = Γ/λ², κ_eff(α,χ) = Γ̂/(α⁸r̃+), C = κ_eff/κ_hyd. In fp3, replace `kap*cbh` by κ_eff (fp3 already multiplies by α¹¹r̃+).
- WM (58) and BGLS (B30) divide the absorbed energy by μ rather than ω_X. This changes the rate by about α²e/2 (5% for 211×211→322 at α = 0.5).

**211×211→322 in detail.**
- m_X = 0 and ω_X ≈ μ(1−7α²/36) < μ, so the forced mode is bound.
- ℓ_X ∈ {0,2,4} by parity. ℓ_X = 2 and 4 are suppressed by α⁸ and α¹⁶ (BGLS, text after (B27)).
- The flux factor is ω_X > 0, so there is no superradiant factor.
- Spin dependence enters four ways:
  - the factor 2r+, which C divides out;
  - the near-horizon transmission. It is large: the 200 decay rate divided by its leading order (r̃+ included) at α = 0.2 is 2.98, 3.26, 3.65, 4.58, 5.09, 5.91 for χ = 0, 0.5, 0.7, 0.9, 0.95, 0.99 (my Leaver check);
  - the bound-state shapes near r+, through σ_i ∝ ω_i − m_iΩ_H;
  - Re ω_i(χ), through the α⁵χ terms.

All 31 processes in §5 have m_X = 0. The 150 bound combinations with m_X ≠ 0 are zero in hyd. For m_X > 0 at high spin they have ω_X−m_XΩ_H < 0, so the process would run backward [inferred; compare BGLS App. E, rule m′ ≤ 0]. I recommend computing a few of them at α = 0.3–0.5 to bound their size.

**Resonances.** W vanishes at the complex quasi-bound frequencies of the (ℓ_X, m_X) family, which for ℓ_X = 0 are the n00 levels. In the hydrogenic limit 1/W ∝ Γ(1−ν), so hyd's α⁷ formula is the pole limit. The Kerr Green function regulates the pole automatically with the true detuning and the n00 width. At α ≈ 0.3 that width is comparable to the hydrogenic detuning.

## 2. Emission, i×j → k×∞

This uses the same Green function, with ω_X > μ and an outgoing R^up:
  R^up ~ r⁻¹ exp{i[kr + ((2ω²−μ²)/k) ln r]},  k = √(ω²−μ²).
Note that WM (60) writes e^{−ikr*} and has no log phase. Then
  Z^∞_ℓ = W⁻¹∫R^in𝒯_ℓ dr,  Ṅ_∞ = 2kΣ|Z^∞_ℓ|² [inferred].
The same solution also has a small horizon part (it can be negative for m_X > 0); the rate for the kinetics is Ṅ_∞ + Ṅ_H.

It is harder than the BH case:
- many partial waves (m_X up to 7, ℓ_X up to 11; in E21, ℓ = 6 gave 20% of 422×544→322);
- R^in must be split into incoming and outgoing Coulomb waves beyond the source;
- processes change class when ω_X crosses μ (at O(α⁴), only 411×511→322 approaches threshold, near α ≈ 0.5);
- WM Fig. 3 (left) shows ±20% wiggles that need checking.

Do the BH processes first.

## 3. Ingredients and hard spots

**(a) Bound states.** `leaver.py` returns only ω. Extend it to store the backward ratios r_n = a_n/a_{n−1} = −γ_n/(β_n+α_n r_{n+1}). Then a_n = Πr_k gives R(r) through WM (19).
- Scratch test: double precision, N = 6×10⁴, 0.4 s per state, max|a_n| ≈ 1 (no cancellation).
- 211 at α = 0.1 and 0.2, and 322 at α = 0.3 (χ = 0.95), match the hydrogenic peak to 0.5–1.6%. At r = 2r+ the relativistic |R| is 0.84, 0.70 and 0.55 of hydrogenic.
- For α ≲ 0.05, use the series only for r ≲ 10. Beyond that, integrate the ODE inward from about 30n²/α² and match the two pieces; WM §II B report that the two methods agree.
- **Root jumping:** at α = 0.3, `bound_state` asked for 200 returned the 100 root for χ ≤ 0.9. Fix it with a complex initial guess (Im part from BCST 2.32), continuation in α and χ, and node counting.

**(b) Angular functions.** |c²| ≤ χ²α⁴/n² ≤ 0.016, which gives about 3×10⁻⁴ admixture of neighbouring Y_ℓ and a rate change of about 10⁻³. Expanding S in 3–5 Y_ℓ (Gaunt integrals with sympy) and including the a²cos²θ term costs almost nothing, so include both.

**(c) Green-function solutions.**
- R^in: start at r+(1+10⁻⁶) with a Frobenius series and integrate outward (stable direction).
- R^up: start from the asymptotic series and integrate inward (stable direction).
- Use r* near the horizon and ln r far away. Cut the horizon end where Δ/(r²+a²) < 10⁻¹².
- Rescale R^in, since e^{κr} reaches about e^{450}.

**(d) Overlaps.** Ordinary quadrature on the R^up grid.

**Hard spots.**
1. **Resonant and near-threshold processes.** There are 13 resonant processes, plus 411×411→322 (ν = 8.49) and 511×511→422 (ν = 7.56), whose ν_X can cross high n00 poles.
   - My O(α⁵) estimate puts 411×411→322 on the 800 pole at α ≈ 0.24–0.26 for χ = 0.95. This is plausibly the ×150 peak and ×0.2 dip in WM Fig. 3 [inferred].
   - Tabulate the smooth pieces separately: F_ℓ = ∫R^up𝒯_ℓ, and W ≈ W′(ω_p)(ω_X−ω_p) with ω_p from Leaver.
   - Refine the α grid where the detuning crosses zero, and flag points that violate the perturbativity condition, BGLS (A14).
2. **Small α.** The radial range spans 1 to 10⁶. C−1 grows like α² with a large coefficient: the 100 decay rate over its leading order is 1.015, 1.056 and 1.32 at α = 0.01, 0.02 and 0.05 (χ = 0.9). The last value matches BGLS Fig. 21.
3. **Large ℓ_X** in emission.
4. **χ → 1.** The near-horizon factor grows like αr+/b, so use a denser χ grid near 1.

## 4. Validation

**(a) α→0 limit.** Compute C at α = 0.003, 0.005, 0.01 and 0.02 for:
- 211×211→322 (κ_hyd = 4.21e-7);
- 322×322→544 (κ_hyd = 1.93e-9);
- 211×411→422, resonant (κ_res = 3.22e-11);
- 411×411→322 (κ_hyd = 1.52e-11).

Fit 1 + c₁α + c₂α² + c₃α³ and require the intercept to equal 1 to 10⁻³. Expect c₁ ≠ 0 only for the resonant process.

**(b) Earlier estimates.**
- E7's read-off for 211×211→322 at ã = 0.95: expect agreement to about 20%.
- BGLS Fig. 21 bottom (a* = 0.9): about 3.9 at α = 0.13, a peak of about 5.2 near 0.18, about 1.9 at 0.30.
- BGLS Fig. 21 top: the n00 decay ratios.

**(c) WM figures** (graphical only, about ±20% read-off):
- Fig. 3 right (ã = 0.95, α = 0.03–0.37). All 8 curves are in our set: 411×411→422, 211×211→422, 411×411→322, 211×211→322, 211×422→433, 322×411→433, 211×411→322, 211×322→433.
- Fig. 4 bottom: 322×322→544, 322×422→544, 422×422→544, 322×411→533.
- Fig. 3 left and Fig. 4 top: emission.
- Fig. 2: wavefunctions at ã = 0.99.
- Fig. 11: the n44 ratios.

**(d) Internal checks.**
- −2Im ω_n00 = Ṅ_H/N for the n00 quasi-bound states. This fixes both the flux formula and the normalization at every α and χ.
- Optical identity: Ṅ_X = 2Σ Im∫R_ℓ*𝒯_ℓ dr (it follows from WM (61)).
- The Wronskian is constant to 10⁻¹⁰.
- No dependence on the matching radius, r_max or the horizon cutoff.
- Convergence in ℓ_X, in N and in the ODE tolerance.
- Symmetry under i↔j.
- mpmath and double precision agree.
- Results are smooth in (α, χ).

## 5. BH processes among the 8 levels

Rule (hyd.py; same loop as `e21_level422.py` around line 237): e > 0 and nonzero angular integral, which means m_X = 0 and ℓ_X = 0. This gives 31 processes.

Non-resonant, γ = κα¹¹(M_pl/f)⁴r̃+μ:

| process | κ | process | κ |
|---|---|---|---|
| 211×211→322 | 4.21e-7 | 211×211→422 | 1.53e-7 |
| 211×422→533 | 1.15e-7 | 211×322→433 | 9.06e-8 |
| 211×322→533 | 3.08e-8 | 211×511→322 | 2.93e-8 |
| 211×411→322 | 2.46e-8 | 322×411→533 | 2.06e-8 |
| 411×433→544 | 5.21e-9 | 322×511→433 | 2.74e-9 |
| 422×422→544 | 2.27e-9 | 322×322→544 | 1.93e-9 |
| 211×433→544 | 1.12e-9 | 411×422→533 | 9.39e-10 |
| 211×511→422 | 2.03e-10 | 322×422→544 | 1.21e-11 |
| **411×411→322** (ν 8.49) | 1.52e-11 | **511×511→422** (ν 7.56) | 4.95e-13 |

Resonant, γ = κ_res α⁷(M_pl/f)⁴r̃+μ:

| n00 hit | processes (κ_res) |
|---|---|
| 200 | 211×411→422 (3.22e-11), 211×422→433 (7.83e-11), 211×533→544 (4.58e-13) |
| 300 | 322×411→433 (3.81e-11), 322×511→533 (6.26e-12) |
| 400 | 411×411→422 (2.19e-11), 411×422→433 (2.28e-11), 411×533→544 (6.57e-12), 422×511→533 (6.38e-13), 433×511→544 (5.96e-12) |
| 500 | 411×511→422 (2.64e-11), 422×511→433 (4.91e-12), 511×533→544 (2.06e-12) |

At α = 0.1, κ_res α⁷ is comparable to κα¹¹ for the leading processes, so the resonant ones are not negligible. There are also 43 emission candidates (e < 0).

## 6. Cost [estimates]

Per (α, χ) point:
- Shared work: Leaver roots for the 8 levels and n00 with n = 1–9, plus 8 wavefunctions, about 10 s.
- BH: about 2 s per process, so about 60 s for all 31.
- Emission: about 5 s per process, so about 215 s for all 43.

Grid: α = 0.01–0.5 (5 log-spaced points up to 0.05, then steps of 0.01, plus 30% refinement), about 65 values; χ ∈ {0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99}. That is about 520 points.

Totals:
- BH: about 10 CPU-h, under 1 h wall on 14 cores.
- Emission: about 33 CPU-h, about 2.5 h wall.
- Small-α validation: 1–2 CPU-h.

Budget one day including reruns. Measure one worker's RSS before running in parallel.

**Riskiest step:** the resonant and near-threshold processes. This means picking the correct quasi-bound roots (n00 widths are comparable to their spacing at α ≳ 0.25), handling the poles, and staying within perturbation theory. Second risk: O(α²) normalization conventions. The n00 identity and the α→0 limit test both risks.

## 7. Open questions

1. **Normalization of relativistic states.** Neither paper states it. WM normalize ∫|Φ|² = 1 only in the hydrogenic section.
2. **μ vs ω_X.** WM (58) and BGLS (B30) divide the absorbed energy by μ; it should be ω_X.
3. **WM (55).** Its first form, 𝒜_H|ω|²|λR|², is 2π times its second form, 4α²r̃+λ²|R|². BGLS (B28) is ambiguous about the Y₀₀ factor. hyd.py sidesteps this by calibrating on the n00 decay rate.
4. **WM (59)–(60).** The 1/Δ is wrong for s = 0. The exponent sign is wrong, and (60) omits the decaying (ω < μ) branch and the log phase.
5. **ℓ_X.** WM do not say which ℓ_X they keep. BGLS keep only ℓ = 0.
6. **Domain.** WM compute only where all three states are superradiant (text after (64)). We also need C outside that region.
7. **Tables.**
   - 411×411→322: BGLS 9.8e-11, WM 1.6e-11; hyd gives 1.52e-11.
   - 411×411→422: WM 2.3e-7 α⁷; hyd gives 2.19e-11.
   - WM Table I omits r̃_p for 211×322→433.
   - E21's unexplained emission mismatches include 533×533→322.
8. **Unexplained features in WM's figures.** WM explain only the downturn in emission (Fig. 4 caption). Not explained:
   - the 211×322→433 cusp at α ≈ 0.15;
   - the 411×411→322 peak (×150) and dip (×0.2);
   - jumps in the n = 5 curves at α ≈ 0.65;
   - the ±20% oscillations in the emission ratios.
9. **One C_BH for all processes (E15, `cbh_all`) does not hold.** The n00 ratio (4.5 at α = 0.2) is close to the forced-oscillation ratio (5.1, BGLS), but WM Fig. 3 shows C ranging from 0.2 to 150 across processes.
