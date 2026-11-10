# Theory rates for self-interacting scalar superradiance (for an n ≤ 3 / n ≤ 5 kinetic code)

Compiled 2026-11-03 from the **full LaTeX sources and arXiv PDFs** (latest arXiv versions as of today) in
`Refs/arxiv/<ID>/` (`src/` = extracted e-print, `<ID>.pdf`, `<ID>.txt`/`raw.txt` = pdftotext).
Equation numbers are those printed in the arXiv PDF (journal versions may renumber).

Legend:
- **[confirmed]** = read in the full text (LaTeX + PDF number checked).
- **[unconfirmed]** = not stated explicitly in the full text, or only visible in a figure.
- **[my check]** = arithmetic/consistency check done here (scripts not saved); not a statement of any paper.

---

## 0. Sources

| arXiv | Authors (as in metadata) | Short title | Notes |
|---|---|---|---|
| 1908.10370 | Baumann, Chia, **Stout, ter Haar** (not Porto) | The Spectra of Gravitational Atoms ("BCST") | analytic Γ_nlm, ω to O(α⁵) |
| 1804.03208 | Baumann, Chia, Porto | Probing Ultralight Bosons with Binary BHs ("BCP18") | PRD 99, 044001 |
| 2211.03845 | Siemonsen, May, East | SuperRad | PRD 107, 104003; code: `www.bitbucket.org/weast/superrad` (python) |
| 0705.2880 | Dolan | Instability of massive KG field on Kerr | numerical (continued fraction) rates |
| 1411.0686 | Brito, Cardoso, Pani | BHs as particle detectors ("BCP14") | Schwarzschild GW flux |
| 1312.2326 | Yoshino, Kodama | GW from an axion cloud: superradiant phase ("YK") | flat-space GW formula (downloaded in addition) |
| 1411.2263 | Arvanitaki, Baryakhtar, Huang | Discovering the QCD axion with BHs and GWs ("ABH15") | PRD 91, 084011 |
| 2011.11646 | Baryakhtar, Galanis, Lasenby, Simon | BH superradiance of self-interacting scalars ("BGLS") | PRD 103, 095019 |
| 2412.03655 | Witte, Mummery | Stepping up superradiance constraints on axions ("WM") | v2 2025-04-07 |
| 2507.21788 | Caputo, Franciolini, Witte | Superradiance constraints from GW231123 ("CFW") | updated 2026-03-27 |
| 2607.01317 | Ning, Safdi, Welch | No evidence for SR axions in GWTC-5 ("NSW") | 2026-07-01 |

---

## 1. Conventions per paper

| Paper | Units | r_g | α | spin | horizon | Ω_H | "Γ" means |
|---|---|---|---|---|---|---|---|
| BCST 1908.10370 | ħ=c=1; often μ=1 | $r_g\equiv GM/c^2$ (eq. 1.1 text) | $\alpha=GM\mu/\hbar c$ (1.1) | $\tilde a=a/M$ | $r_\pm=M\pm\sqrt{M^2-a^2}$, $\tilde r_\pm=r_\pm/M$ (text below 2.3); in Sec. 3 also $\tilde r_\pm=r_\pm/\alpha$ in μ=1 units (same number) | $\Omega_H=a/(2Mr_+)$ | $\omega=E+i\Gamma$ (2.27): **Γ = Im ω (field amplitude)** |
| BCP18 1804.03208 | ħ=c=1 | $GM$ | $GM\mu$ | $\tilde a$ | $\tilde r_+=r_+/M$ | — | Γ = Im ω (2.11) |
| BCP14 1411.0686 | G=c=1 | — | $M\mu$ | $a/M$ | $r_\pm=M\pm\sqrt{M^2-a^2}$ | — | $\omega_I$ = Im ω; $\dot E_S=2M_S\omega_I$ (17) |
| ABH15 1411.2263 | ħ=c=1 | $r_g\equiv G_NM$ (2) | $r_g\mu_a$ (2) | $a_*=a/r_g$ | — | $\omega^+=\frac{1}{2}\frac{a_*}{1+\sqrt{1-a_*^2}}r_g^{-1}$ (5) | $\dot N=\Gamma_{\rm sr}N$ (7): N-growth |
| BGLS 2011.11646 | ħ=c=1, signature −+++, **$M_{\rm pl}\equiv 1/\sqrt G$** (=1.22×10¹⁹ GeV) | $r_g\equiv GM$ | $GM\mu$ | $a_*=J/GM^2$ | $\tilde r_+\equiv r_+/r_g=1+\sqrt{1-a_*^2}$ | $\Omega_H=\frac{1}{2r_g}\frac{a_*}{1+\sqrt{1-a_*^2}}$ | $\Gamma_{\rm SR}\equiv 2\omega_i$ (text below eq. 2): **N-growth** |
| WM 2412.03655 | ħ=c=1, lengths in GM | GM | GM μ | $\tilde a=a/GM$ | $r_\pm=GM(1\pm\sqrt{1-\tilde a^2})$, $\tilde r_\pm=r_\pm/GM$ | $\Omega_H=\tilde a/(2r_+)$ | $\Gamma_{nlm}$ = Im ω (16); $\gamma^{\rm SR}$ in ε-eqs is N-growth (2Γ) **[unconfirmed: not stated in WM text; stated in CFW eq. (1)]** |
| CFW 2507.21788 | natural | — | — | $\tilde a=a/GM$ | $r_+=GM(1+\sqrt{1-\tilde a^2})$ | $\tilde a/(2r_+)$ | $\omega=\omega_r+i\Gamma$; $\tau_{211}\equiv(2\Gamma_{211})^{-1}$; $\dot\epsilon=2\Gamma\epsilon+\dots$ (1) |
| NSW 2607.01317 | — | — | $GMm_a$ (1) | χ | $r_+=M(1+\sqrt{1-\chi^2})$ | $\chi/[2M(1+\sqrt{1-\chi^2})]$ | $\Gamma_{\rm SR}$ = "rate of cloud growth" (4) |
| SuperRad 2211.03845 | G=c=1 | — | GMμ | $a_*$ | $r_+$ | — | $\omega_I$, $\tau_I=1/\omega_I$ (12); cloud mass grows as $e^{2\omega_I t}$ (36) |

Quantum-number conventions: BCST/BGLS/WM/CFW use principal $n\ge\ell+1$ (211 = n=2). BCP14 and NSW use radial node number (211 is "n=0"; NSW: $n_h=n+\ell+1$). SuperRad: $n_S=\ell_S+1+\hat n_S$ (10).

Self-coupling sign: BGLS $\mathcal L=-\frac12(D\varphi)^2-\frac12\mu^2\varphi^2+\frac{\lambda}{4!}\varphi^4$, EoM $(D^2-\mu^2)\varphi=-\frac{\lambda}{6}\varphi^3$ (9); **λ>0 attractive**, $\lambda\equiv\mu^2/f^2$; matches $V=\mu^2f^2(1-\cos\varphi/f)=\frac12\mu^2\varphi^2-\frac{\mu^2}{24f^2}\varphi^4+\dots$. So the user's "λ = −μ²/f²" is the same physics in the $V\supset+\lambda\varphi^4/4!$ convention. WM writes $\mathcal L\supset-\frac{\lambda}{4!}\Phi^4$ with $\lambda=\mu^2/f_a^2$ (text above eq. 34) but uses the same EoM (34) $(\Box-\mu^2)\Phi=-\frac{\lambda}{6}\Phi^3$ as BGLS → the Lagrangian sign in WM is a typo **[my check]**. CFW: $V\simeq\frac12\mu^2\Phi^2-\frac{1}{4!}\lambda\Phi^4$, $\lambda=(\mu/f_\Phi)^2$ (text). Rates ∝ λ² are sign-blind; frequency shifts are not.

---

## 2. Free-field superradiance

### 2.1 Detweiler-type rate (BCST 1908.10370) [confirmed]

$$\Gamma_{n\ell m}=2\tilde r_+\,C_{n\ell}\,g_{\ell m}(\tilde a,\alpha,\omega)\,(m\Omega_H-\omega_{n\ell m})\,\alpha^{4\ell+5}\qquad\text{(2.32)}$$
$$C_{n\ell}\equiv\frac{2^{4\ell+1}(n+\ell)!}{n^{2\ell+4}(n-\ell-1)!}\left[\frac{\ell!}{(2\ell)!(2\ell+1)!}\right]^2\qquad\text{(2.33)}$$
$$g_{\ell m}(\tilde a,\alpha,\omega)\equiv\prod_{k=1}^{\ell}\left(k^2(1-\tilde a^2)+(\tilde a m-2r_+\omega)^2\right)\qquad\text{(2.34)}$$

- Here Γ = Im ω (2.27). Dimensionful units restored: $r_+=M(1+\sqrt{1-\tilde a^2})$ is dimensionful in (2.34), so $2r_+\omega=2\tilde r_+\alpha\,(\omega/\mu)$; $\tilde r_+=r_+/M$; $\Omega_H=a/(2Mr_+)$ (text below 2.3).
- Equivalent: $2\tilde r_+(m\Omega_H-\omega)=(m\tilde a-2r_+\omega)/M$, so $\mathrm{Im}\,\omega=C_{n\ell}\,g_{\ell m}\,\alpha^{4\ell+5}\,(m\tilde a-2r_+\omega)/M$.
- Derivation from the matching result (3.21) $\nu_0=n+iP_+[2\alpha^2(\tilde r_+-\tilde r_-)]^{2\ell+1}\frac{(n+\ell)!}{n^{2\ell+1}(n-\ell-1)!}[\dots]^2\prod_k(k^2+4P_+^2)$ with $P_+=\alpha(\tilde am-2r_+\omega)/(r_+-r_-)$ (2.13) and $\omega=\mu\sqrt{1-\alpha^2/\nu^2}$ (2.27) reproduces (2.32)–(2.34) exactly **[my check]**.
- Same formula: BCP18 (2.12)–(2.13) [with $\omega\to\mu$ inside g: $(\tilde am-2\tilde r_+\alpha)^2$], WM (16)–(18), SuperRad (13)+(B4) [SuperRad (13) is printed with $(\omega_R-m\Omega_H)$, sign typo relative to (2.32) **[my check]**].
- BCP18 footnote to (2.13): Refs. Pani et al. 2012 (PRD 86, 104017) and Rosa & Dolan 2012 "point out a missing factor of 1/2 in Detweiler's original computation". Dolan (3) quotes Detweiler's original $\tau\approx24(a/M)^{-1}(\mu M)^{-9}GM/c^3$, i.e. Im ω = (a/M)μ(Mμ)⁸/24 — **factor 2 too large**.

### 2.2 Small-α limit for 211 — exact statement

n=2, ℓ=m=1: $C_{21}=1/48$; $g_{11}=(1-\tilde a^2)+(\tilde a-2\tilde r_+\alpha\,\omega/\mu)^2=1-4\tilde a\tilde r_+\alpha+4\tilde r_+^2\alpha^2+O(\alpha^3)$. Hence

$$\boxed{\mathrm{Im}\,\omega_{211}=\frac{\mu\,\alpha^8}{48}\,(\tilde a-2\alpha\tilde r_+)\,g_{11}\;\simeq\;\frac{\mu\alpha^8}{48}(\tilde a-2\alpha\tilde r_+)},\qquad \boxed{\Gamma^{\rm SR}_{211}\equiv\frac{\dot N}{N}=2\,\mathrm{Im}\,\omega_{211}\simeq\frac{\mu\alpha^8}{24}(\tilde a-2\alpha\tilde r_+)}$$

with $\tilde r_+=1+\sqrt{1-\tilde a^2}$ (r_+ in units of GM, r_g = GM not 2GM).

- So "Γ_211 ≈ μα⁸(χ − 2αr̃₊)/24" is **correct as the occupation-number (cloud-mass) growth rate**, and **a factor 2 too large as Im ω**.
- Confirmations: BGLS text below (2): "211 … has $\omega_i=\frac{a_*}{48}\alpha^8\mu$ … $\Gamma_{\rm SR}\equiv2\omega_i$"; BGLS Table I: $\Gamma^{\rm SR}_{211}\simeq4\times10^{-2}\alpha^8(a_*-2\alpha(1+\sqrt{1-a_*^2}))\mu$ (1/24 = 0.0417); BCP14 (4) γ₁ = 48 and text below (17) "$M\omega_I=\frac{1}{48}(a/M-2\mu r_+)(M\mu)^9$"; CFW (1) $\dot\epsilon=2\Gamma\epsilon$.
- Strict leading order at fixed ã is $\tilde a\mu\alpha^8/48$; the $-2\alpha\tilde r_+$ term is kept to put the zero at the right threshold (BGLS comment after eq. 2). Note $g_{11}$ is **not** ≈1 at moderate α: g₁₁ = 0.57 at α = 0.1, ã = 0.9, and at threshold ($\tilde am=2r_+\omega$) $g_{\ell m}\to(\ell!)^2(1-\tilde a^2)^\ell$ **[my check]**.

### 2.3 Level energies (needed for transition/emission frequencies) [confirmed]

BCST (2.28)–(2.29):
$$E_{n\ell m}=\mu\left(1-\frac{\alpha^2}{2n^2}-\frac{\alpha^4}{8n^4}+\frac{f_{n\ell}}{n^3}\alpha^4+\frac{h_\ell}{n^3}\tilde a m\alpha^5+\cdots\right),\quad \mathrm{Re}\,\nu_{n\ell m}=n+f_{n\ell}\alpha^2+h_\ell\tilde am\alpha^3$$
$$f_{n\ell}\equiv-\frac{6}{2\ell+1}+\frac{2}{n},\qquad h_\ell\equiv\frac{16}{2\ell(2\ell+1)(2\ell+2)}$$
- ℓ=0: $E_{n00}=\mu(1-\frac{\alpha^2}{2n^2}-\frac{\alpha^4}{8n^4}+\frac{(2-6n)\alpha^4}{n^4})$ (3.25); no hyperfine term.
- Equivalent forms: BCP18 (2.10) $\Delta\omega=\mu(-\frac{\alpha^4}{8n^4}+\frac{(2\ell-3n+1)\alpha^4}{n^4(\ell+1/2)}+\frac{2\tilde am\alpha^5}{n^3\ell(\ell+1/2)(\ell+1)})$ (identical, checked); BGLS (App. "Perturbations from BH companion", unnumbered-in-text eq. `baumannenergies`); SuperRad (9),(B2),(B3); WM (15) (no α⁵ term).
- BCST: (2.28) "is an excellent approximation for α ≲ 0.2" (Fig. 6 caption). SuperRad: analytic ω_R error ∝ α⁶ at small α and ≲10⁻² in the relativistic regime (Sec. 3.2, Fig. "wr_model").

Explicit coefficients **[my check]** ($E/\mu=1+c_2\alpha^2+c_4\alpha^4+c_5\tilde a\alpha^5$):

| level | c₂ | c₄ | c₅ |
|---|---|---|---|
| 211 | −1/8 | −17/128 | 1/12 |
| 311 | −1/18 | −11/216 | 2/81 |
| 322 | −1/18 | −23/1080 | 4/405 |
| 411 | −1/32 | −49/2048 | 1/96 |
| 422 | −1/32 | −117/10240 | 1/240 |
| 433 | −1/32 | −87/14336 | 1/448 |
| 544 | −1/50 | −7/3000 | 4/5625 |

$\omega_{322}-\omega_{211}=\mu[\frac{5}{72}\alpha^2+\frac{1927}{17280}\alpha^4-\frac{119}{1620}\tilde a\alpha^5]$; $2\omega_{322}-\omega_{211}-\mu=+\frac{\alpha^2}{72}\mu$ (>0 → emission to ∞, BGLS Sec. III B); $2\omega_{211}-\omega_{322}-\mu=-\frac{7\alpha^2}{36}\mu$ (<0 → bound, BH-absorbed forced oscillation, BGLS Sec. III C).

Cloud self-energy shifts (time-dependent frequency):
- Self-interaction (attractive): BGLS (B15)
$\Delta\omega_{211}\simeq-\alpha^5\mu(M_{\rm pl}/f)^2(1.2\times10^{-4}\varepsilon_{211}+3.5\times10^{-5}\varepsilon_{322})$, $\Delta\omega_{322}\simeq-\alpha^5\mu(M_{\rm pl}/f)^2(3.5\times10^{-5}\varepsilon_{211}+1.4\times10^{-5}\varepsilon_{322})$; WM (40)–(43) give the same 211/322 coefficients plus 433, 544 entries.
- Self-gravity: BGLS (G11)–(G12) $\Delta\omega_{211}\simeq-\frac{\alpha^3\mu}{GM^2}(0.19N_{211}+0.11N_{322})$, $\Delta\omega_{322}\simeq-\frac{\alpha^3\mu}{GM^2}(0.11N_{211}+0.09N_{322})$; combined in (H1a,b).

### 2.4 Accuracy vs numerical relativistic rates

- Dolan 0705.2880 [confirmed]: 211 maximum growth $M\,\mathrm{Im}\,\omega\approx1.5\times10^{-7}$ at $M\mu\approx0.42$, $a=0.99$; Table I maxima: (a, Mμ, M Im ω) = (0.7, 0.187, 3.33e-10), (0.8, 0.231, 2.16e-9), (0.9, 0.293, 1.55e-8), (0.95, 0.343, 4.88e-8), (0.98, 0.393, 1.11e-7), (0.99, 0.421, 1.50e-7).
- **[my check]** BCST (2.32) at those points (with ω = μ(1−α²/2n²)) gives 1.85e-10, 9.3e-10, 4.4e-9, 8.8e-9, 1.1e-8, 9.8e-9: the analytic formula **underestimates the numerical peak rate by ×1.8 (α≈0.19) up to ×15 (α≈0.42, a=0.99)**.
- SuperRad (Sec. 3.3, Fig. "wi_model"): relative error of analytic ω_I is "∼O(1)" in the relativistic regime; error ∝ α at small α. SuperRad fits $\omega_I$ with ansatz (15): $\frac{\omega_IM\alpha^{-\kappa}G^{-1}}{2r_+(\omega_R-m\Omega_H)}-1=\sum_{p,q}\alpha^p[\hat b_{pq}a_*^{q+1}+\hat c_{pq}(1-a_*^2)^{q/2}]$ (κ_S = 4m+5, (14)); coefficients are inside the code, not printed **[unconfirmed values]**. Only m = 1, 2 scalar modes are tabulated in SuperRad.
- BCP18 below (2.13): "reasonable agreement with numerical results for α < 0.5" (citing Brito–Cardoso–Pani review Fig. 28) — loose; my Dolan comparison shows O(1–10) differences near the peak.
- WM Sec. II B: non-relativistic wavefunctions deviate from relativistic ones "by many orders of magnitude" near the horizon for α ≳ 0.2 (Fig. 2); WM compute all Γ_nlm with Leaver continued fractions (N ~ 10³–10⁴ terms, 256-bit BigFloat).
- Practical: use numerical Im ω (SuperRad for 211/322; Leaver/Dolan continued fraction or WM code for 311, 411, 422, 433, 5ℓm) for α ≳ 0.1; keep (2.32) only for α ≲ 0.1 and as a check.

---

## 3. GW emission from free clouds (flag: most analytic results are flat/Schwarzschild approximations)

Dimensionless convention: $P_{\rm GW}=\tilde P\,(M_c/M)^2\,c^5/G$, $M_c=\mu N$.

### 3.1 Annihilation (frequency $2\omega_R$)

| Source | Formula | Approximation |
|---|---|---|
| BCP14 (13), (14) | $\dot E_{\rm GW}=\frac{484+9\pi^2}{23040}\frac{M_S^2}{M^2}(M\mu)^{14}=0.0249\,(M_S/M)^2\alpha^{14}$; $\dot J_{\rm GW}=\dot E_{\rm GW}/\omega_R$ | **Schwarzschild** background, 211 (n_r=0, ℓ=m=1); stated to be an **upper bound** |
| YK (57a–c), Table 1 | $\dot E_{\rm GW}=C_{n\ell}(E_a/M)^2(\mu M)^{Q_\ell}$, $Q_\ell=4\ell+10$, $C_{n\ell}=\frac{16^{\ell+1}\ell(2\ell-1)\Gamma(2\ell-1)^2\Gamma(\ell+n+1)^2}{n^{4\ell+8}(\ell+1)\Gamma(\ell+1)^4\Gamma(4\ell+3)\Gamma(n-\ell)^2}$; C₂₁ = 1/640 = 1.56e-3, C₃₂ = 1.89e-7, C₄₃ = 2.89e-11, C₅₄ = 3.17e-15 | **flat**; YK Sec. 3.4: gauge-dependent at leading order, "C_nℓ … not reliable and may be changed by a factor"; radiation-gauge result differs by $(\ell-1)^2/4$ |
| SuperRad (31)–(34),(E9) | $\tilde P=H\alpha^{\eta}$, $\eta_S=4m+10$; $H^{\rm Schw}=\frac{484+9\pi^2}{23040}$, $H^{\rm flat}=\frac1{640}$ for m=1; numerical Teukolsky values lie in between | flat underestimates, Schw overestimates (~10× spread) |
| BCP18 (2.15) | $dE_{\rm gw}/dt\simeq0.01(M_c/M)^2\alpha^{14}$ (fit to Yoshino–Kodama/Brito et al. numerics) | fit; good for α ≲ 0.1, upper bound beyond (footnote) |
| BGLS Table IV / Table I | $\gamma^{\rm GW}_{211\times211}\simeq10^{-2}\alpha^{14}\mu$ entering $\dot\varepsilon_{211}\supset-2\gamma\varepsilon_{211}^2$ ⇒ $\tilde P=2\times10^{-2}\alpha^{14}$ **[my check]**; $\gamma^{\rm GW}_{322\times322}\simeq3\times10^{-8}\alpha^{18}\mu$ ⇒ $\tilde P=6\times10^{-8}\alpha^{18}$ | leading order; origin of coefficients not stated **[unconfirmed]** |
| ABH15 Table VI, (A12) | $dP/d\Omega$ per $N^2$ for 2p, 3d, 4f; $\Gamma_a\equiv\int dP/d\Omega/(2\omega N^2)$ | flat, hydrogenic; authors: 2p leading term cancels → "lower bound" |

**[my check]** on ABH15 Table VI: 3d integrates to $\tilde P=4.72\times10^{-8}\alpha^{18}$ = exactly YK's $C_{32}\times(\ell-1)^2/4$ (radiation gauge); 2p ∝ $(\ell-1)^2=0$ at leading order (their "cancellation"); 4f entry is $2^{20}$ larger than YK's $C_{43}$ (and YK×(ℓ−1)²/4 = C₄₃ for ℓ=3) → typo in ABH15 Table VI. BGLS's 322 coefficient (3e-8 → $\tilde P$ = 6e-8) matches ABH15's, i.e. ~3× below YK flat.

Recommendation: 211/322 annihilation from SuperRad (numerical Teukolsky, Kerr). Cross-annihilation 211×322→GW: $P\propto\alpha^{16+2(\ell+\ell')}$ per BGLS Sec. VII A text, "numerical GW power calculations … have not yet been performed" **[note: BGLS's α-power statements in Sec. II and VII A are not consistent with their Table I scaling (α¹² at fixed N, μ for 211); trust Table I / YK $Q_\ell$]**.

Single-level dynamics (ABH15 (23)): $\dot N=\Gamma_{\rm sr}N-\Gamma_aN^2$; after saturation $N(t)=N_{\max}/(1+\Gamma_aN_{\max}t)$ (ABH15 eq. annihtime). In BGLS normalization $\dot N_{211}\supset-2\Gamma^{\rm GW}_{211\times211}N_{211}^2$ (31) and $\tau_{\rm ann}\approx[2\Gamma^{\rm GW}_{211\times211}N_{211,\max}]^{-1}\simeq3\times10^9\,{\rm yr}\,(M/10M_\odot)(0.05/\alpha)^{15}$ for α ≲ 0.2 (39). Beware ABH15's Γ_a ↔ BGLS's Γ^GW differ by the multiplicity convention.

Strain (ABH15 (12)): $h=(4G P/r^2\omega^2)^{1/2}$, ω the GW angular frequency.

### 3.2 Transitions 322 → 211 + GW [confirmed, flat-space hydrogenic]

BGLS (99)–(101) (total from Arvanitaki–Dubovsky 2010):
$$\frac{dP}{d\Omega}=\frac{GN_{322}N_{211}}{\pi r_g^4}\alpha^{14}\left(\frac{2^5}{3^65^8}(1-\cos^4\theta)+\frac{(27+28\cos2\theta+9\cos4\theta)\sin^2\theta}{2^23^65^{10}7^2}\right)\quad(99)$$
$$P=\frac{2^8\times5717}{3^55^{11}7^3}\frac{GN_{322}N_{211}}{r_g^4}\alpha^{14}\quad(100),\qquad \omega=\omega_{322}-\omega_{211}\simeq\frac{5}{72}\alpha^2\mu$$
$$\dot\varepsilon_{322}\simeq-5\times10^{-6}\alpha^{10}\varepsilon_{211}\varepsilon_{322}\quad(101)$$
**[my check]**: (100)/ω gives 5.18e-6 — consistent with (101) and Table IV; BGLS text before (50) uses 3e-6 (inconsistent). ABH15 two-level transition dynamics (15), quadrupole estimate (16) $\Gamma_t\sim\frac{2G\omega^5}{5}\mu^2r_c^4$, Table VII (Δn=1 transitions with Δℓ=0, e.g. 6g→5g; flat space).

---

## 4. Self-interaction kinetics — BGLS 2011.11646

### 4.1 Building blocks [confirmed]
Non-relativistic reduction $\varphi=\frac{1}{\sqrt{2\mu}}(\psi e^{-i\mu t}+{\rm c.c.})$, $\int dV|\psi|^2\simeq N$ (10); GP equation $(i\partial_t+\frac{\nabla^2}{2\mu}+\frac{\alpha}{r})\psi\simeq\frac{-3}{24\mu^2}\lambda\psi^2\psi^*$ (12). Parametric treatment (20) and the second-order frequency shift (A11):
$$\delta\hat\omega=\tilde\lambda\!\int\!\Psi_i^*|\Psi_{211}|^2\Psi_i-\tilde\lambda^2\sum_k\frac{|\int\Psi_k^*\Psi_{211}^2\Psi_i^*|^2}{\hat\omega_i+\hat\omega_k^*}+\tilde\lambda^2\sum_k\frac{|\int\Psi_k^*|\Psi_{211}|^2\Psi_i^*|^2}{\hat\omega_i-\hat\omega_k},\quad\tilde\lambda\equiv-\frac{3\lambda}{24\mu^2}$$
(mixing with decaying $n00$ modes gives Im δω̂ > 0 → growth of 322). Damping rate definition: $\Gamma_{\rm damp}=P_{\rm abs}/(\mu N_{211}^2N_{322})$ (App. B3).

### 4.2 (a) Leading processes, powers and rates

Power formulas (leading order in α):
- $211\times211\times211\to\infty$ (relativistic 3→1): $P\simeq1.5\times10^{-8}\alpha^{17}\mu^2\lambda^2N_{211}^3$ (15).
- $322\times322\to211\times\infty$ (NR emission into ℓ=3,5; m=3): $P\simeq10^{-8}\alpha^4\lambda^2\mu^2N_{322}^2N_{211}$ (18) (classical; quantum: $N_{211}\to N_{211}+1$; a factor 4 below Gruzinov 2016 because of a √2 normalization error there).
- $211\times211\to322\times{\rm BH}$ (m=0 forced oscillation absorbed): $P\simeq4\times10^{-7}\alpha^7\lambda^2(1+\sqrt{1-a_*^2})\mu^2N_{211}^2N_{322}$ (21).
- Resonant $211\times311\to322\times{\rm BH}$ (200 forced oscillation nearly on shell): $P\simeq3\times10^{-10}\alpha^3\lambda^2(1+\sqrt{1-a_*^2})\mu^2N_{211}N_{311}N_{322}$ (23).
- Cubic $g\varphi^3/3!$, $g=C\mu^2/f$: $211\times211\to\infty$, $P\simeq10^{-4}\alpha^{12}C^2(\mu^4/f^2)N_{211}^2$ (28); effective quartic $1/f^2\to(1+\frac53C^2)/f^2$ for NR processes (29).

Rate constants (event rate = P/μ for BH/NR-∞ processes; Table I; ε ≡ N/(GM²), λ = μ²/f², $\tilde r_+=1+\sqrt{1-a_*^2}$):

| process | N-units | ε-units (γ) |
|---|---|---|
| 211 SR | $\Gamma^{\rm SR}_{211}\simeq4\times10^{-2}\alpha^8(a_*-2\alpha\tilde r_+)\mu$ | same |
| 211×211→GW | $10^{-2}\alpha^{12}(\mu/M_{\rm pl})^2\mu$ | $10^{-2}\alpha^{14}\mu$ |
| 211×211→322×BH | $4\times10^{-7}\alpha^7\lambda^2\tilde r_+\mu$ | $4\times10^{-7}\alpha^{11}(M_{\rm pl}/f)^4\tilde r_+\mu$ (Table V: 4.3e-7) |
| 322×322→211×∞ | $10^{-8}\alpha^4\lambda^2\mu$ | $10^{-8}\alpha^8(M_{\rm pl}/f)^4\mu$ (Table VI: 1.1e-8) |
| 322 SR | $8\times10^{-5}\alpha^{12}(a_*-\alpha\tilde r_+)\mu$ | (Table I prints "$\gamma^{\rm SR}_{322}=\Gamma^{\rm SR}_{211}$" — typo for $\Gamma^{\rm SR}_{322}$) |

Further tables: Table IV (gravitational: SR 211/322/433/544, GW ann 211/322, GW tr 322→211), Table V (13 BH-absorbed processes, e.g. $\Gamma^{422\times BH}_{211\times211}=1.5\times10^{-7}$, $\Gamma^{433\times BH}_{211\times322}=9.1\times10^{-8}$, $\Gamma^{544\times BH}_{322\times322}=1.9\times10^{-9}$, all ×$\alpha^{11}(M_{\rm pl}/f)^4\tilde r_+\mu$; resonant $\Gamma^{433\times200}_{211\times422}=1.1\times10^{-9}\alpha^7$), Table VI (21 NR-emission processes, ×$\alpha^8(M_{\rm pl}/f)^4\mu$), Table VII (relativistic: $\Gamma^{3\to1}_{211}=5\times10^{-9}\alpha^{21}(M_{\rm pl}/f)^4$, $\Gamma^{2\to1}_{211}({\rm cubic})=1.9\times10^{-4}\alpha^{14}|C|^2(M_{\rm pl}/f)^2$).

Scalings of γ (Sec. IV A list): SR ∝ α^{4ℓ+4}; NR emission to ∞ ∝ α⁸(M_pl/f)⁴; BH-damped forced oscillation of angular momentum ℓ ∝ α^{11+4ℓ}(M_pl/f)⁴ (except resonant); 3→1 ∝ α^{2(ℓ+ℓ'+ℓ'')+15}; GW annihilation ∝ α^{10+2(ℓ+ℓ')}.

**[my check]** BGLS Table IV SR coefficients for ℓ ≥ 2 (8e-5, 2e-8, 2e-12) equal $2mC_{n\ell}\,g_{\ell m}(\tilde a{=}0)$; the leading-order spin dependence of $g_{\ell m}$ ($(\tilde am)^2$ terms) is dropped, so at $a_*=0.9$ they are low by ×3.4 (322), ×15 (433), ×74 (544) relative to BCST (2.32). 211's 4e-2 is exact at leading order. Use (2.32) or numerics instead.

Finite-α: "rate a factor few larger than the leading-order value at α ∼ 0.2" for 211×211→322×BH, tracking the $n00$ horizon-flux enhancement (Fig. 21); BGLS's plots use these semi-analytic/numerical rates. The correction function is only shown graphically **[unconfirmed numbers]**.

### 4.3 (b) Coupled evolution equations [confirmed]

General (30), (32): $\dot\varepsilon_j=\gamma^{\rm SR}_j\varepsilon_j+\sum_{j'}(-c\gamma^{\rm GW}_{j\times j'}+\gamma^{\rm GW}_{j'\to j}-\gamma^{\rm GW}_{j\to j'})\varepsilon_j\varepsilon_{j'}+\sum_{j',j''}(\gamma^{j\times k}_{j'\times j''}-c\gamma^{j''\times k}_{j\times j'}-c\gamma^{j'\times k}_{j\times j''}-c\gamma^\infty_{j\times j'\times j''})\varepsilon_j\varepsilon_{j'}\varepsilon_{j''}$ ("cΓ encodes multiplicity factors"). N-form for 211 (31):
$$\dot N_{211}=\Gamma^{\rm SR}_{211}N_{211}-2\Gamma^{\rm GW}_{211\times211}N_{211}^2-\Gamma^{\rm GW}_{211\to322}N_{211}N_{322}-2\Gamma^{322\times{\rm BH}}_{211\times211}N_{211}^2N_{322}+\Gamma^{211\times\infty}_{322\times322}N_{211}N_{322}^2-3\Gamma^\infty_{(211)^3}N_{211}^3-2\Gamma^\infty_{(211)^2\times322}N_{211}^2N_{322}+\dots$$

Two-level system, leading order in α (34)–(37) ($\gamma=\kappa\times$[leading α, f, a* dependence]$\times\mu$):
$$\begin{aligned}\frac{\dot\varepsilon_{211}}{\mu}={}&\kappa^{\rm SR}_{211}\alpha^8(a_*-2\alpha\tilde r_+)\varepsilon_{211}-2\kappa^{322\times{\rm BH}}_{211\times211}\alpha^{11}(M_{\rm pl}/f)^4\tilde r_+\varepsilon_{211}^2\varepsilon_{322}+\kappa^{211\times\infty}_{322\times322}\alpha^8(M_{\rm pl}/f)^4\varepsilon_{322}^2\varepsilon_{211}-2\kappa^{\rm GW}_{211\times211}\alpha^{14}\varepsilon_{211}^2\\&+\Big(-\kappa^{\rm GW}_{211\times322}\alpha^{16}\varepsilon_{211}\varepsilon_{322}+\kappa^{\rm GW}_{322\to211}\alpha^{10}\varepsilon_{211}\varepsilon_{322}-3\kappa^\infty_{(211)^3}\alpha^{21}(M_{\rm pl}/f)^4\varepsilon_{211}^3\\&\quad-2\kappa^\infty_{(211)^2\times322}\alpha^{23}(M_{\rm pl}/f)^4\varepsilon_{211}^2\varepsilon_{322}-\kappa^\infty_{211\times(322)^2}\alpha^{25}(M_{\rm pl}/f)^4\varepsilon_{211}\varepsilon_{322}^2\Big)\end{aligned}\quad(34)$$
$$\begin{aligned}\frac{\dot\varepsilon_{322}}{\mu}={}&\kappa^{\rm SR}_{322}\alpha^{12}(a_*-\alpha\tilde r_+)\varepsilon_{322}+\kappa^{322\times{\rm BH}}_{211\times211}\alpha^{11}(M_{\rm pl}/f)^4\tilde r_+\varepsilon_{211}^2\varepsilon_{322}-2\kappa^{211\times\infty}_{322\times322}\alpha^8(M_{\rm pl}/f)^4\varepsilon_{322}^2\varepsilon_{211}\\&+\Big(-2\kappa^{\rm GW}_{322\times322}\alpha^{18}\varepsilon_{322}^2-\kappa^{\rm GW}_{211\times322}\alpha^{16}\varepsilon_{211}\varepsilon_{322}-\kappa^{\rm GW}_{322\to211}\alpha^{10}\varepsilon_{211}\varepsilon_{322}\\&\quad-3\kappa^\infty_{(322)^3}\alpha^{27}(M_{\rm pl}/f)^4\varepsilon_{322}^3-\kappa^\infty_{(211)^2\times322}\alpha^{23}(M_{\rm pl}/f)^4\varepsilon_{211}^2\varepsilon_{322}-2\kappa^\infty_{211\times(322)^2}\alpha^{25}(M_{\rm pl}/f)^4\varepsilon_{211}\varepsilon_{322}^2\Big)\end{aligned}\quad(35)$$
$$\frac{\dot a_*}{\mu}=-\kappa^{\rm SR}_{211}\alpha^8(a_*-2\alpha\tilde r_+)\varepsilon_{211}-2\kappa^{\rm SR}_{322}\alpha^{12}(a_*-\alpha\tilde r_+)\varepsilon_{322}\quad(36)$$
$$\frac{\dot M}{\mu^2GM^2}\simeq-\kappa^{\rm SR}_{211}\alpha^8(a_*-2\alpha\tilde r_+)\varepsilon_{211}-\kappa^{\rm SR}_{322}\alpha^{12}(a_*-\alpha\tilde r_+)\varepsilon_{322}+\kappa^{322\times{\rm BH}}_{211\times211}\alpha^{11}(M_{\rm pl}/f)^4\tilde r_+\varepsilon_{211}^2\varepsilon_{322}\quad(37)$$
Notation: $\kappa^{\rm SR}\equiv\kappa^{\rm SR}_{211}\approx4\times10^{-2}$ (=1/24), $\kappa^{\rm BH}\equiv\kappa^{322\times BH}_{211\times211}=4.3\times10^{-7}$, $\kappa^\infty\equiv\kappa^{211\times\infty}_{322\times322}=1.1\times10^{-8}$ (Sec. IV B 3). The m=0 BH leg carries no angular momentum (absent from (36)) but carries energy ≈ μ into the BH (+ term in (37)). Common simplification: neglect $\dot M$ ⇒ $\varepsilon^{\max}_{211}=a_*(t_0)-4\alpha/(1+4\alpha^2)$ (text after (37); (F1)).

Late-time simplified system (46a–c): drop SR and GW; $a_*\approx4\alpha/(1+4\alpha^2)$; exact attractive ratio $\eta^B=\frac{1}{\gamma^\infty}\big(\gamma^{\rm BH}-\gamma^\infty+\sqrt{(\gamma^{\rm BH})^2-\gamma^{\rm BH}\gamma^\infty+(\gamma^\infty)^2}\big)$ (App. D), small-α $\eta^B\simeq\frac12\frac{\kappa^{\rm BH}\alpha^3\tilde r_+}{\kappa^\infty}\simeq4\times10^{-5}(\alpha/0.01)^3$ (47); decay $\varepsilon_{211}(t)\simeq\varepsilon_{211}(t_*)/\sqrt{1+2\varepsilon_{211}^2(t_*)(t-t_*)/\tau_{\rm scalar}}$ (48), $\tau_{\rm scalar}\equiv\frac{4}{3\mu}\frac{\kappa^\infty}{(\kappa^{\rm BH}\tilde r_+)^2\alpha^{14}}(f/M_{\rm pl})^4\approx0.1\,{\rm yr}(0.1/\alpha)^{14}(10^{-12}{\rm eV}/\mu)(f/10^{17}{\rm GeV})^4$ (49); onset time $t_*\simeq GM\log(GM^2)/[\kappa^{\rm BH}\tilde r_+\alpha^{12}(M_{\rm pl}/f)^4(\varepsilon^{\max}_{211})^2]$ (45).

### 4.4 (c) Quasi-equilibrium (fixed $a_*$; from Gruzinov 2016) [confirmed]
$$\varepsilon^{\rm eq}_{211}\approx\frac{2}{\sqrt3}\frac{\sqrt{\kappa^\infty\kappa^{\rm SR}(a_*-2\alpha\tilde r_+)}}{\alpha^3\kappa^{\rm BH}\tilde r_+}\Big(\frac{f}{M_{\rm pl}}\Big)^2=0.25\Big(\frac{0.01}{\alpha}\Big)^3\Big(\frac{a_*}{0.9}\Big)^{1/2}\Big(\frac{f}{10^{15}{\rm GeV}}\Big)^2\quad(55a)$$
$$\varepsilon^{\rm eq}_{322}\approx\sqrt{\frac13\frac{\kappa^{\rm SR}(a_*-2\alpha\tilde r_+)}{\kappa^\infty}}\Big(\frac{f}{M_{\rm pl}}\Big)^2=6.9\times10^{-6}\Big(\frac{a_*}{0.9}\Big)^{1/2}\Big(\frac{f}{10^{15}{\rm GeV}}\Big)^2\quad(55b)$$
$\eta^{\rm eq}=\gamma^{\rm BH}/(2\gamma^\infty)\approx(\eta^B)_{\rm small\,\alpha}$ (57). Both ∝ f²; numerical values reproduced from the κ's **[my check]**. Equivalently $\varepsilon^{\rm eq}_{211}=(f/f_{\rm eq})^2\varepsilon^{\max}_{211}$ with
$$f_{\rm eq}\approx M_{\rm pl}\Big(\frac{\sqrt3}{2}\frac{\alpha^3\kappa^{\rm BH}\tilde r_+\varepsilon^{\max}_{211}}{\sqrt{\kappa^{\rm SR}\kappa^\infty(a_*-2\alpha\tilde r_+)}}\Big)^{1/2}\approx2\times10^{15}{\rm GeV}\Big(\frac{\alpha}{0.01}\Big)^{3/2}\Big(\frac{a_*(t_0)}{0.9}\Big)^{1/4}\quad(56)$$
Overshoot: $\varepsilon^{\rm thresh}_{211}\approx\sqrt{2\gamma^{\rm SR}_{211}\log(GM^2)/\gamma^{322\times BH}_{211\times211}}\approx(f/f_{\rm thresh})^2\varepsilon^{\max}_{211}$ (58); overshoot occurs for α ≳ 0.04. Exact time to reach SR-halting (E5) with Lambert W; boundary (E7).

### 4.5 (d) Spin-down in the self-interacting regime and the f-criterion [confirmed unless marked]

Regime boundaries (Table II; Fig. 3 for M = 10 M☉, a* = 0.9, T = 10¹⁰ yr):
- A|B: 322 grows within $T_{\rm BH}$: $\gamma^{322\times BH}_{211\times211}(\varepsilon^{\max}_{211})^2\gtrsim\log(GM^2)/T_{\rm BH}$ (41) and not pre-empted by GW annihilation (42)–(43); $f_{\rm AB}\approx\min[3\times10^{16}{\rm GeV}(T_{\rm BH}/10^{10}{\rm yr})^{1/4}(\mu/10^{-13}{\rm eV})^{1/4}(\alpha/0.01)^{11/4},\;8\times10^{18}{\rm GeV}(0.01/\alpha)^{3/4}(a_*/0.9)^{1/4}]$.
- B|C (211 SR halted before saturation): $\gamma^{322\times BH}_{211\times211}(\varepsilon^{\max}_{211})^2\gtrsim2\log(GM^2)\gamma^{\rm SR}_{211}$ (53);
$f_{\rm thresh}\approx M_{\rm pl}\big(\frac{\alpha^3}{2\log(GM^2)}\frac{\kappa^{\rm BH}\tilde r_+(\varepsilon^{\max}_{211})^2}{\kappa^{\rm SR}a_*(t_0)}\big)^{1/4}\approx6\times10^{15}{\rm GeV}(\alpha/0.01)^{3/4}(a_*(t_0)/0.9)^{1/4}$ (54); $f_{\rm BC}\equiv\min[f_{\rm thresh},f_{\rm eq}]$ (59) ≈ $2\times10^{16}{\rm GeV}(a_*(t_0)/0.9)^{1/4}\min[(\alpha/0.04)^{3/4},(\alpha/0.04)^{3/2}]$ (Table II).
- Spin-down timescale in quasi-equilibrium (from (36) with $\varepsilon_{211}=\varepsilon^{\rm eq}_{211}(a_*)$):
$$\tau_{\rm sd}(a_*)\approx\frac{\sqrt3}{2\alpha^5\mu}\frac{\kappa^{\rm BH}\tilde r_+(M_{\rm pl}/f)^2}{\sqrt{\kappa^\infty}\,(\kappa^{\rm SR}(a_*-2\alpha\tilde r_+))^{3/2}}\approx10^7{\rm yr}\Big(\frac{0.01}{\alpha}\Big)^5\Big(\frac{10^{-12}{\rm eV}}{\mu}\Big)\Big(\frac{0.9}{a_*}\Big)^{3/2}\Big(\frac{10^{15}{\rm GeV}}{f}\Big)^2\quad(60)$$
(τ_sd = 1/|ȧ*|.)
- C|D (no full spin-down within age): $\tau_{\rm sd}=T_{\rm BH}$ ⇒
$$f_{\rm CD}\approx3\times10^{14}{\rm GeV}\Big(\frac{10^{10}{\rm yr}}{T_{\rm BH}}\Big)^{1/2}\Big(\frac{10^{-13}{\rm eV}}{\mu}\Big)^{1/2}\Big(\frac{0.01}{\alpha}\Big)^{5/2}\Big(\frac{0.9}{a_*(t_0)}\Big)^{3/4}\quad(62)$$
and for $f\ll f_{\rm CD}$: $|\Delta a_*|\simeq T_{\rm BH}/\tau_{\rm sd}(a_*(t_0))$ (63). Exact (unexpanded) form: $f_{\rm CD}=M_{\rm pl}\big[\frac{\sqrt3}{2}\frac{\kappa^{\rm BH}\tilde r_+}{\alpha^5\mu T_{\rm BH}\sqrt{\kappa^\infty}(\kappa^{\rm SR}(a_*-2\alpha\tilde r_+))^{3/2}}\big]^{1/2}$ **[my algebra from (60)]**.
- **[my check — numerical inconsistency]** With κ^SR = 0.04, κ^BH = 4.3e-7, κ^∞ = 1.1e-8, a* = 0.9, α = 0.01, μ = 10⁻¹² eV, f = 10¹⁵ GeV, (60)'s first line gives τ_sd ≈ 2.4×10⁷ yr (printed: 10⁷ yr). Setting τ_sd = 10¹⁰ yr at μ = 10⁻¹³ eV gives f_CD ≈ 1.5×10¹⁴ GeV from the κ form, 1.0×10¹⁴ GeV from (60)'s numeric line, but (62) prints 3×10¹⁴ GeV. Use the κ form, not the printed prefactors.
- Closed-form forward map for regime C/D **[derived here, not in BGLS]**: with $x\equiv a_*-2\alpha\tilde r_+$ (treat $\tilde r_+$, M fixed), (36)+(55a) give $\dot x=-Kx^{3/2}$,
$K=\frac{2}{\sqrt3}\frac{(\kappa^{\rm SR})^{3/2}\sqrt{\kappa^\infty}\,\alpha^5\mu}{\kappa^{\rm BH}\tilde r_+}\big(\frac{f}{M_{\rm pl}}\big)^2$, so
$$x(t)=\big[x_0^{-1/2}+\tfrac{1}{2}Kt\big]^{-2},\qquad t(a_i\to a_f)=\frac{2}{K}\big[(a_f-2\alpha\tilde r_+)^{-1/2}-(a_i-2\alpha\tilde r_+)^{-1/2}\big].$$
Valid only for $f\ll f_{\rm BC}$ (equilibrium reached at negligible spin-down), α ≲ 0.15–0.2 (two-level validity, see §5), leading-order κ's. Check: $\tau_{\rm sd}=1/(Kx^{3/2})$ reproduces (60).

Regime map summary: $f>f_{\rm AB}$ free; $f_{\rm AB}>f>f_{\rm BC}$ spin-down as free (322 earlier); $f_{\rm BC}>f>f_{\rm CD}$ slow spin-down; $f<f_{\rm CD}$ no spin-down within $T_{\rm BH}$.

Bounds-plot conventions (Sec. VI, Fig. 10): X-ray binaries Cyg X-1, GRS 1915+105, M33 X-7, GRO J1655-40, LMC X-1 (figure legend). Solid contour = **intersection** of m=1 spin-down regions over BH parameters (mass, spin, lifetime, binary period, companion mass) within the observational error intervals; light shading = union. For m=2 only regimes A and B are included. Low-μ edge set by companion tidal mixing, eqs. (condition1) $\frac{M_c}{M}\frac{144\pi^2\sqrt3}{a_*\alpha^7(1+M_c/M)(\mu T)^2}\gtrsim1$ and resonance (condition2) $\frac16a_*\alpha^5\mu T\simeq1$ (App. J; T = orbital period). The actual lifetime values per BH are not tabulated in BGLS **[unconfirmed]**; WM states (citing private communication) that BGLS set limits by equating the BH age with the quasi-equilibrium spin-depletion time and truncated m=2 for f ≲ 10¹⁷ GeV. ABH15 Table III (Cyg X-1 4.8×10⁶ yr, GRS 1915+105 4×10⁹ yr → τ_Edd/10 used, GRO J1655-40 3.4×10⁸ yr, …) is the likely input **[unconfirmed]**. Fig. 3 uses M = 10 M☉, a* = 0.9, T = 10¹⁰ yr.

### 4.6 (e) GW signals in the strong-coupling regime [confirmed]
- Moderate (B): peak annihilation strain unchanged, duration shortened; $h\propto\sqrt{\tau_{\rm scalar}/t}$ at late times (98) instead of $h_{\rm peak}/(1+t/\tau_{\rm ann})$ (97).
- Large (C, D): peak annihilation strain drops ∝ $(f/f_{\rm BC})^2$; "self-interactions strong enough to avoid BH spin-down constraints also render GW annihilation signals undetectable" with current detectors (since $f_{\rm BC}>f_{\rm CD}$); for f ≲ 10¹⁶ GeV expected blind-search events ≲ 10⁻³ (Sec. VII A). Quasi-equilibrium strain O(1–5) below the overshoot peak but lasts ∝ $(f_{\rm BC}/f)^2$.
- Transitions 322→211 appear for $f<f_{\rm AB}$ (eqs. 99–101; Fig. 15): Advanced LIGO only for μ ≳ 10⁻¹¹ eV around ~3 M☉ BHs; better for MAGIS/DECIGO.
- Frequency drifts: (102)–(106), e.g. $\dot\nu_{\lambda}\simeq10^{-10}\,{\rm Hz/s}\,(10^{17}{\rm GeV}/f)^6(\mu/10^{-12}{\rm eV})^2(\alpha/0.1)^{19}$ (106). Note (103) prints $\Delta\omega_\lambda\simeq-3.5\times10^{-5}\mu\alpha^5\varepsilon_{211}(M_{\rm pl}/f)^2$ while App. (B15) (and WM) give 1.2×10⁻⁴ for the 211–211 term **[my check: inconsistency]**.
- Main signal in C/D is non-relativistic axion waves from 322×322→211×∞ (Sec. VIII; requires SM couplings).

### 4.7 (f) Bosenova [confirmed]
Variational 211 ansatz (83), effective potential (86), critical occupation
$$\varepsilon_{\rm crit}=\frac{1024\pi f^2}{9\alpha^3M_{\rm pl}^2}\quad(88),\qquad \varepsilon_{\rm crit}=\frac{32}{711\alpha^2}\sqrt{75840\pi(f/M_{\rm pl})^2+225\alpha^2}-\frac{160}{237\alpha}\ \text{(with self-gravity)}\quad(89)$$
Field amplitude $|\theta_{\max}|\approx\alpha^{5/2}\sqrt{\varepsilon_{211}}(M_{\rm pl}/f)\sqrt{1/8\pi}\,e^{-1}$ (79); $|\theta^{\rm crit}_{\max}|=\frac{8\sqrt2}{3e}\alpha\simeq0.42(\alpha/0.3)$ (text, matches Yoshino–Kodama simulations 0.40–0.45 at α = 0.3). **Not reached**: max $\varepsilon_{211}/\varepsilon_{\rm crit}\lesssim0.3$ for α ≲ 0.2 (Fig. 8); $\varepsilon^{\rm eq}_{211}/\varepsilon^{\rm crit}_{211}\sim0.1$ at small α, small f (Sec. VI). Compare: ABH15 $N_{\rm bosenova}\simeq10^{78}c_0\frac{n^4}{\alpha^3}(M/10M_\odot)^2(f_a/M_{\rm pl})^2$, $c_0\sim5$, $M_{\rm pl}=2\times10^{18}$ GeV (9) [reduced Planck mass], bosenova before full spin extraction if $f_a\lesssim2\times10^{16}{\rm GeV}\,n^{-1/2}((\alpha/n)/0.4)^{3/2}(\Delta a_*/0.1)^{1/2}(5/c_0)^{1/2}$ (10); Arvanitaki–Dubovsky criterion quoted in BGLS (95): $N\gtrsim16\pi\frac{\ell^4}{\alpha}\frac{f^2}{\mu^2}$; WM (44)–(45): $N_{bn}\sim GM^2(n/2)^4\frac{1024\pi f_a^2}{9\alpha^3M_{\rm Pl}^2}$, $\epsilon_{bn}\simeq\frac{357}{\alpha^3}(n/2)^4(f_a/M_{\rm Pl})^2$ (WM use it as a growth cap in plots; also conclude it is not reached).

Beyond two levels (BGLS Sec. IV C): for α ≲ 0.2 only 211, 322 populated; for α ≳ 0.2 411 grows after equilibrium with effective rate $\gamma^{\rm SR\text{-}eff}_{411}=\gamma^{\rm SR}_{411}-(\gamma^{322\times BH}_{211\times411}+\gamma^{211\times\infty}_{411\times322})\varepsilon^{\rm eq}_{211}\varepsilon^{\rm eq}_{322}$ (65); if 322 is the first SR level, 544 is the next to be pumped. Selection rules for pumping via mixing with a damped state (App. E1): m' ≤ 0, m'' = m'+m, |ℓ−ℓ''| ≤ ℓ' ≤ ℓ+ℓ'', ℓ+ℓ'+ℓ'' even.

---

## 5. Witte & Mummery 2412.03655 (relativistic, n ≤ 5)

What they add/correct vs BGLS [confirmed]:
1. **Levels**: n ≤ 3 (211, 311, 322), n ≤ 4 (+411, 422, 433), n ≤ 5 (+511, 522, 533, 544); 321, 421, 431, 432, 532, 542, 543 argued negligible.
2. **Rates**: hydrogen-like tables extended — Table I (n ≤ 4; ~37 processes) and Table II (n = 5 subset; "all scattering permutations between active states are included" in practice). New processes include e.g. $\gamma^{422\times BH}_{411\times411}=2.3\times10^{-7}\alpha^7$, $\gamma^{522\times BH}_{211\times211}=7.5\times10^{-8}\alpha^{11}$, $\gamma^{533\times BH}_{322\times411}=2\times10^{-8}\alpha^{11}$, $\gamma^{211\times\infty}_{311\times311}=5.1\times10^{-8}\alpha^8$, $\gamma^{211\times\infty}_{422\times322}=1.6\times10^{-8}\alpha^8$ (all ×(M_p/f_a)⁴[r̃_p]μ).
3. **Relativistic rates** via Green's functions on Kerr (in/up solutions, (59)–(64)): $R=\frac{1}{W}[R^{\rm up}\int_{r_+}^rR^{\rm in}\mathcal T+R^{\rm in}\int_r^\infty R^{\rm up}\mathcal T]$, $\mathcal T_{\ell m\omega}=-\frac{\lambda}{6}\int d\phi\,d\cos\theta\,\rho^2(\Phi^{(0)}_\omega)^3S^*_{\ell m}$; energy flux (46) through horizon or sphere at ∞. Ratios γ_R/γ_NR shown only in Figs. 3–4 (ã = 0.95) **[magnitudes unconfirmed; graphical only]**. They compute relativistic rates only while all three participating states are superradiant.
4. Hydrogen-like horizon loss (55): $\frac{dE}{dt}=\mathcal A_H|\omega|^2|\lambda R_1(r_+)|^2\simeq4\alpha^2(1+\sqrt{1-\tilde a^2})\lambda^2|R^{(1)}(r_+)|^2$, $\mathcal A_H=8\pi r_gr_+$; emission to ∞ (57): $\frac{dE}{dt}\simeq\frac{\omega\lambda^2\xi^2}{2k}|\int dr'\psi_1\psi_2\psi_3^*\psi_{k\ell m}|^2$ ("equivalent to BGLS"); projection coefficients (53) with ξ = ½ (identical initial states) or 1; rate definition (58) $\dot N\simeq-\frac1\mu\frac{dE}{dt}\equiv\Gamma\,NN'N''$, γ from Γ by $\lambda^2\to\alpha^4(M_p/f_a)^4$.
5. Self-energy shifts (40)–(43) (agree with BGLS (B15)); relativistic corrections to these O(1).

Evolution equations [confirmed]: n ≤ 3, (65)–(69):
$$\begin{aligned}\dot\epsilon_{211}&=\gamma^{\rm SR}_{211}\epsilon_{211}-2\gamma^{322\times BH}_{211\times211}\epsilon_{211}^2\epsilon_{322}-\gamma^{322\times BH}_{211\times311}\epsilon_{211}\epsilon_{311}\epsilon_{322}+\gamma^{211\times\infty}_{311\times311}\epsilon_{311}^2\epsilon_{211}+\gamma^{211\times\infty}_{311\times322}\epsilon_{311}\epsilon_{322}\epsilon_{211}\\&\quad+\gamma^{211\times\infty}_{322\times322}\epsilon_{322}^2\epsilon_{211}-2\gamma^{\rm GW}_{211\times211}\epsilon_{211}^2-3\gamma^\infty_{211^3}\epsilon_{211}^3+\gamma^{\rm GW}_{322\to211}\epsilon_{211}\epsilon_{322}\quad(65)\\
\dot\epsilon_{311}&=\gamma^{\rm SR}_{311}\epsilon_{311}-\gamma^{322\times BH}_{211\times311}\epsilon_{211}\epsilon_{311}\epsilon_{322}-2\gamma^{211\times\infty}_{311\times311}\epsilon_{311}^2\epsilon_{211}-\gamma^{211\times\infty}_{311\times322}\epsilon_{311}\epsilon_{322}\epsilon_{211}-2\gamma^{322\times BH}_{311\times311}\epsilon_{311}^2\epsilon_{322}\quad(66)\\
\dot\epsilon_{322}&=\gamma^{\rm SR}_{322}\epsilon_{322}+\gamma^{322\times BH}_{211\times211}\epsilon_{211}^2\epsilon_{322}+\gamma^{322\times BH}_{211\times311}\epsilon_{211}\epsilon_{311}\epsilon_{322}-\gamma^{211\times\infty}_{311\times322}\epsilon_{311}\epsilon_{322}\epsilon_{211}-2\gamma^{211\times\infty}_{322\times322}\epsilon_{322}^2\epsilon_{211}\\&\quad+\gamma^{322\times BH}_{311\times311}\epsilon_{311}^2\epsilon_{322}-2\gamma^{\rm GW}_{322\times322}\epsilon_{322}^2-\gamma^{\rm GW}_{322\to211}\epsilon_{211}\epsilon_{322}\quad(67)\\
\dot{\tilde a}&=-\gamma^{\rm SR}_{211}\epsilon_{211}-\gamma^{\rm SR}_{311}\epsilon_{311}-2\gamma^{\rm SR}_{322}\epsilon_{322}\quad(68)\\
\frac{\dot M}{\mu GM^2}&=-\gamma^{\rm SR}_{211}\epsilon_{211}-\gamma^{\rm SR}_{322}\epsilon_{322}+\gamma^{322\times BH}_{211\times211}\epsilon_{211}^2\epsilon_{322}+\gamma^{322\times BH}_{311\times311}\epsilon_{311}^2\epsilon_{322}+\gamma^{322\times BH}_{211\times311}\epsilon_{211}\epsilon_{311}\epsilon_{322}\quad(69)\end{aligned}$$
(WM's (69) omits $-\gamma^{\rm SR}_{311}\epsilon_{311}$ though (68) includes it.) n ≤ 4: (70) ε̇211, unnumbered ε̇322, (71) ε̇411, (72) ε̇422, (73) ε̇433; general $\dot{\tilde a}=-\sum_{n\ell m}m\gamma^{\rm SR}_{n\ell m}\epsilon_{n\ell m}$ (74), $\frac{\dot M}{\mu GM^2}=-\sum\gamma^{\rm SR}\epsilon+\sum_{ijk}\gamma^{k\times BH}_{i\times j}\epsilon_i\epsilon_j\epsilon_k$ (75). n = 5 equations not written out ("straightforward to infer"). CFW (1)–(3) give the compact generic form $\dot\epsilon_{n\ell m}=2\Gamma_{n\ell m}\epsilon_{n\ell m}+\sum\zeta_{ij}\gamma^{k,L}_{i,j}\epsilon_i\epsilon_j\epsilon_k$, $\dot{\tilde a}=-\sum2m\Gamma\epsilon$, $\dot M/(\mu GM^2)=-\sum2\Gamma\epsilon+\sum\gamma^{k,BH}_{i,j}\epsilon_i\epsilon_j\epsilon_k$ (ζ = 2 for i = j; signs implicit).

Size of corrections / physics conclusions [confirmed]:
- 311 never reaches non-negligible occupation (n ≤ 3 is effectively two-level).
- Relativistic corrections to $\gamma^{433\times BH}_{211\times322}$ ≫ those to $\gamma^{211\times\infty}_{322\times433}$ ⇒ 433 grows; two-level system accurate only for **α ≲ 0.15** (vs ≲ 0.2–0.22 in BGLS), in agreement with Omiya et al.
- At large α the 322/422/433 grow to the bosenova cap unless the n = 5 sink (544) is included: "for the m = 2 states, this leading order contribution will arise at |544⟩, while for the m = 3 states … at |766⟩". Trust n ≤ 5 results for m ≤ 2 spin-down, α ≲ 0.66.
- m = 1 limits: relativistic corrections "slightly suppress" the constraint; n = 4, 5 have "minimal impact on the m = 1 spin down region"; m = 2 limit suppressed by "a factor of nearly 10 in f_a at its peak" when n = 5 is included; BN-cut-only treatment extends limits to ~10× smaller f_a than level mixing.
- **f-sensitivity change**: "prior limits had significantly underestimated the sensitivity to f_a (by roughly an order of magnitude)" — attributed to the **limit-setting procedure** (MCMC over (f_a, m_a, M_i, a_i), evolving from $t_0-\tau_{\max}$ to $t_0$, 2σ marginal posterior) versus BGLS's "age = quasi-equilibrium spin-depletion time" criterion: ~10× for Cyg X-1 (high-precision spin), ~2× for M33 X-7; including m = 2 extends to ~2× heavier axions.
- Conventions: $\tau_{\max}={\rm Min}[\tau_{bh},\tau_{\rm Sal}]$, $\tau_{\rm Sal}\sim4.5\times10^7$ yr, lower limit of τ_bh used; priors (eqs. a_prior, mass_prior); BH table (Cyg X-1 fid. 21.2 ± 2.2 M☉, ã > 0.9985, τ 4.8–7.6 Myr; cnsrv ã = 0.92⁺⁰·⁰⁵₋₀.₀₇; GRS 1915+105 12.4⁺²·⁰₋₁.₀ M☉, ã > 0.95, 3–5 Gyr; GRO J1655-40; LMC X-1; M33 X-7); companion criterion quoted from BGLS $\alpha\gtrsim0.013(\frac{M}{M_\odot}\frac{1\,{\rm day}}{P_{\rm orb}})^{2/9}\tilde a^{-1/9}(\frac{M}{M_{\rm comp}+M})^{1/9}$. Fig. 1: combined 2σ limit from Cyg X-1 (conservative) + GRS 1915+105, fiducial m ≤ 2.
- WM report "good agreement with BGLS … except for two processes" (not named). **[my check]** Table comparison: $\gamma^{322\times BH}_{411\times411}$ = 1.6e-11 (WM) vs 9.8e-11 (BGLS Table V); $\gamma^\infty_{211^3}$ = 1.5e-8 (WM) vs 5e-9 (BGLS Table VII; BGLS divides the 3→1 power by 3).
- **Code**: footnote "All code … can be found at" `https://github.com/SamWitte/Axion_SR_TDE` (biblio.bib). As of 2026-11-03 that URL redirects (HTTP 301) to the public repo **`https://github.com/SamWitte/Axion_SR`** (Julia/Jupyter; last push 2026-11-03) **[checked via GitHub API]**. NSW used this codebase for numerical Γ_SR and a self-interaction test.

---

## 6. Compact f-criteria in NSW 2607.01317 and CFW 2507.21788

- **Neither paper gives a compact f_crit(α, M, τ) formula** [confirmed by full-text search].
- NSW (QCD axion only): rate (4) $\Gamma^{[n,\ell,m]}_{\rm SR}=2\alpha^{4\ell+5}M^{-1}(m\chi-2m_ar_+)C_{n\ell m}\mathcal K(\chi)$, $\mathcal K=\prod_{k=1}^{\ell}[k^2(1-\chi^2)+(m\chi-2m_ar_+)^2]$, "$C_{0,1,1}=1/24$"; threshold (5) $\frac{\chi}{2(1+\sqrt{1-\chi^2})}=\frac{\alpha}{m}$; saturation spin (6) $\Gamma_{\rm SR}(m_a,M,\chi_{\rm sat})=N_{\rm efolds}/\tau_m$, $N_{\rm efolds}=180$, valid if $\tau_m>\tau_{\rm grow}\equiv N_{\rm efolds}\tau_{\rm SR}(\chi_I)$; $\tau_m$ log-uniform in $[10^7,\tau_H]$ yr; modes m = 1–6. On self-interactions (SM "Superradiance for higher α"): with WM's code "they have a generally negligible effect on the final spins when assuming the QCD axion decay constant … final spin χ may be altered by up to ∼0.05"; "a more intensive analysis for other f_a values could be an avenue for future work".
  **[my check — discrepancy]** Taken literally, (4) with C₀₁₁ = 1/24 gives Γ_SR,211 = μα⁸(χ−2αr̃₊)/12: 2× the occupation-number rate (1/24) and 4× Im ω (1/48). Consistent only if C₀₁₁ = 1/48 (BCST C₂₁) with the explicit 2 being Γ_N = 2 Im ω. Effect on χ_sat is logarithmic-ish (shifts the curve), but should be fixed in any reimplementation.
- CFW (GW231123): evolve (1)–(3) with WM rates, n ≤ 5, $\tau_{\max}=10^5$–$10^7$ yr; result $f_\Phi\gtrsim10^{14}$ GeV for $\mu\sim[0.6,5]\times10^{-13}$ eV (abstract; 90% lower limit from marginal posterior).
- Practical compact criterion available: BGLS (60)/(62) (two-level, leading order) or the closed-form map of §4.5 [derived].

---

## 7. Implementation notes

### Which formulas
**n ≤ 3 code (two-level; trust for α ≲ 0.15, cf. WM):**
1. State vector $(\varepsilon_{211},\varepsilon_{322},a_*,M)$; optionally 311 (WM (66), negligible in practice).
2. $\gamma^{\rm SR}_{n\ell m}=2\,\mathrm{Im}\,\omega_{n\ell m}(\alpha,a_*)$ evaluated **at the current (α, a*)**; use numerical Im ω (SuperRad interpolation for 211/322, or Leaver continued fraction), BCST (2.32)–(2.34) only as a small-α cross-check. Recompute r̃₊ and Ω_H each step; α = GMμ changes with M.
3. Interaction rates: BGLS Table I/V/VI κ's with leading-α scalings as in (34)–(35); for α ≳ 0.1 multiply by the finite-α correction (BGLS Fig. 21: κ^BH up "a factor few" at α ∼ 0.2; WM Figs. 3–4), or compute them yourself with WM's Green's-function method / WM code. Keep 211×211→GW (κ ≈ 1e-2, or SuperRad P_GW) and 322→211 GW (5e-6 α¹⁰) for the f ≳ f_AB regime.
4. Spin/mass: (36)/(37) [or WM (68)/(69)].
5. Initial conditions: ε ~ 1/(GM²) (one quantum); ≈180 e-folds to saturation for 10 M☉ (log GM² ≈ 179 for 10 M☉). Use a stiff integrator in log ε (rates span >10 orders).
6. Shortcut for f ≪ f_BC: quasi-equilibrium (55a,b) + spin-down ODE (36), or the closed form in §4.5.

**n ≤ 5 code (to extend to m = 2 region, α ≲ 0.66):** WM (65)–(75) with Table I/II rates (hydrogen-like) times relativistic correction factors; include every permutation allowed by selection rules (m conservation with an m = 0 BH leg or an emitted wave; energy ω_i+ω_j−ω_k < μ → BH leg, > μ → emission to ∞; BGLS App. E1 rules) — WM's code is the reference implementation. Include 544 whenever 322 can grow (sink for m = 2), and keep in mind m ≥ 3 spin-down needs n ≥ 7 (766) sinks.

### Pitfalls
- **Im ω vs N-growth**: BCST/WM Γ_nlm and SuperRad ω_I are amplitude rates; occupation and cloud mass grow at 2 Im ω (BGLS Γ_SR ≡ 2ω_i; CFW eq. 1; SuperRad (36) $\dot M_c=2\omega_IM_c+\dots$).
- **Detweiler factor 2**: 1/24 in Detweiler's original Im ω is wrong by 2 (BCP18 footnote); correct Im ω coefficient is C₂₁ = 1/48. NSW's printed (4) with C₀₁₁ = 1/24 is 2×(N-growth) if read literally.
- **r_g = GM** (not 2GM) in all papers here; $\tilde r_+=1+\sqrt{1-a_*^2}$; in BCST (2.34) $r_+\omega$ is dimensionless product with dimensionful r₊ ($=\tilde r_+\alpha\,\omega/\mu$). BCP18 replaces ω→μ inside g.
- **Planck mass**: BGLS M_pl = 1/√G = 1.22×10¹⁹ GeV in all (M_pl/f) factors; ABH15 uses M_pl = 2×10¹⁸ GeV (reduced) in (9)–(10); NSW quote reduced m_pl for the QCD relation. Mixing M_pl with the reduced Planck mass (M_pl/√(8π) = 2.4×10¹⁸ GeV) is a factor √(8π) ≈ 5 in f, i.e. (8π)² ≈ 630 in f⁴ (corrected 2026-11-10; the earlier "2.4 in f, 35 in f⁴" was wrong).
- **N counting**: occupation number in units of ħ, $N_{\max}\simeq GM^2\Delta a_*/m$ (ABH15 (8)); ε ≡ N/(GM²) with GM² = (M/M_pl)² ≈ 8×10⁷⁷ (M/10M☉)²; BGLS App. F normalizes ε by the **initial** BH mass. ε_max = Δa* only if Ṁ is neglected (cloud can reach ~7% of M, App. F).
- **Real vs complex field**: BGLS/ABH15/WM: real field, $\varphi=\frac{1}{\sqrt{2\mu}}(\psi e^{-i\mu t}+{\rm c.c.})$, $\int|\psi|^2=N$ (BGLS (10); ABH15 (A5) $\sqrt{N/2\mu}\Psi$; WM text above (36)). Complex field: pure-particle configurations have static T_μν and no annihilation GWs (BCP14 below (8); BGLS Sec. VII A; BCP18 Sec. 2.2). Gruzinov 2016 hydrogenic normalization √2 too large ⇒ rates ×4 too large (BGLS after (18)). BCP14's action is for a complex Ψ with ½ normalization but the cloud model (8) is real.
- **Multiplicities**: for $j\times j\to k\times{\rm BH}$: $\dot N_j=-2\Gamma N_j^2N_k$, $\dot N_k=+\Gamma N_j^2N_k$; for $j\times j\to k\times\infty$ same structure; 3→1: −3Γ (BGLS Γ = P/(3μN³); WM tabulate the power coefficient 1.5e-8 instead of 5e-9). BH-absorbed processes add mass (+μ per event) but no angular momentum when the BH leg has m = 0.
- **Spin equation**: BGLS (36) and WM (68) set $\dot a_*=\dot J/GM^2$ and omit $-2a_*\dot M/M$; include it (BGLS App. F (F2)–(F4)) for α ≳ 0.2 or when the cloud mass is non-negligible.
- **GW angular momentum**: $\dot J_{\rm GW}=(m/\omega_R)\dot E_{\rm GW}$ per level (BCP14 (14) for m = 1).
- **Leading-order tables**: BGLS Table I/IV SR coefficients for ℓ ≥ 2 use g at ã = 0 (§4.2); several BGLS printed numerical prefactors are mutually inconsistent at the factor 2–3 level (τ_sd (60) vs κ form; f_CD (62); GW transition 3e-6 vs 5e-6; Δω_λ (103) vs (B15)) — always rebuild numbers from κ's.
- **Flat-space GW coefficients are gauge-ambiguous** (YK Sec. 3.4): flat vs Schwarzschild differ ~16× for 211; use SuperRad numerical P_GW.
- **Numerics**: continued-fraction Im ω at small α needs N ≳ 1000–5000 terms (BCST Fig. "contFracResults": ~1000 at α ∼ 0.02; SuperRad N ≥ 5000), high precision (WM: 256-bit BigFloat).
- **Sign of λ**: only matters for frequency shifts/bosenova; rates ∝ λ².

---

## 8. Discrepancies found (summary)

1. Detweiler 1/24 vs correct 1/48 for Im ω₂₁₁ (BCP18 footnote; BGLS; BCP14). "μα⁸(χ−2αr̃₊)/24" is right only as the N-growth rate.
2. NSW eq. (4): prefactor 2 together with C₀₁₁ = 1/24 ⇒ 2× (vs N-growth) / 4× (vs Im ω) if literal.
3. BGLS Table I: $\gamma^{\rm SR}_{322}=\Gamma^{\rm SR}_{211}$ typo; Table IV ℓ ≥ 2 SR coefficients omit the spin dependence of g_ℓm.
4. BGLS printed τ_sd (60) prefactor ~2.4× low vs its own κ expression; f_CD (62) ~2–3× high vs (60).
5. BGLS GW transition coefficient 3e-6 (near (50)) vs 5e-6 ((101), Table IV; recomputed 5.18e-6).
6. BGLS Δω_λ (103) coefficient 3.5e-5 vs 1.2e-4 in (B15) and WM (40).
7. BGLS text scalings $P\propto\alpha^{16+4\ell}$ (Sec. II) / $\alpha^{16+2(\ell+\ell')}$ (Sec. VII A) inconsistent with Table I (α¹² for 211 at fixed N, μ).
8. WM vs BGLS rates: $\gamma^{322\times BH}_{411\times411}$ (1.6e-11 vs 9.8e-11); 3→1 normalization (factor 3); WM note two unexplained discrepancies.
9. ABH15 Table VI: 2p and 3d correspond to YK's radiation-gauge flat result [(ℓ−1)²/4 factor]; 4f entry off by 2²⁰ (typo); BGLS's 322 annihilation coefficient follows ABH15 (≈3× below YK flat).
10. Analytic Im ω₂₁₁ (BCST 2.32) underestimates Dolan's numerical peak rates by ×1.8–15 for α = 0.19–0.42 (my check).
11. WM Lagrangian sign ($\mathcal L\supset-\lambda\Phi^4/4!$ with λ = μ²/f²) contradicts their own EoM (34) and the cosine potential; SuperRad (13) and (36) have sign typos ($\omega_R-m\Omega_H$; $+P_{\rm GW}$).
12. BGLS: two-level valid for α ≲ 0.2; WM (relativistic): α ≲ 0.15.
