# Manuscript source: model deviations & validation notes

Working draft for the paper's model/deviations discussion and validation subsection.
Derived from RESULTS.md + session findings (2026-09-19). Numbers are final unless
marked otherwise. Companion file: RESULTS.md (check IDs referenced below as V0–V4).

---

## 1. Deviations from the reference implementation (Dirkx & Mooij 2017)

| # | Item | D&M reference | Present implementation | Quantified impact | Disposition |
|---|------|---------------|------------------------|-------------------|-------------|
| 1 | Bank-modulation law | Eq 2.41, cos σ = (m/L)·[(g−V²/r)cos γ − 2ωV cos δ **sin χ**]; active only while 0 < cos σ < 1 | Identical form (sin χ), rate-limited: σ is a propagated 7th state, dσ/dt clamped at 5°/s | Kills 0↔90° chatter of the discontinuous command | Adopted (see §3 sensitivity) |
| 2 | Eq 2.40 vs 2.41 | Book prints sin γ in Eq 2.40, sin χ in Eq 2.41 | sin χ used (EOM-consistent); sin γ tested as variant only | — | Eq 2.40's sin γ treated as a typo in the book |
| 3 | Lee-side pressure law | Table 3.2: ACM (Eq 3.58) on blunt lee surfaces, PM on rounded | Prandtl–Meyer in both regimes | ≤ 0.05 on CD(α = 0) | Documented simplification, now quantified |
| 4 | Tangent-cone beyond shock detachment | Fig 3.5 correlation continues to rise past θ_detach | Exact Taylor–Maccoll, saturating + Newtonian fill past detach | V1b core ≤ 0.0605 Cp (band 0.065); worst at M = 3 where θ_detach = 50° | Documented method difference |
| 5 | Base pressure | ACM base term available | Not modelled | — | Documented simplification |
| 6 | Shadowing | Available for non-convex bodies | Not required (capsule is convex) | — | Not applicable |
| 7 | Cm sign/reference | Cm about c.g., nose-up positive, ref length 3.9116 m | Cm_DM = −Cm_ours · L_total / 3.9116 | Conversion exact (magnitudes ≤ 0.005 after flip, same trim α and slope) | Convention mapping stated in paper |
| 8 | Mach grid | 3–25 | 3–30 (10 nodes; entry starts at M ≈ 27.6) | PM at vacuum limit above M = 25: ΔCp ≈ 7×10⁻⁴ | Negligible extension |
| 9 | Load-factor definition | Eq 6.13: n = √(D²+L²)/(mg) (total aero force) | Code metric is D-only (thesis convention); n_tot computed post-hoc for the paper | n_D 5.06 vs n_tot 5.33 g | Paper constraint uses Eq 6.13; both reported |

**Not resurrected:** the thesis's modified-Newtonian aero is not re-run for the paper;
thesis Ch-4 MN numbers are cited as the deterministic baseline, and all paper
validation is against the local-inclination implementation.

## 2. Validation targets and data provenance

| Layer | Reference data | Result |
|-------|----------------|--------|
| Atmosphere (US-76 + 86–120 km extension) | US-76 anchor values | V0 6/6 PASS (anchors ≤ 2.1×10⁻⁵ rel; 86-km join ≤ 4.2×10⁻⁴; a₀ = 340.294 m/s) |
| Surface-pressure laws (MN, TC, PM, blend) | D&M Figs 3.5/3.6 correlations (digitized) | V1a–V1d all PASS (MN ≤ 0.026 Cp; TC core ≤ 0.061; PM anchors ≤ 0.005°; blend w(5)=0, w(12)=1) |
| Apollo force coefficients | Apollo wind-tunnel data **as re-plotted by Dirkx & Mooij (Fig 7.7)**, M = 3 and 10, digitized (PlotDigitizer, per-subpanel calibration) | V2 9/9 PASS: L/D absolute ≤ 0.029; CD/CL/Cm shape ≤ 0.015 except CD-M3 0.061 (see #4 above) |
| M = 10 CD over-prediction (inherent LIM trait) | D&M WT curve | Present 1.610 vs WT 1.523 (+0.087); D&M's own LIM 1.704 (+0.181) — present model sits between, closer to WT |
| Convective heating (Reff) | Zoby–Sullivan chart (digitized) | V3a ≤ 0.0143 over the design space |
| Shoulder heating ratio | Marvin & Sinclair; Jones; Wadham (primary WT, via thesis Table 4.1) | V3b ≤ 10.2%, trends reproduced (rises with α, falls with Rs) |
| Trim | D&M Fig 7.17a | V4a PASS (0.25° max dev) |
| Trajectory / end-to-end | D&M Fig 7.17 + Table 7.1, exact Apollo case (V 7.83 km/s, γ = −4°, m = 4532 kg; case identity verified against book text pp. 152–153) | See §3 |

Provenance wording for the paper (coefficients): *"validated against Apollo
wind-tunnel measurements as re-plotted by Dirkx & Mooij"* — a primary-source
comparison (TN D-3890 Fig 8) remains available as a revision-stage upgrade.

Digitization protocol note: figures digitized in PlotDigitizer with per-subpanel
calibration; tracing artifacts identified (mis-traced return branches, stray
points) were removed with a documented automated filter (rolling median in
log-space + monotonicity truncation) and all verdicts re-verified stable.

## 3. Guidance-law sensitivity & reference-data inconsistency (headline finding)

**3.1 Reference-data inconsistency.** Panel (c) of D&M Fig 7.17 plots q̇(M) and
n(M) for one simulation; these are algebraically coupled independent of
trajectory model:

    n = q̇² · R_eff · C_D · A / (2 m g k² V⁴)

(identity validated to ~2% on our own trajectory). Applied to the digitized
curves, panel/implied ratio swings 0.17–0.65 over M = 6.5–22.6; no single
parameter (mass, R_eff, C_D, heating constant) reconciles them (required q̇
rescale varies 2.4×–1.24×). The n-axis itself is correctly digitized (0–8, peak
≈ 4.4–5 g, confirmed by direct reading). The book's own q̇ curve implies
n_peak ≈ 11 g; the panel prints 4.35 g. Table 7.1 (Qs = 53.06 MJ/m²) is below
every member of our variant family while the book's q̇ curve runs 1.5–2.7× ours
below M ≈ 13 — table and figure disagree as well. **Consequence: the n panel is
treated as qualitative only; peak-load-factor is not a valid validation target.**

Independent anchor: Apollo CM flights peaked at ~3 g (LEO-class) to ~7 g (lunar
return); Apollo 11 ≈ 6.6 g (Hillje TN D-5399; Moseley & Wells TN D-5514). For
this shallower 7.83 km/s case the 5–7 g class is expected — consistent with the
present BASE result (n_D 5.06, n_tot 5.33 g), inconsistent with the panel's 4.35 g.

**3.2 Guidance-family sensitivity.** Variant matrix (D&M: Qs 53.06 MJ/m²,
Sg 1430.2 km, q̇_pk 0.616 MW/m²):

| variant | Qs | Sg | q̇_pk | n_D | n_tot | bank onset |
|---|---|---|---|---|---|---|
| BASE (Eq 2.41, sin χ) | 70.12 (+32%) | 1657 (+15.9%) | **0.621 (0.7%)** | 5.06 | 5.33 | M 24.4 |
| GATE ε=0.15 + cap 82° | 65.67 (+24%) | 1519 (+6.2%) | 0.663 | 6.56 | 6.91 | M 25.8 |
| GATE ε=0.30 + cap 82° | 64.15 (+21%) | 1477 (+3.3%) | 0.683 | 7.13 | 7.52 | M 26.6 |
| Eq 2.40 as printed (sin γ) + cap 82° | 61.56 (+16%) | **1414.4 (−1.1%)** | 0.722 (+17%) | 7.99 | 8.41 | entry |

Sg and q̇ move in opposite directions across the family: **no member reproduces
the printed Table 7.1 and the figure panels simultaneously**; the family floor
for Qs is +16%. Per-metric agreement is achievable (q̇ to 0.7%, Sg to −1.1%),
which is the framing used for V4b/V4c reporting.

## 4. Corrections carried from the thesis into the paper

1. **"Rm = 0.196" typo** → Rm = 1.956 m (0.196 is Rs).
2. **Peak-q̇ inconsistency (0.76 vs 0.4 MW/m²)** → 0.76 MW/m² confirmed by
   independent reproduction of the thesis pipeline (the 0.4 value is unsupported).
3. **Peak n = 9.5 g (thesis Table 4.2)** → was an artifact of the clipped bank
   law (σ saturating 150–180°, lift effectively down); with the corrected law
   n_tot = 5.33 g, in the Apollo flight-data class (§3.1).
4. **Validation geometry**: paper validation uses the exact Apollo geometry
   (Rm = 1.956 m, Lc = 2.662 m, design vars 0.4235/0.4632/0.5091 m, c.g.
   inverted through the c.g. model, z_cg sign set for stable +α trim, trim
   ≈ 21.9°/25.3° at M = 10/3); the thesis validated with Rm = 2, Lc = 2 without
   stating so.
5. **Thesis Table 6.2 Sg 2930–4695 km vs validation 1330 km**: different entry
   angle (γ = −2° optimization cases vs −4° Apollo validation case) — pre-emptive
   sentence added.
6. **Thesis Ch-4 MN-based validation numbers are cited, not re-run** (§1, last
   paragraph).
