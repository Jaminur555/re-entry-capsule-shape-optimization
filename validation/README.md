# Validation Plan — against Dirkx (2017), the motivation paper

All validation data come from the motivation paper (Dirkx 2017, `1.dirkx2017.pdf`,
book-page references below) and were digitized by the user in PlotDigitizer
(`../Extracted Data.txt` at repo parent). Every check compares **our
`capsule_opt` implementation against the paper's own published curves**, so a
pass certifies the code as *D&M-faithful* — the credibility claim Phase 1 needs.

## 1. Provenance — where each dataset comes from in the paper

| Our dataset | Paper source | Content |
|---|---|---|
| `reference_data/cp_methods_dm2017_digitized.csv` | Fig 3.5 & 3.6 (pp. 66–71) | Cp vs local inclination θ: Modified Newtonian, Tangent Wedge, Tangent Cone (correlation); M = 3, 5, 10, ∞ |
| `reference_data/cone_tables_reference.csv` | Analytic Taylor–Maccoll (no digitization) | Exact-TC benchmark values, e.g. M=10: θ=5/15/30/40/55° → Cp = 0.019/0.144/0.529/0.875/1.465 |
| `reference_data/prandtl_meyer_reference.csv` | Analytic (Eq. 3.55) | Closed-form PM values for the expansion branch |
| `reference_data/apollo_flight_wind_tunnel_data.csv` | Fig 7.7 (pp. 144–146) | Apollo CD, CL, Cm, L/D vs α at M = 3 & 10; **both** D&M's LIM curves and the NAA (1965) wind-tunnel curves |
| `reference_data/zoby_sullivan_chart_digitized.csv` | Fig 3.9 (p. 81) | Rm/Reff vs Rm/Rn for Rs/Rm = 0, 0.1, 0.2, 0.25, 0.3 |
| `reference_data/apollo_entry_reference_dm2017.csv` | Fig 7.17 (pp. 152–153) | Entry profiles vs M (α_trim, L/D, bank, q̇_cs, n) and h vs t, for both WT-aero and LIM-aero runs |
| `reference_data/dm2017_table7_1.csv` | Table 7.1 (p. 153) | Scalars: heat load (WT 49.232 / LIM 53.06 MJ/m²), ground track (1425.6 / 1430.2 km) |

Paper facts verified 2026-09-03 against the PDF: Apollo geometry
**RN = 4.694 m, Rm = 1.956 m, RS = 0.196 m, θc = 33°, Lc = 2.662 m**
(Hirschel & Weiland 2009); reference quantities **cref = 3.9116 m** (= 2·Rm),
**rcom = (1.0367, 0, 0.1369) m**; validation-entry ICs **h = 120 km,
τ = 225.5°, δ = −23.75°, V_R = 7.83 km/s, γ = −4°, χ = 49.6°, m = 4532 kg**,
integrated by D&M with fixed-step RK4 at 0.1 s. D&M's own stated accuracy for
the LIM-vs-WT comparison: α_trim off ≈ 1° (M=10) / 0.1° (M=3); L/D error ≈ 10%.
Table 7.1 percentages (+7.7754 / +0.3227) reproduced exactly in the extraction.

## 2. Conventions (must be applied before any overlay)

| Item | Paper (D&M) | Us | Mapping |
|---|---|---|---|
| AoA sign | trim at α ≈ −25° (negative) | trim at α ≈ +25° | **α_DM = −α_ours** |
| Lift sign | CL > 0 at α < 0 | CL > 0 at α > 0 | consistent under the flip |
| Moment ref | rcom, normalized by Sref·cref | CG, normalized by A_ref·L_ref | place CG at rcom, then **Cm_DM = −Cm_ours · L_total/3.9116** (sign flipped — the book's positive Cm is our negative; established empirically in V2: magnitudes match to ≤0.005 after the flip) |
| Reference area | printed "Sref = 39.441 m²" | A_ref = π·Rm² | use π·Rm² = 12.02 m² (the printed 39.441 is inconsistent with CD ≈ 1.5 in Fig 7.7 — treat as a book typo; note in manuscript) |
| Mach validity | DB spans full entry (M ≈ 0.28–27.6 in figures) | aero valid M ≥ 3, DB grid 3–25 | **restrict trajectory comparison to M ≥ 3**; entry starts at M ≈ 27.6 → clipped at 25 (consider extending `mach_nodes` to 30) |
| Shadowing | none for capsules (§7.3.1) | none | identical (our earlier "simplification" is actually D&M-faithful — upgrade the manuscript wording) |

**Geometry setup for the Apollo validation case** (all aero/trajectory checks):
override `config.Rm_fixed = 1.956`, `config.Lc_fixed = 2.662` in-script (both
are read at call time), design params from the inverse map
`rn = (4.694−3)/4`, `rs = (0.196−0.02)/0.38`, `r_theta = (33−5)/55`.
Solve `cg_params = (dx, dz)` so that (x_cg, z_cg) equals rcom mapped into our
axes — the D&M z-offset is **+0.1369 m (up)** while our `z_cg = −dz·h_local`,
so dz < 0 is expected; verify the sign empirically by requiring trim near
+25° on the stable branch. The thesis-fixed geometry (Rm = 2, Lc = 2) is *not*
used for validation — deltas vs D&M would be geometry, not model, error.

## 3. The checks

**V1 — Cp(θ) law level** (`aerodynamics/validate_newtonian.py`,
`validate_tangent_cone.py`, `validate_prandtl_meyer.py`)
- Modified Newtonian: our `cp_newtonian(θ; Cpmax(M))` vs digitized curves at
  M = 3, 5, 10, ∞. Pass: |ΔCp| ≤ 0.03 or ≤ 3% per point (digitization noise
  dominates; digitized Cp_max(∞) ≈ 1.852 vs our limit ✓).
- Tangent cone: our cached exact Taylor–Maccoll table vs the digitized
  **correlation** curves (M = 3, 5, 10). The paper's panel is explicitly the
  correlation (§3.3.1: correlations "give good results below the critical
  value and high Mach numbers"), ours is the exact solution → expect small
  systematic offset, growing toward low M / high θ. Pass: |ΔCp| ≤ 0.05 for
  θ ≤ 55°; document the offset as exact-vs-correlation, not error. Cross-check
  the table itself against `cone_tables_reference.csv` (exact analytic values,
  must match to 1e-3).
- Prandtl–Meyer: analytic self-check vs closed form (no D&M curve exists).
- Blend continuity (`validate_blend_continuity.py`): Eq. 3.48–3.49 continuity
  across M ∈ [5, 12] (our continuous-form implementation, already verified).

**V2 — Apollo force & moment coefficients** (`aerodynamics/` — the core check;
data = Fig 7.7). α-sweep 0–30° at M = 3 and 10 on the exact-Apollo setup
above; overlay on both D&M-LIM and NAA-WT curves with the §2 mappings.
Pass bands (anchored to D&M's own stated accuracy):
- CD: our-LIM within 0.05 of D&M-LIM over α ∈ [0, 30]; reproduces the known
  M=10 over-prediction vs WT near α = 0 (CD_LIM ≈ 1.70 vs CD_WT ≈ 1.51) —
  this *method behaviour* must appear, it is not a failure;
- CL: within 0.03 where digitized data exist (CL M=10 LIM was never digitized
  — compare only M = 3 LIM and M = 10 WT);
- Cm: within 0.01 after the cref conversion; dCm/dα of the correct sign;
- α_trim within ~1° of −25.0° (DM convention); L/D(trim) ≈ 0.30 ± 0.03.

**V3 — Heating, two scripts (thesis §4.3.1/§4.3.2 kept separate):**
- `heating/validate_chapman.py` — (a) Zoby–Sullivan: our
  `effective_nose_radius()` interpolates an in-code digitization of Fig 3.9;
  cross-validate against the user's *independent* digitization
  (`zoby_sullivan_chart_digitized.csv`). Pass: Apollo point (Rm/RN = 0.417,
  RS/Rm = 0.100) and a 5×5 grid over Rm/RN ∈ [0, 2.5], RS/Rm ∈ [0, 0.3]
  agree within 3% on Reff. (b) Stagnation-flux sanity on the Apollo entry:
  peak q̇_cs within ~5% of the digitized Fig 7.17c LIM peak.
- `heating/validate_shoulder_heating.py` — reproduce thesis Table 4.1:
  shoulder/stagnation ratio vs the Marvin&Sinclair / Jones / Wadham source
  points (`shoulder_ratio_reference.csv`). Pass: within ~10% (thesis level),
  correct trends vs α and RS.

**V4 — Trim + entry trajectory, end-to-end** (`trim/validate_trim_solver.py`,
`trajectory/validate_trajectory.py`, `system/validate_apollo_end_to_end.py`;
data = Fig 7.17 + Table 7.1). Exact-Apollo setup, rcom CG, m = 4532 kg, ICs
with **γ = −4°** (`propagate_trajectory` default is −2° — the optimization
case; pass −4 explicitly). Our RK45 adaptive vs their RK4 0.1 s — compare via
interpolated overlays on the M ≥ 3 portion. Two comparison layers:
- vs D&M's **LIM** curves (same method → code-vs-book fidelity; tight bands):
  α_trim(M) within ~1° incl. the slope change below M = 10 (regime blending);
  L/D plateau ≈ 0.332 ± 0.015 then decline; bank-angle modulation onset/peak/
  end within a few degrees; peak q̇_cs within 5%; peak n within 0.2 g;
  h(t) within ~1 km through mid-entry;
- Table 7.1 scalars: heat load 53.06 MJ/m² ± 8%, ground track 1430.2 km ± 2%.
  (The WT curves provide context only; the LIM-vs-WT gap ≈ 0.03 in L/D should
  already be established by V2.)

**V0 — atmosphere** (`atmosphere/validate_atmosphere.py`): US-76 values vs
`US76_reference_table.csv`, incl. the 86–120 km extension join (continuity,
no step in ρ, T, p). Minor, run first as a smoke test.

## 4. Known deviations from D&M (state explicitly in the manuscript)

1. Gravity: D&M use oblate Earth + J2 (EGM96); we use the US-76 geopotential
   spherical Earth (r = 6356.767 km), no J2. Small extra trajectory delta.
2. Integrator: RK45 adaptive (D&M validation: RK4 fixed 0.1 s).
3. Atmosphere 86–120 km: tabulated extension (variable composition) rather
   than plain US-76 extrapolation.
4. Tangent-cone: exact Taylor–Maccoll vs D&M's correlation (V1 quantifies it).
5. Lee-side: Prandtl–Meyer everywhere incl. blunt base (D&M use an ACM/
   base-pressure empirical for the lee of blunt parts; rear cap is small).
6. Aero DB Mach grid tops at 25 vs entry-start M ≈ 27.6 (clipped; optional
   extension to 30).
7. Trajectory compared on M ≥ 3 only (our model validity).

## 5. Digitization hygiene (applied during CSV conversion)

- Sort by θ where PlotDigitizer exported out of order (TW M=10, TC M=5 are
  descending; TW M=3 has two appended out-of-sequence points).
- Trim endpoints grazing θ > 90° (90.25–90.31°) to 90°.
- One tab-separated row in the original → parser accepts , and \t.
- Provenance column on every row (`lim` | `wind_tunnel`, mach) — "CL M=10"
  is labeled wind-tunnel only (the LIM CL panel was not digitized); compare
  CL only where data exist.
- Keep full float precision as extracted; do not smooth.

## 6. Scope & order

Run order: **V0 → V1 → V2 → V3 → V4** (each layer certifies the next).
Out of scope: Space Shuttle (Fig 7.8–7.19, Table 7.2) — winged vehicle, not
in this thesis; lateral/directional coefficients (axisymmetric, β = 0);
viscous effects (Appendix A of the paper). Figures for the manuscript are
produced by the scripts into `../results/figures/`.
