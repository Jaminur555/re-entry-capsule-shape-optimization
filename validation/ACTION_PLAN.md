# Action Plan — implement the D&M-based validation (V0–V4)

Companion to `README.md` (which defines *what* checks pass/fail). This file
defines *what to do, how, and in what order*, given the current repo state
and the thesis's existing Chapter-4 validation. Last updated 2026-09-03.

## 0. Baseline: what the thesis (Ch. 4) already validated

| Thesis § | Check | Result (thesis) |
|---|---|---|
| 4.2 | Aero vs Apollo NAA wind-tunnel, **modified-Newtonian only**, M = 3,6,9,10,25 (Fig 4.1) | L/D within ~10%, trim within ~2° (≈20°); drag high at high M (single-method artifact) |
| 4.3 | Heating: Chapman/Sutton-Graves + Zoby–Sullivan Reff; shoulder-ratio vs Marvin&Sinclair / Jones / Wadham (Table 4.1) | ratios within ~10%; stag. peak 0.76 MW/m² |
| 4.4 | Trajectory vs D&M's Apollo validation entry (γ = −4°, m = 4532 kg) (Fig 4.2, Table 4.2) | Qs = 58.3 MJ/m² (vs 53.06 LIM / 49.23 WT); Sg = 1330 km (vs 1430.2 / 1425.6) |
| 4.5 | Atmosphere vs US-76 incl. 86–120 km extension (Fig 3.4) | exact at layer bases; plotting-accuracy elsewhere |

The paper's validation keeps this skeleton and upgrades it: full LIM aero +
direct overlays on D&M's own digitized curves. **The MN-only "present
implementation" column of thesis Table 4.2 is cited from the thesis — no need
to resurrect modified-Newtonian-only code.**

## 1. Required code/config changes (small, surgical)

1. **DONE 2026-09-03** — `heating/heating.py`: three def-time
   `Rm=config.Rm_fixed` defaults → `Rm=None` + call-time guards
   (`effective_nose_radius`, `stagnation_heat_flux`, `shoulder_heat_flux`).
   Verified: regression exact (L/D = 0.3075 / Cps = 1.8317), override test
   passes (Reff follows config.Rm_fixed at call time).
2. **DONE 2026-09-03** — `config.py`: `mach_nodes` extended to 30
   (10 nodes); DB rebuild verified.
3. **New `validation/apollo_case.py`** — single source of truth for the
   exact-Apollo validation case (see step 3 below).
4. **New `validation/aerodynamics/validate_apollo_coefficients.py`** — the
   core V2 script; the scaffold has no file for the Fig-7.7 comparison.
5. **New `validation/heating/validate_shoulder_heating.py`** +
   `reference_data/shoulder_ratio_reference.csv` — thesis §4.3.2 / Table 4.1
   shoulder-ratio check (kept separate from `validate_chapman.py` = §4.3.1,
   mirroring the thesis subsection split).
6. **New `validation/build_reference_data.py` + `validation/utils.py`** —
   data converter (Extracted Data.txt → CSVs) and shared report helpers.
7. Fill the empty scaffold: remaining `validation/**/validate_*.py`,
   `reference_data/*.csv`, plus repo `README.md` and `requirements.txt`
   (both currently 0 bytes — content needed eventually for the paper's
   code-availability statement; not blocking).

Structure decision (settled 2026-09-03): the module-based tree
(`aerodynamics/ heating/ atmosphere/ trajectory/ trim/ system/` +
`reference_data/`) is KEPT — it mirrors both thesis Ch. 4 and D&M §7.2
organization. Only the additions above; no renaming, no section-numbered
layout. MN-only aero is NOT resurrected (thesis Table 4.2 cited instead).

## 2. Step-by-step

**Step 1 — code changes (§1 items 1–2).** Smoke-test after: default-config
Apollo-like case still L/D = 0.3075 @ M=10/α=20; with override Rm=1.956 the
`effective_nose_radius` result changes accordingly.

**Step 2 — `validation/build_reference_data.py` + `validation/utils.py`.**
**DONE 2026-09-03.** Ran converter; 6 CSVs written & spot-checked:
apollo_flight_wind_tunnel_data.csv (123 rows, all 16 curves: 4 coeffs × M∈{3,10} × {lim,WT}, sorted, monotone),
cp_methods_dm2017_digitized.csv (128 rows, 12 curves θ-sorted/clip≤90, MN(∞)=1.852, Cp ordering TW>TC>MN ✓),
zoby_sullivan_chart_digitized.csv (55 rows, 5 curves; interp at Apollo Rm/Rn=0.4168 → 0.510 vs in-code 0.502),
dm2017_table7_1.csv (2 rows), apollo_entry_reference_dm2017.csv (358 rows, 6 profiles × 2 sources + ICs parsed),
shoulder_ratio_reference.csv (9 rows, thesis Table 4.1 manual). Digitization self-check |L/D−CL/CD|max=0.010
(interior points; endpoint α-grid mismatch trimmed). Mach column writes '3' not '3.0' (%.12g) — consumers use float().
Parse `../Extracted Data.txt` → long-format CSVs in `reference_data/`:
- `cp_methods_dm2017_digitized.csv` (method, mach, theta_deg, cp)
- `apollo_flight_wind_tunnel_data.csv` (coeff, mach, source, alpha_deg, value)
- `zoby_sullivan_chart_digitized.csv` (rs_over_rm, rm_over_rn, rm_over_reff)
- `apollo_entry_reference_dm2017.csv` (profile, aero_source, indep, value)
- `dm2017_table7_1.csv` (scalar, wind_tunnel, local_inclination)
- `shoulder_ratio_reference.csv` — thesis Table 4.1 source points
  (Marvin&Sinclair M=10.5 α=0 RS∈{0.05,0.15,0.25}; Jones M=8
  α∈{0,15,30}; Wadham M=10 α∈{0,20,28}) — entered manually, not digitized.
`utils.py` holds the shared helpers every V-script uses: reference-CSV path
resolution, pass-band evaluation, and the RESULTS.md/JSON appender.
Hygiene per README §5: sort out-of-order blocks by θ, clip θ>90→90, accept
`,` and `\t`, keep full precision. Spot-check counts/ranges after the run.

**Step 3 — `validation/apollo_case.py`.** **DONE 2026-09-03.** Self-test
PASS: CG lands exactly on rcom (our z sign is FLIPPED vs D&M: z_cg = −0.1369
trims stably at +α; documented in module); trim M=10 = 21.88° (D&M LIM −22.2°),
M=3 = 25.31° (−24.8°); L/D(trim, M=10) = 0.333 (D&M ~0.335); A_ref = 12.0195.
CG offsets land INSIDE the parameterized bounds (x_centroid = 1.080, so
dx_over_L = −0.012) → standard cg_params path, NO package edits needed.
Run scripts as `python -m validation.<module>` from repo root (namespace pkg).
Gotcha fixed: default cg_params=(0.5,0.7) → dx_over_L = MIN+0.5·(MAX−MIN) =
0.025, not 0.5·(MAX−MIN) — the CG inversion must use the true default.
Contains:
- config overrides: `Rm_fixed=1.956`, `Lc_fixed=2.662`, `mach_nodes→30`;
- inverse design vars: `rn=(4.694−3)/4`, `rs=(0.196−0.02)/0.38`,
  `r_theta=(33−5)/55`;
- CG at rcom: **direct inversion** (both mappings are linear):
  `dx=(1.0367−x_centroid)/L_total`, `dz=∓0.1369/h_local(x_cg)` — sign chosen
  empirically so trim lands at +25° on the stable branch (assert in module);
- `apollo_props()`, `apollo_database()` (wraps `build_aero_database`),
  `ENTRY_ICS` (γ=−4°, χ=49.6°, m=4532 kg, h=120 km, τ=225.5°, δ=−23.75°,
  V=7.83 km/s), and the Cm conversion factor `L_total/3.9116`.
- Module self-test (`python -m validation.apollo_case`): trim ≈ 25° ± 2°,
  L/D(trim) ≈ 0.30 ± 0.03, A_ref ≈ 12.02 m².

**Steps 4–8 — run V0→V4** (each script: load CSVs, run ours, save overlay to
`../results/figures/`, append a pass/fail row to `../results/logs/` JSON and
to `validation/RESULTS.md`):
- V0 `atmosphere/validate_atmosphere.py`
- V1 `aerodynamics/validate_newtonian.py`, `validate_tangent_cone.py`,
  `validate_prandtl_meyer.py`, `validate_blend_continuity.py`
- V2 `aerodynamics/validate_apollo_coefficients.py` ← core; α 0–30° at
  M=3,10; α_DM = −α; Cm×(L_total/3.9116); overlay LIM+WT digitized.
- V3 `heating/validate_chapman.py` — Fig 3.9 cross-check of the in-code
  ZS table vs the independent digitization + reproduce thesis Table 4.1
  shoulder ratios (source values: M&S 10.5/0° rs∈{0.05,0.15,0.25} → 1.35/
  1.28/1.20; Jones 8°/α∈{0,15,30} → 1.35/1.55/2.0; Wadham 10/α∈{0,20,28} →
  1.25/1.50/1.70) + stagnation-peak sanity on the Apollo entry.
- V4 `trim/validate_trim_solver.py` (α_trim(M) vs Fig 7.17a LIM),
  `trajectory/validate_trajectory.py` (profiles, M ≥ 3, both sources),
  `system/validate_apollo_end_to_end.py` (Table 7.1 scalars + peak q̇, n).

**Step 9 — consolidate:** `validation/RESULTS.md` summary table + figures =
the skeleton of the paper's model-validation section.

## 3. Thesis issues found while reading (fix in the manuscript, not code)

1. §4.2 para: "Rn = 4.694 m, **Rm = 0.196 m**, θc = 33°" — typo: 0.196 is
   **RS**; Rm (1.956 m) is missing entirely.
2. §4.3.1 says stagnation peak **0.76 MW/m²**, §4.4 says "about 0.4 MW/m²"
   (D&M: ≈0.66). Reconcile with the V4 run and use one number.
3. §4.4 peak load factor **9.5 g** vs D&M's ≈ 4.4–4.5 g — large; likely the
   MN-only drag artifact or a units slip. V4 must resolve this before the
   paper repeats it.
4. Thesis's own "Apollo" validation ran with fixed Rm = 2, Lc = 2 (config),
   not the stated Apollo values — the paper validation uses the exact
   geometry via `apollo_case.py`; state the upgrade explicitly.
5. Sanity note: thesis Table 6.2 optimization Sg = 2930–4695 km vs validation
   Sg = 1330 km — different ICs (γ = −2° vs −4°) explains direction; one
   sentence in the paper to preempt reviewer confusion.

## 4. Open decisions (defaults recommended)

1. Exact-Apollo override (recommended) vs thesis-fixed Rm = 2 / Lc = 2.
2. `mach_nodes` → 30 (recommended yes).
3. Optional extra digitization: CL @ M=10 LIM panel (missing), and/or
   M = 6 / 9 / 25 WT panels to enrich the paper figure beyond M = 3, 10.
