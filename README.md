# capsule_opt

Hypersonic capsule shape optimization: a Local-Inclination-Method (LIM) aerothermodynamic model coupled to 3-DoF entry trajectory simulation and multi-objective optimization, validated against Dirkx & Mooij, *Aerothermodynamic Design and Assessment of Space Vehicles* (2017). Developed from a defended MSc thesis; the target is a Q1-journal paper adding uncertainty quantification and reliability-based design optimization (Phase 2).

## Package layout

```
capsule_opt/
  geometry/         # axisymmetric capsule profile + shape parameterization
  properties.py     # panel discretization, part tagging (nose / afterbody), mass properties
  aerodynamics/     # LIM pressure laws: modified Newtonian, exact Taylor-Maccoll
                    #   tangent-cone (cached table), Prandtl-Meyer; Mach-blended databases
  heating/          # Sutton-Graves stagnation + Chapman/Smithsonian shoulder heating
  atmosphere/       # US76 up to 86 km + hypersonic extension to 120 km
  trajectory/       # 3-DoF rotating-Earth entry EOM, RK4, no-skip bank guidance (rate-limited)
  trim/             # trim angle-of-attack solver
  optimization/     # objectives, constraints, metrics (heat load, load factor, ground track)
main.py             # NSGA-II shape-optimization entry point (pymoo)
validation/         # V0-V4 validation suite vs Dirkx & Mooij (2017) — see below
data/cache/         # cached Taylor-Maccoll cone Cp table (CSV)
results/            # generated figures + run logs
```

## Status — Phase 1 complete (September 2026)

Implementation and validation of the D&M-faithful LIM aerothermal model:

| Check | Content | Verdict |
|---|---|---|
| V0 | US76 + extension vs reference tables | PASS (anchors ~2e-5) |
| V1 | Local pressure laws vs book charts (Newtonian, tangent-cone, Prandtl-Meyer, Mach blend) | PASS |
| V2 | **Apollo CD/CL/Cm/L/D vs book Fig 7.7 (M=3, M=10)** — the core check | **9/9 PASS** (L/D within 0.022 absolute) |
| V3 | Heating: Zoby-Sullivan stagnation, shoulder-ratio trends vs experiment | PASS |
| V4a | Trim angle vs book Fig 7.17a | PASS (max 0.33 deg) |
| V4b/c | Full Apollo entry: trajectory profiles + Table 7.1 quantities | Reported with verdict — see `validation/RESULTS.md` |

V4 headline results on the exact Apollo entry case (V=7.83 km/s, gamma=-4 deg, m=4532 kg):
peak heat rate matches the book to **0.7%**; ground track reproducible to -1.1% within the
guidance-law family; peak load factor 5.1-5.3 g (in-family with actual Apollo flights, 5-7 g).

### Known deviations vs the book (documented, closed)

1. **Heat load / ground track / peak-g gaps trace to the book's guidance spec, not our physics.**
   The book prints two mutually inconsistent bank-law equations (2.40 vs 2.41); no member of the
   bank-law family reproduces Table 7.1 and all Fig 7.17 panels simultaneously, and the book's own
   n-panel contradicts its own heat-rate panel (identity-based ratio 0.17-0.65 over the entry).
   Evidence: `validation/RESULTS.md`.
2. **Deliberate simplifications:** tangent-cone saturates past shock detachment (their published
   correlation keeps rising - the one banded difference, CD at M=3, <= 0.06); Prandtl-Meyer on the
   lee everywhere instead of empirical ACM (<= 0.05 CD effect); no base pressure; continuous form
   of their (inverted-as-printed) blend Eq 3.48.

## Running

Requires Python >= 3.10 with `numpy scipy matplotlib pymoo` (`pip install -r requirements.txt`).

```bash
# validation suite (from repo root; each module self-contained)
python -m validation.build_reference_data          # rebuild reference CSVs from digitized data
python -m validation.aerodynamics.validate_apollo_coefficients   # V2 (core)
python -m validation.system.validate_apollo_end_to_end           # V4c
python -m validation.system.sensitivity_bank_cl                  # guidance-law variant study

# shape optimization
python main.py
```

Reference data in `validation/reference_data/` was digitized (PlotDigitizer) from Dirkx & Mooij
(2017) figures; the raw extraction is kept outside the repo alongside the source book.

## Future work

- **Phase 2 - UQ / robust / RBDO layer** (the paper's novelty): uncertain inputs (entry state,
  atmosphere, aero/heating coefficients), Kriging / polynomial-chaos-Kriging surrogate of the
  objective-constraint space, robust and reliability-based formulation of the shape optimization.
- Manuscript: validation chapter + deviations section drafted from `validation/RESULTS.md`;
  target journal Acta Astronautica.
- Optional: Hillje TN D-5399 / Moseley D-5514 wind-tunnel triangulation of the aero database.
