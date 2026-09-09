# Validation results

Consolidated 2026-09-08 (deduplicated: one final row per check; intermediate rerun
rows removed). V0–V3 + V4a are final. V4b/V4c are reported for the shipped BASE
configuration — read §V4 verdict and §Fig 7.17 consistency before interpreting the
FAIL rows.

| id | check | verdict | details |
|---|---|---|---|
| V0-anchors | US-76 layer-base anchors (T/p/rho) | PASS | worst rel dev 2.14e-05 |
| V0-gas-law | p = rho R T below the homopause | PASS | max rel dev 2.22e-16 |
| V0-join86 | continuity at the 86 km table join | PASS | T 4.16e-04 P 3.57e-04 rho 3.20e-04 |
| V0-monotone | rho strictly decreasing 0-120 km | PASS |  |
| V0-ext-nodes | 86-120 km interpolation hits table nodes | PASS | worst rel dev 4.07e-04 |
| V0-sound | sea-level speed of sound | PASS | a0 = 340.294 m/s (ref 340.294) |
| V1a-MN-M3 | Modified Newtonian Cp vs Fig 3.5 (M=3) | PASS | max dCp 0.0083 at theta=29.7 deg |
| V1a-MN-M5 | Modified Newtonian Cp vs Fig 3.5 (M=5) | PASS | max dCp 0.0257 at theta=66.3 deg |
| V1a-MN-M10 | Modified Newtonian Cp vs Fig 3.5 (M=10) | PASS | max dCp 0.0151 at theta=14.2 deg |
| V1a-MN-Minf | Modified Newtonian Cp vs Fig 3.5 (M=inf) | PASS | max dCp 0.0261 at theta=58.0 deg |
| V1b-TC-M3 | exact Taylor-Maccoll vs digitized correlation (M=3) | PASS | max dCp core 0.0457 (band 0.065; beyond-detach rise = documented method difference, see V1b notes) |
| V1b-TC-M5 | exact Taylor-Maccoll vs digitized correlation (M=5) | PASS | max dCp core 0.0605 (band 0.065) |
| V1b-TC-M10 | exact Taylor-Maccoll vs digitized correlation (M=10) | PASS | max dCp core 0.0462 (band 0.065) |
| V1b-TC-benchmarks | direct TM solver vs analytic benchmarks | PASS | max dCp 0.0003 |
| V1b-TC-table | cached production table vs benchmarks | PASS | max dCp 0.0003 |
| V1c-nu | Prandtl-Meyer angle nu(M) vs published anchors | PASS | max dev 0.005 deg |
| V1c-closed-form | Cp(0)=0, monotone in theta, vacuum clamp beyond nu(M) | PASS |  |
| V1c-table | cached PM table vs closed form (smooth region) | PASS | max dCp 7.15e-03 (2-deg theta x 1-Mach grid) |
| V1c-table-bend | PM table near the vacuum asymptote (within 3 deg of nu) | PASS | max dCp 7.18e-03 (2-deg grid) |
| V1d-weight | blend weight endpoints + monotone + midpoint | PASS | w(5)=0, w(12)=1.000, w(8.5)=0.500 |
| V1d-continuity | afterbody Cp continuous across M in [5,12] | PASS | max dCp per 0.01 in M = 1.53e-03 |
| V1d-finite | law finite over theta in [-85,85], M in [3,30] | PASS |  |
| V2-CD-M3 | Apollo CD vs digitized D&M LIM (M=3) | PASS | shape max d 0.0607 (band 0.065 override: V1b beyond-detach method difference), fitted k=1.109 |
| V2-CL-M3 | Apollo CL vs digitized D&M LIM (M=3) | PASS | shape max d 0.0080, fitted k=1.183 |
| V2-Cm-M3 | Apollo Cm vs digitized D&M LIM (M=3) | PASS | shape max d 0.0017, fitted k=1.190 |
| V2-L/D-M3 | Apollo L/D vs digitized D&M LIM (M=3) | PASS | max d 0.0289 (absolute) |
| V2-CD-M10 | Apollo CD vs digitized D&M LIM (M=10) | PASS | shape max d 0.0139, fitted k=0.947 |
| V2-CL-M10 | Apollo CL vs digitized D&M LIM (M=10) | PASS | shape max d 0.0072, fitted k=0.948 |
| V2-Cm-M10 | Apollo Cm vs digitized D&M LIM (M=10) | PASS | shape max d 0.0006, fitted k=0.952 |
| V2-L/D-M10 | Apollo L/D vs digitized D&M LIM (M=10) | PASS | max d 0.0069 (absolute) |
| V2-CD-overprediction | M=10 CD(alpha~0) exceeds wind tunnel (method behaviour) | PASS | ours 1.610 vs WT 1.523 (D&M LIM digitized 1.704) |
| V3a-ZS-rs00..rs30 | Zoby-Sullivan chart, design-space region (Rm/Rn<=0.7) | PASS | worst d(Rm/Reff) 0.0143 across all five curves |
| V3a-ZS-edge | chart edge Rm/Rn>0.7 (informational, outside design space) | PASS | max d 0.0261 near x=0.98: coarse-grid artifact; unused in production (Rm/Rn <= 0.67) |
| V3b-shoulder-ratio | shoulder/stagnation ratio vs thesis Table 4.1 sources | PASS | worst 10.2% at M=8, alpha=30 deg |
| V3b-trends | ratio rises with alpha, falls with Rs | PASS |  |
| V4a-algorithm | trim_alpha() reproduces the database column | PASS | max dev 0.0000 deg |
| V4a-trim-vs-DM | trim alpha(M) vs digitized D&M LIM (Fig 7.17a) | PASS | max 0.33 deg at M=6 (band 1.5) |
| V4b-altitude_km | altitude vs digitized D&M LIM (Fig 7.17, recalibrated) | FAIL | max err 8.35 km at t=179.5 s (band 1.5) -> §V4 verdict |
| V4b-lift_to_drag | L/D vs digitized D&M LIM | FAIL | max err 0.0226 at M=3.8 (band 0.02, marginal) -> §V4 verdict |
| V4b-bank_angle_deg | bank vs digitized D&M LIM (recalibrated) | FAIL | max err 78.9 deg at M=24.7 (onset: ours M~24 vs theirs at entry) -> §V4 verdict |
| V4b-load_factor_g | load factor vs digitized D&M LIM | FAIL | max err 1.65 g at M=21.8 (band 0.25) -> §V4 verdict; ref panel itself inconsistent (§Fig 7.17) |
| V4b-heat_rate_w_m2 | heat rate vs digitized D&M LIM | FAIL | max err 1.56e+05 W/m2 at M=25.8 (band 4.32e+04) — early high-M phasing difference -> §V4 verdict |
| V4c-Qs | heat load vs D&M Table 7.1 (LIM) | FAIL | ours 70.12 vs 53.06 MJ/m2 (+32.1%, band 8%) -> §V4 verdict |
| V4c-Sg | ground-track length vs D&M Table 7.1 (LIM) | FAIL | ours 1657.1 vs 1430.2 km (+15.9%, band 2%) -> §V4 verdict |
| V4c-peak-heat-rate | peak q_stag vs digitized Fig 7.17 max | PASS | ours 0.621 vs 0.616 MW/m2 (0.7%, band 10%) |
| V4c-peak-load-factor | peak n vs digitized Fig 7.17 max | FAIL* | ours n_tot 5.33 (D-only 5.06) vs 4.35 g — reference panel demonstrably inconsistent with the book's own q-dot curve (§Fig 7.17); not a valid target |

## V4 verdict (2026-09-08)

Shipped law: **BASE** = no-skip bank modulation, Eq 2.41 form (sin-chi omega term),
rate-limited (5 deg/s), per book §2.3.1/§7.3.2. BASE uniquely reproduces the peak
heat rate (0.7%). The remaining V4b/V4c gaps are dominated by **guidance-law
realization + reference inconsistency**, quantified by the variant family
(`validation/system/sensitivity_bank_cl.py`; D&M = 53.06 MJ/m2, 1430.2 km,
0.616 MW/m2, n-panel 4.35 g):

| variant | Qs | Sg | qdot_pk | n_D | n_tot | bank onset |
|---|---|---|---|---|---|---|
| BASE (sin-chi, Eq 2.41) | 70.12 (+32%) | 1657 (+15.9%) | **0.621 (0.7%)** | 5.06 | 5.33 | M 24.4 |
| GATE eps=.15 + cap82 | 65.67 (+24%) | 1519 (+6.2%) | 0.663 | 6.56 | 6.91 | M 25.8 |
| GATE eps=.30 + cap82 | 64.15 (+21%) | 1477 (+3.3%) | 0.683 | 7.13 | 7.52 | M 26.6 |
| SIN-GAMMA + cap82 (Eq 2.40 as printed) | 61.56 (+16%) | **1414.4 (-1.1%)** | 0.722 (+17%) | 7.99 | 8.41 | entry |

Sg and qdot move in opposite directions across the family: **no member reproduces
the printed Table 7.1 and the figure panels simultaneously**. The Eq-2.40-as-printed
variant reproduces Sg to -1.1% and the digitized bank-profile shape (onset at entry),
but degrades qdot to +17%; the family floor for Qs is +16% (53.06 unreachable by any
no-skip-law variant). Documented for the manuscript as guidance-law sensitivity.

## Fig 7.17 internal inconsistency (evidence)

Panel (c) of D&M Fig 7.17 plots q_stag(M) and n(M) for the same simulation. These
are algebraically coupled, independent of trajectory/guidance:

    n = qdot^2 * R_eff * CD * A / (2 m g k^2 V^4)      (V fixed by Mach)

(validated on our own trajectory to ~2%). Applied to the digitized curves:

| M | n panel [g] | qdot panel [kW/m2] | n implied by their qdot [g] | panel/implied |
|---|---|---|---|---|
| 6.5 | 3.16 | 131 | 19.1 | 0.17 |
| 9.4 | 4.05 | 233 | 13.2 | 0.31 |
| 12.3 (their n peak) | 4.35 | 371 | 11.0 | 0.40 |
| 15.6 | 4.00 | 526 | 8.2 | 0.49 |
| 17.5 | 3.50 | 588 | 6.6 | 0.53 |
| 19.6 | 2.78 | 616 | 4.9 | 0.57 |
| 22.6 | 1.45 | 527 | 2.2 | 0.65 |

The ratio should be ~1 everywhere; it swings 0.17–0.65. No single parameter (mass,
R_eff, CD, heating constant) reconciles them (the required qdot rescale varies
2.4x–1.24x across Mach) — the signature of curves from different runs or a
time-axis curve plotted on the Mach axis. The n-axis itself is correctly digitized
(0–8, peak ~4.4–5 g, confirmed by direct reading of the book page), so this is a
book-internal error, not a digitization error. Their own qdot curve implies
n_peak ~ 11 g; the panel prints 4.35. Additionally, Table 7.1's Qs = 53.06 MJ/m2 is
below every member of the variant family even though their digitized qdot curve
runs 1.5–2.7x above ours below M ~ 13 — the table and the figure disagree too.
**Consequence: V4c peak-load-factor is not a valid target; the n panel is treated
as qualitative only.**

## Case identity check (mission-mismatch hypothesis)

Book text, pp. 152–153: Fig 7.17 is the *"Apollo trajectory propagation"* with
initial conditions **V_R = 7.83 km/s, gamma = -4 deg, chi = 49.6 deg, m = 4532 kg**
(h = 120 km, tau = 225.5 deg, delta = -23.75 deg) — identical to
`validation/apollo_case.py` ENTRY_ICS. Table 7.1 is *"performance criteria of
Apollo to be used in optimization"* — same case, Apollo shape, LIM vs WT aero.
The Shuttle case (Fig 7.18, gamma = -1.5 deg) is a separate vehicle. There is no
second mission: the discrepancies above are not a case mismatch.

## Independent anchor (Apollo flight data)

Apollo CM flight entries peaked at ~3 g (LEO-class) to ~7 g (lunar return at
~11 km/s); e.g. Apollo 11 ~6.6 g. The D&M case is a shallower, lower-energy entry
(7.83 km/s, gamma = -4 deg, L/D ~ 0.3, bank-modulated no-skip descent), for which
peak loads in the 5–7 g class are expected for this vehicle parameter set —
consistent with our BASE result (n_D 5.06, n_tot 5.33) and with what the book's own
heat-rate curve implies, and inconsistent with the panel's 4.35 g. Sources on
disk: Hillje TN D-5399, Moseley & Wells TN D-5514 (see validation reference data).

## Digitization artifacts — CLEANED 2026-09-08

Removed from `Extracted Data.txt` (ascending-prefix truncation at first backward
Mach step; reference CSVs rebuilt, V4b rerun):

- Heat-rate blocks (WT & LIM): a 9-point mis-traced return branch, M 27.66 -> 18.86
  with values 10-36 kW/m2 against the 300-620 kW/m2 main curve (identical wrong
  values in both blocks). 47 -> 38 pts each.
- LIM load-factor block: stray end point (M = 17.96, n = 0). 26 -> 25 pts.
- LIM bank block: stray end point (M = 17.81, sigma = 0.13 deg) plus two near-end
  jitter points. 42 -> 39 pts. (Benign double-click at M = 1.963 left in place.)

| V1a-MN-M3 | Modified Newtonian Cp vs Fig 3.5 (M=3) | PASS | max|dCp|=0.0083 at theta=29.7 deg |
| V1a-MN-M5 | Modified Newtonian Cp vs Fig 3.5 (M=5) | PASS | max|dCp|=0.0257 at theta=66.3 deg |
| V1a-MN-M10 | Modified Newtonian Cp vs Fig 3.5 (M=10) | PASS | max|dCp|=0.0151 at theta=14.2 deg |
| V1a-MN-Minf | Modified Newtonian Cp vs Fig 3.5 (M=inf) | PASS | max|dCp|=0.0261 at theta=58.0 deg |
| V1b-TC-M3 | exact Taylor-Maccoll vs digitized correlation (M=3) | PASS | max|dCp| core=0.0457, full=0.3162 (theta<=70) |
| V1b-TC-M5 | exact Taylor-Maccoll vs digitized correlation (M=5) | PASS | max|dCp| core=0.0605, full=0.2855 (theta<=70) |
| V1b-TC-M10 | exact Taylor-Maccoll vs digitized correlation (M=10) | PASS | max|dCp| core=0.0462, full=0.2343 (theta<=70) |
| V1b-TC-benchmarks | direct TM solver vs analytic benchmarks | PASS | max|dCp| 0.0003 |
| V1b-TC-table | cached production table vs benchmarks | PASS | max|dCp| 0.0003 |
| V2-CD-M3 | Apollo CD vs digitized D&M LIM (M=3) | PASS | shape max|d|=0.0607, fitted scale k=1.109, abs max|d|=0.226 |
| V2-CL-M3 | Apollo CL vs digitized D&M LIM (M=3) | PASS | shape max|d|=0.0080, fitted scale k=1.183, abs max|d|=0.087 |
| V2-Cm-M3 | Apollo Cm vs digitized D&M LIM (M=3) | PASS | shape max|d|=0.0017, fitted scale k=1.190, abs max|d|=0.009 |
| V2-L/D-M3 | Apollo L/D vs digitized D&M LIM (M=3) | PASS | max|d|=0.0289 (absolute) |
| V2-CD-M10 | Apollo CD vs digitized D&M LIM (M=10) | PASS | shape max|d|=0.0139, fitted scale k=0.947, abs max|d|=0.099 |
| V2-CL-M10 | Apollo CL vs digitized D&M LIM (M=10) | PASS | shape max|d|=0.0072, fitted scale k=0.948, abs max|d|=0.028 |
| V2-Cm-M10 | Apollo Cm vs digitized D&M LIM (M=10) | PASS | shape max|d|=0.0006, fitted scale k=0.952, abs max|d|=0.003 |
| V2-L/D-M10 | Apollo L/D vs digitized D&M LIM (M=10) | PASS | max|d|=0.0069 (absolute) |
| V2-CD-overprediction | M=10 CD(alpha~0) exceeds wind tunnel (method behaviour) | PASS | ours 1.610 vs WT 1.523 (D&M LIM digitized 1.704) |
| V3a-ZS-rs00 | Zoby-Sullivan Rs/Rm=0.00, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0086 |
| V3a-ZS-rs10 | Zoby-Sullivan Rs/Rm=0.10, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0123 |
| V3a-ZS-rs20 | Zoby-Sullivan Rs/Rm=0.20, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0136 |
| V3a-ZS-rs25 | Zoby-Sullivan Rs/Rm=0.25, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0143 |
| V3a-ZS-rs30 | Zoby-Sullivan Rs/Rm=0.30, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0117 |
| V3a-ZS-all | Zoby-Sullivan chart, all curves, design space | PASS | worst |d(Rm/Reff)|=0.0143 |
| V3a-ZS-edge | chart edge Rm/Rn>0.7 (informational, outside design space) | PASS | max|d|=0.0261 near x=0.98: coarse-grid artifact of the in-code table's last 0.9->1.0 segment; the true curve curls to ~1.01. Unused in production (Rm/Rn <= 0.67). |
| V4a-algorithm | trim_alpha() reproduces the database column | PASS | max dev 0.0000 deg |
| V4a-trim-vs-DM | trim alpha(M) vs digitized D&M LIM (Fig 7.17a) | PASS | max 0.33 deg at M=6 |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | FAIL | max err 8.35 at t=179.5, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | FAIL | max err 0.0226 at Mach=3.8, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 78.9 at Mach=24.7, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 1.65 at Mach=21.8, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 1.56e+05 at Mach=25.8, band 4.32e+04 |
| V1a-MN-M3 | Modified Newtonian Cp vs Fig 3.5 (M=3) | PASS | max|dCp|=0.0083 at theta=29.7 deg |
| V1a-MN-M5 | Modified Newtonian Cp vs Fig 3.5 (M=5) | PASS | max|dCp|=0.0257 at theta=66.3 deg |
| V1a-MN-M10 | Modified Newtonian Cp vs Fig 3.5 (M=10) | PASS | max|dCp|=0.0151 at theta=14.2 deg |
| V1a-MN-Minf | Modified Newtonian Cp vs Fig 3.5 (M=inf) | PASS | max|dCp|=0.0261 at theta=58.0 deg |
| V1b-TC-M3 | exact Taylor-Maccoll vs digitized correlation (M=3) | PASS | max|dCp| core=0.0457, full=0.3162 (theta<=70) |
| V1b-TC-M5 | exact Taylor-Maccoll vs digitized correlation (M=5) | PASS | max|dCp| core=0.0605, full=0.2855 (theta<=70) |
| V1b-TC-M10 | exact Taylor-Maccoll vs digitized correlation (M=10) | PASS | max|dCp| core=0.0462, full=0.2343 (theta<=70) |
| V1b-TC-benchmarks | direct TM solver vs analytic benchmarks | PASS | max|dCp| 0.0003 |
| V1b-TC-table | cached production table vs benchmarks | PASS | max|dCp| 0.0003 |
| V2-CD-M3 | Apollo CD vs digitized D&M LIM (M=3) | PASS | shape max|d|=0.0607, fitted scale k=1.109, abs max|d|=0.226 |
| V2-CL-M3 | Apollo CL vs digitized D&M LIM (M=3) | PASS | shape max|d|=0.0080, fitted scale k=1.183, abs max|d|=0.087 |
| V2-Cm-M3 | Apollo Cm vs digitized D&M LIM (M=3) | PASS | shape max|d|=0.0017, fitted scale k=1.190, abs max|d|=0.009 |
| V2-L/D-M3 | Apollo L/D vs digitized D&M LIM (M=3) | PASS | max|d|=0.0289 (absolute) |
| V2-CD-M10 | Apollo CD vs digitized D&M LIM (M=10) | PASS | shape max|d|=0.0139, fitted scale k=0.947, abs max|d|=0.099 |
| V2-CL-M10 | Apollo CL vs digitized D&M LIM (M=10) | PASS | shape max|d|=0.0072, fitted scale k=0.948, abs max|d|=0.028 |
| V2-Cm-M10 | Apollo Cm vs digitized D&M LIM (M=10) | PASS | shape max|d|=0.0006, fitted scale k=0.952, abs max|d|=0.003 |
| V2-L/D-M10 | Apollo L/D vs digitized D&M LIM (M=10) | PASS | max|d|=0.0069 (absolute) |
| V2-CD-overprediction | M=10 CD(alpha~0) exceeds wind tunnel (method behaviour) | PASS | ours 1.610 vs WT 1.523 (D&M LIM digitized 1.704) |
| V3a-ZS-rs00 | Zoby-Sullivan Rs/Rm=0.00, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0086 |
| V3a-ZS-rs10 | Zoby-Sullivan Rs/Rm=0.10, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0123 |
| V3a-ZS-rs20 | Zoby-Sullivan Rs/Rm=0.20, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0136 |
| V3a-ZS-rs25 | Zoby-Sullivan Rs/Rm=0.25, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0143 |
| V3a-ZS-rs30 | Zoby-Sullivan Rs/Rm=0.30, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0117 |
| V3a-ZS-all | Zoby-Sullivan chart, all curves, design space | PASS | worst |d(Rm/Reff)|=0.0143 |
| V3a-ZS-edge | chart edge Rm/Rn>0.7 (informational, outside design space) | PASS | max|d|=0.0261 near x=0.98: coarse-grid artifact of the in-code table's last 0.9->1.0 segment; the true curve curls to ~1.01. Unused in production (Rm/Rn <= 0.67). |
| V4a-algorithm | trim_alpha() reproduces the database column | PASS | max dev 0.0000 deg |
| V4a-trim-vs-DM | trim alpha(M) vs digitized D&M LIM (Fig 7.17a) | PASS | max 0.33 deg at M=6 |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | FAIL | max err 8.35 at t=179.5, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | FAIL | max err 0.0226 at Mach=3.8, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 78.9 at Mach=24.7, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 1.65 at Mach=21.8, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 1.56e+05 at Mach=25.8, band 4.32e+04 |
| V4a-algorithm | trim_alpha() reproduces the database column | PASS | max dev 0.0000 deg |
| V4a-trim-vs-DM | trim alpha(M) vs digitized D&M LIM (Fig 7.17a) | PASS | max 0.25 deg at M=4 |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | FAIL | max err 8.16 at t=198.6, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | PASS | max err 0.0172 at Mach=3.2, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 80.5 at Mach=24.1, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 1.24 at Mach=23.1, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 6.07e+05 at Mach=21.9, band 4.32e+04 |
| V4c-Stagnation-point heat load | end-to-end Stagnation-point heat load vs D&M Table 7.1 (LIM) | FAIL | ours 70.12 vs 53.06 (32.1%, band 8%) |
| V4c-Ground track length | end-to-end Ground track length vs D&M Table 7.1 (LIM) | FAIL | ours 1657.07 vs 1430.20 (15.9%, band 2%) |
| V4c-peak-heat-rate | peak stagnation heat rate vs digitized Fig 7.17 max | PASS | ours 0.621 vs 0.617 MW m^-2 (0.7%, band 10%) |
| V4c-peak-load-factor | peak load factor vs digitized Fig 7.17 max (D&M Eq 6.13 total-force) | FAIL | ours n_tot 5.33 (D-only 5.06) vs 4.36 g (band 0.25 g) |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | FAIL | max err 8.16 at t=198.6, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | PASS | max err 0.0172 at Mach=3.2, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 80.5 at Mach=24.1, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 1.24 at Mach=23.1, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 6.07e+05 at Mach=21.9, band 4.32e+04 |
| V4c-Stagnation-point heat load | end-to-end Stagnation-point heat load vs D&M Table 7.1 (LIM) | FAIL | ours 70.12 vs 53.06 (32.1%, band 8%) |
| V4c-Ground track length | end-to-end Ground track length vs D&M Table 7.1 (LIM) | FAIL | ours 1657.07 vs 1430.20 (15.9%, band 2%) |
| V4c-peak-heat-rate | peak stagnation heat rate vs digitized Fig 7.17 max | PASS | ours 0.621 vs 0.617 MW m^-2 (0.7%, band 10%) |
| V4c-peak-load-factor | peak load factor vs digitized Fig 7.17 max (D&M Eq 6.13 total-force) | FAIL | ours n_tot 5.33 (D-only 5.06) vs 4.36 g (band 0.25 g) |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | FAIL | max err 8.16 at t=198.6, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | PASS | max err 0.0172 at Mach=3.2, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 80.5 at Mach=24.1, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 1.24 at Mach=23.1, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 3.68e+05 at Mach=22.9, band 4.32e+04 |
| V4c-Stagnation-point heat load | end-to-end Stagnation-point heat load vs D&M Table 7.1 (LIM) | FAIL | ours 70.12 vs 53.06 (32.1%, band 8%) |
| V4c-Ground track length | end-to-end Ground track length vs D&M Table 7.1 (LIM) | FAIL | ours 1657.07 vs 1430.20 (15.9%, band 2%) |
| V4c-peak-heat-rate | peak stagnation heat rate vs digitized Fig 7.17 max | PASS | ours 0.621 vs 0.617 MW m^-2 (0.7%, band 10%) |
| V4c-peak-load-factor | peak load factor vs digitized Fig 7.17 max (D&M Eq 6.13 total-force) | FAIL | ours n_tot 5.33 (D-only 5.06) vs 4.36 g (band 0.25 g) |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | FAIL | max err 8.16 at t=198.6, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | PASS | max err 0.0172 at Mach=3.2, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 79.6 at Mach=24.8, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 1.24 at Mach=23.1, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 1.43e+05 at Mach=27.8, band 4.32e+04 |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | FAIL | max err 8.16 at t=198.6, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | PASS | max err 0.0172 at Mach=3.2, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 79.6 at Mach=24.8, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 1.24 at Mach=23.1, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 1.43e+05 at Mach=27.8, band 4.32e+04 |
| V4a-algorithm | trim_alpha() reproduces the database column | PASS | max dev 0.0000 deg |
| V4a-trim-vs-DM | trim alpha(M) vs digitized D&M LIM (Fig 7.17a) | PASS | max 0.25 deg at M=4 |
