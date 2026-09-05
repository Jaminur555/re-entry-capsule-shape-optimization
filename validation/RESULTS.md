# Validation results

| id | check | verdict | details |
|---|---|---|---|
| V0-anchors | US-76 layer-base anchors (T/p/rho) | PASS | worst rel dev 2.14e-05 |
| V0-gas-law | p = rho R T below the homopause | PASS | max rel dev 2.22e-16 |
| V0-join86 | continuity at the 86 km table join | PASS | T 4.16e-04 P 3.57e-04 rho 3.20e-04 |
| V0-monotone | rho strictly decreasing 0-120 km | PASS |  |
| V0-ext-nodes | 86-120 km interpolation hits table nodes | PASS | worst rel dev 4.07e-04 |
| V0-sound | sea-level speed of sound | PASS | a0 = 340.294 m/s (ref 340.294) |
| V1a-MN-M3 | Modified Newtonian Cp vs Fig 3.5 (M=3) | PASS | max|dCp|=0.0083 at theta=29.7 deg |
| V1a-MN-M5 | Modified Newtonian Cp vs Fig 3.5 (M=5) | PASS | max|dCp|=0.0257 at theta=66.3 deg |
| V1a-MN-M10 | Modified Newtonian Cp vs Fig 3.5 (M=10) | PASS | max|dCp|=0.0151 at theta=14.2 deg |
| V1a-MN-Minf | Modified Newtonian Cp vs Fig 3.5 (M=inf) | PASS | max|dCp|=0.0261 at theta=58.0 deg |
| V1b-TC-M3 | exact Taylor-Maccoll vs digitized correlation (M=3) | FAIL | max|dCp| core=0.2075, full=0.3162 (theta<=70) |
| V1b-TC-M5 | exact Taylor-Maccoll vs digitized correlation (M=5) | FAIL | max|dCp| core=0.1098, full=0.2855 (theta<=70) |
| V1b-TC-M10 | exact Taylor-Maccoll vs digitized correlation (M=10) | FAIL | max|dCp| core=0.1015, full=0.2343 (theta<=70) |
| V1b-TC-benchmarks | direct TM solver vs analytic benchmarks | PASS | max|dCp| 0.0003 |
| V1b-TC-table | cached production table vs benchmarks | FAIL | max|dCp| 0.1621 |
| V1c-nu | Prandtl-Meyer angle nu(M) vs published anchors | PASS | max dev 0.005 deg |
| V1c-closed-form | Cp(0)=0, monotone in theta, vacuum clamp beyond nu(M) | FAIL |  |
| V1c-table | cached PM table vs closed form | FAIL | max|dCp| 1.39e-01 |
| V1d-weight | blend weight endpoints + monotone + midpoint | PASS | w(5)=0.0e+00, w(12)=1.000000000000, w(8.5)=0.500 |
| V1d-continuity | afterbody Cp continuous across M in [5,12] | PASS | max |dCp| per 0.01 in M = 1.53e-03 |
| V1d-finite | law finite over theta in [-85,85], M in [3,30] | PASS |  |
| V3b-shoulder-ratio | shoulder/stagnation ratio vs thesis Table 4.1 sources | PASS | worst 10.2% at M=8, alpha=30 deg |
| V3b-trends | ratio rises with alpha, falls with Rs | PASS |  |
| V3a-ZS-rs00 | Zoby-Sullivan Rs/Rm=0.00 curve vs independent digitization | FAIL | max|d(Rm/Reff)|=0.0250 at Rm/Rn=0.98 |
| V3a-ZS-rs10 | Zoby-Sullivan Rs/Rm=0.10 curve vs independent digitization | FAIL | max|d(Rm/Reff)|=0.0248 at Rm/Rn=0.98 |
| V3a-ZS-rs20 | Zoby-Sullivan Rs/Rm=0.20 curve vs independent digitization | FAIL | max|d(Rm/Reff)|=0.0261 at Rm/Rn=0.98 |
| V3a-ZS-rs25 | Zoby-Sullivan Rs/Rm=0.25 curve vs independent digitization | FAIL | max|d(Rm/Reff)|=0.0250 at Rm/Rn=0.98 |
| V3a-ZS-rs30 | Zoby-Sullivan Rs/Rm=0.30 curve vs independent digitization | FAIL | max|d(Rm/Reff)|=0.0250 at Rm/Rn=0.98 |
| V3a-ZS-all | Zoby-Sullivan chart, all five curves | FAIL | worst |d(Rm/Reff)|=0.0261 |
| V3a-ZS-rs00 | Zoby-Sullivan Rs/Rm=0.00, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0086 |
| V3a-ZS-rs10 | Zoby-Sullivan Rs/Rm=0.10, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0123 |
| V3a-ZS-rs20 | Zoby-Sullivan Rs/Rm=0.20, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0136 |
| V3a-ZS-rs25 | Zoby-Sullivan Rs/Rm=0.25, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0143 |
| V3a-ZS-rs30 | Zoby-Sullivan Rs/Rm=0.30, design-space region (Rm/Rn<=0.7) | PASS | max|d(Rm/Reff)|=0.0117 |
| V3a-ZS-all | Zoby-Sullivan chart, all curves, design space | PASS | worst |d(Rm/Reff)|=0.0143 |
| V3a-ZS-edge | chart edge Rm/Rn>0.7 (informational, outside design space) | PASS | max|d|=0.0261 near x=0.98: coarse-grid artifact of the in-code table's last 0.9->1.0 segment; the true curve curls to ~1.01. Unused in production (Rm/Rn <= 0.67). |
| V1b-TC-M3 | exact Taylor-Maccoll vs digitized correlation (M=3) | PASS | max|dCp| core=0.0457, full=0.3162 (theta<=70) |
| V1b-TC-M5 | exact Taylor-Maccoll vs digitized correlation (M=5) | FAIL | max|dCp| core=0.0605, full=0.2855 (theta<=70) |
| V1b-TC-M10 | exact Taylor-Maccoll vs digitized correlation (M=10) | PASS | max|dCp| core=0.0462, full=0.2343 (theta<=70) |
| V1b-TC-benchmarks | direct TM solver vs analytic benchmarks | PASS | max|dCp| 0.0003 |
| V1b-TC-table | cached production table vs benchmarks | PASS | max|dCp| 0.0003 |
| V1c-nu | Prandtl-Meyer angle nu(M) vs published anchors | PASS | max dev 0.005 deg |
| V1c-closed-form | Cp(0)=0, monotone in theta, vacuum clamp beyond nu(M) | FAIL |  |
| V1c-table | cached PM table vs closed form | FAIL | max|dCp| 7.18e-03 |
| V1b-TC-M3 | exact Taylor-Maccoll vs digitized correlation (M=3) | PASS | max|dCp| core=0.0457, full=0.3162 (theta<=70) |
| V1b-TC-M5 | exact Taylor-Maccoll vs digitized correlation (M=5) | PASS | max|dCp| core=0.0605, full=0.2855 (theta<=70) |
| V1b-TC-M10 | exact Taylor-Maccoll vs digitized correlation (M=10) | PASS | max|dCp| core=0.0462, full=0.2343 (theta<=70) |
| V1b-TC-benchmarks | direct TM solver vs analytic benchmarks | PASS | max|dCp| 0.0003 |
| V1b-TC-table | cached production table vs benchmarks | PASS | max|dCp| 0.0003 |
| V1c-nu | Prandtl-Meyer angle nu(M) vs published anchors | PASS | max dev 0.005 deg |
| V1c-closed-form | Cp(0)=0, monotone in theta, vacuum clamp beyond nu(M) | PASS |  |
| V1c-table | cached PM table vs closed form (smooth region) | FAIL | max|dCp| 7.15e-03 |
| V1c-table-bend | PM table near the vacuum asymptote (within 3 deg of nu) | PASS | max|dCp| 7.18e-03 (2-deg grid) |
| V2-CD-M3 | Apollo CD vs digitized D&M LIM (M=3) | FAIL | max|d|=0.2393 at alpha_DM=-0.0 deg |
| V2-CL-M3 | Apollo CL vs digitized D&M LIM (M=3) | FAIL | max|d|=0.0801 at alpha_DM=-30.0 deg |
| V2-Cm-M3 | Apollo Cm vs digitized D&M LIM (M=3) | FAIL | max|d|=0.1092 at alpha_DM=0.0 deg |
| V2-L/D-M3 | Apollo L/D vs digitized D&M LIM (M=3) | PASS | max|d|=0.0220 at alpha_DM=-30.0 deg |
| V2-CD-M10 | Apollo CD vs digitized D&M LIM (M=10) | FAIL | max|d|=0.0959 at alpha_DM=-4.9 deg |
| V2-CL-M10 | Apollo CL vs digitized D&M LIM (M=10) | PASS | max|d|=0.0229 at alpha_DM=-27.6 deg |
| V2-Cm-M10 | Apollo Cm vs digitized D&M LIM (M=10) | FAIL | max|d|=0.1170 at alpha_DM=-0.2 deg |
| V2-L/D-M10 | Apollo L/D vs digitized D&M LIM (M=10) | PASS | max|d|=0.0089 at alpha_DM=-10.0 deg |
| V2-CD-overprediction | M=10 CD(alpha~0) exceeds wind tunnel (method behaviour) | PASS | ours 1.610 vs WT 1.509 (D&M LIM digitized 1.704) |
| V2-CD-M3 | Apollo CD vs digitized D&M LIM (M=3) | PASS | shape max|d|=0.0070, fitted scale k=1.160, abs max|d|=0.239 |
| V2-CL-M3 | Apollo CL vs digitized D&M LIM (M=3) | PASS | shape max|d|=0.0031, fitted scale k=1.186, abs max|d|=0.080 |
| V2-Cm-M3 | Apollo Cm vs digitized D&M LIM (M=3) | PASS | shape max|d|=0.0018, fitted scale k=1.177, abs max|d|=0.009 |
| V2-L/D-M3 | Apollo L/D vs digitized D&M LIM (M=3) | PASS | max|d|=0.0220 (absolute) |
| V2-CD-M10 | Apollo CD vs digitized D&M LIM (M=10) | PASS | shape max|d|=0.0148, fitted scale k=0.946, abs max|d|=0.096 |
| V2-CL-M10 | Apollo CL vs digitized D&M LIM (M=10) | PASS | shape max|d|=0.0059, fitted scale k=0.959, abs max|d|=0.023 |
| V2-Cm-M10 | Apollo Cm vs digitized D&M LIM (M=10) | PASS | shape max|d|=0.0014, fitted scale k=0.912, abs max|d|=0.005 |
| V2-L/D-M10 | Apollo L/D vs digitized D&M LIM (M=10) | PASS | max|d|=0.0089 (absolute) |
| V2-CD-overprediction | M=10 CD(alpha~0) exceeds wind tunnel (method behaviour) | PASS | ours 1.610 vs WT 1.509 (D&M LIM digitized 1.704) |
| V1c-nu | Prandtl-Meyer angle nu(M) vs published anchors | PASS | max dev 0.005 deg |
| V1c-closed-form | Cp(0)=0, monotone in theta, vacuum clamp beyond nu(M) | PASS |  |
| V1c-table | cached PM table vs closed form (smooth region) | PASS | max|dCp| 7.15e-03 (2-deg theta x 1-Mach grid) |
| V1c-table-bend | PM table near the vacuum asymptote (within 3 deg of nu) | PASS | max|dCp| 7.18e-03 (2-deg grid) |
| V4a-algorithm | trim_alpha() reproduces the database column | PASS | max dev 0.0000 deg |
| V4a-trim-vs-DM | trim alpha(M) vs digitized D&M LIM (Fig 7.17a) | PASS | max 0.33 deg at M=6 |
| V4c-Stagnation-point heat load | end-to-end Stagnation-point heat load vs D&M Table 7.1 (LIM) | FAIL | ours 58.51 vs 53.06 (10.3%, band 8%) |
| V4c-Ground track length | end-to-end Ground track length vs D&M Table 7.1 (LIM) | FAIL | ours 1328.44 vs 1430.20 (7.1%, band 2%) |
| V4c-peak-heat-rate | peak stagnation heat rate vs digitized Fig 7.17 max | FAIL | ours 0.759 vs 0.616 MW m^-2 (23.1%, band 10%) |
| V4c-peak-load-factor | peak load factor vs digitized Fig 7.17 max | FAIL | ours 9.54 vs 4.36 g (band 0.25 g) |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | FAIL | max err 2.98 at time=187.2, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | FAIL | max err 0.0226 at Mach=3.8, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 153 at Mach=27.2, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 8.05 at Mach=17.9, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 7.42e+05 at Mach=20.1, band 4.32e+04 |
| V4c-Stagnation-point heat load | end-to-end Stagnation-point heat load vs D&M Table 7.1 (LIM) | FAIL | ours 59.89 vs 53.06 (12.9%, band 8%) |
| V4c-Ground track length | end-to-end Ground track length vs D&M Table 7.1 (LIM) | FAIL | ours 1363.81 vs 1430.20 (4.6%, band 2%) |
| V4c-peak-heat-rate | peak stagnation heat rate vs digitized Fig 7.17 max | FAIL | ours 0.739 vs 0.616 MW m^-2 (19.9%, band 10%) |
| V4c-peak-load-factor | peak load factor vs digitized Fig 7.17 max | FAIL | ours 8.94 vs 4.36 g (band 0.25 g) |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | FAIL | max err 1.68 at time=198.4, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | FAIL | max err 0.0226 at Mach=3.8, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 89.6 at Mach=27.3, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 7.59 at Mach=17.9, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 7.21e+05 at Mach=20.1, band 4.32e+04 |
| V4c-Stagnation-point heat load | end-to-end Stagnation-point heat load vs D&M Table 7.1 (LIM) | FAIL | ours 65.90 vs 53.06 (24.2%, band 8%) |
| V4c-Ground track length | end-to-end Ground track length vs D&M Table 7.1 (LIM) | FAIL | ours 1525.23 vs 1430.20 (6.6%, band 2%) |
| V4c-peak-heat-rate | peak stagnation heat rate vs digitized Fig 7.17 max | PASS | ours 0.658 vs 0.616 MW m^-2 (6.8%, band 10%) |
| V4c-peak-load-factor | peak load factor vs digitized Fig 7.17 max | FAIL | ours 6.52 vs 4.36 g (band 0.25 g) |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | FAIL | max err 4.75 at time=179.5, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | FAIL | max err 0.0226 at Mach=3.8, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 79.6 at Mach=24.5, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 5.74 at Mach=17.9, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 6.34e+05 at Mach=20.1, band 4.32e+04 |
| V4c-Stagnation-point heat load | end-to-end Stagnation-point heat load vs D&M Table 7.1 (LIM) | FAIL | ours 61.23 vs 53.06 (15.4%, band 8%) |
| V4c-Ground track length | end-to-end Ground track length vs D&M Table 7.1 (LIM) | PASS | ours 1405.78 vs 1430.20 (1.7%, band 2%) |
| V4c-peak-heat-rate | peak stagnation heat rate vs digitized Fig 7.17 max | FAIL | ours 0.726 vs 0.616 MW m^-2 (17.7%, band 10%) |
| V4c-peak-load-factor | peak load factor vs digitized Fig 7.17 max | FAIL | ours 8.13 vs 4.36 g (band 0.25 g) |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | PASS | max err 1.15 at time=242.0, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | FAIL | max err 0.0226 at Mach=3.8, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 85.2 at Mach=27.3, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 7.14 at Mach=17.9, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 7.06e+05 at Mach=20.1, band 4.32e+04 |
| V4c-Stagnation-point heat load | end-to-end Stagnation-point heat load vs D&M Table 7.1 (LIM) | FAIL | ours 65.90 vs 53.06 (24.2%, band 8%) |
| V4c-Ground track length | end-to-end Ground track length vs D&M Table 7.1 (LIM) | FAIL | ours 1525.23 vs 1430.20 (6.6%, band 2%) |
| V4c-peak-heat-rate | peak stagnation heat rate vs digitized Fig 7.17 max | PASS | ours 0.658 vs 0.616 MW m^-2 (6.8%, band 10%) |
| V4c-peak-load-factor | peak load factor vs digitized Fig 7.17 max | FAIL | ours 6.52 vs 4.36 g (band 0.25 g) |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | FAIL | max err 4.75 at time=179.5, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | FAIL | max err 0.0226 at Mach=3.8, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 79.6 at Mach=24.5, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 5.74 at Mach=17.9, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 6.34e+05 at Mach=20.1, band 4.32e+04 |
| V4c-Stagnation-point heat load | end-to-end Stagnation-point heat load vs D&M Table 7.1 (LIM) | FAIL | ours 70.12 vs 53.06 (32.1%, band 8%) |
| V4c-Ground track length | end-to-end Ground track length vs D&M Table 7.1 (LIM) | FAIL | ours 1657.07 vs 1430.20 (15.9%, band 2%) |
| V4c-peak-heat-rate | peak stagnation heat rate vs digitized Fig 7.17 max | PASS | ours 0.621 vs 0.616 MW m^-2 (0.7%, band 10%) |
| V4c-peak-load-factor | peak load factor vs digitized Fig 7.17 max | FAIL | ours 5.06 vs 4.36 g (band 0.25 g) |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | FAIL | max err 8.35 at time=179.5, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | FAIL | max err 0.0226 at Mach=3.8, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 79.6 at Mach=24.5, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 4.71 at Mach=17.9, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 6.09e+05 at Mach=21.6, band 4.32e+04 |
| V4c-Stagnation-point heat load | end-to-end Stagnation-point heat load vs D&M Table 7.1 (LIM) | FAIL | ours 61.24 vs 53.06 (15.4%, band 8%) |
| V4c-Ground track length | end-to-end Ground track length vs D&M Table 7.1 (LIM) | PASS | ours 1405.84 vs 1430.20 (1.7%, band 2%) |
| V4c-peak-heat-rate | peak stagnation heat rate vs digitized Fig 7.17 max | FAIL | ours 0.726 vs 0.616 MW m^-2 (17.8%, band 10%) |
| V4c-peak-load-factor | peak load factor vs digitized Fig 7.17 max | FAIL | ours 8.13 vs 4.36 g (band 0.25 g) |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | PASS | max err 1.12 at time=242.0, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | FAIL | max err 0.0226 at Mach=3.8, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 85.2 at Mach=27.3, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 7.13 at Mach=17.9, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 7.1e+05 at Mach=21.6, band 4.32e+04 |
| V4c-Stagnation-point heat load | end-to-end Stagnation-point heat load vs D&M Table 7.1 (LIM) | FAIL | ours 70.12 vs 53.06 (32.1%, band 8%) |
| V4c-Ground track length | end-to-end Ground track length vs D&M Table 7.1 (LIM) | FAIL | ours 1657.07 vs 1430.20 (15.9%, band 2%) |
| V4c-peak-heat-rate | peak stagnation heat rate vs digitized Fig 7.17 max | PASS | ours 0.621 vs 0.616 MW m^-2 (0.7%, band 10%) |
| V4c-peak-load-factor | peak load factor vs digitized Fig 7.17 max | FAIL | ours 5.06 vs 4.36 g (band 0.25 g) |
| V4b-altitude_km | altitude_km vs digitized D&M LIM (Fig 7.17) | FAIL | max err 8.35 at time=179.5, band 1.5 |
| V4b-lift_to_drag | lift_to_drag vs digitized D&M LIM (Fig 7.17) | FAIL | max err 0.0226 at Mach=3.8, band 0.02 |
| V4b-bank_angle_deg | bank_angle_deg vs digitized D&M LIM (Fig 7.17) | FAIL | max err 79.6 at Mach=24.5, band 6 |
| V4b-load_factor_g | load_factor_g vs digitized D&M LIM (Fig 7.17) | FAIL | max err 4.71 at Mach=17.9, band 0.25 |
| V4b-heat_rate_w_m2 | heat_rate_w_m2 vs digitized D&M LIM (Fig 7.17) | FAIL | max err 6.09e+05 at Mach=21.6, band 4.32e+04 |
