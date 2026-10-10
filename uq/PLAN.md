# Phase 2 — Uncertainty module + surrogate (advisor-mode plan)

Scope decision (2026-09-19): **robust-only paper** (deterministic vs robust
comparison under uncertainty). RBDO = paper #2. Work mode: Jaminur implements,
Claude specifies & verifies. Companion to `validation/` (Phase 1, closed).

Measured basis (2026-09-19): one `evaluate_shape` call ≈ **13 s** (3 mid-design
points: 13.1/12.1/14.6 s). DOE of 500 pts ≈ 1.8 h — run overnight, rerun freely.

---

## Uncertainty model (v1 — keep it to 4 vars)

| Var | Symbol | Distribution | Rationale (cite in paper) |
|---|---|---|---|
| Entry flight-path angle | gamma0 | N(−2°, 0.5°) 1σ | Typical EI FPA dispersion, entry-corridor studies |
| Entry velocity | V0 | N(7830, 50) m/s | ~0.6% navigation/deorbit budget |
| Drag model-form scale | k_CD | N(1.00, 0.05) | V2 validation: M=10 CD 1.610 vs WT 1.523 = +5.7%; M=3 shape 6.1% |
| Lift model-form scale | k_CL | N(1.00, 0.05) | V2 absolute L/D error ≤ 0.029 ≈ 10% of L/D ≈ 0.3 → 5% per coefficient |

Deferred (do not build now): mass dispersion, atmospheric density family,
correlated k_CD/k_CL. Each is a one-line addition to the table if week 1 is fast.

Surrogate input space: **x = (rn, rs, r_theta, gamma0_deg, V0_km_s, k_CD, k_CL)** — 7 dims.
Outputs (one GP per scalar): **Qs, sg, q_stag_pk, n_tot_pk** (eta_V is pure
geometry = f(design vars) — compute exactly, never surrogate it).

## Step 1 — `capsule_opt/uncertainty/` package

New package, two files:

- `__init__.py` — empty.
- `variables.py`:

```python
from scipy import stats

# name, scipy frozen distribution, hook
UNCERTAIN_VARS = [
    ("gamma0_deg", stats.norm(-2.0, 0.5)),      # entry FPA [deg]
    ("V0",         stats.norm(7830.0, 50.0)),   # entry speed [m/s]
    ("k_CD",       stats.norm(1.00, 0.05)),     # drag model-form scale
    ("k_CL",       stats.norm(1.00, 0.05)),     # lift model-form scale
]

def sample(n, seed=42):
    """(n, 4) array of samples via independent inverse-CDF of a unit LHS."""
    from scipy.stats import qmc
    u = qmc.LatinHypercube(d=len(UNCERTAIN_VARS), seed=seed).random(n)
    import numpy as np
    return np.column_stack([dist.ppf(u[:, j]) for j, (_, dist) in enumerate(UNCERTAIN_VARS)])
```

Self-test (`python -m capsule_opt.uncertainty.variables`): sample(10000) →
means within 2% of nominal, stds within 3%, min/max sane (no |gamma0| > 5σ).

## Step 2 — hook aero scales into `evaluate_shape`

`capsule_opt/optimization/objectives.py`: add kwarg `aero_scales=(k_CD, k_CL)=None`
and wrap right after `build_aero_database`:

```python
if aero_scales is not None:
    k_CD, k_CL = aero_scales
    def _s(interp, k):
        def w(x): return k * interp(x)
        return w
    interp_CD, interp_CL = _s(interp_CD, k_CD), _s(interp_CL, k_CL)
```

Call sites are `interp_CD(np.array([[M, alpha]]))[0]` (trajectory.py:62,
metrics.py:131) — the closure matches that signature as-is. Note: alpha_trim
comes from CM (unscaled) → trim is unaffected by k_CD/k_CL, as intended.
`entry_conditions` kwarg already exists → gamma0/V0 need no code.

**Regression (must pass before anything else):** `aero_scales=(1.0, 1.0)` at
(0.5,0.5,0.5), (0.4,0.3,0.6), (0.6,0.4,0.4) reproduces today's objectives to
float equality: eta_V 0.8347657, Qs 162660012.1, sg 8386376.6 (point 1).

## Step 3 — DOE runner `uq/doe.py` (repo-root run: `python -m uq.doe`)

- Inputs: `N=500` LHS over the 7-dim space: design vars uniform [0,1]³ (own
  qmc LHS, seed 1), uncertain vars from `variables.sample(N, seed=2)`.
- Loop `evaluate_shape(rn, rs, r_theta, aero_scales=(k_CD, k_CL),
  entry_conditions={'gamma0_deg': g, 'Vo': V})`; record the 4 outputs +
  `feasible` + all constraint values (for later robust constraints).
- **Checkpoint/resume: write `uq/doe_results.csv` after every point**
  (append + flush); skip already-done rows on restart. A 2 h run must survive
  a crash. Print progress every 10 pts with ETA.
- Infeasible/non-converged points: record NaN outputs, keep the row, count them
  (report at end; >5% = tell Claude before proceeding).

## Step 4 — surrogate `uq/build_surrogate.py`

- sklearn `GaussianProcessRegressor`, one per output. Normalize inputs
  (gamma0 in deg, V0 in km/s to keep O(1); or a StandardScaler) and outputs.
  Kernel: ConstantKernel * Matern(nu=2.5, length_scale=[...], ARD) + WhiteKernel.
  `n_restarts_optimizer=10`, `random_state=0`.
- Split 450 train / 50 test (fixed seed). Report per output: RMSE_test,
  RMSE / range, R²_test (target: **R² > 0.95 and RMSE < 2% of range** for Qs,
  sg; < 3% for peaks). Save models + scaler to `uq/gp_<output>.joblib`.
- Print a Q²-style table. This is the **go/no-go checkpoint**: if any output
  misses target, do NOT tune blindly — bring the table to Claude (options:
  more DOE pts, log-transform Qs, separate GP per regime).

## Step 5 — sanity checks (before declaring week 1 done)

- Deterministic corner check: surrogate at nominal z (gamma0=−2, V0=7830,
  k=1) vs direct `evaluate_shape` at ~10 random design points → within the
  Step-4 error bands.
- Objective ranges vs thesis Table 6.2 (sg 2930–4695 km was the thesis Pareto
  envelope at gamma=−2°): our full design-space DOE will exceed it (probe point
  (0.5,0.5,0.5) gave sg ≈ 8386 km, feasible) — confirm the *envelope* contains
  the thesis range and note where it doesn't.

## Step 6 — DECIDED 2026-10-08: Option A (statistical moments)

Option B (chance constraints) rejected: reserved territory of paper #2 and
the 2025 RESS capsule RBMDO paper, and its optima sit in the GP tails.

- Objectives (4): max eta_V, min E[Qs], max E[sg], min std[Qs].
- Constraints: pitch-stability + trim-on-nose exact (z-free geometry);
  peaks as E[q] + 2·sigma_q <= limit (config: 700 kW/m², 1000 kW/m², 5 g);
  P(slow) <= 0.05 via the family classifier — carries the range-dispersion
  story, since std[sg] is not pointwise-learnable (sg-slow irreducible).
- Machinery: det and robust fronts BOTH on the surrogate (det = surrogate at
  nominal z0 = (-2°, 7830, 1, 1)); one fixed 2000-pt z-LHS (seed 42) reused
  for every candidate (common random numbers → smooth fitness); post-hoc
  fragility audit of the det front through the same engine; verification MC
  at selected designs (simulator, overnight, resumable).
- Probe-0 FIRST (`uq/probe_dispersion.py`, 5 designs × 60 sims, t_max=14400,
  stages to uq/data/probe_dispersion.csv, resume-safe): true sigma_z at fixed
  geometry vs final-gate RMSE. Verdict rule: sigma/RMSE > ~2 → plain std[Qs]
  objective; ~1 → inflate with GP predictive variance. Smoke (2 rows, D0):
  Qs 3.10, q_stag 2.43, q_shldr 1.94, n_max 0.26 — n_max dispersion may be
  BELOW surrogate error → its E+2sigma constraint would be error-dominated;
  decide after the full probe.

### Step-6 machinery + pre-flight findings (2026-10-08/09, built while probe ran)

- `uq/engine.py` — MomentEngine: CRN z-LHS (fixed seed) → per-candidate
  E/sd (plain + GP-variance-inflated via law of total variance), P(peak >
  limit), P(slow), P(capture); det point predictions at nominal z0; eta_V
  exact. Geometry constraints surrogated (exact = aero-db build per
  candidate): trim-existence classifier 98.5% OOF (in the DOE pitch-stable
  <=> trim exists — 473/2800 no-trim sentinels are exactly the unstable
  rows) + within-trim nose HGBR 98.8% sign-acc → artifacts/geom_margins.
  Self-test (`python -m uq.engine`): DOE-row point predictions 0–8% vs CSV
  (gate-consistent); mixture sg residuals larger on slow rows = known
  irreducible branch.
- `uq/optimize.py` — det/robust pymoo problems (batched), robust = 4 obj +
  E+2·sd ≤ limit peaks (sd infl|plain) + geometry [+ P_slow ≤ budget];
  two-tier CRN: n_z=500 in-loop (~0.65 s/cand under load), final front
  re-evaluated at n_z=2000. Usage in module docstring; smoke-passed both
  modes; full robust run ≈ 4–6 h at pop 150 × 200 gens.
- `uq/diagnostics/diag_robust_feasible.py` (300-design LHS scan, n_z=500):
  peak constraints NON-binding in the interior (98–100% feasible even at
  E+2σ_infl — they bind only in the sharp corners the det front will hug);
  geometry binds (~51%); **P_slow = 0.34–0.59 over the WHOLE design space →
  any small budget is infeasible (0/300 at 0.05)** — loft probability is set
  by the γ0 dispersion, not by geometry; corr(P_slow, E_sg) = +0.72,
  corr(P_slow, σ_Qs_infl) = +0.74 (loft buys range with heat-load
  dispersion). **USER DECISION 2026-10-09: P_slow is REPORT-ONLY** (no
  constraint; carried as a coloring variable + paper finding). CLI default
  'off' reflects this.

### Both fronts done + audit (2026-10-09)

Robust run: 10.5 h overnight, pop 150 × 200 gens, 100% feasible →
`uq/data/opt_robust.csv` (n_z=2000 report columns included). Det audit
(`uq/audit.py`, 426 s): σ_Qs median 30.3 MJ/m² (CV 22.8%), P_exc(q_stag)
>1% on 25% of det designs (max 6.8%), P_exc(q_shldr) max 10%, P(n_max)
≈ 0, robust-feasible 92% of det front. **Headline**: matched-performance
pairs (E_Qs +2%, sg −2%): 5268/5371 = 98% favor robust on σ_Qs, median
σ cut 43.5%; robust front worst-member P_exc halves (4.8%/6.6% vs
6.8%/10%); σ floor identical on both fronts — **16.4 MJ/m² / CV 19.9% =
irreducible dispersion floor set by the uncertainty model, not shape**.
Strict 4-way dominance rare (30/30000) as expected. Next: figures.py +
verification MC at selected designs.

### Det baseline done (2026-10-09)

`python -m uq.optimize det 200 300` — 203 s, seed 1, final pop 200 all
feasible → `uq/data/opt_det.csv` (designs + F/G + point predictions incl.
P_slow for the post-hoc fragility audit). Remaining: probe verdict (σ form
for sd_Qs objective / n_max constraint) → robust run (user-launched,
~4–6 h at pop 150 × 200 gens) → det-front audit → figures.

### Knee designs + strict fronts + det sg-inflation finding (2026-10-09, `uq/knee.py`)

Both final pops already strictly non-dominated at report precision (det
200/200 on (−η_V, Qs, −sg); robust 150/150 on (−η_V, E[Qs], −E[sg], σ_Qs,infl)
at n_z=2000, not the loop's n_z=500 F). Knees: det = row 124 (anchor-hyperplane
knee); robust = row 51 (**pseudo-weights** — anchors degenerate: one design
anchors both E[Qs] and σ_Qs; pw_dev 0.0059). Third column: robust design
nearest the det knee in the (E[Qs], E[sg]) plane (fig2's ±2% window misses the
knee by 0.05% on E[sg] — brittle threshold, nearest is parameter-free).
Table (simulator at z0 + robust lens n_z=2000): `uq/data/knee_designs.csv` +
`knee_sim_z0.csv` cache. Sim at z0: det knee (Rn 5.60 m, θc 37.0°, η_V 0.745)
Qs 128.8 MJ/m², sg 7,465 km; robust knee (Rn 4.84 m, θc 21.6°, η_V 0.857) Qs
183.9, sg 11,321 km — **the robust knee dominates the det knee nominally on
range** despite its higher E[Qs]; matched robust design: E[Qs] −5.8%,
E[sg] −2.1%, σ cut only **3.7%** — the det knee itself sits near the σ floor,
so the robust advantage concentrates in the higher-performance region (fig2).

**FINDING (paper material): det-front sg(z0) inflated by slow-GP
extrapolation.** 54/200 det designs have surrogate sg(z0) > 20,015 km
(physical great-circle max); at the det knee the slow-family GP extrapolates
(μ_s → m_s = 101,226 km, log-sd 0.39 at z0 / up to 1.2 over the LHS) giving
sg(z0) = 26,004 km vs simulator 7,465 km = **3.5× overprediction**; NSGA-II
maximized sg straight into this epistemic artifact. Qs unaffected (1.7% at
z0). Robust front clean: E_sg spans 4,117–6,393 km, no inflated designs.
All headline comparisons (audit, fig2) already use the same CRN robust lens
on both fronts — fair; disclose in the surrogate-limitations paragraph as
"point-prediction optimization exploits epistemic extrapolation; the moment
formulation self-regulates". Verification MC should include the det knee.

### `uq/verify_mc.py` written + smoke-passed (2026-10-09; run = USER)

5 designs x 300 true-simulator z-samples (z-LHS fixed at 500, seed 2026,
sliced; t_max=14400): 0 det_knee, 1 rob_knee, 2 rob_matched (knee table) +
3 det_fragile (det front max audit Pexc q_shldr — the audit headline) +
4 rob_floor (robust front min sd_Qs_infl — the 16.4 MJ/m² floor claim).
Resume-safe staging to uq/data/verify_mc.csv; summary compares MC vs
MomentEngine n_z=2000 CRN: E/sd[Qs], E[sg], P_slow, P_cap, peak E+2σ and
exceedance counts (rule of three). Smoke (2 sims) passed. Launch (from
capsule_opt, foreground+Tee, ~10–12 h):
`C:\Users\xamii\anaconda3\python.exe -m uq.verify_mc | Tee-Object -FilePath uq\logs\verify_mc_log.txt`
(300 3 = knee-only ~6–7 h; summary auto-prints at end; re-run to re-print).

### VERIFY-MC DONE (2026-10-10, 9.7 h): Qs / sigma / P_exc VERIFIED; E[sg] biased high

5×300 sims, 0 failures, 1,499/1,500 captured (sur P_cap ~0.99), ~23 s/sim.
Per design vs surrogate (n_z=2000 CRN):

| check | result |
|---|---|
| E[Qs] | −3.2% … +1.4% (≤2.6 SE) across all 5 — **VERIFIED** |
| σ_Qs (MC vs infl) | MC 15.4–65.9 vs infl 16.4–66.7 — **VERIFIED** (floor design MC 15.4 vs claimed 16.4: floor real, slightly conservative) |
| peak E+2σ | within ~1–3% everywhere — **VERIFIED** |
| **P_exc at det_fragile** | q_shldr **29/299 = 9.7% vs sur 10.0%**; q_stag 7.4% vs 6.8%; n_max 0.3% vs 0.2% — **audit headline CONFIRMED by true simulator** |
| P_slow | MC 34.7–48.8% vs sur 39.1–49.0% (within ~6 pp) — OK |
| E[sg] | **OVERPREDICTED on all 5 designs: +8.8…+28.1% (2.7–10.1 SE)** — mixture lognormal-mean inflation (same mechanism as the det-knee sg(z0) finding), now measured under the CRN |

Paper handling: absolute E[sg] carries ~10–30% surrogate bias — disclose +
report verified-true MC E[sg] for the 5 designs; matched-pair comparison
uses E[sg] on BOTH fronts with uniform bias direction → ranking robust,
absolute values conservative-high. The load-bearing claims (σ_Qs objective,
E+2σ constraints, audit P_exc) are simulator-verified. Files:
`uq/data/verify_mc.csv`, `verify_mc_summary.csv`, `uq/logs/verify_mc_log.txt`.

## Explicitly NOT in this cycle

RBDO/probabilistic constraints (paper #2) · density/mass uncertainty families ·
PCK (plain GP only unless it fails targets) · V2-prime D-3890 validation
(deferred) · shape–guidance co-design (paper #3 candidate).

---

## Progress log

**2026-09-20 — Step 4 executed; gate FAILED for Qs/sg/n_max. Step 5 blocked.**

`uq/build_surrogate.py` + `uq/sanity_check.py` written (sanity NOT run — it
needs a passing surrogate). Both fitting spaces tried (raw | log-y; metrics
always computed in original units; 450/50 split, seed 0, Matérn 5/2 ARD +
WhiteKernel, 10 restarts):

| output | raw R² / RMSE%rng | log R² / RMSE%rng | gate 0.95 / 2%(3%) |
|---|---|---|---|
| Qs | 0.818 / 7.40 | 0.860 / 6.49 | FAIL |
| sg | 0.669 / 7.97 | 0.787 / 6.40 | FAIL |
| q_stag_max | 0.952 / 2.82 | 0.960 / 2.58 | **PASS (log)** |
| q_shldr_max | 0.939 / 3.16 | 0.949 / 2.90 | fail by 0.0013 on R² |
| n_max | 0.873 / 6.37 | 0.847 / 6.99 | FAIL (raw better) |

Diagnosis (all cheap tests done, none sufficient): log-transform buys ~1 pt;
poly2 trend in log space weak (CV-R²: sg 0.16, Qs 0.57) — no universal-Kriging
shortcut; no clean regime split (sg median smooth in γ₀; the 7% long-range
pts >10,000 km are scattered γ₀×V₀×geometry combos, mostly infeasible). GP
assigns 7–20% WhiteKernel noise to a noiseless simulator → response wiggles
below 500-pt resolution. Failures are the integral quantities (Qs, sg, n_max);
local peaks pass. **Verdict: data-limited → extend the DOE** (the plan's
remaining sanctioned option).

DECISION PENDING (user stopped here for the day): extend to 1000 total
(+500 pts ≈ 2.7 h) or 1500 (+1000 ≈ 5.3 h, overnight). Command:
`python -m uq.doe 1500` — resumes at row 500 and appends (the changed-N LHS
design yields fresh distinct points; union coverage is fine for GP fitting).

Next `build_surrogate` run on the extended DOE: switch the single 50-pt test
to 5-fold CV (mean ± std R²/RMSE) for a sturdier gate. On-disk `gp_*.joblib`
are the log-space fits from today — diagnostic only, everything refits on the
bigger DOE. Useful for the paper later: peak-heating ARD says k_CD/k_CL nearly
irrelevant (ls ≈ 32–47) — peak q̇ is set by local ρ,V at peak deceleration.
