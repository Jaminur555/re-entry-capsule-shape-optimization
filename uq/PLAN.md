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

## Step 6 — (week 2 preview, do not build yet)

Robust formulation on the GP: E_z[f] and std_z[f] via 2000-sample MCS on the
surrogate per candidate; deterministic-LIM baseline NSGA-II; robust NSGA-II;
direct verification of selected Pareto points. Formulation choice (mean+std as
4th objective vs probabilistic constraint, Ridolfi-style) — decide with Claude
after seeing week-1 spread.

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
