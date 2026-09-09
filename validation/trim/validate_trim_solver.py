"""V4a - Trim angle of attack vs Mach vs digitized Fig 7.17a (LIM source).

Run from the repo root:  python -m validation.trim.validate_trim_solver
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import validation.apollo_case as apollo
from validation import utils
from capsule_opt import config
from capsule_opt.trim import trim_alpha

BAND = 1.5          # deg, vs the digitized D&M LIM curve


def digitized(profile, source):
    rows = utils.load_reference("apollo_entry_reference_dm2017.csv")
    m, v = [], []
    for r in rows:
        if r["profile"] == profile and r["aero_source"] == source \
                and r["indep_name"] == "mach":
            m.append(float(r["indep_value"]))
            v.append(float(r["value"]))
    o = np.argsort(m)
    return np.array(m)[o], np.array(v)[o]


def main():
    _iCD, _iCL, iCM, atr = apollo.apollo_database()

    # trim algorithm vs the database's stored trim (consistency)
    worst = 0.0
    for k, M in enumerate(config.mach_nodes):
        a_alg = trim_alpha(float(M), iCM)
        worst = max(worst, abs(a_alg - atr[k]))
    utils.report("V4a-algorithm", "trim_alpha() reproduces the database column",
                 worst < 0.05, f"max dev {worst:.4f} deg")

    # our trim curve vs the digitized D&M LIM curve
    m_ref, a_ref = digitized("alpha_trim_deg", "lim")     # negative (D&M)
    sel = m_ref >= 3.0
    m_ref, a_ref = m_ref[sel], -a_ref[sel]                # flip to our sign
    machs = np.array([M for M in config.mach_nodes
                      if m_ref.min() <= M <= m_ref.max()])
    ours = np.interp(machs, config.mach_nodes, atr)
    ref = np.interp(machs, m_ref, a_ref)
    err = np.abs(ours - ref)
    utils.report("V4a-trim-vs-DM",
                 "trim alpha(M) vs digitized D&M LIM (Fig 7.17a)",
                 bool(np.all(err <= BAND)),
                 f"max {err.max():.2f} deg at M={machs[np.argmax(err)]:.0f}")

    # figure
    utils.set_paper_style()
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(config.mach_nodes, atr, "o-", label="Present model")
    ax.plot(m_ref, a_ref, "x", color="green", ms=6,
            label="D&M LIM (digitized)")
    ax.set_xlabel("Mach number")
    ax.set_ylabel(r"$\alpha_{trim}$ [deg]")
    ax.set_title("Trim angle of attack (Fig. 7.17a)")
    ax.legend()
    utils.save_fig(fig, "v4a_trim.png")


if __name__ == "__main__":
    main()
