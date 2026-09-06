"""V4c - End-to-end Apollo entry vs D&M Table 7.1 + digitized Fig 7.17 peaks.

Run from the repo root:  python -m validation.system.validate_apollo_end_to_end
"""

import validation.apollo_case as apollo
from validation import utils
from validation.trajectory.validate_trajectory import digitized
import numpy as np

from capsule_opt import config
from capsule_opt.optimization.metrics import heat_load, load_factor

# fractional bands vs the D&M local-inclination column
BAND_FRAC = {
    "Stagnation-point heat load": 0.08,   # Table 7.1: 53.06 MJ m^-2
    "Ground track length": 0.02,          # Table 7.1: 1430.2 km
}
HEAT_PEAK_FRAC = 0.10                     # peak q-dot vs digitized curve max
N_PEAK_ABS = 0.25                         # g, peak load factor vs digitized max


def main():
    traj = apollo.run_entry()
    iCD, _iCL, _iCM, _atr = apollo.apollo_database()
    props = apollo.apollo_props()
    m_kg = apollo.ENTRY_ICS["mass_kg"]

    heat = heat_load(traj, apollo.APOLLO_RN_NORM, apollo.APOLLO_RS_NORM,
                     apollo.APOLLO_RTHETA_NORM)
    nres = load_factor(traj, apollo.APOLLO_RN_NORM, apollo.APOLLO_RS_NORM,
                       apollo.APOLLO_RTHETA_NORM, m_kg,
                       shape_props=props, interp_CD=iCD)

    ours = {
        "Stagnation-point heat load": heat["Qs"] / 1e6,      # MJ m^-2
        "Ground track length": traj["sg"] / 1e3,             # km
    }

    # Table 7.1 rows (D&M local-inclination column)
    for row in utils.load_reference("dm2017_table7_1.csv"):
        name = row["quantity"]
        ref, band_frac = float(row["local_inclination"]), BAND_FRAC[name]
        err = abs(ours[name] - ref)
        utils.report(f"V4c-{name}", f"end-to-end {name} vs D&M Table 7.1 (LIM)",
                     err <= band_frac * ref,
                     f"ours {ours[name]:.2f} vs {ref:.2f} "
                     f"({err / ref:.1%}, band {band_frac:.0%})")

    # peak quantities vs the digitized Fig 7.17 curve maxima
    _, q_ref = digitized("heat_rate_w_m2")
    _, n_ref = digitized("load_factor_g")
    q_pk, n_pk = heat["q_stag_max"], nres["n_max"]

    # D&M define the load factor as the TOTAL aerodynamic force, Eq (6.13):
    # n_tot = sqrt(D^2 + L^2)/(m g) = (D/mg) * sqrt(1 + (L/D)^2).
    # Recompute the peak like-for-like (our package metric is drag-only).
    _Mc = np.clip(np.asarray(traj["Mach"], float),
                  config.mach_nodes[0], config.mach_nodes[-1])
    _Ac = np.clip(np.asarray(traj["alpha_trim"], float),
                  config.aoa_nodes[0], config.aoa_nodes[-1])
    _CDa = np.array([float(iCD(np.array([[_Mc[k], _Ac[k]]]))[0])
                     for k in range(_Mc.size)])
    _CLa = np.array([float(_iCL(np.array([[_Mc[k], _Ac[k]]]))[0])
                     for k in range(_Mc.size)])
    n_tot_pk = float((nres["n_arr"]
                      * np.sqrt(1.0 + (_CLa / _CDa) ** 2)).max())
    utils.report("V4c-peak-heat-rate",
                 "peak stagnation heat rate vs digitized Fig 7.17 max",
                 abs(q_pk - q_ref.max()) <= HEAT_PEAK_FRAC * q_ref.max(),
                 f"ours {q_pk / 1e6:.3f} vs {q_ref.max() / 1e6:.3f} MW m^-2 "
                 f"({abs(q_pk - q_ref.max()) / q_ref.max():.1%}, band 10%)")
    utils.report("V4c-peak-load-factor",
                 "peak load factor vs digitized Fig 7.17 max (D&M Eq 6.13 total-force)",
                 abs(n_tot_pk - n_ref.max()) <= N_PEAK_ABS,
                 f"ours n_tot {n_tot_pk:.2f} (D-only {n_pk:.2f}) "
                 f"vs {n_ref.max():.2f} g (band {N_PEAK_ABS} g)")


if __name__ == "__main__":
    main()
