"""Central configuration for the capsule shape-optimization study.

Single source of truth for every tunable numerical constant: fixed geometry,
physical constants, the US-76 atmosphere layer tables, the aero-database grid,
the constraint limits, and the normalized design-parameter bounds.

Importing this module is side-effect free apart from building the US-76 base
pressure array "P_b" at import time (a cheap, deterministic loop).
"""

import numpy as np

#=======================================================================================
#                            fixed geometry/Numericals
#=======================================================================================


Rm_fixed = 2         # fixed maximum body radius [m]
Lc_fixed = 2         # fixed afterbody cone-frustum length [m]

eps = 1e-10          # guard against divide-by-zero throughout the pipeline

# Capsule mass model: Apollo-based constant density (Dirkx & Mooij 2017, p. 103).
# Apollo command module m = 5470 kg, V = 20.02 m^3  ->  rho = 273.29 kg/m^3.
# Per-shape mass is then m = CAPSULE_DENSITY * V_capsule.

CAPSULE_DENSITY = 273.29    # [kg m^-3]          uniform capsule density assumption           


#=======================================================================================
#       physical bounds of Design-parameter (normalized -> physical via cap_params)
#=======================================================================================

RN_MAX, RN_MIN           = 7.0, 3.0       # nose-radius         Rn      [m]
RS_MAX, RS_MIN           = 0.4, 0.02      # shoulder-radius     Rs      [m]
THETA_C_MAX, THETA_C_MIN = 60.0, 5.0      # afterbody cone half-angle [deg]

# Center-of-gravity placement bounds (offsets as fractions of length / height).
DX_OVER_L_MAX, DX_OVER_L_MIN = 0.10, -0.05   # axial CG offset / L_total
DZ_OVER_H_MAX, DZ_OVER_H_MIN = 0.10, 0.0     # radial CG offset / local height


#=======================================================================================
#                            Physical constants (US-76, Table 2)
#=======================================================================================

R_star  = 8314.32         # [J kmol^-1 K^-1]          universal gas constant  
M_air   = 28.9644         # [kg kmol^-1]              mean molecular weight of air
R_air   = R_star / M_air  # [J kg^-1 K^-1]            specific gas constant of air  
go      = 9.80665         # [m s^-2]                  standard (sea-level) gravity
gamma   = 1.4                                         # ratio of specific heats
r_earth = 6356766.0       # [m]                       Earth radius (geopotential)


#=======================================================================================
#                   US-76 atmosphere layer tables (Table 4, 0-86 km)
#=======================================================================================

# layer index:        0      1       2       3      4        5         6        7
H_b = np.array([   0.0, 11000,  20000,  32000, 47000,   51000,   71000,   84852], dtype=float)
L_b = np.array([-6.5e-3,   0.0, 1.e-3, 2.8e-3,   0.0, -2.8e-3,   -2e-3,     0.0], dtype=float)
T_b = np.array([ 288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946], dtype=float)

# Base pressures P_b[i] derived recursively from the sea-level value.
P_b = np.zeros(len(H_b))
P_b[0] = 101325.0        # sea-level pressure [Pa]

for i in range(1, 8):
    _dH = H_b[i] - H_b[i - 1]
    if abs(L_b[i - 1]) < eps:
        P_b[i] = P_b[i - 1] * np.exp(-go * _dH / (R_air * T_b[i - 1]))
    else:
        P_b[i] = P_b[i - 1] * (T_b[i] / T_b[i - 1]) ** (go / (-L_b[i - 1] * R_air))


#============================================================================================
#       Extended atmosphere table (86-120 km, geometric altitude)
#============================================================================================

# homopause (~86 km) the atmosphere is no longer well mixed: the mean molecular weight M drops
# with altitude (28.96 -> ~26 kg/kmol by 120 km), so the specific gas constant R = R*/M is not
# constant and rho cannot be recovered from P, T via the ideal-gas law with the sea-level R_air.
# The density is therefore interpolated directly from the table (it is ~10 % below P/(R_air*T) at 120 km,
# which is exactly the variable-composition effect).   

H_EXT   = np.array([86000, 86500, 87000, 87500, 88000, 88500, 89000, 89500, 90000,
                    90500, 91000, 91500, 92000, 92500, 93000, 93500, 94000, 94500,
                    95000, 95500, 96000, 96500, 97000, 97500, 98000, 98500, 99000,
                    99500, 100000, 101000, 102000, 103000, 104000, 105000, 106000,
                    107000, 108000, 109000, 110000, 111000, 112000, 113000, 114000,
                    115000, 116000, 117000, 118000, 119000, 120000], dtype=float)

T_EXT   = np.array([186.87, 186.87, 186.87, 186.87, 186.87, 186.87, 186.87, 186.87,
                    186.87, 186.87, 186.87, 186.89, 186.96, 187.08, 187.25, 187.47,
                    187.74, 188.05, 188.42, 188.84, 189.31, 189.83, 190.40, 191.04,
                    191.72, 192.47, 193.28, 194.15, 195.08, 197.16, 199.53, 202.23,
                    205.31, 208.84, 212.89, 217.63, 223.33, 230.33, 240.00, 252.00,
                    264.00, 276.00, 288.00, 300.00, 312.00, 324.00, 336.00, 348.00,
                    360.00], dtype=float)

P_EXT   = np.array([3.733800e-01, 3.416300e-01, 3.125900e-01, 2.860200e-01, 2.617300e-01,
                    2.395100e-01, 2.191900e-01, 2.006000e-01, 1.835900e-01, 1.680400e-01,
                    1.538100e-01, 1.407800e-01, 1.288700e-01, 1.179800e-01, 1.080100e-01,
                    9.889600e-02, 9.056000e-02, 8.293700e-02, 7.596600e-02, 6.959200e-02,
                    6.376500e-02, 5.844617e-02, 5.357100e-02, 4.912200e-02, 4.505700e-02,
                    4.134200e-02, 3.794800e-02, 3.484600e-02, 3.201100e-02, 2.719200e-02,
                    2.314400e-02, 1.974200e-02, 1.688200e-02, 1.447700e-02, 1.245400e-02,
                    1.075100e-02, 9.318800e-03, 8.114200e-03, 7.104200e-03, 6.261400e-03,
                    5.554700e-03, 4.970253e-03, 4.447300e-03, 4.009600e-03, 3.631200e-03,
                    3.302200e-03, 3.014400e-03, 2.761500e-03, 2.538200e-03], dtype=float)

RHO_EXT = np.array([6.958000e-06, 6.365799e-06, 5.824000e-06, 5.328000e-06, 4.875000e-06,
                    4.460000e-06, 4.081000e-06, 3.734000e-06, 3.416000e-06, 3.126000e-06,
                    2.860000e-06, 2.616000e-06, 2.393000e-06, 2.188000e-06, 2.000000e-06,
                    1.828000e-06, 1.670000e-06, 1.526000e-06, 1.393000e-06, 1.273000e-06,
                    1.162000e-06, 1.061000e-06, 9.685000e-07, 8.842000e-07, 8.071000e-07,
                    7.367000e-07, 6.725000e-07, 6.139000e-07, 5.604000e-07, 4.695000e-07,
                    3.695000e-07, 3.300000e-07, 2.769000e-07, 2.325000e-07, 1.954000e-07,
                    1.643000e-07, 1.381000e-07, 1.161000e-07, 9.708000e-08, 8.111000e-08,
                    6.838000e-08, 5.811000e-08, 4.975000e-08, 4.289000e-08, 3.720000e-08,
                    3.246000e-08, 2.847000e-08, 2.509000e-08, 2.222000e-08], dtype=float)

H_ATM_MAX = H_EXT[-1]     # 120000 m -- upper limit of the atmosphere model


#=======================================================================================
#                            Earth rotation
#=======================================================================================

omega_earth = 7.292115e-5   # Earth sidereal rotation rate [rad s^-1]


#=======================================================================================
#                            Aerodynamic database grid
#=======================================================================================

mach_nodes = np.array([3, 4, 5, 6, 8, 10, 15, 20, 25], dtype=float)
aoa_nodes = np.linspace(-5, 35, 41)

#=======================================================================================
#                               Local Inclination Methods specification
#=======================================================================================

mach_law = np.linspace(3, 25, 23)
theta_law = np.linspace(-89, 90, 180)
AERO_METHOD = 'lim'

M_LOW_MAX  = 12.0
M_HIGH_MIN = 5.0

CONE_MACHS  = np.linspace(3.0, 12.0, 19)
CONE_THETAS = np.linspace(0.0, 75, 31) 
CONE_CACHE  = "cp_cone_table.csv"

PRANDTL_MACHS  = np.linspace(3.0, 25.0, 23)   # Prandtl-Meyer table: Mach grid
PRANDTL_THETAS = np.linspace(0.0, 89.0, 45)   # expansion angle |theta| [deg]


#=======================================================================================
#                            Constraint limits
#=======================================================================================

C_Q_STAG_MAX  = 700e3    # max stagnation heat flux      [W m^-2]
C_Q_SHLDR_MAX = 1000e3   # max shoulder heat flux        [W m^-2]
C_N_MAX       = 5.0      # max deceleration load factor  [g]
C_CMA_MAX     = 0.0      # pitch-stability margin dCm/da must be <= 0 at trim

