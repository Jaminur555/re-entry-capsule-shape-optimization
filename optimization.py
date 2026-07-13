"""pymoo multi-objective problem definition for capsule shape optimization."""

import numpy as np
from pymoo.core.problem import ElementwiseProblem

from .objectives import evaluate_shape


class CapsuleOptimization(ElementwiseProblem):
    """NSGA-II elementwise problem: optimize capsule shape over 3 design params.

    Variables (normalized in [0, 1]): "rn" (nose radius), "rs" (shoulder
    radius) and "r_theta" (afterbody cone half-angle).

    Objectives (3): maximize volumetric efficiency "eta_V", minimize heat
    load "Qs" and maximize ground-track range "sg". Maximization targets
    are negated and objectives are scaled to a comparable ~0.1-10 range.

    Inequality constraints (5, "g(x) <= 0" feasible): peak stagnation heat
    flux, peak shoulder heat flux, load factor, pitch stability and
    trim-on-nose margin.

    Parameters
    ----------
    m : float or None, optional
        Vehicle mass [kg]. If ``None`` (default), mass is derived per shape from
        the Apollo-based constant density (Dirkx & Mooij 2017, p. 103).
    cg_params : tuple, optional
        CG placement offsets.
    **kwargs
        Forwarded to :class:'pymoo.core.problem.ElementwiseProblem' (e.g. an
        "elementwise_runner" for multiprocessing).

    Notes
    -----
    "__getstate__"/"__setstate__" drop the (unpicklable)
    "elementwise_runner" so the problem can be sent to "multiprocessing"
    pool workers.
    """

    def __init__(self, m=None, cg_params=(0.5, 0.7), **kwargs):
        super().__init__(
            n_var=3,                                  # rn, rs, r_theta
            n_obj=3,                                  # eta_V, Qs, sg
            n_ieq_constr=5,                           # q_stag, q_shldr, n_max, pitch, nose
            xl=np.array([0.0, 0.0, 0.0]),             # lower bound (avoid exact 0)
            xu=np.array([1.0, 1.0, 1.0]),             # upper bound (avoid exact 1)
            **kwargs
        )
        self.m = m
        self.cg_params = cg_params

    def __getstate__(self):
        state = self.__dict__.copy()
        if 'elementwise_runner' in state:
            del state['elementwise_runner']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        if not hasattr(self, 'elementwise_runner'):
            self.elementwise_runner = None

    def _evaluate(self, x, out, *arg, **kwargs):
        rn, rs, r_theta = x

        try:
            res = evaluate_shape(rn, rs, r_theta,
                                 cg_params=self.cg_params,
                                 m=self.m)

            # Objectives (pymoo minimizes, so negate maximization targets).
            # Scaled so they are roughly the same order of magnitude (~0.1-10).
            f1 = -res['objectives']['eta_V']         # maximize  -> minimize (-eta_V)
            f2 = res['objectives']['Qs'] / 1e7       # minimize Qs (scaled by 1e7)
            f3 = -res['objectives']['sg'] / 1e6      # maximize  -> minimize (-sg scaled by 1e6)

            out['F'] = [f1, f2, f3]

            # ----- constraints (g(x) <= 0 to be feasible) -----
            out['G'] = [
                res['constraints']['g_q_stag'],        # <= 0
                res['constraints']['g_q_shldr'],       # <= 0
                res['constraints']['g_n_max'],         # <= 0
                res['constraints']['g_pitch_stable'],  # <= 0
                res['constraints']['g_trim_on_nose']   # <= 0
            ]

        except Exception as e:
            print(f"Evaluation Error for x = {x} : {e}")
            out['F'] = [10.0, 1e6, -10.0]              # worst possible objectives
            out['G'] = [10.0, 10.0, 10.0, 10.0, 10.0]  # all constraints violated
