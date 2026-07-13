"""Orchestration entry point: set up and run the capsule shape optimization."""

import time
from multiprocessing import Pool

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.lhs import LHS
from pymoo.optimize import minimize

from . import plotting
from .optimization import CapsuleOptimization


def main(n_cores=4, pop_size=100, n_gen=250, seed=42):
    """Run the NSGA-II capsule shape optimization and analyze the Pareto front.

    Sets up a multiprocessing pool, builds the problem and algorithm, runs the
    optimization, then prints a trade-off table and saves the Pareto-front
    figures.

    Parameters
    ----------
    n_cores : int, optional
        Number of parallel worker processes for evaluation.
    pop_size, n_gen : int, optional
        Population size and number of generations.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    pymoo.Result
        The optimization result object.
    """
    # set up the multiprocessing pool
    pool = Pool(n_cores)
    problem = CapsuleOptimization(elementwise_runner=pool.map)

    print(f"Starting optimization on {n_cores} cores...")
    print(f"Population: {pop_size} | Generation: {n_gen} | Total Evals: {pop_size * n_gen}")

    t0 = time.perf_counter()

    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=LHS(),
        crossover=SBX(prob=0.9, eta=15),
        mutation=PM(eta=20),
        eliminate_duplicates=True
    )

    res = minimize(problem, algorithm,
                   ('n_gen', n_gen),
                   seed=seed,
                   verbose=True,
                   save_history=True)

    pool.close()
    pool.join()

    t1 = time.perf_counter()
    print("\nOptimization Finished!")
    print(f"Total time taken to done this Optimization: {t1 - t0}")
    print(f"Number of optimal trade-off solutions found: {len(res.F)}")

    # ---- analyze result ----
    if res.F is not None and len(res.F) > 0:
        plotting.save_results(res)
        p = plotting.unpack_pareto(res)
        rn_arr, rs_arr, r_theta_arr = p['rn_arr'], p['rs_arr'], p['r_theta_arr']
        Qs_phys = p['Qs_phys']

        plotting.print_tradeoff_table(res)

        plotting.plot_pareto_3d(res, fname="Pareto_Front_3D.png")

        plotting.plot_shadow_pareto(res, rn_arr, r"$r_{R_N}$ [-]", cmap='viridis',
                                    title=r"Pareto Front — Nose Radius Parameter $r_{R_N}$ [-]", fname="Pareto_Rn.png")
        plotting.plot_shadow_pareto(res, rs_arr, r"$r_{R_S}$ [-]", cmap='viridis',
                                    title=r"Pareto Front — Shoulder Radius Parameter $r_{R_S}$ [-]", fname="Pareto_Rs.png")
        plotting.plot_shadow_pareto(res, r_theta_arr, r"$r_{\theta_c}$ [-]", cmap='viridis',
                                    title=r"Pareto Front — Cone Angle Parameter $r_{\theta_c}$ [-]", fname="Pareto_r_theta.png")
        plotting.plot_shadow_pareto(res, Qs_phys / 1e7, r"Heat Load $Q_s$ [J/m²] $\times 10^7$", cmap='magma',
                                    title=r"Pareto Front — StagnationHeat Load $Q_s$ [J/m²] $\times 10^7$", fname="Pareto_Heat_Load.png")
    else:
        print("No feasible solutions found")

    return res


if __name__ == '__main__':
    main()
