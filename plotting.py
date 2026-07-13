"""Pareto-front visualization and trade-off reporting for the optimization result."""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def unpack_pareto(res):
    """Derive the physical objective arrays and 3-D axis geometry from "res".

    Converts pymoo's scaled objective matrix back to physical values
    ('eta_V', 'Qs' [J m^-2], 'sg' [m]) and the plotting-scaled values,
    plus the design-variable arrays and the shadow-wall coordinates.

    Parameters
    ----------
    res : pymoo.Result
        Optimization result with 'res.F' and 'res.X'.

    Returns
    -------
    dict
        Physical and scaled objective arrays, design variables, axis limits,
        paddings and wall positions.
    """
    F_phys = res.F.copy()
    eta_V_phys = -F_phys[:, 0]
    Qs_phys = F_phys[:, 1] * 1e7
    sg_phys = -F_phys[:, 2] * 1e6

    # plotting-scaled values
    Qs_plot, sg_plot = Qs_phys / 1e7, sg_phys / 1e6

    rn_arr, rs_arr, r_theta_arr = res.X[:, 0], res.X[:, 1], res.X[:, 2]

    x_min, x_max = Qs_plot.min(), Qs_plot.max()
    y_min, y_max = sg_plot.min(), sg_plot.max()
    z_min, z_max = eta_V_phys.min(), eta_V_phys.max()

    pad_x = (x_max - x_min) * 0.1
    pad_y = (y_max - y_min) * 0.1
    pad_z = (z_max - z_min) * 0.1

    x_wall = x_min - pad_x          # left wall
    y_wall = y_min - pad_y          # back wall
    z_wall = z_min - pad_z          # bottom wall

    return {
        'eta_V_phys': eta_V_phys, 'Qs_phys': Qs_phys, 'sg_phys': sg_phys,
        'Qs_plot': Qs_plot, 'sg_plot': sg_plot,
        'rn_arr': rn_arr, 'rs_arr': rs_arr, 'r_theta_arr': r_theta_arr,
        'x_min': x_min, 'x_max': x_max,
        'y_min': y_min, 'y_max': y_max,
        'z_min': z_min, 'z_max': z_max,
        'pad_x': pad_x, 'pad_y': pad_y, 'pad_z': pad_z,
        'x_wall': x_wall, 'y_wall': y_wall, 'z_wall': z_wall,
    }


def plot_pareto_3d(res, fname="Pareto_Front_3D.png"):
    """3-D Pareto-front scatter with projected shadow walls; saves to "fname"."""
    p = unpack_pareto(res)
    eta_V_phys, Qs_plot, sg_plot = p['eta_V_phys'], p['Qs_plot'], p['sg_plot']
    x_wall, y_wall, z_wall = p['x_wall'], p['y_wall'], p['z_wall']

    sns.set_theme(style="whitegrid", font_scale=1.2)

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')

    sc = ax.scatter(Qs_plot, sg_plot, eta_V_phys,
                    color='black',
                    s=20,
                    alpha=0.9,
                    edgecolors='k',
                    linewidths=0.5,
                    depthshade=True)
    ax.scatter(Qs_plot, sg_plot, zs=z_wall, zdir='z', color='#FFC800', s=20, alpha=0.4, edgecolors=None, depthshade=False)
    ax.scatter(x_wall, sg_plot, eta_V_phys, color='#FFC800', s=20, alpha=0.4, edgecolors=None, depthshade=False)
    ax.scatter(Qs_plot, y_wall, eta_V_phys, color='#FFC800', s=20, alpha=0.4, edgecolors=None, depthshade=False)

    ax.set_xlim(x_wall, p['x_max'] + p['pad_x'])
    ax.set_ylim(y_wall, p['y_max'] + p['pad_y'])
    ax.set_zlim(z_wall, p['z_max'] + p['pad_z'])

    ax.set_xlabel(r"Heat Load $Q_s$ [J/m²] $\times 10^7$", fontsize=10, labelpad=12)
    ax.set_ylabel(r"Ground Track $S_g$ [m] $\times 10^6$", fontsize=10, labelpad=12)
    ax.set_zlabel(r"Volumetric Efficiency ($\eta_V$ [-])", fontsize=10, labelpad=12)

    ax.set_title('Multi-Objective Pareto Front\nCapsule Shape Optimization', fontsize=16, fontweight='bold', pad=20)

    ax.view_init(elev=20, azim=60)
    ax.tick_params(axis='both', which='major', labelsize=10)

    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('lightgray')
    ax.yaxis.pane.set_edgecolor('lightgray')
    ax.zaxis.pane.set_edgecolor('lightgray')
    ax.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.show()
    return sc


def plot_shadow_pareto(res, color_val, cbar_label, cmap='copper', title=None, fname=None):
    """3-D Pareto front colored by "color_val" with projected shadow walls.

    Parameters
    ----------
    res : pymoo.Result
        Optimization result.
    color_val : ndarray
        Per-design value used to color the points.
    cbar_label : str
        Colorbar label.
    cmap : str, optional
        Matplotlib colormap name.
    title : str, optional
        Axes title.
    fname : str, optional
        If given, the figure is saved to this filename.
    """
    p = unpack_pareto(res)
    eta_V_phys, Qs_plot, sg_plot = p['eta_V_phys'], p['Qs_plot'], p['sg_plot']
    x_wall, y_wall, z_wall = p['x_wall'], p['y_wall'], p['z_wall']

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    sc1 = ax.scatter(Qs_plot, sg_plot, zs=z_wall, zdir='z', c=color_val, cmap=cmap,
                     s=20, alpha=0.4, edgecolors=None, depthshade=False)
    ax.scatter(np.full_like(Qs_plot, x_wall), sg_plot, eta_V_phys, c=color_val, cmap=cmap,
               s=20, alpha=0.4, edgecolors=None, depthshade=False)
    ax.scatter(Qs_plot, np.full_like(sg_plot, y_wall), eta_V_phys, c=color_val, cmap=cmap,
               s=20, alpha=0.4, edgecolors=None, depthshade=False)

    ax.set_xlim(x_wall, p['x_max'] + p['pad_x'])
    ax.set_ylim(y_wall, p['y_max'] + p['pad_y'])
    ax.set_zlim(z_wall, p['z_max'] + p['pad_z'])

    ax.set_xlabel(r"Heat Load $Q_s$ [J/m²] $\times 10^7$", fontsize=10, labelpad=12)
    ax.set_ylabel(r"Ground Track $S_g$ [m] $\times 10^6$", fontsize=10, labelpad=12)
    ax.set_zlabel(r"Volumetric Efficiency ($\eta_V$  [-])", fontsize=10, labelpad=12)

    if title:
        ax.set_title(title, fontsize=12, fontweight='bold', pad=20)

    ax.view_init(elev=20, azim=60)
    ax.tick_params(axis='both', which='major', labelsize=9)

    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor('0.85')
        pane.set_linewidth(0.6)

    ax.grid(True, linestyle=':', linewidth=0.6, alpha=0.4)

    cbar = fig.colorbar(sc1, ax=ax, shrink=0.55, pad=0.08, aspect=20)
    cbar.set_label(cbar_label, fontsize=11)
    cbar.ax.tick_params(labelsize=9)
    cbar.outline.set_linewidth(0.5)

    plt.tight_layout()

    if fname:
        plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.show()
    return sc1


def print_tradeoff_table(res):
    """Print the Pareto designs sorted by ascending heat load."""
    p = unpack_pareto(res)
    eta_V_phys, Qs_phys, sg_phys = p['eta_V_phys'], p['Qs_phys'], p['sg_phys']

    print("\nTrade-off Designs (Sorted by lowest Heat Load):")
    print(f"{'rn':>6} {'rs':>6} {'r_theta':>8} | {'eta_V':>7} {'Qs[MJ/m2]':>10} {'sg[km]':>7}")
    print("-" * 55)

    sorted_idx = np.argsort(Qs_phys)
    for i in sorted_idx[:len(res.X)]:
        x_opt = res.X[i]
        print(f"{x_opt[0]:6.3f} {x_opt[1]:6.3f} {x_opt[2]:8.3f} | {eta_V_phys[i]:7.4f} {Qs_phys[i]/1e6:10.2f} {sg_phys[i]/1e3:7.1f}")


def save_results(res, outdir="Final_result_25000"):
    """Persist the Pareto-front designs for the Results chapter / Appendix.

    Writes, under ``outdir``:
      pareto_designs.csv   design vars + physical objectives, sorted by heat load
      pareto_raw.npz       raw ``res.X`` and ``res.F`` (reloadable for re-plotting)
      tradeoff_table.txt   human-readable trade-off table
    """
    import os
    import csv
    os.makedirs(outdir, exist_ok=True)

    p = unpack_pareto(res)
    eta_V_phys, Qs_phys, sg_phys = p['eta_V_phys'], p['Qs_phys'], p['sg_phys']
    rn_arr, rs_arr, r_theta_arr = p['rn_arr'], p['rs_arr'], p['r_theta_arr']
    order = np.argsort(Qs_phys)

    with open(os.path.join(outdir, "pareto_designs.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rn", "rs", "r_theta", "eta_V", "Qs_J_m2", "sg_m"])
        for i in order:
            w.writerow([f"{rn_arr[i]:.6f}", f"{rs_arr[i]:.6f}", f"{r_theta_arr[i]:.6f}",
                        f"{eta_V_phys[i]:.6f}", f"{Qs_phys[i]:.1f}", f"{sg_phys[i]:.1f}"])

    np.savez(os.path.join(outdir, "pareto_raw.npz"), X=res.X, F=res.F)

    with open(os.path.join(outdir, "tradeoff_table.txt"), "w", encoding="utf-8") as f:
        f.write(f"Pareto-optimal designs ({len(res.X)} solutions), sorted by heat load\n")
        f.write(f"{'rn':>7} {'rs':>7} {'r_theta':>9} | {'eta_V':>7} {'Qs[MJ/m2]':>11} {'sg[km]':>8}\n")
        f.write("-" * 60 + "\n")
        for i in order:
            f.write(f"{rn_arr[i]:7.3f} {rs_arr[i]:7.3f} {r_theta_arr[i]:9.3f} | "
                    f"{eta_V_phys[i]:7.4f} {Qs_phys[i]/1e6:11.2f} {sg_phys[i]/1e3:8.1f}\n")

    print(f"saved: {outdir}/pareto_designs.csv, pareto_raw.npz, tradeoff_table.txt")
