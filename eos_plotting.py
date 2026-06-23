import numpy as np
import matplotlib.pyplot as plt


def _symlog_autozoom(ax, *arrays, linthresh=1, pad=0.15):
    """Set ylim and ticks for a symlog axis based on the data range."""
    all_vals = np.concatenate([np.asarray(a).ravel() for a in arrays])
    ymin, ymax = all_vals.min(), all_vals.max()
    margin = abs(ymax - ymin) * pad
    ylim = (ymin - margin, ymax + margin)
    cands = np.unique(
        [-c * 10**e for e in range(-2, 4) for c in (1, 2, 5)]
        + [0]
        + [ c * 10**e for e in range(-2, 4) for c in (1, 2, 5)]
    )
    ticks = cands[(cands >= ylim[0]) & (cands <= ylim[1])]
    ax.set_yscale('symlog', linthresh=linthresh)
    ax.set_ylim(ylim)
    ax.set_yticks(ticks)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:g}'))


def plot_pulay_stress_sweep(a_values, stress_residuals, ylim=(-11, 0.1),
                            yticks=(-10, -1, -0.1, 0, 0.1)):
    """Overlay Pulay stress residuals ΔP = P_QE − P_BM for a cutoff sweep.

    Parameters
    ----------
    a_values : array_like
        Lattice parameter values (Å).
    stress_residuals : dict[int, array_like]
        {ecutwfc (Ry): ΔP array (GPa)}.
    ylim, yticks : passed to set_ylim / set_yticks.
    """
    fig, ax = plt.subplots(figsize=(7, 4))
    for ecut, residuals in stress_residuals.items():
        ax.plot(a_values, residuals, marker='o', label=f'ecutwfc = {ecut} Ry')
    ax.axhline(0, color='gray', lw=0.8, ls=':')
    ax.set_yscale('symlog', linthresh=1)
    ax.set_ylim(ylim)
    ax.set_yticks(list(yticks))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:g}'))
    ax.set_xlabel('Lattice parameter $a$ (Å)')
    ax.set_ylabel('$\\Delta P = P_{\\rm QE} - P_{\\rm BM}$ (GPa)')
    ax.set_title('Pulay stress vs plane-wave cutoff')
    ax.legend()
    fig.tight_layout()
    plt.show()


def plot_pulay_stress_comparison(a_values, dP_std, dP_ecf, ecutwfc):
    """Compare Pulay stress residuals for standard vs ecfixed at one cutoff.

    Parameters
    ----------
    a_values : array_like
        Lattice parameter values (Å).
    dP_std, dP_ecf : array_like
        ΔP = P_QE − P_BM (GPa) for the standard and ecfixed runs.
    ecutwfc : int
        Cutoff label shown in the legend and title.
    """
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(a_values, dP_std, marker='o', label=f'Standard ({ecutwfc} Ry)')
    ax.plot(a_values, dP_ecf, marker='s', label=f'ecfixed ({ecutwfc} Ry)')
    ax.axhline(0, color='gray', lw=0.8, ls=':')
    _symlog_autozoom(ax, dP_std, dP_ecf)
    ax.set_xlabel('Lattice parameter $a$ (Å)')
    ax.set_ylabel('$\\Delta P = P_{\\rm QE} - P_{\\rm BM}$ (GPa)')
    ax.set_title(f'Pulay stress at {ecutwfc} Ry: standard vs ecfixed')
    ax.legend()
    fig.tight_layout()
    plt.show()
