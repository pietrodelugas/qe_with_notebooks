"""Plot utilities for QE EOS and Pulay-stress notebooks.

Each class stores sweep data and exposes:
  - plot_all(save=None)   — the full figure (single panel for these classes)
  - plot(save=None)       — alias for plot_all
Pass a file path to `save` to write the figure to disk instead of displaying it.
"""

import numpy as np
import matplotlib.pyplot as plt


def _show_or_save(fig, save):
    fig.tight_layout()
    if save is not None:
        fig.savefig(save, bbox_inches='tight')
        plt.close(fig)
    else:
        plt.show()


def _symlog_autozoom(ax, *arrays, linthresh=1, pad=0.15):
    """Set ylim and ticks for a symlog axis based on the data range."""
    all_vals = np.concatenate([np.asarray(a).ravel() for a in arrays])
    ymin, ymax = all_vals.min(), all_vals.max()
    margin = abs(ymax - ymin) * pad
    ylim = (ymin - margin, ymax + margin)
    cands = np.unique(
        [-c * 10**e for e in range(-2, 4) for c in (1, 2, 5)]
        + [0]
        + [c * 10**e for e in range(-2, 4) for c in (1, 2, 5)]
    )
    ticks = cands[(cands >= ylim[0]) & (cands <= ylim[1])]
    ax.set_yscale('symlog', linthresh=linthresh)
    ax.set_ylim(ylim)
    ax.set_yticks(ticks)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:g}'))


class PulayStressSweepPlot:
    """Overlay of Pulay stress residuals ΔP = P_QE − P_BM for a cutoff sweep."""

    def __init__(self, a_values, stress_residuals,
                 ylim=(-11, 0.1), yticks=(-10, -1, -0.1, 0, 0.1)):
        self.a_values = a_values
        self.stress_residuals = stress_residuals
        self.ylim = ylim
        self.yticks = yticks

    def _draw(self, ax):
        for ecut, residuals in self.stress_residuals.items():
            ax.plot(self.a_values, residuals, marker='o', label=f'ecutwfc = {ecut} Ry')
        ax.axhline(0, color='gray', lw=0.8, ls=':')
        ax.set_yscale('symlog', linthresh=1)
        ax.set_ylim(self.ylim)
        ax.set_yticks(list(self.yticks))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:g}'))
        ax.set_xlabel('Lattice parameter $a$ (Å)')
        ax.set_ylabel(r'$\Delta P = P_{\rm QE} - P_{\rm BM}$ (GPa)')
        ax.set_title('Pulay stress vs plane-wave cutoff')
        ax.legend()

    def plot_all(self, save=None):
        fig, ax = plt.subplots(figsize=(7, 4))
        self._draw(ax)
        _show_or_save(fig, save)

    plot = plot_all


class PulayStressComparisonPlot:
    """Comparison of Pulay stress residuals for standard vs ecfixed at one cutoff."""

    def __init__(self, a_values, dP_std, dP_ecf, ecutwfc):
        self.a_values = a_values
        self.dP_std = dP_std
        self.dP_ecf = dP_ecf
        self.ecutwfc = ecutwfc

    def _draw(self, ax):
        ax.plot(self.a_values, self.dP_std, marker='o',
                label=f'Standard ({self.ecutwfc} Ry)')
        ax.plot(self.a_values, self.dP_ecf, marker='s',
                label=f'ecfixed ({self.ecutwfc} Ry)')
        ax.axhline(0, color='gray', lw=0.8, ls=':')
        _symlog_autozoom(ax, self.dP_std, self.dP_ecf)
        ax.set_xlabel('Lattice parameter $a$ (Å)')
        ax.set_ylabel(r'$\Delta P = P_{\rm QE} - P_{\rm BM}$ (GPa)')
        ax.set_title(f'Pulay stress at {self.ecutwfc} Ry: standard vs ecfixed')
        ax.legend()

    def plot_all(self, save=None):
        fig, ax = plt.subplots(figsize=(7, 4))
        self._draw(ax)
        _show_or_save(fig, save)

    plot = plot_all


# ─── Backward-compatible wrappers ─────────────────────────────────────────────

def plot_pulay_stress_sweep(a_values, stress_residuals, ylim=(-11, 0.1),
                            yticks=(-10, -1, -0.1, 0, 0.1)):
    PulayStressSweepPlot(a_values, stress_residuals, ylim, yticks).plot_all()


def plot_pulay_stress_comparison(a_values, dP_std, dP_ecf, ecutwfc):
    PulayStressComparisonPlot(a_values, dP_std, dP_ecf, ecutwfc).plot_all()
