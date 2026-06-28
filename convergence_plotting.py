"""Plot utilities for QE convergence notebooks.

Each class stores sweep data and exposes:
  - plot_all(save=None)          — all panels in one row
  - plot_<panel>(save=None)      — one panel at a time
Pass a file path to `save` to write the figure to disk instead of displaying it.
"""

from __future__ import annotations

import matplotlib.pyplot as plt


def _show_or_save(fig, save):
    fig.tight_layout()
    if save is not None:
        fig.savefig(save, bbox_inches='tight')
        plt.close(fig)
    else:
        plt.show()


# ─── Ecut convergence ────────────────────────────────────────────────────────

class EcutConvergencePlot:
    """Cutoff convergence: energy, convergence error, wall time."""

    def __init__(self, ecut_x, ecut_e_ev_atom, dE_mev_atom, ecut_t,
                 thr_mev_atom, idx_conv_ecut=None, ecut_conv=None):
        self.ecut_x = ecut_x
        self.ecut_e_ev_atom = ecut_e_ev_atom
        self.dE_mev_atom = dE_mev_atom
        self.ecut_t = ecut_t
        self.thr_mev_atom = thr_mev_atom
        self.idx_conv_ecut = idx_conv_ecut
        self.ecut_conv = ecut_conv

    def _draw_energy(self, ax):
        ax.plot(self.ecut_x, self.ecut_e_ev_atom, 'o-')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('Total energy (eV/atom)')
        ax.set_title('Energy vs cutoff')

    def _draw_convergence(self, ax):
        ax.semilogy(self.ecut_x, self.dE_mev_atom, 'o-')
        ax.axhline(self.thr_mev_atom, ls='--', color='tab:red',
                   label=f'{self.thr_mev_atom:g} meV/atom')
        if self.idx_conv_ecut is not None and self.ecut_conv is not None:
            ax.axvline(self.ecut_conv, ls='--', color='tab:green',
                       label=f'converged ~ {self.ecut_conv} Ry')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('|ΔE| vs highest cutoff (meV/atom)')
        ax.set_title('Strict convergence criterion')
        ax.legend()

    def _draw_timing(self, ax):
        ax.plot(self.ecut_x, self.ecut_t, 'o-')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('Wall time (s)')
        ax.set_title('Timing vs cutoff')

    def plot_energy(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_energy(ax)
        _show_or_save(fig, save)

    def plot_convergence(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_convergence(ax)
        _show_or_save(fig, save)

    def plot_timing(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_timing(ax)
        _show_or_save(fig, save)

    def plot_all(self, save=None):
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        self._draw_energy(ax[0])
        self._draw_convergence(ax[1])
        self._draw_timing(ax[2])
        _show_or_save(fig, save)


# ─── K-point convergence ─────────────────────────────────────────────────────

class KConvergencePlot:
    """K-mesh convergence: energy, convergence error, wall time, irreducible k-points."""

    def __init__(self, k_x, k_e_ev_atom, dE_k_mev_atom, k_t, k_irr,
                 thr_mev_atom, idx_conv_k=None, nk_conv=None):
        self.k_x = k_x
        self.k_e_ev_atom = k_e_ev_atom
        self.dE_k_mev_atom = dE_k_mev_atom
        self.k_t = k_t
        self.k_irr = k_irr
        self.thr_mev_atom = thr_mev_atom
        self.idx_conv_k = idx_conv_k
        self.nk_conv = nk_conv

    def _draw_energy(self, ax):
        ax.plot(self.k_x, self.k_e_ev_atom, 'o-')
        ax.set_xlabel('k-mesh size n (n×n×n)')
        ax.set_ylabel('Total energy (eV/atom)')
        ax.set_title('Energy vs k-mesh')

    def _draw_convergence(self, ax):
        ax.semilogy(self.k_x, self.dE_k_mev_atom, 'o-')
        ax.axhline(self.thr_mev_atom, ls='--', color='tab:red',
                   label=f'{self.thr_mev_atom:g} meV/atom')
        if self.idx_conv_k is not None and self.nk_conv is not None:
            ax.axvline(self.nk_conv, ls='--', color='tab:green',
                       label=f'converged ~ {self.nk_conv}³')
        ax.set_xlabel('k-mesh size n (n×n×n)')
        ax.set_ylabel('|ΔE| vs densest mesh (meV/atom)')
        ax.set_title('Strict convergence criterion')
        ax.legend()

    def _draw_timing(self, ax):
        ax.plot(self.k_x, self.k_t, 'o-')
        ax.set_xlabel('k-mesh size n (n×n×n)')
        ax.set_ylabel('Wall time (s)')
        ax.set_title('Timing vs k-mesh')

    def _draw_kpoints(self, ax):
        ax.plot(self.k_x, self.k_irr, 'o-')
        ax.set_xlabel('k-mesh size n (n×n×n)')
        ax.set_ylabel('Irreducible k-points')
        ax.set_title('Symmetry-reduced k-points')

    def plot_energy(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_energy(ax)
        _show_or_save(fig, save)

    def plot_convergence(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_convergence(ax)
        _show_or_save(fig, save)

    def plot_timing(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_timing(ax)
        _show_or_save(fig, save)

    def plot_kpoints(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_kpoints(ax)
        _show_or_save(fig, save)

    def plot_all(self, save=None):
        fig, ax = plt.subplots(1, 4, figsize=(19, 4))
        self._draw_energy(ax[0])
        self._draw_convergence(ax[1])
        self._draw_timing(ax[2])
        self._draw_kpoints(ax[3])
        _show_or_save(fig, save)


# ─── Force/ecut convergence ───────────────────────────────────────────────────

class ForceEcutConvergencePlot:
    """Force convergence vs cutoff: force value, convergence error, wall time."""

    def __init__(self, ecut_x, fz, dF_mev_ang, ecut_t, threshold_mev_ang,
                 idx_conv=None, ecut_conv=None):
        self.ecut_x = ecut_x
        self.fz = fz
        self.dF_mev_ang = dF_mev_ang
        self.ecut_t = ecut_t
        self.threshold_mev_ang = threshold_mev_ang
        self.idx_conv = idx_conv
        self.ecut_conv = ecut_conv

    def _draw_force(self, ax):
        ax.plot(self.ecut_x, self.fz, 'o-')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('Fz (eV/Å)')
        ax.set_title('Force vs cutoff')

    def _draw_convergence(self, ax):
        ax.semilogy(self.ecut_x, self.dF_mev_ang, 'o-')
        ax.axhline(self.threshold_mev_ang, ls='--', color='tab:red',
                   label=f'{self.threshold_mev_ang:g} meV/Å')
        if self.idx_conv is not None and self.ecut_conv is not None:
            ax.axvline(self.ecut_conv, ls='--', color='tab:green',
                       label=f'converged ~ {self.ecut_conv} Ry')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('|ΔFz| vs highest cutoff (meV/Å)')
        ax.set_title('Convergence criterion')
        ax.legend()

    def _draw_timing(self, ax):
        ax.plot(self.ecut_x, self.ecut_t, 'o-')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('Wall time (s)')
        ax.set_title('Timing vs cutoff')

    def plot_force(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_force(ax)
        _show_or_save(fig, save)

    def plot_convergence(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_convergence(ax)
        _show_or_save(fig, save)

    def plot_timing(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_timing(ax)
        _show_or_save(fig, save)

    def plot_all(self, save=None):
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        self._draw_force(ax[0])
        self._draw_convergence(ax[1])
        self._draw_timing(ax[2])
        _show_or_save(fig, save)


# ─── Force + stress / ecut convergence ───────────────────────────────────────

class ForceStressEcutConvergencePlot:
    """Force and stress convergence vs cutoff: force, force error, stress, wall time."""

    def __init__(self, ecut_x, force_z_ev_ang, dF_mev_ang, stress_zz_kbar,
                 dS_kbar, ecut_t, force_threshold_mev_ang, stress_threshold_kbar,
                 idx_conv_both=None, ecut_conv_both=None):
        self.ecut_x = ecut_x
        self.force_z_ev_ang = force_z_ev_ang
        self.dF_mev_ang = dF_mev_ang
        self.stress_zz_kbar = stress_zz_kbar
        self.dS_kbar = dS_kbar
        self.ecut_t = ecut_t
        self.force_threshold_mev_ang = force_threshold_mev_ang
        self.stress_threshold_kbar = stress_threshold_kbar
        self.idx_conv_both = idx_conv_both
        self.ecut_conv_both = ecut_conv_both

    def _draw_force(self, ax):
        ax.plot(self.ecut_x, self.force_z_ev_ang, 'o-')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('Fz (eV/Å)')
        ax.set_title('Force vs cutoff')

    def _draw_force_convergence(self, ax):
        ax.semilogy(self.ecut_x, self.dF_mev_ang, 'o-')
        ax.axhline(self.force_threshold_mev_ang, ls='--', color='tab:red',
                   label=f'{self.force_threshold_mev_ang:g} meV/Å')
        if self.idx_conv_both is not None and self.ecut_conv_both is not None:
            ax.axvline(self.ecut_conv_both, ls='--', color='tab:green',
                       label=f'combined conv ~ {self.ecut_conv_both} Ry')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('|ΔFz| (meV/Å)')
        ax.set_title('Force convergence')
        ax.legend()

    def _draw_stress(self, ax):
        ax.plot(self.ecut_x, self.stress_zz_kbar, 'o-')
        ax.plot(self.ecut_x, self.dS_kbar, 's--', alpha=0.8)
        ax.axhline(self.stress_threshold_kbar, ls='--', color='tab:purple')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('Stress zz / |ΔStress zz| (kbar)')
        ax.set_title('Stress and stress error')

    def _draw_timing(self, ax):
        ax.plot(self.ecut_x, self.ecut_t, 'o-')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('Wall time (s)')
        ax.set_title('Timing vs cutoff')

    def plot_force(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_force(ax)
        _show_or_save(fig, save)

    def plot_force_convergence(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_force_convergence(ax)
        _show_or_save(fig, save)

    def plot_stress(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_stress(ax)
        _show_or_save(fig, save)

    def plot_timing(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_timing(ax)
        _show_or_save(fig, save)

    def plot_all(self, save=None):
        fig, ax = plt.subplots(1, 4, figsize=(20, 4))
        self._draw_force(ax[0])
        self._draw_force_convergence(ax[1])
        self._draw_stress(ax[2])
        self._draw_timing(ax[3])
        _show_or_save(fig, save)


# ─── conv_thr sweep ───────────────────────────────────────────────────────────

class ConvThrSweepPlot:
    """conv_thr sweep: force value, error metrics (SCF correction + ΔF), wall time."""

    def __init__(self, conv_thr_values, fz_ev_ang, scf_corr_ev_ang, dF_mev_ang,
                 time_s, force_threshold_mev_ang=None):
        self.conv_thr_values = conv_thr_values
        self.fz_ev_ang = fz_ev_ang
        self.scf_corr_mev = [v * 1000 for v in scf_corr_ev_ang]
        self.dF_mev_ang = dF_mev_ang
        self.time_s = time_s
        self.force_threshold_mev_ang = force_threshold_mev_ang

    def _draw_force(self, ax):
        ax.semilogx(self.conv_thr_values, self.fz_ev_ang, 'o-')
        ax.invert_xaxis()
        ax.set_xlabel('conv_thr (Ry)  [tighter →]')
        ax.set_ylabel('Fz (eV/Å)')
        ax.set_title('Force vs SCF threshold')

    def _draw_errors(self, ax):
        # dF_mev_ang[-1] is zero by construction — drop to keep log scale valid
        ax.loglog(self.conv_thr_values, self.scf_corr_mev, 's-',
                  label='SCF correction (QE estimate)')
        ax.loglog(self.conv_thr_values[:-1], self.dF_mev_ang[:-1], 'o--',
                  label='|ΔFz| vs tightest')
        if self.force_threshold_mev_ang is not None:
            ax.axhline(self.force_threshold_mev_ang, ls='--', color='tab:red',
                       label=f'{self.force_threshold_mev_ang:g} meV/Å')
        ax.invert_xaxis()
        ax.set_xlabel('conv_thr (Ry)  [tighter →]')
        ax.set_ylabel('Force error (meV/Å)')
        ax.set_title('SCF correction vs threshold')
        ax.legend(fontsize=8)

    def _draw_timing(self, ax):
        ax.semilogx(self.conv_thr_values, self.time_s, 'o-')
        ax.invert_xaxis()
        ax.set_xlabel('conv_thr (Ry)  [tighter →]')
        ax.set_ylabel('Wall time (s)')
        ax.set_title('Timing vs SCF threshold')

    def plot_force(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_force(ax)
        _show_or_save(fig, save)

    def plot_errors(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_errors(ax)
        _show_or_save(fig, save)

    def plot_timing(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_timing(ax)
        _show_or_save(fig, save)

    def plot_all(self, save=None):
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        self._draw_force(ax[0])
        self._draw_errors(ax[1])
        self._draw_timing(ax[2])
        _show_or_save(fig, save)


# ─── K-shift convergence ──────────────────────────────────────────────────────

class KShiftConvergencePlot:
    """K-point convergence via shift difference: both energies, shift error, k-point counts."""

    def __init__(self, nk_x, e_gamma_ev_atom, e_shift_ev_atom, dE_shift_mev_atom,
                 nk_irr_gamma, nk_irr_shift, thr_mev_atom,
                 idx_conv=None, nk_conv=None):
        self.nk_x = nk_x
        self.e_gamma_ev_atom = e_gamma_ev_atom
        self.e_shift_ev_atom = e_shift_ev_atom
        self.dE_shift_mev_atom = dE_shift_mev_atom
        self.nk_irr_gamma = nk_irr_gamma
        self.nk_irr_shift = nk_irr_shift
        self.thr_mev_atom = thr_mev_atom
        self.idx_conv = idx_conv
        self.nk_conv = nk_conv

    def _draw_energy(self, ax):
        ax.plot(self.nk_x, self.e_gamma_ev_atom, 'o-', label='sk=0 (Γ-centred)')
        ax.plot(self.nk_x, self.e_shift_ev_atom, 's--', label='sk=1 (shifted)')
        ax.set_xlabel('k-mesh size n (n×n×n)')
        ax.set_ylabel('Total energy (eV/atom)')
        ax.set_title('Energy vs k-mesh')
        ax.legend()

    def _draw_shift_difference(self, ax):
        ax.semilogy(self.nk_x, self.dE_shift_mev_atom, 'o-')
        ax.axhline(self.thr_mev_atom, ls='--', color='tab:red',
                   label=f'{self.thr_mev_atom:g} meV/atom')
        if self.idx_conv is not None and self.nk_conv is not None:
            ax.axvline(self.nk_conv, ls='--', color='tab:green',
                       label=f'converged ~ {self.nk_conv}³')
        ax.set_xlabel('k-mesh size n (n×n×n)')
        ax.set_ylabel('|E(sk=1) − E(sk=0)| (meV/atom)')
        ax.set_title('Shift difference (error gauge)')
        ax.legend()

    def _draw_kpoints(self, ax):
        ax.plot(self.nk_x, self.nk_irr_gamma, 'o-', label='sk=0')
        ax.plot(self.nk_x, self.nk_irr_shift, 's--', label='sk=1')
        ax.set_xlabel('k-mesh size n (n×n×n)')
        ax.set_ylabel('Irreducible k-points')
        ax.set_title('Symmetry-reduced k-points')
        ax.legend()

    def plot_energy(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_energy(ax)
        _show_or_save(fig, save)

    def plot_shift_difference(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_shift_difference(ax)
        _show_or_save(fig, save)

    def plot_kpoints(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_kpoints(ax)
        _show_or_save(fig, save)

    def plot_all(self, save=None):
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        self._draw_energy(ax[0])
        self._draw_shift_difference(ax[1])
        self._draw_kpoints(ax[2])
        _show_or_save(fig, save)


# ─── Force + stress / k-point convergence ────────────────────────────────────

class ForceStressKConvergencePlot:
    """Force and stress convergence vs k-mesh: force, force error, stress, timing, k-points."""

    def __init__(self, k_x, force_z_ev_ang, dF_mev_ang, stress_zz_kbar,
                 dS_kbar, k_t, k_irr, force_threshold_mev_ang, stress_threshold_kbar,
                 idx_conv_both=None, nk_conv_both=None):
        self.k_x = k_x
        self.force_z_ev_ang = force_z_ev_ang
        self.dF_mev_ang = dF_mev_ang
        self.stress_zz_kbar = stress_zz_kbar
        self.dS_kbar = dS_kbar
        self.k_t = k_t
        self.k_irr = k_irr
        self.force_threshold_mev_ang = force_threshold_mev_ang
        self.stress_threshold_kbar = stress_threshold_kbar
        self.idx_conv_both = idx_conv_both
        self.nk_conv_both = nk_conv_both

    def _draw_force(self, ax):
        ax.plot(self.k_x, self.force_z_ev_ang, 'o-')
        ax.set_xlabel('k-mesh size n (n×n×n)')
        ax.set_ylabel('Fz (eV/Å)')
        ax.set_title('Force vs k-mesh')

    def _draw_force_convergence(self, ax):
        ax.semilogy(self.k_x, self.dF_mev_ang, 'o-')
        ax.axhline(self.force_threshold_mev_ang, ls='--', color='tab:red',
                   label=f'{self.force_threshold_mev_ang:g} meV/Å')
        if self.idx_conv_both is not None and self.nk_conv_both is not None:
            ax.axvline(self.nk_conv_both, ls='--', color='tab:green',
                       label=f'combined conv ~ {self.nk_conv_both}³')
        ax.set_xlabel('k-mesh size n (n×n×n)')
        ax.set_ylabel('|ΔFz| (meV/Å)')
        ax.set_title('Force convergence')
        ax.legend()

    def _draw_stress(self, ax):
        ax.plot(self.k_x, self.stress_zz_kbar, 'o-')
        ax.plot(self.k_x, self.dS_kbar, 's--', alpha=0.8)
        ax.axhline(self.stress_threshold_kbar, ls='--', color='tab:purple')
        ax.set_xlabel('k-mesh size n (n×n×n)')
        ax.set_ylabel('Stress zz / |ΔStress zz| (kbar)')
        ax.set_title('Stress and stress error')

    def _draw_timing(self, ax):
        ax.plot(self.k_x, self.k_t, 'o-')
        ax.set_xlabel('k-mesh size n (n×n×n)')
        ax.set_ylabel('Wall time (s)')
        ax.set_title('Timing vs k-mesh')

    def _draw_kpoints(self, ax):
        ax.plot(self.k_x, self.k_irr, 'o-')
        ax.set_xlabel('k-mesh size n (n×n×n)')
        ax.set_ylabel('Irreducible k-points')
        ax.set_title('Symmetry-reduced k-points')

    def plot_force(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_force(ax)
        _show_or_save(fig, save)

    def plot_force_convergence(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_force_convergence(ax)
        _show_or_save(fig, save)

    def plot_stress(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_stress(ax)
        _show_or_save(fig, save)

    def plot_timing(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_timing(ax)
        _show_or_save(fig, save)

    def plot_kpoints(self, save=None):
        fig, ax = plt.subplots(figsize=(6, 4))
        self._draw_kpoints(ax)
        _show_or_save(fig, save)

    def plot_all(self, save=None):
        fig, ax = plt.subplots(1, 5, figsize=(24, 4))
        self._draw_force(ax[0])
        self._draw_force_convergence(ax[1])
        self._draw_stress(ax[2])
        self._draw_timing(ax[3])
        self._draw_kpoints(ax[4])
        _show_or_save(fig, save)


# ─── Backward-compatible wrappers ─────────────────────────────────────────────

def plot_ecut_convergence(ecut_x, ecut_e_ev_atom, dE_mev_atom, ecut_t,
                          thr_mev_atom, idx_conv_ecut=None, ecut_conv=None):
    EcutConvergencePlot(ecut_x, ecut_e_ev_atom, dE_mev_atom, ecut_t,
                        thr_mev_atom, idx_conv_ecut, ecut_conv).plot_all()


def plot_k_convergence(k_x, k_e_ev_atom, dE_k_mev_atom, k_t, k_irr,
                       thr_mev_atom, idx_conv_k=None, nk_conv=None):
    KConvergencePlot(k_x, k_e_ev_atom, dE_k_mev_atom, k_t, k_irr,
                     thr_mev_atom, idx_conv_k, nk_conv).plot_all()


def plot_force_ecut_convergence(ecut_x, fz, dF_mev_ang, ecut_t,
                                threshold_mev_ang, idx_conv=None, ecut_conv=None):
    ForceEcutConvergencePlot(ecut_x, fz, dF_mev_ang, ecut_t,
                             threshold_mev_ang, idx_conv, ecut_conv).plot_all()


def plot_force_stress_ecut_convergence(ecut_x, force_z_ev_ang, dF_mev_ang,
                                       stress_zz_kbar, dS_kbar, ecut_t,
                                       force_threshold_mev_ang, stress_threshold_kbar,
                                       idx_conv_both=None, ecut_conv_both=None):
    ForceStressEcutConvergencePlot(ecut_x, force_z_ev_ang, dF_mev_ang,
                                   stress_zz_kbar, dS_kbar, ecut_t,
                                   force_threshold_mev_ang, stress_threshold_kbar,
                                   idx_conv_both, ecut_conv_both).plot_all()


def plot_conv_thr_sweep(conv_thr_values, fz_ev_ang, scf_corr_ev_ang, dF_mev_ang,
                        time_s, force_threshold_mev_ang=None):
    ConvThrSweepPlot(conv_thr_values, fz_ev_ang, scf_corr_ev_ang, dF_mev_ang,
                     time_s, force_threshold_mev_ang).plot_all()


def plot_k_shift_convergence(nk_x, e_gamma_ev_atom, e_shift_ev_atom,
                              dE_shift_mev_atom, nk_irr_gamma, nk_irr_shift,
                              thr_mev_atom, idx_conv=None, nk_conv=None):
    KShiftConvergencePlot(nk_x, e_gamma_ev_atom, e_shift_ev_atom,
                          dE_shift_mev_atom, nk_irr_gamma, nk_irr_shift,
                          thr_mev_atom, idx_conv, nk_conv).plot_all()


def plot_force_stress_k_convergence(k_x, force_z_ev_ang, dF_mev_ang,
                                     stress_zz_kbar, dS_kbar, k_t, k_irr,
                                     force_threshold_mev_ang, stress_threshold_kbar,
                                     idx_conv_both=None, nk_conv_both=None):
    ForceStressKConvergencePlot(k_x, force_z_ev_ang, dF_mev_ang,
                                stress_zz_kbar, dS_kbar, k_t, k_irr,
                                force_threshold_mev_ang, stress_threshold_kbar,
                                idx_conv_both, nk_conv_both).plot_all()
