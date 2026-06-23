"""Plot utilities for QE convergence notebooks.

These helpers keep plotting code outside notebooks so notebook cells can stay
focused on data preparation and interpretation.
"""

from __future__ import annotations

import matplotlib.pyplot as plt


def plot_ecut_convergence(
    ecut_x,
    ecut_e_ev_atom,
    dE_mev_atom,
    ecut_t,
    thr_mev_atom,
    idx_conv_ecut=None,
    ecut_conv=None,
):
    """Plot energy, convergence error, and timing for cutoff sweep."""
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))

    ax[0].plot(ecut_x, ecut_e_ev_atom, "o-")
    ax[0].set_xlabel("ecutwfc (Ry)")
    ax[0].set_ylabel("Total energy (eV/atom)")
    ax[0].set_title("Energy vs cutoff")

    ax[1].semilogy(ecut_x, dE_mev_atom, "o-")
    ax[1].axhline(
        thr_mev_atom,
        ls="--",
        color="tab:red",
        label=f"{thr_mev_atom:g} meV/atom",
    )
    if idx_conv_ecut is not None and ecut_conv is not None:
        ax[1].axvline(
            ecut_conv,
            ls="--",
            color="tab:green",
            label=f"converged ~ {ecut_conv} Ry",
        )
    ax[1].set_xlabel("ecutwfc (Ry)")
    ax[1].set_ylabel("|Delta E| vs highest cutoff (meV/atom)")
    ax[1].set_title("Strict convergence criterion")
    ax[1].legend()

    ax[2].plot(ecut_x, ecut_t, "o-")
    ax[2].set_xlabel("ecutwfc (Ry)")
    ax[2].set_ylabel("Wall time (s)")
    ax[2].set_title("Timing vs cutoff")

    plt.tight_layout()
    plt.show()


def plot_k_convergence(
    k_x,
    k_e_ev_atom,
    dE_k_mev_atom,
    k_t,
    k_irr,
    thr_mev_atom,
    idx_conv_k=None,
    nk_conv=None,
):
    """Plot energy, convergence error, timing, and irreducible k-points."""
    fig, ax = plt.subplots(1, 4, figsize=(19, 4))

    ax[0].plot(k_x, k_e_ev_atom, "o-")
    ax[0].set_xlabel("k-mesh size n (n x n x n)")
    ax[0].set_ylabel("Total energy (eV/atom)")
    ax[0].set_title("Energy vs k-mesh")

    ax[1].semilogy(k_x, dE_k_mev_atom, "o-")
    ax[1].axhline(
        thr_mev_atom,
        ls="--",
        color="tab:red",
        label=f"{thr_mev_atom:g} meV/atom",
    )
    if idx_conv_k is not None and nk_conv is not None:
        ax[1].axvline(
            nk_conv,
            ls="--",
            color="tab:green",
            label=f"converged ~ {nk_conv}^3",
        )
    ax[1].set_xlabel("k-mesh size n (n x n x n)")
    ax[1].set_ylabel("|Delta E| vs densest mesh (meV/atom)")
    ax[1].set_title("Strict convergence criterion")
    ax[1].legend()

    ax[2].plot(k_x, k_t, "o-")
    ax[2].set_xlabel("k-mesh size n (n x n x n)")
    ax[2].set_ylabel("Wall time (s)")
    ax[2].set_title("Timing vs k-mesh")

    ax[3].plot(k_x, k_irr, "o-")
    ax[3].set_xlabel("k-mesh size n (n x n x n)")
    ax[3].set_ylabel("Irreducible k-points")
    ax[3].set_title("Symmetry-reduced k-points")

    plt.tight_layout()
    plt.show()


def plot_force_ecut_convergence(
    ecut_x,
    fz,
    dF_mev_ang,
    ecut_t,
    threshold_mev_ang,
    idx_conv=None,
    ecut_conv=None,
):
    """Plot force, convergence error, and timing for a cutoff sweep."""
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))

    ax[0].plot(ecut_x, fz, "o-")
    ax[0].set_xlabel("ecutwfc (Ry)")
    ax[0].set_ylabel("Fz (eV/Ang)")
    ax[0].set_title("Force vs cutoff")

    ax[1].semilogy(ecut_x, dF_mev_ang, "o-")
    ax[1].axhline(threshold_mev_ang, ls="--", color="tab:red",
                  label=f"{threshold_mev_ang:g} meV/Ang")
    if idx_conv is not None and ecut_conv is not None:
        ax[1].axvline(ecut_conv, ls="--", color="tab:green",
                      label=f"converged ~ {ecut_conv} Ry")
    ax[1].set_xlabel("ecutwfc (Ry)")
    ax[1].set_ylabel("|Delta Fz| vs highest cutoff (meV/Ang)")
    ax[1].set_title("Convergence criterion")
    ax[1].legend()

    ax[2].plot(ecut_x, ecut_t, "o-")
    ax[2].set_xlabel("ecutwfc (Ry)")
    ax[2].set_ylabel("Wall time (s)")
    ax[2].set_title("Timing vs cutoff")

    plt.tight_layout()
    plt.show()


def plot_force_stress_ecut_convergence(
    ecut_x,
    force_z_ev_ang,
    dF_mev_ang,
    stress_zz_kbar,
    dS_kbar,
    ecut_t,
    force_threshold_mev_ang,
    stress_threshold_kbar,
    idx_conv_both=None,
    ecut_conv_both=None,
):
    """Plot force/stress convergence and timing for cutoff sweep."""
    fig, ax = plt.subplots(1, 4, figsize=(20, 4))

    ax[0].plot(ecut_x, force_z_ev_ang, "o-")
    ax[0].set_xlabel("ecutwfc (Ry)")
    ax[0].set_ylabel("Fz (eV/Ang)")
    ax[0].set_title("Force vs cutoff")

    ax[1].semilogy(ecut_x, dF_mev_ang, "o-")
    ax[1].axhline(
        force_threshold_mev_ang,
        ls="--",
        color="tab:red",
        label=f"{force_threshold_mev_ang:g} meV/Ang",
    )
    if idx_conv_both is not None and ecut_conv_both is not None:
        ax[1].axvline(
            ecut_conv_both,
            ls="--",
            color="tab:green",
            label=f"combined conv ~ {ecut_conv_both} Ry",
        )
    ax[1].set_xlabel("ecutwfc (Ry)")
    ax[1].set_ylabel("|Delta Fz| (meV/Ang)")
    ax[1].set_title("Force convergence")
    ax[1].legend()

    ax[2].plot(ecut_x, stress_zz_kbar, "o-")
    ax[2].plot(ecut_x, dS_kbar, "s--", alpha=0.8)
    ax[2].axhline(stress_threshold_kbar, ls="--", color="tab:purple")
    ax[2].set_xlabel("ecutwfc (Ry)")
    ax[2].set_ylabel("Stress zz / |Delta Stress zz| (kbar)")
    ax[2].set_title("Stress and stress error")

    ax[3].plot(ecut_x, ecut_t, "o-")
    ax[3].set_xlabel("ecutwfc (Ry)")
    ax[3].set_ylabel("Wall time (s)")
    ax[3].set_title("Timing vs cutoff")

    plt.tight_layout()
    plt.show()


def plot_conv_thr_sweep(
    conv_thr_values,
    fz_ev_ang,
    scf_corr_ev_ang,
    dF_mev_ang,
    time_s,
    force_threshold_mev_ang=None,
):
    """Plot results of a conv_thr sweep: force value, error metrics, and timing.

    Parameters
    ----------
    conv_thr_values : sequence of float
        conv_thr values used (Ry), e.g. [1e-6, 1e-7, ..., 1e-10].
    fz_ev_ang : sequence of float
        z-force on the displaced atom (eV/Ang) for each conv_thr.
    scf_corr_ev_ang : sequence of float
        Mean absolute per-component SCF correction (eV/Ang) from the verbose block.
    dF_mev_ang : sequence of float
        |ΔFz| relative to the tightest conv_thr run (meV/Ang).
    time_s : sequence of float
        Wall time (s) for each run.
    force_threshold_mev_ang : float or None
        If given, draws a horizontal threshold line on the error panel.
    """
    scf_corr_mev = [v * 1000 for v in scf_corr_ev_ang]

    fig, ax = plt.subplots(1, 3, figsize=(15, 4))

    # --- Force value ---
    ax[0].semilogx(conv_thr_values, fz_ev_ang, "o-")
    ax[0].invert_xaxis()
    ax[0].set_xlabel("conv_thr (Ry)  [tighter →]")
    ax[0].set_ylabel("Fz (eV/Å)")
    ax[0].set_title("Force vs SCF threshold")

    # --- Error metrics ---
    # dF_mev_ang[-1] is zero by construction (reference point) — drop it to keep log scale valid
    ax[1].loglog(conv_thr_values, scf_corr_mev, "s-", label="SCF correction (QE estimate)")
    ax[1].loglog(conv_thr_values[:-1], dF_mev_ang[:-1], "o--", label="|ΔFz| vs tightest")
    if force_threshold_mev_ang is not None:
        ax[1].axhline(force_threshold_mev_ang, ls="--", color="tab:red",
                      label=f"{force_threshold_mev_ang:g} meV/Å")
    ax[1].invert_xaxis()
    ax[1].set_xlabel("conv_thr (Ry)  [tighter →]")
    ax[1].set_ylabel("Force error (meV/Å)")
    ax[1].set_title("SCF correction vs threshold")
    ax[1].legend(fontsize=8)

    # --- Timing ---
    ax[2].semilogx(conv_thr_values, time_s, "o-")
    ax[2].invert_xaxis()
    ax[2].set_xlabel("conv_thr (Ry)  [tighter →]")
    ax[2].set_ylabel("Wall time (s)")
    ax[2].set_title("Timing vs SCF threshold")

    plt.tight_layout()
    plt.show()


def plot_force_stress_k_convergence(
    k_x,
    force_z_ev_ang,
    dF_mev_ang,
    stress_zz_kbar,
    dS_kbar,
    k_t,
    k_irr,
    force_threshold_mev_ang,
    stress_threshold_kbar,
    idx_conv_both=None,
    nk_conv_both=None,
):
    """Plot force/stress convergence and timing for k-point sweep."""
    fig, ax = plt.subplots(1, 5, figsize=(24, 4))

    ax[0].plot(k_x, force_z_ev_ang, "o-")
    ax[0].set_xlabel("k-mesh size n (n x n x n)")
    ax[0].set_ylabel("Fz (eV/Ang)")
    ax[0].set_title("Force vs k-mesh")

    ax[1].semilogy(k_x, dF_mev_ang, "o-")
    ax[1].axhline(
        force_threshold_mev_ang,
        ls="--",
        color="tab:red",
        label=f"{force_threshold_mev_ang:g} meV/Ang",
    )
    if idx_conv_both is not None and nk_conv_both is not None:
        ax[1].axvline(
            nk_conv_both,
            ls="--",
            color="tab:green",
            label=f"combined conv ~ {nk_conv_both}^3",
        )
    ax[1].set_xlabel("k-mesh size n (n x n x n)")
    ax[1].set_ylabel("|Delta Fz| (meV/Ang)")
    ax[1].set_title("Force convergence")
    ax[1].legend()

    ax[2].plot(k_x, stress_zz_kbar, "o-")
    ax[2].plot(k_x, dS_kbar, "s--", alpha=0.8)
    ax[2].axhline(stress_threshold_kbar, ls="--", color="tab:purple")
    ax[2].set_xlabel("k-mesh size n (n x n x n)")
    ax[2].set_ylabel("Stress zz / |Delta Stress zz| (kbar)")
    ax[2].set_title("Stress and stress error")

    ax[3].plot(k_x, k_t, "o-")
    ax[3].set_xlabel("k-mesh size n (n x n x n)")
    ax[3].set_ylabel("Wall time (s)")
    ax[3].set_title("Timing vs k-mesh")

    ax[4].plot(k_x, k_irr, "o-")
    ax[4].set_xlabel("k-mesh size n (n x n x n)")
    ax[4].set_ylabel("Irreducible k-points")
    ax[4].set_title("Symmetry-reduced k-points")

    plt.tight_layout()
    plt.show()
