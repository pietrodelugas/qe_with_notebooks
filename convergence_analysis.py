"""Analysis utilities for QE convergence notebooks.

These helpers convert raw sweep results into NumPy arrays and convergence
summaries, so notebook cells can stay compact and easier to read.
"""

from __future__ import annotations

import numpy as np


def _pick_float_array(results, key):
    """Return a float array with values extracted from a list of dicts."""
    return np.array([row[key] for row in results], dtype=float)


def _pick_optional_float_array(results, key):
    """Return a float array; missing values become np.nan."""
    return np.array([np.nan if row[key] is None else row[key] for row in results], dtype=float)


def summarize_ecut_results(
    ecut_results,
    ry_to_ev,
    nat,
    threshold_mev_atom,
    first_globally_converged_index,
):
    """Build arrays and convergence summary for an ecut sweep."""
    ecut_x = _pick_float_array(ecut_results, "ecutwfc")
    ecut_e_ry = _pick_float_array(ecut_results, "energy_ry")
    ecut_t = _pick_float_array(ecut_results, "wall_s")

    ecut_e_ev_atom = ecut_e_ry * ry_to_ev / nat
    dE_mev_atom = np.abs(ecut_e_ev_atom - ecut_e_ev_atom[-1]) * 1000.0

    idx_conv_ecut = first_globally_converged_index(dE_mev_atom, threshold_mev_atom)
    ecut_conv = int(ecut_x[idx_conv_ecut]) if idx_conv_ecut is not None else int(ecut_x[-1])

    return {
        "ecut_x": ecut_x,
        "ecut_e_ev_atom": ecut_e_ev_atom,
        "dE_mev_atom": dE_mev_atom,
        "ecut_t": ecut_t,
        "idx_conv_ecut": idx_conv_ecut,
        "ecut_conv": ecut_conv,
    }


def summarize_k_results(
    k_results,
    ry_to_ev,
    nat,
    threshold_mev_atom,
    first_globally_converged_index,
):
    """Build arrays and convergence summary for a k-mesh sweep."""
    k_x = _pick_float_array(k_results, "nk")
    k_e_ry = _pick_float_array(k_results, "energy_ry")
    k_t = _pick_float_array(k_results, "wall_s")
    k_irr = _pick_optional_float_array(k_results, "nk_irr")

    k_e_ev_atom = k_e_ry * ry_to_ev / nat
    dE_k_mev_atom = np.abs(k_e_ev_atom - k_e_ev_atom[-1]) * 1000.0

    idx_conv_k = first_globally_converged_index(dE_k_mev_atom, threshold_mev_atom)
    nk_conv = int(k_x[idx_conv_k]) if idx_conv_k is not None else int(k_x[-1])

    return {
        "k_x": k_x,
        "k_e_ev_atom": k_e_ev_atom,
        "dE_k_mev_atom": dE_k_mev_atom,
        "k_t": k_t,
        "k_irr": k_irr,
        "idx_conv_k": idx_conv_k,
        "nk_conv": nk_conv,
    }


def summarize_force_stress_ecut_results(
    ecut_results,
    force_threshold_mev_ang,
    stress_threshold_kbar,
    first_globally_converged_index,
):
    """Build arrays and convergence summary for force/stress vs ecut."""
    ecut_x = _pick_float_array(ecut_results, "ecutwfc")
    force_z_ev_ang = _pick_float_array(ecut_results, "force_z_ev_ang")
    stress_zz_kbar = _pick_float_array(ecut_results, "stress_zz_kbar")
    ecut_t = _pick_float_array(ecut_results, "wall_s")

    dF_mev_ang = np.abs(force_z_ev_ang - force_z_ev_ang[-1]) * 1000.0
    dS_kbar = np.abs(stress_zz_kbar - stress_zz_kbar[-1])

    idx_conv_force = first_globally_converged_index(dF_mev_ang, force_threshold_mev_ang)
    idx_conv_stress = first_globally_converged_index(dS_kbar, stress_threshold_kbar)

    max_norm_error = np.maximum(
        dF_mev_ang / force_threshold_mev_ang,
        dS_kbar / stress_threshold_kbar,
    )
    idx_conv_both = first_globally_converged_index(max_norm_error, 1.0)
    ecut_conv_both = int(ecut_x[idx_conv_both]) if idx_conv_both is not None else int(ecut_x[-1])

    return {
        "ecut_x": ecut_x,
        "force_z_ev_ang": force_z_ev_ang,
        "stress_zz_kbar": stress_zz_kbar,
        "dF_mev_ang": dF_mev_ang,
        "dS_kbar": dS_kbar,
        "ecut_t": ecut_t,
        "idx_conv_force": idx_conv_force,
        "idx_conv_stress": idx_conv_stress,
        "idx_conv_both": idx_conv_both,
        "ecut_conv_both": ecut_conv_both,
    }


def summarize_force_stress_k_results(
    k_results,
    force_threshold_mev_ang,
    stress_threshold_kbar,
    first_globally_converged_index,
):
    """Build arrays and convergence summary for force/stress vs k-mesh."""
    k_x = _pick_float_array(k_results, "nk")
    force_z_ev_ang = _pick_float_array(k_results, "force_z_ev_ang")
    stress_zz_kbar = _pick_float_array(k_results, "stress_zz_kbar")
    k_t = _pick_float_array(k_results, "wall_s")
    k_irr = _pick_optional_float_array(k_results, "nk_irr")

    dF_mev_ang = np.abs(force_z_ev_ang - force_z_ev_ang[-1]) * 1000.0
    dS_kbar = np.abs(stress_zz_kbar - stress_zz_kbar[-1])

    idx_conv_force = first_globally_converged_index(dF_mev_ang, force_threshold_mev_ang)
    idx_conv_stress = first_globally_converged_index(dS_kbar, stress_threshold_kbar)

    max_norm_error = np.maximum(
        dF_mev_ang / force_threshold_mev_ang,
        dS_kbar / stress_threshold_kbar,
    )
    idx_conv_both = first_globally_converged_index(max_norm_error, 1.0)
    nk_conv_both = int(k_x[idx_conv_both]) if idx_conv_both is not None else int(k_x[-1])

    return {
        "k_x": k_x,
        "force_z_ev_ang": force_z_ev_ang,
        "stress_zz_kbar": stress_zz_kbar,
        "dF_mev_ang": dF_mev_ang,
        "dS_kbar": dS_kbar,
        "k_t": k_t,
        "k_irr": k_irr,
        "idx_conv_force": idx_conv_force,
        "idx_conv_stress": idx_conv_stress,
        "idx_conv_both": idx_conv_both,
        "nk_conv_both": nk_conv_both,
    }
