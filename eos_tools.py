import numpy as np
from scipy.optimize import curve_fit

EV_ANG3_TO_GPA = 160.2176   # 1 eV/Å³ in GPa


def fit_parabola(x, y):
    """Fit y = c2*x^2 + c1*x + c0; return (poly1d, x0, y0) at the vertex."""
    coeffs = np.polyfit(x, y, 2)           # [c2, c1, c0]
    poly   = np.poly1d(coeffs)
    x0     = -coeffs[1] / (2 * coeffs[0])  # vertex
    return poly, x0, poly(x0)


def birch_murnaghan(V, E0, V0, B0, B0p):
    """Third-order Birch–Murnaghan energy EOS (B0 in eV/Å³)."""
    eta = (V0 / V) ** (2 / 3)
    return E0 + (9 * V0 * B0 / 16) * (
        B0p * (eta - 1) ** 3 + (eta - 1) ** 2 * (6 - 4 * eta)
    )


def bm_pressure_gpa(V, V0, B0_gpa, B0p):
    """3rd-order Birch–Murnaghan pressure P = -dE/dV in GPa.

    P = (3B₀/2)·(η⁷ − η⁵)·[1 + (3/4)(B₀′−4)(η²−1)],  η = (V₀/V)^(1/3)
    """
    eta = (V0 / V) ** (1 / 3)
    return 1.5 * B0_gpa * (eta**7 - eta**5) * (1 + 0.75 * (B0p - 4) * (eta**2 - 1))


def parabola_pressure_gpa(V, poly_v):
    """Pressure from a parabolic E(V) fit: P = −dE/dV in GPa.

    Parameters
    ----------
    poly_v : np.poly1d
        Result of ``fit_parabola(V, E)[0]`` with E in eV and V in Å³.
    """
    return -poly_v.deriv()(V) * EV_ANG3_TO_GPA


def fit_bm_eos(V, E):
    """
    Fit a 3rd-order Birch–Murnaghan EOS to E(V) data.

    Returns a dict with:
      V0      – equilibrium volume (Å³)
      B0_gpa  – bulk modulus (GPa)
      B0p     – pressure derivative B₀′ (dimensionless)
      E0      – equilibrium energy (eV)
      popt    – raw (E0, V0, B0_eV_A3, B0p) for use with birch_murnaghan()
    """
    poly_v, V0_par, E0_par = fit_parabola(V, E)
    B0_par = 2 * poly_v.coeffs[0] * V0_par   # eV/Å³
    p0 = [E0_par, V0_par, B0_par, 4.0]
    popt, _ = curve_fit(birch_murnaghan, V, E, p0=p0)
    E0_fit, V0_fit, B0_fit, B0p_fit = popt
    return {
        'E0':     E0_fit,
        'V0':     V0_fit,
        'B0_gpa': B0_fit * EV_ANG3_TO_GPA,
        'B0p':    B0p_fit,
        'popt':   popt,
    }
