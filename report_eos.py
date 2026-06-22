#!/usr/bin/env python3
"""
Parse pw.x output files, fit Murnaghan EOS from E(V),
and print a table of computed vs fitted pressures at each lattice constant.

Usage:
    report_eos.py [GLOB_PATTERN]

    GLOB_PATTERN  shell glob for output files (default: si_r2scan_eos_*.out)
"""

import glob
import re
import sys
import numpy as np
from scipy.optimize import curve_fit

RY_TO_EV       = 13.605703976
BOHR_TO_ANG    = 0.529177210903
RY_BOHR3_TO_GPA = 14710.507
KBAR_TO_GPA    = 0.1


def parse_pw_output(filename):
    with open(filename) as f:
        txt = f.read()
    m = re.search(r'celldm\(1\)\s*=\s*([\d.]+)', txt)
    celldm = float(m.group(1)) if m else None
    # parse unit-cell volume directly so the script works for any ibrav;
    # use last occurrence so vc-relax output gives the final relaxed volume
    vhits = re.findall(r'unit-cell volume\s*=\s*([\d.]+)', txt)
    vol = float(vhits[-1]) if vhits else None
    hits = re.findall(r'!\s+total energy\s+=\s+([-\d.]+)\s+Ry', txt)
    energy = float(hits[-1]) if hits else None
    phits = re.findall(r'P=\s*([-\d.]+)', txt)
    pressure = float(phits[-1]) if phits else None
    return celldm, vol, energy, pressure


def murnaghan(V, E0, V0, B0, B0p):
    return E0 + B0 * V0 / B0p * (
        (1.0 / (B0p - 1.0)) * (V0 / V) ** (B0p - 1.0)
        + V / V0 - B0p / (B0p - 1.0)
    )


def murnaghan_P(V, V0, B0, B0p):
    """Murnaghan pressure in same units as B0."""
    return (B0 / B0p) * ((V0 / V) ** B0p - 1.0)


def fit_murnaghan(volumes, energies):
    i0 = np.argmin(energies)
    p0 = [energies[i0], volumes[i0], 0.005, 4.0]
    popt, _ = curve_fit(murnaghan, volumes, energies, p0=p0, maxfev=20000)
    rms = np.sqrt(np.mean((murnaghan(volumes, *popt) - energies) ** 2))
    return popt, rms


def main():
    pattern = sys.argv[1] if len(sys.argv) > 1 else "si_r2scan_eos_*.out"
    outfiles = sorted(glob.glob(pattern))
    if not outfiles:
        sys.exit(f"No files matching '{pattern}' found.")

    rows = []
    skipped = []
    for fn in outfiles:
        celldm, vol, energy, pressure = parse_pw_output(fn)
        if celldm is None or vol is None or energy is None:
            skipped.append(fn)
            continue
        rows.append((celldm, vol, energy, pressure))

    if skipped:
        print(f"Skipped (not converged): {skipped}\n")

    if len(rows) < 4:
        sys.exit("Need at least 4 converged points.")

    rows.sort()
    celldms, volumes, energies, pressures_raw = map(np.array, zip(*rows))

    # --- Murnaghan fit from E(V) ---
    popt, rms = fit_murnaghan(volumes, energies)
    E0, V0, B0, B0p = popt
    # infer the celldm→volume factor from the data (works for any ibrav)
    vol_factor = np.mean(volumes / celldms ** 3)
    a0    = (V0 / vol_factor) ** (1.0 / 3.0)
    B0_gpa = B0 * RY_BOHR3_TO_GPA

    print("Murnaghan EOS  [from E(V)]")
    print("=" * 40)
    print(f"  E0   = {E0:.6f} Ry  ({E0 * RY_TO_EV:.4f} eV)")
    print(f"  V0   = {V0:.4f} Bohr^3  ({V0 * BOHR_TO_ANG**3:.4f} Ang^3)")
    print(f"  a0   = {a0:.6f} Bohr  ({a0 * BOHR_TO_ANG:.6f} Ang)")
    print(f"  B0   = {B0_gpa:.2f} GPa")
    print(f"  B0'  = {B0p:.4f}")
    print(f"  RMS  = {rms * RY_TO_EV * 1e3:.4f} meV/cell")

    # Murnaghan P(V) in kbar (B0 in Ry/Bohr^3 → convert to kbar)
    B0_kbar = B0_gpa / KBAR_TO_GPA
    P_murn_kbar = murnaghan_P(volumes, V0, B0_kbar, B0p)

    # --- Pressure table ---
    print()
    print("Pressure table")
    print("=" * 75)
    hdr = f"{'a (Bohr)':>10}  {'V (Bohr^3)':>11}  {'E (Ry)':>13}  "
    hdr += f"{'P_QE (kbar)':>12}  {'P_Murn (kbar)':>14}  {'dP (kbar)':>10}"
    print(hdr)
    print("-" * 75)
    for i, (a, V, E, P_qe) in enumerate(zip(celldms, volumes, energies, pressures_raw)):
        P_m = P_murn_kbar[i]
        if P_qe is not None:
            dP = P_qe - P_m
            print(f"{a:10.4f}  {V:11.3f}  {E:13.8f}  "
                  f"{P_qe:12.2f}  {P_m:14.2f}  {dP:10.2f}")
        else:
            print(f"{a:10.4f}  {V:11.3f}  {E:13.8f}  "
                  f"{'N/A':>12}  {P_m:14.2f}  {'N/A':>10}")
    print("-" * 75)
    print(f"  dP = P_QE - P_Murnaghan (should be ~0 if stress is correct)")


if __name__ == "__main__":
    main()
