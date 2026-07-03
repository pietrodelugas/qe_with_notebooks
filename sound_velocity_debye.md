# Sound velocity from phonon dispersion and DOS

Notes on the unit conversions and the Debye model fitting used in
`qe_phonon_mgo_demo.ipynb`, section 11.

---

## 1. Unit conversion: dispersion slope → velocity

matdyn.x writes q-points in **Cartesian reciprocal coordinates in units of
2π/alat**, where alat is the lattice parameter used in the pw.x run (for
`ibrav=2` this is the conventional cubic cell parameter a₀).  The cumulative
path stored in the `.freq.gp` file is in the same units.

So a path value of 1.0 corresponds to an actual wavevector magnitude:

    |q| = 1.0 × (2π / alat)   [rad/m]

The phonon frequency in cm⁻¹ is related to angular frequency by:

    ω [rad/s] = 2π × c [m/s] × ν̃ [cm⁻¹] × 100 [m⁻¹ / cm⁻¹]

The sound velocity is dω/d|q|:

    v_s = dω / d|q|
        = (2π × c × 100 × Δν̃) / (Δpath × 2π / alat)
        = c [m/s] × 100 [m/cm] × alat [m] × (Δν̃ / Δpath)

In the notebook the slope is estimated as the mean of ν̃/path over a few
near-Γ q-points (i.e. a slope through the origin).  A more robust estimate
uses least squares forced through the origin:

    slope = Σ(path_i · ν̃_i) / Σ(path_i²)    [minimises Σ(ν̃ - slope·path)²]

### Getting alat

For `ibrav=2`, `celldm(1)` = a₀ (conventional cubic cell, in bohr).
ASE's `bulk('MgO', 'rocksalt', a=4.211)` returns the **primitive** FCC cell
whose vectors have length a₀/√2, so:

    a_conv_ang = mgo.cell.lengths()[0] * np.sqrt(2)   # Å
    alat_m     = a_conv_ang * 1e-10                    # m

---

## 2. Debye density of states formula

### Starting point

In the Debye (long-wavelength) approximation the density of states for a
single acoustic branch with velocity v_s in a crystal of primitive cell
volume V_cell is:

    g_s(ω) = V_cell / (2π²)  ×  ω² / v_s³

Summing over all 3 acoustic branches and defining the **Debye velocity** v_D
through

    3 / v_D³  ≡  1/v_LA³ + 1/v_TA1³ + 1/v_TA2³

gives:

    g(ω) = V_cell × 3ω² / (2π² v_D³)

### Converting to cm⁻¹

The matdyn DOS uses wavenumber ν̃ [cm⁻¹] as the frequency axis.
Substituting ω = 2πcν̃ (where c is in cm/s) and using dω = 2πc dν̃:

    g(ν̃) = g(ω) × |dω/dν̃|
          = [V_cell × 3(2πcν̃)² / (2π² v_D³)] × 2πc
          = V_cell × 3 × 4π²c²ν̃² × 2πc / (2π² v_D³)
          = V_cell × 12π c³ / v_D³  ×  ν̃²

So the low-frequency DOS fits **g(ν̃) = A ν̃²** with coefficient:

    A = 12π c³ V_cell / v_D³

Inverting for v_D (all in CGS: c in cm/s, V_cell in cm³, A in cm³):

    v_D [cm/s] = c × (12π V_cell / A)^(1/3)

Or equivalently in SI (c in m/s, V_cell in m³, A converted to m³):

    v_D [m/s] = c × (12π V_cell [m³] × (100 [m/cm])³ / A [states·m²])^(1/3)

The factor of 100³ converts c³ from (cm/s)³ to (m/s)³ while keeping V_cell
in m³ (1 cm³ = 10⁻⁶ m³, and A in states·cm³ → states·m³ = A × 10⁻⁶).

### Where does 12π come from?

The chain of factors is:

    3          — three acoustic branches
    ×  (2π)²  — from ω² = (2πcν̃)²
    ×  2π      — from |dω/dν̃| = 2πc
    ÷  2π²     — from the DOS prefactor 1/(2π²)
    =  3 × (4π²) × (2πc) / (2π²c²) × c²  →  12π

---

## 3. Robust fitting

### Dispersion: least squares through origin

For N points (path_i, ν̃_i) near Γ, minimise Σ(ν̃_i − slope·path_i)²:

```python
slope = np.dot(path[1:N], freq[1:N, branch]) / np.dot(path[1:N], path[1:N])
```

This avoids bias from the Γ point (which has both path=0 and freq=0 by ASR,
contributing nothing to the sum but regularising a polyfit intercept).

### DOS: least squares through origin

For M points (ν̃_j, g_j) in the Debye region, minimise Σ(g_j − A·ν̃_j²)²:

```python
# Substitute X = ν̃², Y = g, then A = Σ(X·Y)/Σ(X²)
A = np.dot(E[mask]**2, dos[mask]) / np.dot(E[mask]**2, E[mask]**2)
```

Equivalently, fitting `dos` as a linear function of `E**2` (polyfit degree 1
forced through the origin) gives the same result.

### Choosing the fit range

- **Dispersion**: use only the first 5–10 q-points along Γ→X.  Beyond ~10%
  of the BZ the relationship becomes nonlinear and the slope overestimates v_s.
- **DOS**: use ν̃ < ~80 cm⁻¹ for MgO (optical branches start ~400 cm⁻¹;
  the quadratic regime is well-established up to ~100 cm⁻¹).  Exclude the
  first bin (ν̃ < 5 cm⁻¹) which can be noisy due to Gaussian smearing
  artifacts.

---

## 4. Consistency check

The Debye velocity from the DOS fit and from the dispersion slopes should
agree to within a few percent for a well-converged calculation.  Differences
arise from:

- **Direction average**: the DOS integrates over all q-directions while the
  dispersion slope uses only Γ→X.  For a cubic crystal the acoustic branches
  are isotropic at long wavelength, so the effect is small.
- **IFC quality**: a coarse q-grid (e.g. 2×2×2) makes the IFC short-ranged,
  which can slightly alter the near-Γ curvature.
- **Gaussian broadening**: the `degauss` parameter in the DOS calculation
  smears spectral weight to lower frequencies, inflating A and reducing v_D.
  Use a small `degauss` (≤ 5 cm⁻¹) for the Debye fit.
