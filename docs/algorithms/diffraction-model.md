# Diffraction model: profiles, radiation and single crystals

This page records the Phase 14 model shared by PXRD, Le Bail, Pawley,
Rietveld and the single-crystal helpers.

## Peak shape

The Gaussian full width at half maximum follows Caglioti's polynomial in
`tan(theta)`:

```text
FWHM_G(2 theta) = sqrt(U tan^2 theta + V tan theta + W)
```

A Lorentzian width is either constant (`FWHM_L`) or follows the same
polynomial with its own `U, V, W`. The pseudo-Voigt is the area-preserving
mixture

```text
pV(x) = eta L(x) + (1 - eta) G(x),      0 <= eta <= 1
```

with `x = 2 theta - 2 theta_Bragg`. TCH splits the Lorentzian into a size term
(`1 / cos theta`) and a strain term (`tan theta`). Every profile is normalised
to unit area; the default remains the Gaussian used in `1.6.0`.

## Asymmetry and polarisation

The Finger-Cox-Jephcoat correction models axial divergence with the sample
height, detector height and diffractometer radius. It is applied as a
multiplicative correction on each side of the peak and is only active for
`BraggBrentano` geometry. The monochromator polarisation enters the
Lorentz-polarisation factor:

```text
LP = (1 + cos^2 2 theta) / (sin^2 theta cos theta)      (unpolarised)
LP = (1 + cos^2 2 theta cos^2 2 alpha) / ...            (monochromator)
```

## Radiation and scattering lengths

`Radiation` keeps the laboratory wavelengths (Cu, Mo, Co, Fe, Cr, Ag), an
arbitrary X-ray wavelength, a neutron wavelength and an electron wavelength.
X-rays use the International Tables Gaussian form factors; neutrons use the
bound coherent scattering lengths `b_c` (including the negative values for
Mn, H and Ti) and the incoherent cross-section is ignored in the coherent
pattern. Electron diffraction uses the electron scattering factor at the
requested accelerating voltage.

## Anisotropic displacement

The Debye-Waller factor is

```text
T(hkl) = exp(-2 pi^2 sum_ij h_i h_j a_i^* a_j^* U_ij)
```

with `U_ij` in Angstrom squared. A site with only `U_iso` keeps the isotropic
form `exp(-8 pi^2 U_iso sin^2 theta / lambda^2)`. Rietveld and the
single-crystal structure factors read the same tensor.

## Single crystal and HKL

`structure_factors` returns `h, k, l, |F|, F_real, F_imag, phase, d` for every
reflection in a `d` range. `write_hkl` emits the SHELX HKLF 4 layout.
`electron_diffraction` builds the same reflections with electron form factors
and returns the kinematic intensities.

## Refinable parameters

Le Bail, Pawley and Rietveld option structs carry `eta`, Caglioti `U, V, W`
and the asymmetry parameters. A parameter is refined only when its flag is set;
the old option structs remain valid and keep the `1.6.0` defaults.

## Validation

- Synthetic patterns with a known `U, V, W` recover the widths within 1e-3.
- A pseudo-Voigt with `eta = 1` equals the Lorentzian and `eta = 0` the
  Gaussian within 1e-9.
- Anisotropic and isotropic Debye-Waller agree for a diagonal `U_ij` equal to
  `U_iso`.
- HKL round trips keep the reflection list and phases.
- Neutron and electron patterns are checked against the bundled constants and
  the analytic two-atom interference.
