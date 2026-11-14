# Neutron and electron diffraction

## Nuclear neutron scattering

The coherent bound scattering length `b` is a constant per isotope/element
from the embedded table (femtometre converted to Angstrom). Unlike X-rays, `b`
does not fall with angle and can be negative (H, Ti, Mn, V), which changes the
phase of the structure factor. The pattern uses the Phase 7 profile model:

- constant-wavelength mode: the 2-theta axis, with the instrument profile;
- time-of-flight mode: `t = 252.78 * L * d` microseconds for path length `L`
  in metres and d-spacing in Angstrom, with the same reflection list mapped to
  the TOF axis.

`NeutronOptions` sets radiation, geometry, path length, probe and the pattern
range. `NeutronPattern` carries the sampled points and the reflections.

## Magnetic neutron scattering

The magnetic structure factor uses the component of the moment perpendicular
to the scattering vector:

```text
F_mag(Q) = r_e gamma/2 * sum_j f_j(Q) (m_j - (m_j . Qhat) Qhat) exp(i Q . r_j)
```

with `r_e gamma/2 = 0.2695` in the conventional units and `f_j(Q)` the
magnetic form factor. The implementation uses the dipole approximation for the
form factor with a documented radial table and falls back to a smooth
`exp(-a s^2)` for elements without one. `magnetic_powder` returns the magnetic
pattern; `with_nuclear(&NeutronPattern)` adds it to the nuclear pattern. A
structure without moments is refused instead of returning zero intensity.

## Electron diffraction

`zone_axis_pattern` is kinematic SAED:

- the zone axis defines the reciprocal plane;
- electron scattering factors use the Mott-Bethe relation from the X-ray
  form factor, `f_e(s) = 0.023934 (Z - f_x(s)) / s^2`, capped at small `s`;
- spot positions are the in-plane reciprocal vectors, intensities include the
  multiplicity and the Lorentz factor;
- `SaedPattern::image()` rasterises the spots into `SaedImage`, with
  `to_u8` for an 8-bit display buffer.

Dynamical scattering, thickness fringes and Kikuchi lines are out of scope.

## Determinism and data

All three paths sort reflections deterministically and use the shared P1
expansion. The neutron scattering-length table and the magnetic form-factor
table are embedded data with their sources in `THIRD_PARTY.md`; the fixtures
used by the tests are hand-written representative patterns, not measured data.

## Limits

- No absorption, extinction or multiple-scattering correction.
- The magnetic form factor is an approximation; it is not a full
  symmetry-adapted magnetic representation analysis.
- SAED is kinematic and uses a flat Ewald sphere.
