# Electronic structure

`xtalforge-engines` analyses the eigenvalues and density of states that an
engine parser extracted. It never runs an engine by itself.

## Bands and gaps

`ElectronicStructure` stores `energies[band][kpoint]` in eV, optional
occupations, an optional spin-down channel and an optional Fermi/reference
energy. The parser result `BandData` uses the transposed
`energies[kpoint][band]`; the surface adapters transpose it before calling the
analysis.

`band_gap` classifies a state as occupied when its occupation is at least 0.5,
or, without occupations, when its energy is below the Fermi level. The
valence-band maximum and conduction-band minimum give the gap; the gap is
direct when both edges share a k-point. `spin_gaps` repeats the calculation per
spin channel. A metal reports a zero or negative gap with `metallic = true`.

## Effective masses

`effective_mass` fits `E(k) = E0 + a (k - k0) + b (k - k0)^2` over a window
around a band edge by least squares. The mass is
`m* = ELECTRON_MASS_FACTOR / (2 b)` in electron masses, with the factor
`7.619964222971922 eV Angstrom^2`. A hole band has a negative mass; the
magnitude is reported separately. `band_edge_masses` applies the fit to both
edges of a gap and `spin_band_edge_masses` to both spin channels. When no
Cartesian reciprocal path is supplied, the fractional path distance is used
and the result is a plotting-scale mass rather than a physical one; callers
that need physical masses must set `distances` in inverse Angstrom.

## Density of states

`smear_levels` places unit-area Gaussian levels on the energy grid with the
requested sigma and sums them, so the smeared curve integrates to the number of
states that were smeared. `density_of_states` builds the total curve and
groups the projections by element and orbital (`s`, `p`, `d`, `f`); the
`DosBuffer` keeps the total and projected arrays. `d_band_centre` is
`integral E * DOS_d(E) dE / integral DOS_d(E) dE` over the supplied grid.

## Limits

- No native band calculation: the bands come from VASP, Quantum ESPRESSO or
  another parser.
- The gap classification uses a 0.5 occupation threshold, not a temperature
  dependent Fermi-Dirac integral.
- Effective masses use a parabolic window; non-parabolic bands are reported
  with the fit residual through the returned curvature.
- Projections are whatever the engine prints (PROCAR/projwfc); no symmetry
  expansion or site mapping beyond the file order is attempted.
