# Vibrational spectroscopy

## IR intensities

`ir_spectrum` keeps the fixed-charge dipole model. `BornCharges` adds a
Born-effective-charge path:

- `BornCharges::from_tensors(Vec<[[f64; 3]; 3]>)` takes one tensor per atom;
- `BornCharges::from_text` reads a simple `atom`/nine-value table;
- `BornCharges::from_charges(&Charges)` builds the isotropic diagonal;
- `BornCharges::qeq` / `eqeq` fill the diagonal from the native charge models;
- `neutralised` subtracts the cell-average charge so the acoustic sum rule
  holds.

`ir_spectrum_with_born(&structure, &modes, &BornCharges, &IrOptions)` returns
the same `IrSpectrum { sticks, points, frequencies, intensities }` as the
fixed-charge form, with Lorentzian broadening.

## Raman

`RamanTensorSet` reads Raman tensors from a text/JSON file
(`from_text`) or from an explicit list (`from_tensors`). `raman_spectrum`
computes the powder average
`(1/15)( (Tr T)^2 + 2 sum_ij T_ij^2 )` per mode and broadens it with the
`RamanOptions` Lorentzian. The result `RamanSpectrum { peaks, points,
frequencies, activities }` is deterministic. No native polarizability model is
claimed: the tensors must come from an engine output.

## Phonon DOS and inelastic neutron scattering

`phonon_dos(&PhononResult, &DosOptions)` builds the Gaussian phonon DOS.
`NeutronWeights` supplies the per-atom coherent scattering length (from a
table, a callback or an explicit list), and
`inelastic_neutron_dos(&result, &options, &weights)` weights every mode by
`|b|^2 / mass` before broadening. `mode_neutron_weights` returns the
per-mode weights for inspection. The normalization is documented on the
algorithm page; the weights are squared scattering lengths, not a full
dynamical structure factor.

Parity keys: `ir_born`, `raman`, `phonon_dos`, `phonon_dos_neutron`.
