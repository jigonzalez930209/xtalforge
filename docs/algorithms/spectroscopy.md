# Vibrational spectroscopy

## IR from Born charges

The Gamma-point modes give a mass-weighted eigenvector `e_i` per atom and a
frequency. The dipole derivative of a mode is

```text
d mu / dQ = sum_i Z_i^* e_i
```

with `Z_i^*` the Born effective-charge tensor. The intensity is proportional
to `|d mu/dQ|^2`; the existing `IrOptions` Lorentzian broadening turns the
sticks into the spectrum. `BornCharges::neutralised` subtracts the cell-average
charge to enforce the acoustic sum rule. The diagonal `qeq`/`eqeq` constructors
are a documented approximation: they use the scalar charge on the tensor
diagonal, not a finite-difference Born tensor.

## Raman

Raman activities require the derivative of the polarizability. xtalforge does
not ship a polarizability model, so `RamanTensorSet` reads the tensors from an
engine output or a text file. For a mode with tensor `T`, the powder average is

```text
I = (1/15) [ (Tr T)^2 + 2 sum_ij T_ij^2 ]
```

and `raman_spectrum` broadens the sticks with a Lorentzian. An isotropic tensor
gives `Tr^2/15`; a traceless symmetric tensor gives the documented factor. The
parser accepts a JSON list of 3x3 matrices or a simple text table.

## Phonon DOS

`phonon_dos` samples the q-mesh frequencies and sums unit-area Gaussians on the
requested energy grid. The DOS integrates to the number of modes when the grid
covers the full frequency range. `DosOptions` controls the grid, sigma and
normalisation.

## Inelastic neutron scattering

The neutron weight of a mode is `sum_i |b_i|^2 |e_i|^2 / m_i`, with `b_i` the
coherent scattering length and `m_i` the mass. `NeutronWeights` supplies the
lengths; `inelastic_neutron_dos` applies the weight before broadening. This is
the incoherent-limit weighting, not the full coherent dynamical structure
factor; the difference is documented rather than hidden.

## Limits

- No ab-initio Born or Raman tensors; they come from an engine or the QEq/EQeq
  diagonal.
- The powder average assumes randomly oriented crystallites.
- Anharmonicity and mode mixing are not modelled.
