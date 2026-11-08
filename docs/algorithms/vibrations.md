# Vibrations, phonons and thermodynamics

## Gamma-point Hessian

`hessian(calculator, structure, &HessianOptions)` displaces each Cartesian
coordinate of the P1 cell by `+step` and `-step`, reads the analytic forces and
forms the central-difference Hessian, symmetrised as `(H + H^T)/2`. The cost is
`6N + 1` force evaluations. The structure must be P1; a dummy atom is rejected
because it has no force-field parameters.

## Modes

`vibrational_modes(&hessian, structure, options)` mass-weights the Hessian with
the element masses, projects out the rigid translations (and the molecular
rotations when `DisplacementMode::Molecule` is selected) and diagonalises the
symmetric matrix with a deterministic Jacobi eigensolver. Frequencies are in
cm^-1, ascending. A mode within `ZERO_MODE_TOLERANCE_CM` is a zero mode; a
negative frequency below `IMAGINARY_BELOW_CM` is imaginary and is reported,
not hidden.

`VibrationalModes { frequencies, modes, zero_modes, imaginary }` carries the
normalised Cartesian displacement fields.

## IR spectrum

`ir_spectrum(structure, modes, &IrOptions)` uses fixed QEq charges and a
finite-difference dipole derivative with respect to each mode. The intensity is
the squared dipole derivative in km/mol. The module states the approximation:
Born effective charges are replaced by the fixed QEq charges, which is the
same approximation Crystal-Builder documents. The sticks are broadened with a
Gaussian or Lorentzian into `IrSpectrum { points, peaks }`.

## Phonons

`phonons` builds the Gamma-point dynamical matrix from the same Hessian and
returns `PhononResult`. The full supercell force constants, a q-mesh and the
phonon DOS are documented as pending; a q-grid other than `[1, 1, 1]` returns a
typed error rather than a silent approximation.

## Thermochemistry

`thermochemistry(frequencies, &ThermoOptions)` sums the harmonic oscillator
contributions: zero-point energy, thermal U and H, entropy and Gibbs energy at
a temperature. `ideal_gas_thermochemistry` adds the translational and
rigid-rotor terms at a pressure. Imaginary, zero and non-finite modes do not
enter the partition function; the caller sees them in the mode result.

## Quasi-harmonic free energy

`quasi_harmonic_free_energy(points, options)` evaluates
`F(V, T) = E(V) + F_vib(V, T)` at each supplied volume. With three or more
points the lowest sampled value is refined by a parabola through its
neighbours. The result carries the points and the refined minimum. Thermal
expansion and bulk modulus versus temperature are not yet derived.
