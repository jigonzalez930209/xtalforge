# Vibrations, phonons and thermodynamics

## Gamma-point Hessian

`hessian(calculator, structure, &HessianOptions)` displaces each Cartesian
coordinate of the P1 cell by `+step` and `-step`, reads the analytic forces and
forms the central-difference Hessian, symmetrised as `(H + H^T)/2`. The cost is
`6N + 1` force evaluations. The structure must be P1; a dummy atom is rejected
because it has no force-field parameters.

## Modes

![Hessian to modes](/diagrams/vibrations-pipeline.svg)


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

`phonons` with `q_grid = [1, 1, 1]` and `supercell = [1, 1, 1]` uses the
primitive-cell Hessian. Any other grid builds the supercell, takes the
finite-difference Hessian there, folds the force constants onto the primitive
cell with the shortest periodic image, applies the acoustic sum rule and
diagonalizes the mass-weighted dynamical matrix at every Monkhorst-Pack point.
The grid q-points are `(n + 1/2) / N - 1/2` per axis, in fractional reciprocal
coordinates; an explicit `q_path` replaces the grid, which is how the Task 2.10
band path is sampled. The optional DOS is a Gaussian-smeared histogram normalized to
unit area.

The Gamma path remains the acoustic reference. The q-mesh path is validated
against the analytic dispersion of a nearest-neighbour spring model and
against a phonopy-generated golden file
(`tests/golden/phase-4/phonons/spring-sc.json`) for the same model; the
generator is `scripts/golden/phonons_spring.py`.

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
