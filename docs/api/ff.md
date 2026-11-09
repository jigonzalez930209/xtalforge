# Force fields

`xtalforge-ff` is the energy-engine layer. `Calculator` answers energy and
forces for a structure; `UffCalculator` is the native Universal Force Field
with UFF4MOF types.

```rust
use xtalforge::ff::{CalculatorOptions, UffCalculator};

let calculator = UffCalculator::new(CalculatorOptions::default())?;
let energy = calculator.energy(&structure)?;
let (energy, forces) = calculator.energy_forces(&structure)?;
```

Units at the boundary are eV, eV/Angstrom and GPa; UFF keeps kcal/mol and
Angstrom internally. The calculator builds the P1 topology once per call;
`UffModel` exposes the prepared topology for optimizers that move positions.

`CalculatorRegistry` lists engines and availability through `CalculatorInfo`
(name, label, description, capabilities); `RunRecord { calculator, energy,
n_atoms, converged, steps, notes }` is the serializable run record.
Dummy `X` atoms are excluded from the energy and return a zero force.

The UFF table covers the 218 Rappe/UFF4MOF types; the typer uses element and
bond geometry, with UFF4MOF metal types for framework nodes. Charges, Ewald,
optimizers, vibrations and thermochemistry build on the same `Calculator`
contract and are documented below.

## Charges

`qeq_charges(structure, &ChargeOptions)`, `eqeq_charges` and `eeq_charges`
solve charge equilibration and return `Charges { per_atom, total }`;
`apply_charges` writes them into the sites. `EwaldOptions` and
`ewald_energy`/`ewald_energy_forces` return an `EwaldResult` and provide the
periodic electrostatic term.
EEQ reuses the local charge table (the sci-form fitted table is not vendored);
the module documents the approximation.

## Optimizers

`minimize(calculator, structure, &OptimizeOptions, Minimizer::Fire)` runs FIRE
or L-BFGS over the P1 coordinates (asymmetric mode keeps the asymmetric unit)
and returns `OptimizationResult { structure, energy, max_force, steps,
converged }`. `apply_optimized_geometry` writes the result back. `fit_eos`
fits a Birch-Murnaghan or Murnaghan equation of state from energy-volume
points and returns an `EosFit { v0, e0, b0, b0_prime, rms }`.

## Scans

`ScanCoordinate` and `Constraint` describe a distance/angle/torsion scan; the
driver steps the coordinate and relaxes the rest.

## Vibrations and thermodynamics

`hessian(calculator, structure, &HessianOptions)` builds the Gamma-point
Hessian by central differences of the analytic forces.
`vibrational_modes(&hessian, structure, &HessianOptions)` mass-weights it and
returns `VibrationalModes { frequencies, modes, zero_modes, imaginary }` in
cm^-1.

`ir_spectrum(structure, &modes, &IrOptions)` uses QEq charges and a
finite-difference dipole derivative (documented approximation) and broadens
the sticks into `IrSpectrum { points, peaks }`.

`phonons(...)` uses the primitive-cell Hessian at Gamma. With
`PhononOptions { q_grid, supercell, q_path, dos }` it builds the supercell
force constants, folds them onto the primitive cell, samples a Monkhorst-Pack
mesh or the explicit `q_path` (a band path) and returns `q_points`,
`frequencies` per q-point, the Gamma `modes` and an optional Gaussian-smeared
`PhononDos`. `thermochemistry(&modes, &ThermoOptions)` and
`ideal_gas_thermochemistry` give ZPE, thermal U/H/S/G and Cv; imaginary modes
are excluded and reported. `quasi_harmonic_free_energy(&points, options)`
minimizes `F(V, T) = E(V) + F_vib(V, T)` by a parabola through the best sampled
volume.

## Preparation pipeline

`prepare(&calculator, &structure, &PrepareOptions)` runs the optional cleanup
and force-field stages in order: duplicate merging, graph-aware disorder
choice guided by the declared formula and cell charge, solvent removal,
optional hydrogens, FIRE relaxation, QEq charges,
periodic Ewald, Gamma frequencies, close-contact listing and the charge-balance
check. `PrepareReport` records every step, including skipped ones, the prepared
`structure`, removed-site count, close contacts, energy, total charge, Ewald
energy, frequencies and the optional `ChargeBalance`. Apply the prepared
structure as one undoable edit with `xtalforge_commands::ReplaceStructure`.

## External engines

`probe(EngineKind)` checks a binary by running `--version`; on WASM it reports
`Missing`. `write_dftb_hsd`, `write_xtb_input` and `write_orca_input` prepare
inputs; `parse_dftb_output` and `parse_xtb_output` return an `EngineResult` with the
energy, forces and optional stress. ORCA is input-only and `orca_catalogue` lists the bundled methods.

## Quantum adapter

`quantum::adapter(method)` is the reserved in-process sci-form path. Without
the `sci-form` feature it returns a typed error that names the feature;
`capabilities` records periodic/forces/stress so a future adapter cannot
silently overclaim.
