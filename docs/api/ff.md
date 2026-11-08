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

`CalculatorRegistry` lists engines and availability. `RunRecord { calculator,
energy, n_atoms, converged, steps, notes }` is the serializable run record.
Dummy `X` atoms are excluded from the energy and return a zero force.

The UFF table covers the 218 Rappe/UFF4MOF types; the typer uses element and
bond geometry, with UFF4MOF metal types for framework nodes. Charges, Ewald,
optimizers and thermochemistry are later Phase 4 tasks.

## Charges

`qeq_charges(structure, &ChargeOptions)` and `eqeq_charges` solve charge
equilibration and return `Charges { per_atom, total }`; `apply_charges` writes
them into the sites. `EwaldOptions` and the Ewald energy/forces provide the
periodic electrostatic term.

## Optimizers

`minimize(structure, calculator, &OptimizeOptions)` runs FIRE or L-BFGS over
the P1 coordinates (asymmetric mode keeps the asymmetric unit) and returns
`OptimizationResult { structure, energy, max_force, steps, converged }`.
`apply_optimized_geometry` writes the result back. `fit_eos` fits a
Birch-Murnaghan or Murnaghan equation of state from energy-volume points.

## Scans

`ScanCoordinate` and `Constraint` describe a distance/angle/torsion scan; the
driver steps the coordinate and relaxes the rest.

## Vibrations and thermodynamics

`hessian(structure, calculator, &HessianOptions)` builds the Gamma-point
Hessian by central differences of the analytic forces.
`vibrational_modes(&hessian, masses, options)` mass-weights it and returns
`VibrationalModes { frequencies, modes, zero_modes, imaginary }` in cm^-1.

`ir_spectrum(&modes, structure, &IrOptions)` uses QEq charges and a
finite-difference dipole derivative (documented approximation) and broadens
the sticks into `IrSpectrum { points, peaks }`.

`phonons(...)` returns the Gamma-point dynamical matrix and frequencies;
the q-mesh and DOS are pending. `thermochemistry(&modes, &ThermoOptions)` and
`ideal_gas_thermochemistry` give ZPE, thermal U/H/S/G and Cv; imaginary modes
are excluded and reported. `quasi_harmonic_free_energy(&points, options)`
minimizes `F(V, T) = E(V) + F_vib(V, T)` by a parabola through the best sampled
volume.
