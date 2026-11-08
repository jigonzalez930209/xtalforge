# Optimization, scans and equations of state

## Degrees of freedom

`RelaxMode::Asymmetric` is the default for a non-P1 structure: the variables
are the asymmetric sites, the calculator works on the expanded P1 cell, and
the forces are projected back onto the symmetry-allowed displacements. The
space group is preserved on every step. `RelaxMode::P1` is the ordinary
Cartesian relaxation.

`OptimizeOptions.cell` adds the six symmetric strain components as variables;
holding selected lengths or the volume while the shape relaxes is the
`CellFreedom` idea for flexible frameworks. External pressure enters through
the stress term.

## Minimizers

- **FIRE** (`Minimizer::Fire`): molecular dynamics with the velocity steered
  downhill, an adaptive time step and a force criterion.
- **L-BFGS** (`Minimizer::Lbfgs`): a limited-memory quasi-Newton method with a
  line search and restarts.

Both stop on `force_tolerance` (eV/Angstrom), `energy_tolerance` (eV) or
`max_steps`. `OptimizationResult { structure, energy, max_force, steps,
converged }` is deterministic. `apply_optimized_geometry` writes the result
back onto the original representation.

## Constraints and scans

`ScanCoordinate` is a distance, angle or torsion between atom groups; an
endpoint group of one is an atom and a larger group is its centroid, so a
coordinate never needs a dummy atom. `Constraint` holds a coordinate while the
rest relaxes by projecting its gradient out, rather than freezing the atoms
that define it.

`relaxed_scan` steps the coordinate, constrains it, relaxes everything else,
and records `ScanPoint { value, energy, structure, converged }`. The driver
writes each point before starting the next and reports a typed error if the
scan axis is invalid.

## Equation of state

`fit_eos(points, EosForm)` fits energy-volume points with Birch-Murnaghan or
Murnaghan and returns `EosFit { v0, e0, b0, b0_prime, rms }`. `b0` is in eV per
Angstrom cubed internally and converted to GPa at the boundary. The points
normally come from a constant-volume scan in which the shape and atoms are
free.
