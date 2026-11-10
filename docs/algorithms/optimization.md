# Optimization, scans and equations of state

## Degrees of freedom

`RelaxMode::Asymmetric` is the default for a non-P1 structure: the variables
are the asymmetric sites, the calculator works on the expanded P1 cell, and
the forces are projected back onto the symmetry-allowed displacements. The
space group is preserved on every step. `RelaxMode::P1` is the ordinary
Cartesian relaxation.

`OptimizeOptions.cell` (a boolean) adds the strain components as variables (averaged over the point group: a
cubic cell keeps one free strain, a triclinic cell six), so the cell relaxes with the atoms. Holding selected lengths or the
volume while the shape relaxes (the `CellFreedom` idea of the reference) and an
external pressure term are **not implemented** in this release; they are
recorded as known subsets in the roadmap. `OptimizeOptions.fixed_atoms` freezes
the listed atoms.

## Minimizers

![FIRE and L-BFGS](/diagrams/optimization.svg)


- **FIRE** (`Minimizer::Fire`): molecular dynamics with the velocity steered
  downhill, an adaptive time step and a force criterion.
- **L-BFGS** (`Minimizer::Lbfgs`): a limited-memory quasi-Newton method with a
  line search and restarts.

Both stop on `force_tolerance` (eV/Angstrom), `energy_tolerance` (eV) or
`max_steps`. `OptimizationResult { structure, energy, max_force, steps,
converged }` is deterministic. `apply_optimized_geometry` writes the result
back onto the original representation.

## Constraints and scans

`ScanCoordinate` is a distance, an angle or a torsion between **single atoms**
(P1 site indices, minimum image through the lattice). Groups of atoms with a
centroid endpoint, the angle between planes and the cell parameters or volume as
coordinates are not implemented. `Constraint` holds a coordinate while the rest
relaxes by projecting its gradient out, rather than freezing the atoms that
define it; `OptimizeOptions.fixed_atoms` freezes whole atoms.

`relaxed_scan` steps the coordinate, constrains it, relaxes everything else
starting from the previous point, and returns
`ScanPoint { value, energy, structure, converged }` for every value in memory. It
does not write `scan.csv` or per-point files and cannot resume an interrupted
scan; a caller that wants that writes each returned point itself. An invalid axis
(non-P1 structure, empty values, a degenerate angle) is a typed error.

The optimizers take no cancellation token or progress callback and write no
trajectory; they stop on `force_tolerance`, `energy_tolerance` or `max_steps`.

## Equation of state

`fit_eos(points, EosForm)` fits energy-volume points with Birch-Murnaghan or
Murnaghan and returns `EosFit { v0, e0, b0, b0_prime, rms }`. `b0` is in eV per
Angstrom cubed internally and converted to GPa at the boundary. The points
normally come from a constant-volume scan in which the shape and atoms are
free.
