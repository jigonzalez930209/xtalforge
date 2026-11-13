# Molecular dynamics

`xtalforge-ff` adds deterministic molecular dynamics on top of any
`Calculator` that returns energy and forces (the native UFF calculator is the
default). Every run takes an explicit seed; the same structure, options and
seed produce bit-identical trajectories.

## Ensembles

| Ensemble | Integrator | Thermostat/barostat |
|---|---|---|
| NVE | velocity Verlet | none |
| NVT Langevin | velocity Verlet | friction `gamma` and seeded Gaussian noise |
| NVT Nose-Hoover | velocity Verlet | thermostat mass from the target temperature and coupling time |
| NPT | velocity Verlet + Parrinello-Rahman cell | Nose-Hoover thermostat and barostat |

The timestep is in femtoseconds at the API boundary. NPT uses the same
point-group-averaged strain as the optimizer, so a cubic cell keeps one cell
variable and a triclinic cell six.

## Output

`MdResult` carries the final structure, the per-step energy, temperature,
pressure and the accepted trajectory. `MdOptions.trajectory` selects whether
the trajectory is kept in memory or streamed to a writer. The integration is
symplectic where the ensemble allows it, and every reduction (kinetic energy,
temperature, pressure) uses a fixed summation order.

## Validation

- NVE energy drift over 1000 steps is below 1e-6 eV per atom for a harmonic
  test potential.
- Langevin and Nose-Hoover reach the target temperature within the documented
  statistical tolerance and reproduce the same sequence for the same seed.
- NPT on a cubic UFF cell keeps the cell cubic by symmetry.
- The MD trajectory is deterministic across thread counts.

## Limits

The force field is the native UFF/UFF4MOF implementation; external engines
are not started from the MD loop. The barostat is isotropic on the
point-group-averaged strain, not a full Parrinello-Rahman cell matrix.
