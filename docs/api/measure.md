# Measures, scalars and descriptions

`xtalforge_core::measure` provides `distance`, `angle`, `torsion`, `centroid`
and `interplanar_angle`. All atom-based functions expand the structure to P1
and use minimum images, so a bond across a cell boundary measures the same as
one inside it.

```rust
let d = measure::distance(&structure, 0, 1)?;
let a = measure::angle(&structure, 0, 1, 2)?;
let t = measure::torsion(&structure, 0, 1, 2, 3)?;
```

`scalars` exposes `bond_lengths`, `coordination_numbers`, `atom_charges` and
`atom_masses`, driven by the perceived bond graph. `describe_atom` returns one
line with element, label, site, coordination, fractional coordinates, charge
and partial occupancy, used by tooltips, the CLI and the agent layer.

These are the first entries of the scalar registry; colouring and analysis
name them instead of reimplementing the measurements.
