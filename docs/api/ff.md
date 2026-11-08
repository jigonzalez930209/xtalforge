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
