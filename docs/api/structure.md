# Structures and sites

A `Structure` is a lattice, a space-group identity (the Hall symbol) and the
asymmetric-unit sites. The current Phase 1 model stores the cell contents in
the asymmetric unit and reports P1; P1 expansion with provenance arrives in
Task 1.6, after which the same type carries the asymmetric unit and its derived
cell.

## API

::: code-group

```rust
use xtalforge::core::{Lattice, Site, StructureBuilder};

let cell = Lattice::cubic(5.64)?;
let structure = StructureBuilder::new("NaCl", cell)
    .site(Site::new("Na1", "Na", [0.0, 0.0, 0.0])?)
    .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5])?)
    .build()?;
let summary = structure.summary()?;
```

```bash
xtalforge info structure.json
xtalforge --json info structure.json
```

```python
summary = xtalforge.structure_info(open("structure.json").read())
```

```js
import { structure_info } from "xtalforge-wasm";
const summary = structure_info(jsonText);
```

:::

## Site

`Site::new(label, element, frac)` normalizes the element (`Fe2+` becomes
`Fe`), validates a non-empty label and finite coordinates, and starts fully
occupied. The chainable methods are `with_occupancy`, `with_u_iso`,
`with_u_aniso`, `with_charge`, `with_disorder_group` and `with_prop`; the
readers are `label`, `element`, `frac`, `occupancy`, `u_iso`, `u_aniso`,
`charge`, `disorder_group`, `props` and `is_partial`.

## Structure

`StructureBuilder::new(name, lattice)` starts a P1 structure; `hall`, `site`
and `metadata` add data, and `build` produces the immutable `Structure`.
Reads are `name`, `lattice`, `hall`, `sites`, `metadata`, `is_p1`,
`cart_sites`, `cell_contents`, `formula`, `cell_mass`, `molar_mass`,
`density`, `z`, `total_charge`, `charge_balance` and `summary`.
`with_name`, `with_lattice` and `with_site` return modified copies.
`Structure::molecule(name, points, padding)` builds a non-periodic molecule
in a padded P1 box.

`StructureSummary` carries `name`, `formula`, `hall`, `n_sites`, `volume`,
`density`, `total_charge` and `is_p1` to the four surfaces. The Phase 1 JSON
bridge is `structure_info(text)`; CIF input replaces it in Task 1.15 and keeps
the same summary.

## Conventions

- Formula counts include occupancy and follow first appearance in the sites.
- Density uses the cell contents: `mass * 1.66053906660 / volume` g/cm3.
- `z` is 1 while the asymmetric unit is the cell contents; the P1 expansion
  refines it.
- `total_charge` is the occupancy-weighted sum of formal charges.
