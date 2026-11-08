# Structures and sites

A `Structure` is a lattice, a space-group identity (the Hall symbol) and the
asymmetric-unit sites. `Structure::p1()` derives the full cell with provenance:
every `P1Atom` carries its parent site, the operation that generated it and the
lattice shift. A structure with no operations is already P1 and reports
`is_p1() == true`.

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
occupied. The chainable methods are `with_occupancy`, `with_frac`, `with_u_iso`,
`with_u_aniso`, `with_charge`, `with_disorder_group` and `with_prop`; the
readers are `label`, `element`, `frac`, `occupancy`, `u_iso`, `u_aniso`,
`charge`, `disorder_group`, `props` and `is_partial`.

## Structure

`StructureBuilder::new(name, lattice)` starts a P1 structure; `hall`,
`operation`, `site` and `metadata` add data, and `build` produces the
immutable `Structure`.
Reads are `name`, `lattice`, `hall`, `sites`, `operations`, `metadata`, `is_p1`,
`cart_sites`, `cell_contents`, `formula`, `cell_mass`, `molar_mass`,
`density`, `z`, `total_charge`, `charge_balance` and `summary`.
`with_name`, `with_lattice` and `with_site` return modified copies.
`Structure::molecule(name, points, padding)` builds a non-periodic molecule
in a padded P1 box.

`StructureSummary` carries `name`, `formula`, `hall`, `n_sites`, `volume`,
`density`, `total_charge` and `is_p1` to the four surfaces. `structure_info`
accepts CIF, POSCAR, XYZ and pymatgen JSON through the format registry and
returns the same summary.

## Conventions

- Formula counts include occupancy and follow first appearance in the sites.
- Density uses the cell contents: `mass * 1.66053906660 / volume` g/cm3.
- `z` is the number of formula units in the cell contents; a structure with
  operations derives it from the P1 expansion.
- `total_charge` is the occupancy-weighted sum of formal charges.
