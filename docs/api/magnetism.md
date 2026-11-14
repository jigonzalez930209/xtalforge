# Magnetism and magnetic structures

## Moments on sites

`MagneticMoment` is an optional site property. It is either collinear
(`MagneticAxis::{A, B, C}` with a signed magnitude in Bohr magnetons) or a
three-component vector in Cartesian coordinates. Sites keep the moment through
`StructureBuilder`, CIF and mCIF I/O, the P1 expansion, supercells and slabs.

```rust
use xtalforge::core::{MagneticAxis, MagneticMoment, Site};

let moment = MagneticMoment::collinear(MagneticAxis::C, 2.2)?;
let site = Site::new("Fe1", "Fe", [0.0, 0.0, 0.0])?
    .with_magnetic_moment(moment);
```

Python exposes `Site.magnetic_moment` and the `with_magnetic_moment` helper;
the structure document carries the same additive field, so the schema version
is unchanged.

## mCIF

`xtalforge_io::read_mcif` / `write_mcif` read and write the magnetic CIF loop
fields:

- `_space_group_symop_magn_operation` and `_space_group_symop_magn_centering`
  (the operations are stored in the structure metadata under
  `magnetic_operations`);
- `_atom_site_moment_crystalaxis_x/y/z` and the magnitude form
  `_atom_site_moment_label`/`_atom_site_moment_crystalaxis_*`;
- the ordinary cell, symmetry, site and occupancy fields.

The reader accepts a block selector (`read_mcif_block`) and the writer keeps
the moments in crystal-axis components. The format is registered as `mcif`
when the registry can resolve the extension; otherwise call the functions
directly.

## Magnetic symmetry

`detect_magnetic_symmetry` (and the `_with_options` form) uses moyo's magnetic
cell support to return the UNI number, BNS/OG symbol, operations and the
moment action. `symmetrize_moments`/`symmetrize_structure_moments` apply the
operations to the moments and return a `MagneticSymmetrization` with the
symmetrised array and the number of operations that leave it invariant.
`magmom_values` and `format_magmom` produce the VASP `MAGMOM` vector/string
used by the Phase 18 input sets.

## Orderings

`MagneticOrderingOptions` enumerates ferro (`MagneticOrderingKind::Ferromagnetic`),
antiferro and ferri collinear orderings for a chosen supercell, ranks them by
magnetic symmetry and returns
`MagneticOrdering { kind, structure, moments, net_moment, uni_number,
operation_count, magmom }`. `magnetic_orderings` uses a cubic repeat;
`magnetic_orderings_in_supercell` takes an explicit cell. The moments are ready
for `MAGMOM`. The enumeration is deterministic and refuses a non-collinear
input rather than silently dropping components.

Parity keys: `magnetic_moments`, `mcif`, `magnetic_symmetry`, `magnetic_orderings`.
