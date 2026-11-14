# CIF input and output

`xtalforge_io::cif` reads and writes CIF data blocks. `read_cif(text)` returns
the first block as a `Structure`; `read_cif_block(text, name)` selects one;
`read_cif_all(text)` returns every block. `write_cif(structure)` writes the
cell, Hall symbol, operation loop and atom sites.

```rust
let structure = xtalforge::io::read_cif(&text)?;
let p1 = structure.p1()?;
let output = xtalforge::io::write_cif(&structure)?;
```

`CifDocument::parse` exposes the raw parse as `CifDocument { blocks }`, with
`CifBlock` for a data block and `CifLoop` for its loops; the `read_cif*`
functions are the typed shortcut.

The parser handles data blocks, tags, loops, single and double quoted values,
semicolon-delimited text, comments, standard uncertainties (`1.234(5)`), `?`
and `.` values, atom sites, occupancy, Uiso/Uij, formal charges (including the
`_atom_type_oxidation_number` fallback), disorder group and assembly, and
symmetry operation loops (`_symmetry_equiv_pos_as_xyz` or
`_space_group_symop_operation_xyz`). The operation loop is the setting when
present; without one, the Hall or H-M symbol is resolved through the symmetry
catalogue. Only an unknown symbol, or an operation loop whose count or
operation set disagrees with the named symbol, records a `symmetry_warning`. The anisotropic loop is
matched by `_atom_site_aniso_label`, so a subset of atoms is read correctly.

`decode_text(bytes)` strips a UTF-8 BOM, reads UTF-8 and falls back to Latin-1
byte by byte. `read_text(path)` is the file wrapper.

## Magnetic CIF

`xtalforge_io::read_mcif` / `write_mcif` (and `read_mcif_block`) handle the
magnetic CIF dialect: `_space_group_symop_magn_operation_xyz` with the
`+1`/`-1` time-reversal suffix, `_space_group_symop_magn_centering_xyz`, and
`_atom_site_moment_crystalaxis_*` / `_atom_site_moment_cartn_*` loops. The
operations are stored in `magnetic_operations`; the moments become typed
`Site::magnetic_moment` values. The registry sniffs a magnetic loop before
plain CIF, so `.mcif` and a hint-less `read_structure` select it. See
[Magnetism](./magnetism.md).

## Status

The reader covers the corpus CIFs and malformed input returns a typed error.
The writer emits the cell, Hall symbol, operation loop, occupancy, Uiso, Uij,
disorder group/assembly, charge and the explicit-bond loop. Explicit bonds are
written with the program's `_xtal_bond_kind`/`_geom_bond_*` tags and read back
for P1 structures. A non-P1 row is mapped through its operation and image to
P1 contacts (MIL53 matches the 36 golden contacts); an unmatched contact is
not added to the graph. A non-P1 explicit set is written back as a
`_xtal_bond_p1_count` plus `_xtal_bond_atom_1`/`_xtal_bond_atom_2` index loop,
and the reader checks the declared count against its own expansion. The fuzz target lives in
`crates/xtalforge-io/fuzz`; inconsistency warnings beyond the unknown-symbol
warning are not implemented.
