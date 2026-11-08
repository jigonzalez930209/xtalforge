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
and `.` values, atom sites, occupancy, Uiso/Uij, charge, and symmetry
operation loops (`_symmetry_equiv_pos_as_xyz` or
`_space_group_symop_operation_xyz`). The symmetry priority is Hall symbol,
then the operation loop, then the H-M symbol/number; the current reader stores
the string but does not warn about inconsistent statements yet.

`decode_text(bytes)` strips a UTF-8 BOM, reads UTF-8 and falls back to Latin-1
byte by byte. `read_text(path)` is the file wrapper.

## Status

The reader covers the corpus CIFs and malformed input returns a typed error.
The writer emits cell, Hall symbol, operation loop, occupancy, Uiso and charge;
it does not write Uij loops or explicit bond loops yet, so those parts are not
round-tripped. The fuzz target lives in `crates/xtalforge-io/fuzz`; the exact
Crystal-Builder inconsistency warnings beyond the missing-operation warning are
not implemented.
