# POSCAR and XYZ

`read_poscar(text)` and `write_poscar(structure)` cover VASP POSCAR/CONTCAR:
comment, positive scale or negative target volume, lattice rows, species and
counts, `Selective dynamics`, and Direct or Cartesian coordinates. A file
without a species line is refused explicitly (VASP 4 carries species only in
POTCAR).

`read_xyz(text)` and `write_xyz(structure)` cover plain and extended XYZ. An
extended file with `Lattice="..."` keeps its cell; a plain file becomes a P1
box padded around the atoms. The writer always emits `Lattice="..."` and
`Properties=species:S:1:pos:R:3`, so a round trip preserves the cell.

`read_structure(text, hint)` and `write_structure(structure, format)` dispatch
by name or extension (`cif`, `poscar`, `vasp`, `xyz`). The registry of Phase 3
replaces the dispatch without changing these signatures.

```rust
let structure = xtalforge::io::read_structure(&text, "poscar")?;
let xyz = xtalforge::io::write_structure(&structure, "xyz")?;
```
