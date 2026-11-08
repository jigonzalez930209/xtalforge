# POSCAR and XYZ

POSCAR stores the lattice rows, a scale and fractional or Cartesian
coordinates. A negative scale is a target volume: the raw lattice is scaled by
`cbrt(target / |det|)`. A missing species line is refused because VASP 4 keeps
species only in POTCAR.

Extended XYZ stores Cartesian coordinates and a `Lattice="..."` attribute. The
reader converts positions through the lattice inverse; a plain XYZ without a
cell becomes a P1 molecule box padded by 5 Angstrom. The writer always emits
the extended form, so a round trip preserves the cell.

Format dispatch accepts names and extensions and is replaced by the registry
in Phase 3 without changing the signatures.
