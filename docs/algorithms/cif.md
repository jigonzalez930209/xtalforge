# CIF parsing

## Pipeline

1. `decode_text` strips a UTF-8 BOM, reads UTF-8 and falls back to Latin-1.
2. The tokenizer walks lines, removes comments outside quotes, keeps quoted
   values together and collects semicolon-delimited text blocks.
3. The parser builds data blocks, single items and loops. Loop values must be
   a multiple of the tag count.
4. `structure_from_block` reads the cell, Hall/H-M symbol, symmetry-operation
   loop, atom sites, occupancy, Uiso/Uij, charge, disorder group/assembly and
   the program's explicit-bond loop, and builds a `Structure`.

## Symmetry priority

The operation loop is the setting when present. Without a loop, the reader
resolves the Hall symbol (or H-M symbol) through
`xtalforge_symmetry::SpaceGroup::parse` and installs its operations, recording
`symmetry_source = "Hall symbol"`. Only an unrecognised symbol falls back to
the identity and records a `symmetry_warning`. This closes the common CIF case
where the space group is named but the operation loop is omitted (graphene is
the corpus example). When both an operation loop and a symbol are present, the
reader compares the full operation set (rotation and translation modulo one)
with the named group. A count or set disagreement records a
`symmetry_warning`.

## Numbers

`1.234(5)` is read as `1.234`; `?` and `.` become missing values. Malformed
numbers or truncated loops return `Error::InvalidInput`, never a panic. A
fuzz-style regression test feeds empty, truncated and mismatched inputs.

## Disorder and anisotropic parameters

`_atom_site_disorder_group` becomes `Site::disorder_group`; a non-numeric group
and `_atom_site_disorder_assembly` are kept in `Site::props` so a save does not
drop them. The anisotropic loop is matched by `_atom_site_aniso_label`, not by
row position, so a file that refines only some atoms still reads correctly.

## Explicit bonds

A `_geom_bond` loop is read only when it carries `_xtal_bond_kind`, the tag
this program writes. Published `_geom_bond` loops are usually distance tables
from a refinement, not a bond graph, and reading them silently would override
perception. The row names sites by label, carries the lattice image and order,
and `suppressed` removes a perceived bond. The writer emits the loop for P1
structures.

For a non-P1 structure the row's `_geom_bond_site_symmetry_2` code is the
operation number plus a three-digit translation offset, and the reader maps
label + operation + image through every symmetry operation to P1 contacts.
The duplicate images of a special-position pair are collapsed by their
minimum-image contact vector. MIL53 maps its two records to the 36 golden
contacts. An unmatched mapped record does not invent a graph edge: only a
contact perception already found is annotated, because the P1 representative
is the one the graph is drawn with.

The writer round-trips a non-P1 explicit set with
`_xtal_bond_p1_count` and the `_xtal_bond_atom_1`/`_xtal_bond_atom_2` index
loop (plus image, order and kind). The reader validates the declared P1 count
against its own expansion before accepting the indices. The label/operation
loop remains the readable form for P1; both are understood on input.

## Writer

`write_cif` writes the cell, Hall symbol, operation loop, sites, occupancy,
Uiso, Uij, disorder group/assembly, charge and the explicit-bond loop. The
round-trip tests read the written structure back.
