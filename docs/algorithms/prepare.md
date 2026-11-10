# Prepare for simulation

A refined crystal describes an average over all cells. Several things true of
that average are not true of any one cell, so preparation is an explicit,
ordered pipeline. Every step is recorded in `PrepareReport::steps`, even when
it is disabled, with the atoms or molecules it changed.

## Order

![Prepare for simulation](/diagrams/prepare-pipeline.svg)


1. **duplicates** — merge sites with the same element at the same periodic
   position; the keeper's occupancy is the sum, clamped to one.
2. **disorder** — group partial sites by `_atom_site_disorder_group`, or by
   their connected bond-graph component when no group is given. Overlapping
   groups are reduced to the compatible subset that best matches the declared
   `_chemical_formula_sum` (per metal when the formula has metals) and a
   neutral cell, with higher occupancy as the tie-break; survivors are set to
   full occupancy. The formula is kept in structure metadata by the CIF
   reader.
3. **solvent** — remove isolated molecules whose heavy-atom formula is a known
   solvent (water, methanol, DMF, DMSO, THF, dichloromethane, ...). A component
   containing a metal is framework and stays; anything unrecognised is left in
   place.
4. **hydrogens** — add missing hydrogens from the perceived valence
   (`xtalforge_build::add_hydrogens`), opt-in.
5. **minimize**, **charges**, **ewald**, **vibrations** — the existing force
   field stages.
6. **close_contacts** — list P1 pairs below the reference clash thresholds
   (H-H 1.4 A, O-O 2.0 A, halogen-halogen 1.9 A, other pairs `0.8 (r_i+r_j)`).
7. **charge_balance** — sum the P1 site charges and report whether the cell is
   neutral within `1e-6 e`.

The cleanup steps are opt-in (`merge_duplicates`, `order_disorder`,
`remove_solvent`, `hydrogens`, `check_close_contacts`,
`check_charge_balance`) so a caller can prepare only what the file needs. The
result is in `PrepareReport::structure`; `removed_sites` counts every site
dropped by the first three steps, and `close_contacts` carries the clash list.
A caller applies the prepared result as one undoable edit with
`xtalforge_commands::ReplaceStructure`; `Document` then exposes the usual
undo/redo.

## Validation

`crates/xtalforge-ff/tests/prepare.rs` covers duplicate merging, a synthetic
two-component disorder, solvent removal and a charged cell. The prepared
samples in `tests/corpus/prepared/` are read back and checked to be fully
occupied with a reported charge balance; this is the golden sanity pass for
the cleanup steps. The disorder choice now consults the declared formula and
the cell charge. The remaining difference from Crystal-Builder is its full
special-position/clash graph and shared-oxygen bookkeeping, which are
documented rather than silently approximated.
