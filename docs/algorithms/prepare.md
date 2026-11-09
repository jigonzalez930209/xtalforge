# Prepare for simulation

A refined crystal describes an average over all cells. Several things true of
that average are not true of any one cell, so preparation is an explicit,
ordered pipeline. Every step is recorded in `PrepareReport::steps`, even when
it is disabled, with the atoms or molecules it changed.

## Order

1. **duplicates** — merge sites with the same element at the same periodic
   position; the keeper's occupancy is the sum, clamped to one.
2. **disorder** — group partial sites by `_atom_site_disorder_group`; an
   ungrouped partial site is its own component unless its periodic position is
   within 2 A of another. Overlapping components are reduced to the one with
   the highest total occupancy (lowest group number on a tie), and survivors
   are set to full occupancy.
3. **solvent** — remove isolated molecules whose heavy-atom formula is a known
   solvent (water, methanol, DMF, DMSO, THF, dichloromethane, ...). A component
   containing a metal is framework and stays; anything unrecognised is left in
   place.
4. **minimize**, **charges**, **ewald**, **vibrations** — the existing force
   field stages.
5. **charge_balance** — sum the P1 site charges and report whether the cell is
   neutral within `1e-6 e`.

The cleanup steps are opt-in (`merge_duplicates`, `order_disorder`,
`remove_solvent`, `check_charge_balance`) so a caller can prepare only what the
file needs. The result is in `PrepareReport::structure`; `removed_sites` counts
every site dropped by the first three steps.

## Validation

`crates/xtalforge-ff/tests/prepare.rs` covers duplicate merging, a synthetic
two-component disorder, solvent removal and a charged cell. The prepared
samples in `tests/corpus/prepared/` are read back and checked to be fully
occupied with a reported charge balance; this is the golden sanity pass for
the cleanup steps. The exact Crystal-Builder disorder planner is richer than
this first pass: it considers graph fragments, disorder across symmetry
elements and charge-driven choices, which are documented rather than silently
approximated.
