# Extinction screening

## What it decides

An indexed pattern fixes a cell and a metric family, but not the space group:
the same cell is compatible with several groups whose general positions make
different reflections systematically absent. `space_group_candidates(system,
hkl)` rejects every group that would make an observed indexed reflection
impossible and returns the survivors.

## Structure-factor test

For a group with operations `(R, t)` in the conventional cell, the structure
factor of a single atom on a general position `r` is

```text
F(hkl) = sum over operations of exp(2 pi i h . (R r + t))
```

This sum is zero for **every** `r` exactly when the reflection is
systematically absent. Two generic fractional points are evaluated; their
normalised `|F|` are below `1e-6` when the group cancels the reflection. The
operations include the centering translations, so the Bravais absences (`F`,
`I`, `C`, `R`) come out of the same sum as the glide and screw absences.

`is_systematically_absent(operations, hkl)` is the public form of that test.
`space_group_candidates` enumerates the catalogue for the cell's crystal
system (trigonal groups are included for a hexagonal cell), keeps the first
Hall setting of each ITA number, counts the indexed reflections each group
forbids and returns the groups with the fewest violations (zero when a group
is perfect). `MAX_CANDIDATES` limits the list; the ranking is by number of
operations, so the highest-symmetry compatible group comes first.

## Why the ranking is not a fit

The reference `rietx` screen fits a Le Bail class model for every extinction
class and ranks the classes by fit quality. This port does not: it uses only
the reflection conditions, which are exact for the indexed list, and orders
the survivors by group order. Two groups with identical absences therefore
keep their catalogue order, and a group that fits the intensities better can
appear below one that fits them worse. The result is a candidate list, not a
decision.

The screen also cannot fix a bad indexing ranking: the automatic table only
screens the candidates the indexing supplied. If the true cell is not in the
leading rows, its space group is not in the table either. The golden rutile
pattern is an example: the indexing ranking puts hexagonal supercells first,
so the tetragonal rutile cell only appears when it is among the enumerated
candidates.

## Validation

`crates/xtalforge-analysis/tests/extinction.rs` checks:

- rutile `P4_2/mnm` (136) is among the candidates of a real rutile pattern,
  and an impossible `001` removes it;
- diamond `Fd-3m` (227) is among the candidates of a diamond pattern, and an
  impossible `200` removes it;
- the reflection conditions of 136 and 227 match the values computed from
  their operations: `00l` even, even index sums for the rutile glide planes,
  the `Fd-3m` diamond-glide conditions, and the fact that 222 is allowed by
  the `Fd-3m` general position (its absence in diamond is a special-position
  effect of the 8a site).

`tests/automatic.rs` runs the whole automatic pipeline on a rutile pattern
and requires the table to contain the space-group candidates.
