# Extinction screening

## What it decides

An indexed pattern fixes a cell and a metric family, but not the space group:
the same cell is compatible with several groups whose general positions make
different reflections systematically absent. `space_group_candidates(system,
hkl)` rejects every group that would make an observed indexed reflection
impossible; `extinction_classes(system, hkl)` groups the survivors by their
full reflection-condition signature so the automatic run can fit one
representative per class.

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

## Extinction classes

A powder pattern shows the extinction symbol; the groups that share it differ
only by elements that produce no absence. `extinction_classes(system, hkl)`
enumerates the catalogue groups of the cell's crystal system, computes each
group's reflection-condition signature over a fixed `hkl` grid, and merges the
groups with the same signature. The representative is the group with the most
operations (rutile's class is represented by `P 4_2/m n m`, not by one of its
subgroups). A class is **refuted** when any observed indexed reflection is
systematically absent in it.

## Ranking

The automatic run ranks the surviving classes by a Le Bail fit of a structure
with one general-position atom in the class's representative group:

```text
BIC = n ln(RSS / n) + k ln(n),   delta_BIC = BIC_class - BIC_lattice,
```

with `RSS` the weighted residual sum of the Le Bail profile, `k` the number of
extracted reflections plus the scale and background terms, and the lattice
class the absence-free class of the same cell. The leading classes are then
Pawley fitted; the table is ordered by Pawley Rwp, with fits within 1 % of each
other treated as a tie and ordered by the indexing and class ranks.

The class symbol in the row is the representative's Hermann-Mauguin symbol (not
the standard extinction symbol); its `space_groups` list names every group in
the class.

## Validation

`crates/xtalforge-analysis/tests/extinction.rs` checks:

- rutile `P4_2/mnm` (136) is the first row of the golden rutile table, and
  an impossible `001` removes it from the candidates;
- diamond `Fd-3m` (227) is among the candidates of a diamond pattern, and an
  impossible `200` removes it;
- the class enumeration of rutile puts 136 in the `P 4_2/m n m` class and
  the golden automatic table selects it first;
- the reflection conditions of 136 and 227 match the values computed from
  their operations: `00l` even, even index sums for the rutile glide planes,
  the `Fd-3m` diamond-glide conditions, and the fact that 222 is allowed by
  the `Fd-3m` general position (its absence in diamond is a special-position
  effect of the 8a site).

`tests/automatic.rs` runs the whole automatic pipeline on a rutile pattern
and requires the table to contain the space-group candidates.
