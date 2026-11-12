# Structure comparison algorithms

**Tasks 9.1.1-9.1.3.** References: pymatgen `StructureMatcher` and
`analysis.structure_matcher`; Niggli (1928); Kabsch (1976) and Horn (1987).

## Matcher

1. Both structures are expanded to P1 and their lattices are Niggli reduced
   (`moyo` through `reduce_lattice`). The fractional coordinates are moved by
   the Cartesian mapping, so the reduction never changes geometry.
2. An integer matrix `M` is searched with `M * L2 ~ L1`: every row of `L1` is
   expressed in the basis of `L2`, rounded to the nearest integer vector, and
   accepted when the fractional error is at most `ltol`. All combinations with
   the required determinant are kept (up to 32, ordered by a length and angle
   score).
3. For each candidate, the second cell is mapped into the first and a
   translation is searched by anchoring the first four atoms of the first
   structure on every species-compatible atom of the second. Atoms are then
   assigned greedily to the nearest unmatched partner within `stol`; the
   assignment must cover every atom.
4. The best RMS displacement over all candidates is returned. Supercell
   matching uses a determinant equal to the atom-count ratio (and the volume
   ratio must agree), and the smaller structure is expanded before matching.

Tolerances default to `ltol = 0.2`, `stol = 0.3` Angstrom and
`angle_tol = 5` degrees. Species are compared exactly, by canonical element,
or ignored.

## RMSD and alignment

The periodic RMSD uses the minimum image of the first cell for the fractional
difference of each mapped pair. Kabsch alignment uses Horn's quaternion
method: the 4x4 matrix built from the cross-covariance is diagonalized with
Jacobi rotations and its largest eigenvector gives the rotation. The result is
deterministic and needs no external linear-algebra routine.

## Equality classes

Structures are visited in input order and joined to the first class whose
representative matches; this gives a deterministic partition independent of
hash iteration order.
