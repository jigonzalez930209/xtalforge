# Defects, orderings and interfaces

**Tasks 9.5.1-9.5.4.** References: pymatgen `analysis.defects` and
`transformations`; Hart and Forcade (2009); Zunger, Wei, Ferreira and
Bernard (1990); Zur and McGill (1984).

## Point defects

A symmetry search of the P1 cell gives the orbits; one representative per
orbit is enough because the other members differ by a symmetry operation.
A vacancy removes the representative, a substitution replaces its element and
an interstitial inserts a species at a candidate position.

Interstitial candidates are the Voronoi vertices of Task 9.3.1, relative to
the central atom, filtered by a minimum clearance to every existing atom and
deduplicated under the detected symmetry operations. The defect supercell
repeats each lattice vector by `ceil(min_image / d_i)`, where `d_i` is the
interplanar spacing, so the shortest distance between defect images reaches
the requested value.

## Orderings of partial occupancies

A disordered site is a group of atoms at the same fractional position with
different elements and partial occupancies (or an explicit
`disorder_group`). The group is replicated into the requested supercell and
every combination of species is enumerated up to `max_orderings`. The
structure matcher of Task 9.1 removes decorations related by a lattice
translation, rotation or origin shift. Each remaining ordered cell is ranked
by the Ewald energy:

```
E = 1/2 sum_ij q_i q_j erfc(alpha r_ij) / r_ij
  + 2 pi / V sum_k exp(-k^2 / 4 alpha^2) / k^2 |S(k)|^2
  - alpha / sqrt(pi) sum_i q_i^2
```

with `S(k) = sum_i q_i exp(i k . r_i)` and the Coulomb constant
`14.399645 eV Angstrom`. The cell must be charge neutral; a non-neutral
decoration reports `energy: null` instead of a meaningless number.

## Special quasirandom structures

The binary disorder is decorated with a seeded SplitMix64 generator and fixed
species counts. Pair correlations are `xi_s = <sigma_i sigma_j>` with
`sigma = +/-1` for each distance shell; the triplet correlation averages
`sigma_i sigma_j sigma_k` over nearest-neighbor triangles. Simulated annealing
swaps two opposite sites, accepts improvements and, with a decreasing
temperature, some uphill moves. The objective is
`w_pair sum xi_s^2 + w_triplet xi_triplet^2`, whose target is zero.

## Interfaces

Zur-McGill matching enumerates integer 2x2 supercells of the two in-plane
cells up to `max_index`, rejects determinants that leave the area bound,
and compares the transformed lengths and angles. A match requires the maximum
relative length mismatch and the angle mismatch to stay inside the tolerances;
matches are ranked by strain and then area.

Stacking maps the film Cartesian coordinates onto the substrate in-plane
cell (which applies the film strain), places the film `gap` Angstrom above the
substrate and adds the requested vacuum. The in-plane shift is a translation
in substrate fractional coordinates.
