# Generation and global search algorithms

This page records what the Phase 20 generators do, which approximations they
make and how the acceptance checks work.

## Wyckoff-constrained crystals

Given a composition and an ITA number, the generator:

1. Lists the Wyckoff positions of the group (`moyo` tables through
   `xtalforge_symmetry::WyckoffPosition`).
2. Finds combinations of at most four distinct positions whose multiplicities
   add up to each element count. The composition may be scaled by a common
   factor up to six to match the conventional cell (for example one formula
   unit of NaCl in a face-centred cell).
3. Evaluates each coordinate template at random free parameters. The template
   grammar is a sum of rational constants and at most one variable per
   component: `x`, `1/2-x`, `-y+1/4`, `2x`, `x-y`.
4. Builds one representative site per orbit and keeps the Hall symbol and the
   explicit operations, so the normal P1 expansion produces the cell.
5. Rejects a candidate when an orbit collapses to fewer images than its
   multiplicity, the periodic minimum distance is below `min_distance`, or a
   symmetry search returns a different ITA number.

A generic free parameter gives the requested group; accidental extra symmetry
is rejected rather than repaired. The lattice follows the crystal system and a
target volume `volume_per_atom * N`; axial ratios are drawn from `[0.85, 1.18]`
and the monoclinic angle from `[96°, 124°]`. Cubic cells have no free ratio, so
two seeds differ only through the free parameters of the Wyckoff orbits (a
structure with only fixed positions is the same for every seed).

### Seed statistics

`generation_statistics` runs the generator over a seed range and counts the
accepted crystals and the mean attempts per success. The repository test
`wyckoff_generation.rs` covers rock salt, hexagonal Mg, a range of groups and
the determinism of one seed; the algorithm page does not claim a 10 000-seed
run in CI. The recommended check is:

```bash
cargo test -p xtalforge-generation --test wyckoff_generation
```

## Molecular crystals

The molecule is kept rigid. Each independent copy is placed at a Wyckoff
centre and rotated by a uniform random quaternion. The generator expands all
operations, groups images by centre, and requires:

- every operation that maps a centre onto itself maps the atom set onto itself
  within `tolerance` (greedy element-aware matching);
- no two images from different centres are closer than `min_distance`.

If either check fails the attempt is redrawn. The output keeps the requested
operations, so the P1 expansion reproduces the images.

## Layer and rod groups

Layers use moyo's 80 layer groups. The cell is `a`, `b` in the layer plane and
`c = 20 Å`; `a`/`b` follow the two-dimensional Bravais class (`mp`, `op`, `oc`,
`tp`, `hp`) and the target area per atom. Acceptance requires the minimum
distance, the operation set preserving the atoms, and
`detect_layer_symmetry` returning the requested layer-group number.

moyo has no rod tables. The rod catalogue lists fourteen operation sets
(identity, inversion, 2, 3, 4 and 6-fold rotations about the chain axis and
their mirrors). A rod is generated in a 20 Å × 20 Å × repeat cell; the
operation set is verified directly on the expanded atoms. This is a symmetry
check, not a 3D space-group approximation.

## AIRSS-style random search

For each candidate `i` the search calls `generate_crystal` with
`seed = base + i`, relaxes with the selected calculator (UFF4MOF by default;
FIRE or L-BFGS), optionally relaxes the cell, builds a `Dataset`, deduplicates
with the Phase 16 matcher and sorts by energy per atom. Ties break by
candidate index, so the ranking is deterministic.

## Evolutionary search

1. The initial population comes from the random search or from a checkpoint.
2. The best `ceil(population * survival)` survive.
3. Each child is cut-and-splice crossover in the first parent's P1 cell: a
   random plane splits both parents and the child takes half the atoms from
   each (mapped into the first cell).
4. With probability `mutation_rate`, each coordinate is displaced by
   `±mutation_sigma` and each lattice vector is strained by
   `±strain_sigma`.
5. The child is relaxed, rejected if it matches a population member, and the
   population is truncated back to `population`.
6. The loop stops after `convergence_generations` without a `1e-4 eV/atom`
   improvement.

Because the child is a P1 assembly, the search does not impose a space group
on it. A caller that needs a symmetric child should symmetrize after the
search.

## Prototype identification

The signature is the ITA number, the sorted `(multiplicity, Wyckoff letter)`
pairs and the sorted stoichiometry counts. Candidates with the same signature
are built from the library, their species are remapped to the input in
first-appearance order, and the structure matcher reports the RMS. The
smallest RMS wins. Element identity is deliberately ignored, so KBr identifies
as the B1 rock-salt prototype.

## MOF enumeration

Combinations are shuffled with the seed and built in order until `limit`. Each
framework records topology, node, linker, structure, optional UFF energy,
framework descriptors, pore descriptors and an error string. Duplicates are
marked by MOFkey. The orientation fit behind the builder is the global Kabsch
fit from Phase 6: every assignment of connection points to incident slot
directions is tried (up to six, 720 permutations) and the smallest residual
wins.

## Determinism and limits

All draws come from a local SplitMix64 stream. Fixed reduction orders and
tie-breaks make the output independent of thread count. The parser is
bounded: 400 atoms per conventional cell, four Wyckoff orbits per element,
six composition scale factors, 64 attempts by default. The full PORMAKE
bundle is a feature (`mof-data`) because of package size.

## Validation

| Check | Test |
|---|---|
| Every generated crystal re-detects its group | `tests/wyckoff_generation.rs` |
| One seed is reproducible | `tests/wyckoff_generation.rs` |
| Molecular placement and P-1 detection | `tests/molecular_crystal.rs` |
| Layer and rod acceptance | `tests/lowdim_generation.rs` |
| Every prototype builds and re-detects | `tests/prototypes.rs` |
| Substitution rules | `tests/substitution.rs` |
| Ranking, deduplication and generations | `tests/search.rs` |
| Checkpoint resume and workflow document | `tests/search_workflow.rs` |
| MOF dataset, orientation and sweeps | `tests/mof_generation.rs` |
| Layer catalogue | `crates/xtalforge-symmetry/tests/layer_groups.rs` |
