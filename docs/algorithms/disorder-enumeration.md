# Disorder enumeration

## Partial occupancies

An occupancy below one and a shared `disorder_group` describe alternative
occupants of one position. The P1 expansion keeps both, so the same structure
can be analysed without first choosing a component. The composition is the
occupancy-weighted sum over the P1 atoms; the density uses the same cell mass
and volume as the ordered path. Diffraction weights are the occupancy-weighted
scattering amplitudes, so a partially occupied site contributes its expected
amplitude rather than one full atom.

## Supercells

The Hermite normal form classifies every supercell of a given index up to a
unimodular change of basis. `hnf_supercells(index)` returns one integer matrix
per class, in a deterministic order. The index limit (128) keeps the
enumeration bounded; the caller can request a smaller set.

## Decorations

`derivative_structures` enumerates the ways to place the requested species on
the supercell sites subject to the concentration window. It:

1. expands the parent cell to the HNF supercell;
2. enumerates the combinations deterministically (lexicographic by site
   index, not by random sampling);
3. rejects a combination outside every concentration interval;
4. removes symmetry-equivalent decorations with the standardized matcher;
5. stops at `max_structures`.

The result stores the supercell matrix, the decoration (site index -> species)
and the number of symmetry operations that leave it invariant. A duplicate
removal that removes a structure lowers `symmetry_count`; the list reports the
survivors, not a post-hoc symmetry count.

## SQS

The binary SQS search keeps its pair/triplet correlation objective.
Multicomponent SQS adds a target for each pair, triplet and quadruplet shell
and an explicit seed. The greedy search proposes swaps, keeps the swap that
reduces the summed squared deviation and never accepts a worse target; the
reported correlations are the achieved values, so a caller can see the
residual. Determinism comes from the seed and the fixed site order.

## Cluster expansion

A `ClusterLattice` stores the neighbour shells (pair, triplet, quadruplet) of
a lattice. A decoration maps each site to a species; the correlation of a
cluster is the product of the species signs (+1/-1 for a binary alloy, an
orthogonal code for a multicomponent one). `cluster_correlation_vectors`
returns one vector per decoration, with the shell and cluster order preserved.
The vectors are the input to an external fit; no coefficients are estimated
here.

## Limits

- HNF index <= 128 and cluster order <= 4.
- The SQS search is greedy, not a guaranteed global optimum.
- Duplicate removal uses the matcher tolerance, so a very loose tolerance can
  collapse genuinely different decorations.
