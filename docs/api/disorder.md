# Disorder, derivative structures and cluster expansion

## Partial occupancies

`Site::occupancy` is in `[0, 1]`; `Site::disorder_group` identifies the
alternative occupants of one position. The P1 expansion carries both fields,
so every operation that expands, supercells or slices a structure keeps the
disorder information.

`xtalforge_analysis::disorder_groups(structure)` returns typed
`DisorderGroup { position, occupants, total_occupancy }` records.
`disorder_composition` aggregates the element occupancies and returns a
`DisorderComposition` with `formula()`, `cell_mass()` and `density(volume)`.
`disorder_formula` and `disorder_density` are the one-call helpers.
`disorder_scattering_weights` returns the per-element X-ray/neutron weights
used by the diffraction code.

## Derivative structures

`hnf_supercells(index)` enumerates the Hermite-normal-form supercells of a
given index in a deterministic order (index up to 128). `derivative_structures`
applies a species decoration under concentration limits and removes
symmetry-equivalent results; `DerivativeOptions` carries the species,
concentrations, supercell index and `max_structures`, and
`DerivativeStructure { supercell, decoration, symmetry_count }` is the result.
No random seed is taken, so the same request gives the same list.

## SQS

`special_quasirandom` keeps its binary API and adds multicomponent
correlation targets: pair, triplet and quadruplet shells plus an explicit
seed. `SqsResult` reports the structure, the achieved correlations and the
target list. The search is deterministic for a seed and reports the deviation
from each target rather than hiding it.

## Cluster-expansion features

`ClusterLattice::new` builds the neighbour shells of a decorated lattice and
`cluster_correlation_vectors(&decorations, &options)` returns a
`ClusterCorrelationSet` with one `ClusterCorrelationVector` per decoration.
Each vector stores the shell label, cluster order and correlation value for
external fitting. No fitting engine is shipped; this is the feature export.

Parity keys: `disorder`, `derivative_structures`, `sqs`, `cluster_expansion`.
