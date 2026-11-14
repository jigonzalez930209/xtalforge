# Machine-learning representations and inference

This page records the Phase 21 algorithms, their conventions and the
validation evidence.

## Crystal graphs

The graph uses the P1 cell in canonical order: atomic number, then fractional
`x`, `y`, `z`, then label. The neighbour list is the exact linked-cell search
of Task 1.10. Every pair contributes two directed edges; each edge stores the
lattice image of the target, the Cartesian vector and the distance. With
`max_neighbors = k` the shortest `k` outgoing edges per node survive and
`truncated` is set. `edge_index` is `[source, target]` in the same order as
`edge_vectors` and `distances`.

Python converts the arrays with `crystal_graph_numpy`; WASM returns JSON
arrays for typed-array construction in JavaScript.

## SOAP

The power spectrum follows the DScribe/`SOAPGTO` structure:

1. `a_k = linspace(1, r_cut, n_max)`.
2. `alpha_lk = -ln(1e-3 / a_k^l) / a_k^2`.
3. Overlap `S_l[k,k'] = 0.5 Gamma(l + 3/2) (alpha_lk + alpha_lk')^(-l-3/2)`.
4. Löwdin coefficients `beta_l = S_l^{-1/2}` via a Jacobi eigendecomposition.
5. Gaussian density width `eta = 1/(2 sigma^2)`; the convolution forms
   `aOa = -alpha / (1 + (1/eta) alpha)` and
   `bOa = (1/eta)^{3/2} beta (1 + (1/eta) alpha)^{-(l+3/2)}`.
6. Real spherical harmonics through degree five use the DScribe polynomial
   normalisation (`1.5708` for `l=0`, `2.7207` for `l=1`, and the published
   `preCoef` constants for `l=2..5`).
7. Per centre, `p_l[k,k'] = pi^4 sqrt(8/(2l+1)) sum_m c_klm c_k'lm` (the
   `pi^3` factor applies for `l > 1`).

Feature order is `l`, `k`, `k'` with `k <= k'`. `average = true` averages the
power spectrum over centres.

### Validation and deviation

`crates/xtalforge-ml/data/soap_dscribe.json` records a DScribe 2.x run
(`r_cut = 5.0`, `n_max = 3`, `l_max = 2`, `sigma = 0.5`, `periodic = true`,
`average = "off"`, `compression = "off"`) on two carbon atoms in a 5 A cubic
cell. The Rust test checks the layout and feature count against that fixture.
**The absolute normalization differs from DScribe in this build**: the
recorded DScribe rows are not reproduced numerically. The implementation is
therefore a GTO power spectrum with the published formula, not a bit-identical
DScribe port. The deviation is deliberate and documented here; the fast
approximation in `xtalforge_analysis::soap_like_descriptor` remains the
default for large screens.

## ACSF

G2 and G4 follow Behler and Parrinello with the cosine cutoff

```text
fc(r) = 0.5 (cos(pi r / r_cut) + 1) for r < r_cut, else 0
G2_i  = sum_j exp(-eta (r_ij - r_s)^2) fc(r_ij)
G4_i  = 2^(1-zeta) sum_{j,k != i} (1 + lambda cos theta_jik)^zeta
        exp(-eta (r_ij^2 + r_ik^2 + r_jk^2)) fc(r_ij) fc(r_ik) fc(r_jk)
```

Every element pair gets a G2 block and every unordered element pair a G4
block in sorted element order. The descriptor is invariant to rotation,
translation and atom order. The tests check shape, finiteness and averaging;
DScribe ACSF comparison is not bit-identical for the same reason as SOAP.

## MBTR

`k=1` broadens pair distances with weight `1/r`; `k=2` broadens the angle at
every neighbour triple with weight `1/(r_ij r_jk)`. Both use a Gaussian of
width `sigma` on a fixed grid. The vector concatenates the element-pair
distance blocks and the element-triple angle blocks; `angular = false` drops
the angle block.

## MOF features

RACs use the bond graph from Phase 1. For property `p` and path length `k`,

```text
RAC(p,k) = sum_{i<j, d_ij=k} (p_i - <p>)(p_j - <p>) / N_pairs
```

with `p` in atomic number, Pauling electronegativity, covalent radius and
mass, and `k` in one to three. The pore block is the Phase 10 distance-grid
descriptor vector. The concatenated vector is the MOF screening representation.

## Splits

Random splits shuffle record indices. Composition splits group by the sorted
element counts. Cluster splits use the Phase 9 matcher with the supplied
tolerances. Prototype splits use the space group plus the sorted
`(multiplicity, Wyckoff letter)` signature. Groups are shuffled with the seed
and assigned whole; largest-remainder rounding makes the requested fractions
sum to the dataset size.

## Inference

`PropertyModel` is a linear predictor with optional mean/scale normalization.
`Ensemble` reports mean, standard deviation and member values. With the
`onnx` feature, `OnnxModel` runs a single-input graph with tract and
`OnnxCalculator` differentiates the energy in Cartesian coordinates and in
the six symmetric strain components. Message-passing potentials stay in the
Python bridge.

## Active learning

`select_batch` supports uncertainty, farthest-point and diversity. The loop
runs predict, select, compute and append, renaming appended records so the
identifier stays unique. `active_learning_workflow` emits the same pipeline as
a Phase 16 workflow.

## Tests

| Check | Test |
|---|---|
| SOAP layout and feature count vs the DScribe fixture | `tests/soap.rs` |
| Graph symmetry, images and kNN truncation | `tests/representations.rs` |
| Split determinism and group integrity | `tests/representations.rs` |
| extxyz, JSON and OPTIMADE imports | `tests/representations.rs` |
| ACSF/MBTR shape and MOF features | `tests/descriptors.rs` |
| Linear, ensemble, selection and active loop | `tests/inference.rs` |
| Parquet round trip and ONNX error paths | `tests/optional_features.rs` |
