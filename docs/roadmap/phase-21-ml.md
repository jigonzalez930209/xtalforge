# Phase 21 — ML-ready data and inference

**Goal:** make xtalforge the preprocessing and inference layer for machine
learning on crystals: graphs and exact descriptors for GNNs and kernel
models, columnar datasets with leakage-safe splits, native inference of ML
potentials, and active-learning loops inside workflows.

**Release gate:** `v1.14.0`.

---

## Task 21.1 — Representations

- [ ] 21.1.1 Crystal graphs
  Description: periodic graphs (nodes, edges, image offsets, distances, unit
  vectors) with radius or k-nearest cutoffs, exported as NumPy arrays and
  optionally as PyTorch Geometric / DGL objects in the Python adapters
  (Phase 13 pattern, optional imports). Typed arrays in WASM.
  Acceptance: identical edges to pymatgen/ASE neighbour lists on the golden
  set.
- [ ] 21.1.2 Exact SOAP
  Description: the full SOAP power spectrum with GTO radial basis, validated
  against dscribe within documented tolerance; the 16.3.3 approximation
  stays as a fast option.
- [ ] 21.1.3 ACSF and MBTR descriptors, validated against dscribe.
- [ ] 21.1.4 MOF-specific features
  Description: RACs (revised autocorrelations) and pore geometry vectors
  combined into one feature table for MOF screening models.

## Task 21.2 — Datasets

- [ ] 21.2.1 Parquet/Arrow container
  Description: columnar datasets through a pure-Rust Arrow/Parquet crate
  behind a feature (known subset 16.3.1); JSON Lines stays the default and
  both convert losslessly.
- [ ] 21.2.2 Leakage-safe splits
  Description: random, composition-grouped, structure-cluster (matcher and
  fingerprint) and prototype-grouped splits, all seeded.
- [ ] 21.2.3 Dataset import
  Description: readers for extxyz training sets, the Materials Project and
  OPTIMADE dumps of Phase 13 and the CoRE MOF/QMOF sets into the dataset
  container with their properties.

## Task 21.3 — Native inference

- [ ] 21.3.1 ONNX inference
  Description: pure-Rust ONNX runtime (`tract`) behind the `onnx` feature, as
  a `Calculator` (energy, forces, stress) for exported MLIPs; the Python ML
  bridge of Phase 8 remains for PyTorch-only models.
  Acceptance: energies and forces equal to the reference PyTorch model within
  documented tolerance on the fixture model.
- [ ] 21.3.2 Property models
  Description: a generic `PropertyModel` (graph or descriptor → scalar) for
  ONNX classifiers and regressors, usable as a workflow step and in WASM.
- [ ] 21.3.3 Uncertainty
  Description: ensemble and committee disagreement reported with each
  prediction.

## Task 21.4 — Active learning

- [ ] 21.4.1 Selection strategies
  Description: uncertainty, farthest-point in descriptor space and
  diversity-batched selection over a dataset.
- [ ] 21.4.2 Loop as a workflow
  Description: generate (Phase 20) → predict (21.3) → select (21.4.1) →
  compute with an engine (Task 18.4) → append to the dataset, with
  checkpoints and provenance.

## Task 21.5 — Notebook experience

- [ ] 21.5.1 Jupyter widget
  Description: an `anywidget`-based viewer reusing the WASM scene buffers
  (Task 12.3) for structures, trajectories and datasets inside notebooks.

## Task 21.6 — Surfaces and docs

- [ ] 21.6.1 Parity matrix; API pages `ml.md` and `datasets.md`; algorithm
  pages for SOAP/ACSF/MBTR and graph construction; tutorial "train and run a
  MOF property model".
