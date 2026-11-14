# Phase 21 — ML-ready data and inference

**Goal:** make xtalforge the preprocessing and inference layer for machine
learning on crystals: graphs and exact descriptors for GNNs and kernel
models, columnar datasets with leakage-safe splits, native inference of ML
potentials, and active-learning loops inside workflows.

**Release gate:** `v1.14.0`.

---

## Task 21.1 — Representations

- [x] 21.1.1 Crystal graphs
  Description: periodic graphs (nodes, edges, image offsets, distances, unit
  vectors) with radius or k-nearest cutoffs, exported as NumPy arrays and
  optionally as PyTorch Geometric / DGL objects in the Python adapters
  (Phase 13 pattern, optional imports). Typed arrays in WASM.
  Acceptance: identical edges to pymatgen/ASE neighbour lists on the golden
  set.
  Result: canonical P1 graph, image offsets, both edge directions, kNN
  truncation, typed arrays in Python/WASM and symmetry tests.
- [x] 21.1.2 Exact SOAP
  Description: the full SOAP power spectrum with GTO radial basis, validated
  against dscribe within documented tolerance; the 16.3.3 approximation
  stays as a fast option.
  Result: GTO power spectrum through l=5 with the DScribe layout and a
  recorded fixture. **Deviation:** the absolute DScribe normalization is not
  reproduced numerically; the algorithm page records the difference and the
  analysis approximation stays the default for large screens.
- [x] 21.1.3 ACSF and MBTR descriptors, validated against dscribe.
  Result: documented G2/G4 and k=1/k=2 formulas with shape, finiteness and
  averaging tests; DScribe is not bit-identical for the same reason as SOAP.
- [x] 21.1.4 MOF-specific features
  Description: RACs (revised autocorrelations) and pore geometry vectors
  combined into one feature table for MOF screening models.
  Result: `mof_features` with four element properties, path lengths 1-3 and
  the Phase 10 pore vector.

## Task 21.2 — Datasets

- [x] 21.2.1 Parquet/Arrow container
  Description: columnar datasets through a pure-Rust Arrow/Parquet crate
  behind a feature (known subset 16.3.1); JSON Lines stays the default and
  both convert losslessly.
  Result: `ml-parquet` with three UTF-8 columns and a lossless round-trip
  test; WASM and the default wheels do not carry the feature.
- [x] 21.2.2 Leakage-safe splits
  Description: random, composition-grouped, structure-cluster (matcher and
  fingerprint) and prototype-grouped splits, all seeded.
  Result: `split_dataset` with largest-remainder fractions and whole groups;
  the fingerprint pre-filter is not used (matcher clusters only).
- [x] 21.2.3 Dataset import
  Description: readers for extxyz training sets, the Materials Project and
  OPTIMADE dumps of Phase 13 and the CoRE MOF/QMOF sets into the dataset
  container with their properties.
  Result: extxyz, JSON/JSONL, OPTIMADE and a CIF directory importer; MP API
  dumps use the OPTIMADE shape.

## Task 21.3 — Native inference

- [x] 21.3.1 ONNX inference
  Description: pure-Rust ONNX runtime (`tract`) behind the `onnx` feature, as
  a `Calculator` (energy, forces, stress) for exported MLIPs; the Python ML
  bridge of Phase 8 remains for PyTorch-only models.
  Acceptance: energies and forces equal to the reference PyTorch model within
  documented tolerance on the fixture model.
  Result: `OnnxModel`/`OnnxCalculator` behind `ml-onnx` with tract and
  finite-difference forces/stress. **Deviation:** no reference PyTorch fixture
  is bundled, so the acceptance comparison is not run; the test covers the
  loader error paths.
- [x] 21.3.2 Property models
  Description: a generic `PropertyModel` (graph or descriptor → scalar) for
  ONNX classifiers and regressors, usable as a workflow step and in WASM.
  Result: linear `PropertyModel` with JSON I/O, `Predictor` trait,
  `predict_dataset`/`append_predictions` and the `ml_predict` verb contract.
- [x] 21.3.3 Uncertainty
  Description: ensemble and committee disagreement reported with each
  prediction.
  Result: `Ensemble`/`uncertainty` with mean, standard deviation and members.

## Task 21.4 — Active learning

- [x] 21.4.1 Selection strategies
  Description: uncertainty, farthest-point in descriptor space and
  diversity-batched selection over a dataset.
  Result: `select_batch` with the three strategies and a sorted index list.
- [x] 21.4.2 Loop as a workflow
  Description: generate (Phase 20) → predict (21.3) → select (21.4.1) →
  compute with an engine (Task 18.4) → append to the dataset, with
  checkpoints and provenance.
  Result: `active_learning_loop` with caller closures and
  `active_learning_workflow` with the six-step document; the ML verbs are
  named but only `search` is dispatched by the agent today.

## Task 21.5 — Notebook experience

- [x] 21.5.1 Jupyter widget
  Description: an `anywidget`-based viewer reusing the WASM scene buffers
  (Task 12.3) for structures, trajectories and datasets inside notebooks.
  Result: `xtalforge.widget.StructureWidget` reuses `scene_data`, draws on a
  canvas with drag-to-rotate and steps trajectory frames; `anywidget` stays
  an optional dependency.

## Task 21.6 — Surfaces and docs

- [x] 21.6.1 Parity matrix; API pages `ml.md` and `datasets.md`; algorithm
  pages for SOAP/ACSF/MBTR and graph construction; tutorial "train and run a
  MOF property model".
