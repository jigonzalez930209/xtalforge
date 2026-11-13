# Phase 16 — Workflows, provenance and high throughput

**Goal:** run thousands of structures reproducibly, record how every result
was produced and let people (and agents) describe a workflow once and run it
on any surface.

**Release gate:** `v1.9.0`.

---

## Task 16.1 — Declarative workflows

- [x] 16.1.1 Workflow file
  Description: a versioned YAML/JSON workflow (`workflow_version`) of steps
  (read, transform, analyse, compute, write) with parameters from the preset
  system (Task 3.8); validated against a JSON Schema published in the docs.
  Note: the shipped parser is JSON; `.yaml`/`.yml` files with non-JSON content
  return a typed parse error because no YAML parser is a dependency in 1.9.
- [x] 16.1.2 Runner
  Description: `xtalforge run wf.yaml --inputs "*.cif" --jobs 8`, Python
  `Workflow.run`, WASM `runWorkflow`; deterministic output order, resume after
  failure, per-item error isolation.
  Note: CLI, agent and Python are complete; WASM can validate metadata but has
  no filesystem for inputs/checkpoints, recorded as a parity exception.
- [x] 16.1.3 Agent integration
  Description: the agent can emit and run workflows; MCP tools `workflow_*`.

## Task 16.2 — Provenance

- [x] 16.2.1 Every result carries a provenance record (xtalforge version,
  options, input hash, data bundle versions, seed).
- [x] 16.2.2 Provenance embedded in CIF comments, extxyz info and project
  files; `xtalforge provenance file.cif` prints it.

## Task 16.3 — Datasets

- [x] 16.3.1 Dataset container: many structures plus properties in one file
  (Parquet/Arrow via a pure-Rust crate, or JSON Lines), with a Python
  `pandas`/`polars` view.
  Note: the shipped container is JSON Lines; Parquet/Arrow remains a future
  option and the Python view uses pandas/polars when installed.
- [x] 16.3.2 Deduplication of structures across a dataset, at scale, using
  the structure matcher of Task 9.1.1 with a fingerprint pre-filter (9.3.4).
- [x] 16.3.3 Descriptors for ML (SOAP-like, RDF, Coulomb matrix, composition
  features) as arrays.

## Task 16.4 — Surfaces and docs

- [x] 16.4.1 High-throughput tutorial on a CoRE MOF subset; parity matrix.
