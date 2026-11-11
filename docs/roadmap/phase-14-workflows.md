# Phase 14 — Workflows, provenance and high throughput

**Goal:** run thousands of structures reproducibly, record how every result
was produced and let people (and agents) describe a workflow once and run it
on any surface.

**Release gate:** `v1.7.0`.

---

## Task 14.1 — Declarative workflows

- [ ] 14.1.1 Workflow file
  Description: a versioned YAML/JSON workflow (`workflow_version`) of steps
  (read, transform, analyse, compute, write) with parameters from the preset
  system (Task 3.8); validated against a JSON Schema published in the docs.
- [ ] 14.1.2 Runner
  Description: `xtalforge run wf.yaml --inputs "*.cif" --jobs 8`, Python
  `Workflow.run`, WASM `runWorkflow`; deterministic output order, resume after
  failure, per-item error isolation.
- [ ] 14.1.3 Agent integration
  Description: the agent can emit and run workflows; MCP tools `workflow_*`.

## Task 14.2 — Provenance

- [ ] 14.2.1 Every result carries a provenance record (xtalforge version,
  options, input hash, data bundle versions, seed).
- [ ] 14.2.2 Provenance embedded in CIF comments, extxyz info and project
  files; `xtalforge provenance file.cif` prints it.

## Task 14.3 — Datasets

- [ ] 14.3.1 Dataset container: many structures plus properties in one file
  (Parquet/Arrow via a pure-Rust crate, or JSON Lines), with a Python
  `pandas`/`polars` view.
- [ ] 14.3.2 Deduplication of structures across a dataset (structure matcher
  with tolerances, compatible with pymatgen `StructureMatcher` results).
- [ ] 14.3.3 Descriptors for ML (SOAP-like, RDF, Coulomb matrix, composition
  features) as arrays.

## Task 14.4 — Surfaces and docs

- [ ] 14.4.1 High-throughput tutorial on a CoRE MOF subset; parity matrix.
