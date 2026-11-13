# Phase 17 — Consolidation and trust

**Goal:** close every open subtask and every documented "known subset" that
blocks production use, so the 1.x line has no silent gaps before the scope
grows. Nothing new is invented here: each subtask finishes work that an
earlier phase started and recorded as open.

**Release gate:** `v1.10.0`.

---

## Task 17.1 — Close the open subtasks of Phases 9-14

- [ ] 17.1.1 Rust-only builders on every surface
  Description: SMILES parsing and embedding, sketches, the fragment library,
  substitution and guest placement (Phases 6 and 11) get CLI subcommands,
  Python functions and `Structure` methods, and WASM exports. Hydrogens are
  inferred from aromatic bond orders (known subset 11.2-11.3).
  Acceptance: the "Open gaps: Rust-only capabilities" table of
  `docs/api/parity.md` is empty; one test per surface per operation,
  including the error paths.
- [ ] 17.1.2 Porosity completeness
  Description: channel percolation, probe-occupiable volume and per-mass
  surface areas on the porosity grid (known subset 5.3-5.4, struck subtask
  14.4.1). Inputs: structure, probe radius, grid spacing; outputs: percolating
  directions, POV in A^3 and cm^3/g, ASA in m^2/g.
  Acceptance: agreement with Zeo++ on ten CoRE MOF structures within the
  tolerance documented in the algorithm page.
- [ ] 17.1.3 Complete PORMAKE and RCSR bundles
  Description: finishes Task 10.5.7 and the RCSR part of 14.4.3 behind the
  `mof-data` feature, loaded lazily (Task 15.3.2), with the license recorded
  in `THIRD_PARTY.md`.
- [ ] 17.1.4 On-demand database fetch
  Description: finishes Task 10.6.4 (CoRE MOF/QMOF) with a cache directory,
  checksums and an offline mode; WASM uses `fetch` and stores nothing on
  disk (parity exception).
- [ ] 17.1.5 Executable recipes and error snapshots
  Description: finishes Task 11.4.5: every recipe of `docs/guide/recipes.md`
  runs in the test suite on every surface and the error messages are
  snapshot-tested.
- [ ] 17.1.6 Python context managers, progress and optional WASM threads
  Description: finishes Tasks 12.1.4 and 12.3.5.
- [ ] 17.1.7 Remaining powder items
  Description: finishes Task 14.1.5 (indexing-engine agreement and the
  powder demo) and Task 14.2.6 (real DFTB+/xTB runs recorded as fixtures).
- [ ] 17.1.8 Playground and e2e references
  Description: finishes Tasks 9.7.2, 10.7.3, 10.7.4 and 10.7.5 (Phase 9
  playground views, MOF tab, Zeo++/MOFid/RASPA e2e references).

## Task 17.2 — Close the known subsets of Phases 15-16

- [ ] 17.2.1 YAML workflows
  Description: a pure-Rust YAML parser (audited, permissive license) so
  `.yaml` workflows accept the full syntax used in the docs; JSON stays the
  canonical form and both validate against the same JSON Schema.
- [ ] 17.2.2 Streaming writers and extxyz columns
  Description: finishes known subset 15.2.2: streaming writers for XYZ,
  extxyz (with forces, energies and arbitrary per-atom columns), LAMMPS dump
  and CIF; memory stays bounded on a one-million-atom trajectory.
- [ ] 17.2.3 Fused SoA neighbour search
  Description: fuse the SoA columns into the linked-cell search (known subset
  15.2.1) and let lazy expansion serve the neighbour list without
  materialising P1; the regression gate of Task 15.1 must show no slowdown.
- [ ] 17.2.4 WASM data on first use and npm tiers
  Description: automatic download-on-first-use of data bundles in WASM and
  publication of the split npm tier names (known subsets 15.3.1 and 15.3.2).
- [ ] 17.2.5 In-memory WASM workflows
  Description: `runWorkflow` in WASM over an in-memory virtual filesystem
  (inputs as strings/bytes, outputs as a map), which narrows the parity
  exception of 16.1.2 to real-disk access only.

## Task 17.3 — Scientific debts

- [ ] 17.3.1 Vinet equation of state in `EosForm` (known subset 9.6).
- [ ] 17.3.2 Terminations with intra-layer cuts and adsorption sites below
  the top layer (known subset 9.4), behind an option that keeps the current
  default.
- [ ] 17.3.3 Matcher on the standardized primitive cell once the primitive
  path audit is finished (known subset 9.1).
- [ ] 17.3.4 CrystalNN bond-order model closer to pymatgen, validated on the
  pymatgen test set (known subset 9.3).

## Task 17.4 — Local quality gate

- [ ] 17.4.1 `scripts/gate.sh`
  Description: one local command that runs format, clippy, rustdoc, tests,
  the rules and parity checks, the Python and WASM tests, the docs build and
  the benchmark regression gate, with a summary table. It replaces the
  removed CI workflows at no runner cost and is required before every tag
  (AGENTS.md section 5.4).
- [ ] 17.4.2 Golden-data refresh: regenerate the oracle data (Task 0.9) and
  record any drift in the algorithm pages.

## Task 17.5 — Surfaces and docs

- [ ] 17.5.1 Parity matrix with no open Rust-only gaps; the "Known subsets"
  list of `index.md` shrinks to what is deliberately out of scope; changelog,
  API pages and AGENTS.md section 6 updated.
