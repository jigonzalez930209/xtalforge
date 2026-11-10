# Phase 8 — Agent tools, ML potentials, distribution and 1.0

**Goal:** the automation layer an AI assistant drives (session verbs,
inspection, diagnostics, change tracking, MCP server, skill documents), the
plugin story, optional machine-learning potentials, polished distribution on
crates.io, PyPI and npm with prebuilt CLI binaries, performance and robustness
work, complete documentation, and the API freeze for **v1.0.0**.

**Reference:** `Crystal-Builder/xtal/agent/` (session 1834 lines, tools 495,
capabilities, answers, inspect, diagnostics, changes, context, render,
serve, proxy, discovery, skill/ with `SKILL.md`, `mcp.md`, prompts and
references), `xtal/plugins.py`, `xtal/cli.py` (1058 lines),
`xtal/ff/{ase_engine.py,mace,orb,mattersim}`.

Crate: `xtalforge-agent` (plus work across all crates).

---

## Task 8.1 — Session, verbs and answers

- [x] 8.1.1 Verb catalogue
  Port `agent/capabilities.py`: the list of verbs (`VERBS`) an assistant can
  call, each with name, parameters (Task 3.8), result type and examples,
  generated from the same definitions as the CLI, Python and WASM surfaces so
  there is no second list to drift. `capabilities` also reports what **this**
  installation can do (engines available, features compiled, data bundles
  present), asked at runtime rather than remembered.
- [x] 8.1.2 Session
  Port `agent/session.py` (1834 lines; split by verb group): a stateful
  session owning documents (Phase 3), running verbs by name with JSON
  arguments and returning typed answers; every mutation goes through commands
  with `Origin::Agent`.
- [x] 8.1.3 Inspection and diagnostics
  Port `agent/inspect.py` and `agent/diagnostics.py`: one call that says what
  a structure is and what is wrong with it (the Info panel, the site table and
  the warnings in one answer): overlapping atoms, missing hydrogens, unusual
  coordination, partial occupancy, charge imbalance, symmetry problems
  (Task 2.9.2), with severity, atoms involved and a suggested verb.
- [x] 8.1.4 Change tracking
  Port `agent/changes.py`: which steps of a document's history an assistant
  made, grouped in bursts, using the `Origin` on each command (Task 3.9.2), so
  a user can review or undo the assistant's work as a unit.
- [x] 8.1.5 Answers
  Port `agent/answers.py`: typed, printable answers that serialize to JSON in
  one call.
- [x] 8.1.6 Pictures for the agent
  Port the intent of `agent/render.py` (a picture of a structure for an agent
  to look at; Crystal-Builder renders with VTK). Implement a small pure-Rust
  orthographic renderer of the scene data (Task 3.14.4) to SVG and PNG
  (`tiny-skia` or similar pure-Rust rasterizer), good enough to see a folded
  linker or a misplaced cluster. Also used by the CLI `render` command.

## Task 8.2 — MCP server and skill

- [x] 8.2.1 MCP tools over stdio
  Port `agent/tools.py` and `agent/serve.py`: one MCP tool per verb, named as
  the verb and taking its parameters, served over stdio by `xtalforge mcp`
  (for clients that start the command themselves). Use a maintained pure-Rust
  MCP implementation or a minimal JSON-RPC layer; record the choice.
- [x] 8.2.2 Context for chat clients
  Port `agent/context.py`: what an MCP client is told (instructions and
  references) as against what it can call, for clients that cannot read the
  skill files from disk.
- [x] 8.2.3 Skill documents
  Port `agent/skill/` (`SKILL.md`, `mcp.md`, prompts `build-mof.md`,
  `porosity-report.md`, `prepare-deposited-cif.md`, references `api`,
  `calculations`, `carbon`, `diagnostics`, `mof`, `polymer`, `porosity`,
  `prepare`, `workspace`) rewritten for xtalforge names; shipped with the CLI
  (`xtalforge skill path`, `xtalforge skill install`) and published in the
  docs as an "Automation guide".
- [x] 8.2.4 Out of scope: window proxy
  `agent/proxy.py` and `agent/discovery.py` forward MCP calls to a running GUI
  window through an `mcp.json` discovery file. There is no window here; the
  decision is recorded in `docs/algorithms/agent.md` section 11 and the
  discovery format is documented in the automation guide for a future front
  end.
- [x] 8.2.5 Safety limits
  Every tool call runs under the standard `Limits` profile (Task 1.3) and the
  supercell, grid and render estimates are checked before allocation. Verbs
  are synchronous; cancellation is applied between calls and a client that
  must preempt a running verb terminates the server process, as recorded in
  `docs/algorithms/agent.md` section 11.

## Task 8.3 — Plugins and extension points

- [x] 8.3.1 Rust extension points
  Document and test registering external formats, calculators, modules,
  scalars and draw styles (the registries of Phases 3 to 5) from a separate
  crate, with an example plugin crate in `examples/`.
- [x] 8.3.2 Python-side plugins
  Python callables registered as formats, calculators (for example wrapping an
  ASE calculator, the replacement for `ff/ase_engine.py`) or modules, with the
  GIL handled correctly and errors mapped.

## Task 8.4 — Machine-learning potentials

- [x] 8.4.1 Backend implementation
  Implement the decision of Task 4.11.1 behind feature `ml`: MACE, ORB-v3 and
  MatterSim calculators, or the subset the chosen backend supports, with
  forces and stress.
- [x] 8.4.2 Model management
  Weights are never bundled: download on demand into a cache with checksum and
  license recorded per model; CLI `models {list,fetch,verify,remove}`.
- [x] 8.4.3 Accuracy and speed
  Compare energies and forces with the Python implementations on a test set
  recorded as golden data; publish timings.

## Task 8.5 — Distribution

- [x] 8.5.1 crates.io
  Publish all crates in dependency order with metadata, `include` allow-lists
  (fixtures, golden data and scripts never reach the tarball; data bundles
  only in the crates that embed them), docs.rs metadata with the right
  features, and crate size checks.
- [x] 8.5.2 PyPI
  abi3 wheels for manylinux and musllinux (x86_64, aarch64), macOS (x86_64,
  arm64) and Windows (x86_64), plus sdist; trusted publishing; the `sci-form`
  and `ml` features in separate wheel variants if they are too large for the
  default wheel.
- [x] 8.5.3 npm
  `xtalforge-wasm` with web, node and bundler targets, TypeScript types, size
  budget, provenance; MOF and RCSR bundles fetched on demand.
- [x] 8.5.4 CLI binaries
  Prebuilt `xtalforge` binaries per platform attached to GitHub releases
  (cargo-dist or equivalent), shell completions, man page; `cargo install
  xtalforge-cli` documented.

## Task 8.6 — Performance and robustness

- [x] 8.6.1 Benchmarks and regression tracking
  Process-level benchmarks for neighbor search, symmetry detection, UFF
  energy/forces, PXRD and the distance grid; `scripts/check-bench-budget.sh`
  is the local regression gate (the repository has no CI), and the Python
  timings are recorded on the ML page. Numbers are in
  `docs/algorithms/benchmarks.md`.
- [x] 8.6.2 Parallelism and determinism review
  Verify identical results across thread counts and platforms for every
  parallel path; document the threaded WASM build.
- [x] 8.6.3 Fuzzing and property tests
  Fuzz every parser (CIF, POSCAR, XYZ, SHELX, LAMMPS, CGD, cube, project zip);
  property tests for symmetry, supercells, reduction and optimizers.
- [x] 8.6.4 Memory and limits audit
  Every allocating path respects `Limits` and fails with a `Verdict`; test with
  structures at and beyond each limit.

## Task 8.7 — Documentation completeness

- [x] 8.7.1 API coverage
  `missing_docs = "deny"`; every public item in all four surfaces documented;
  the docs build fails on a missing page.
- [x] 8.7.2 Algorithm coverage
  Every ported algorithm has its page with references, validation data and
  differences from Crystal-Builder.
- [x] 8.7.3 Tutorials
  End to end: CIF to supercell, bonds and symmetry; Prepare a deposited CIF;
  UFF optimization and a relaxed scan; porosity report; PXRD comparison with a
  measured pattern; powder refinement; building a MOF; building a carbon.
- [x] 8.7.4 Migration guide from Crystal-Builder
  Mapping of every `xtal` module, CLI command (`info`, `render`,
  `capabilities`, `skill`, `mcp`, `symmetry`, `convert`, `bonds`, `types`,
  `formats`, `modules`, `engines`, `run`, builders) and agent verb to its
  xtalforge equivalent, with known behavioural differences.

## Task 8.8 — API freeze and 1.0

- [x] 8.8.1 Public API review
  Audit names, types, error codes, CLI flags and JSON shapes across the four
  surfaces; fix inconsistencies with deprecations; write the SemVer policy
  page (from `AGENTS.md` section 5).
- [x] 8.8.2 Parity matrix closure
  No unexplained gaps in `docs/api/parity.md`; every exception justified.
- [ ] 8.8.3 Release candidate cycle
  `v1.0.0-rc.1` with a stabilization period and fixes only, then `v1.0.0`;
  update `AGENTS.md` section 6 and the changelog.

  Note: `v1.0.0-rc.1` was published on 2026-11-10 (13 crates on crates.io,
  7 platform wheels plus an sdist on PyPI, npm, the GitHub release and Pages).
  The remaining work is the fixes-only stabilization period and the `v1.0.0`
  tag; this checkbox stays open until `1.0.0` is published.
