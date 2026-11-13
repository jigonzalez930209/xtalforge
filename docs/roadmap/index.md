# Roadmap

`xtalforge` reproduces everything the headless core of Crystal-Builder does
(the `xtal/` package; the `xtalapp/` GUI is out of scope) as a pure Rust
library, with four surfaces developed at the same pace: **Rust API, CLI, Python,
WebAssembly**.

Each phase has its own file. Inside a phase, work is split into numbered
**tasks** (`N.M`) and each task into **subtasks** (`N.M.K`) with a checkbox.
A subtask is ticked only when implementation, tests, all four surfaces (when
user-visible) and documentation are done. The rules that govern all of this
are in `AGENTS.md` at the repository root.

## Phases

| Phase | File | Goal | Release gate |
|---|---|---|---|
| 0 | [Foundations](./phase-00-foundations.md) | Workspace, rules, CI, docs site, moyo audit, sci-form reuse plan, dependency replacement map, binding architecture, golden-data oracle | none ([remaining work](./phase-00-status.md)) |
| 1 | [Core](./phase-01-core.md) | Elements, lattice, structure, P1, supercell, slab, neighbors, bonding, rings, groups, measures, limits, CIF/POSCAR/XYZ | with Task 2.1: **v0.0.1** |
| 2 | [Symmetry](./phase-02-symmetry.md) | moyo integration, space groups, Wyckoff, standard cells, reduction, subgroups, point groups, Brillouin zone and k-path | v0.0.x |
| 3 | [Formats, commands, view data](./phase-03-formats-commands.md) | All file formats, project files, undo/redo commands, selection, tracking, workspace, view settings and scene data, STL | v0.0.x |
| 4 | [Force fields](./phase-04-force-fields.md) | Calculator API, UFF/UFF4MOF, QEq/EQeq, Ewald, optimizers, scans, vibrations, phonons, thermodynamics, Prepare for simulation, external engines | **v0.1.0** |
| 5 | [Analysis](./phase-05-analysis.md) | PXRD, grids, porosity and accessibility, topology and RCSR, merged nets, interpenetration, overlays, module registry and reports | v0.y |
| 6 | [Builders](./phase-06-builders.md) | Molecules, substituents, fill, carbon, graphene, nanotubes, fullerenes, polymers, MOF builder | v0.y |
| 7 | [Powder refinement](./phase-07-powder.md) | Pattern data, peaks, indexing, Pawley, Rietveld, energy-weighted refinement, automatic run | v0.y |
| 8 | [Agent and distribution](./phase-08-distribution.md) | Session verbs, inspection, MCP server, ML potentials, distribution, performance, API freeze | **v1.0.0** |
| 9 | [Core features](./phase-09-core-features.md) | Structure matcher, strain/rattle/site edits, Voronoi and coordination, BVS/oxidation states, RDF, unique Miller indices and terminations, Wulff, adsorption sites, defects, orderings/SQS, interfaces, elastic tensor, EOS, convex hull | v1.2.0 |
| 10 | [MOF platform](./phase-10-mof.md) | MOF checker and cleaning, SBU/linker deconstruction, topology from structure, MOFid/MOFkey, open metal sites, native PLD/LCD/ASA, Widom/Henry/GCMC, PSM, MOF defects, CoRE MOF/QMOF screening | v1.3.0 |
| 11 | [Usability audit and 1.x polish](./phase-11-usability.md) | Fix the 2026-11-11 usability findings, close the Rust-only parity gaps, deprecation mechanism, recipes | v1.4.0 |
| 12 | [Native bindings](./phase-12-native-api.md) | Object API in Python and JS, NumPy, typed TS results, typed arrays, workers | v1.5.0 |
| 13 | [Ecosystem interoperability](./phase-13-interop.md) | ASE, pymatgen, phonopy, RDKit adapters; QE/CASTEP/VASP outputs/mmCIF/SDF; OPTIMADE, COD, MP | v1.6.0 |
| 14 | [Scientific completeness](./phase-14-completeness.md) | Close every known subset; instrument model, MD, full topology | v1.7.0 |
| 15 | [Performance and scale](./phase-15-performance.md) | Comparative benchmarks, million-atom structures, streaming IO, split WASM | v1.8.0 |
| 16 | [Workflows and provenance](./phase-16-workflows.md) | Declarative workflows, provenance, datasets, deduplication, ML descriptors | v1.9.0 |
| 17 | [2.0 redesign and freeze](./phase-17-v2.md) | Remove deprecations, unified options and errors, schema v2, migration guide | **v2.0.0** |

## Current status (2026-11-13)

- Phases 0-15 are complete; the Rust API, CLI, Python and WASM surfaces build
  and pass their checks. `v0.0.1`, `v0.0.2`, `v0.0.3`, `v1.0.0-rc.1`, the
  stable `v1.0.0`, `v1.1.0`, `v1.2.0`, `v1.3.0`, `v1.4.0`, `v1.5.0`, `v1.5.1`,
  `v1.6.0`, `v1.7.0` and `v1.8.0` are published on crates.io, PyPI
  (`xtalforge`), npm (`xtalforge-wasm`), the GitHub release and GitHub Pages.
  `v1.0.0` (2026-11-10) ships 13 crates, 7 wheels plus an sdist, the npm
  package and the hosted CLI assets; `v1.1.0` (2026-11-11) adds the 3D
  playground and the end-to-end audit battery; `v1.2.0` (2026-11-11) releases
  Phase 9, the core features of a standard crystal library; `v1.3.0`
  (2026-11-12) releases the MOF platform; `v1.4.0` closes the usability
  audit; `v1.5.0` (2026-11-13) releases the object API, NumPy and typed
  arrays, typed result classes, batch helpers and the typed JavaScript
  interface with the same 13-crate, 7-wheel plus sdist and npm payload;
  `v1.5.1` (2026-11-13) is the icon and package-metadata patch; `v1.6.0`
  (2026-11-13) releases Phase 13, the ecosystem adapters, formats and remote
  data; `v1.7.0` (2026-11-13) releases Phase 14, the diffraction, force-field,
  structure-operation, mesh and report work; `v1.8.0` releases Phase 15, the
  benchmark contract, scale paths, split WASM and lazy data.
- Phase 8 is complete: `xtalforge-agent` (verbs, session, diagnostics, MCP,
  skill), plugins, the Python ML bridge, distribution, benchmarks, robustness,
  tutorials and the API freeze. Phase 9 is complete: comparison, edits, local
  environments, surfaces, defects and disorder, interfaces, elasticity, EOS
  fits and convex hulls on all four surfaces, with five algorithm pages.
  Phase 10 is complete: MOF checking and cleaning, deconstruction, MOFid,
  descriptors, native PLD/LCD/ASA, host-guest force fields, Widom, GCMC,
  RASPA export, PSM, defects, variants, databases and screening. Phase 11 is
  complete: the usability audit fixes, naming and deprecation docs, recipes
  and the newly exposed molecule/sketch/fragment/guest operations. Phase 12
  is shipped: the Python `Structure` object API, typed result classes,
  no-`Any` stubs, NumPy input/output, batch helpers, the WASM `Structure`
  handle, typed scene/PXRD/mesh buffers, the typed TS interface and the
  worker/Node helpers. Phase 13 is complete: the ASE, pymatgen,
  phonopy/spglib and RDKit adapters, the extended format registry, OPTIMADE,
  COD and Materials Project fetchers, Python entry-point plugins, the
  migration guides and the interoperability API page. Phase 14 is complete:
  profiles, anisotropic ADP, single-crystal factors, cell freedom, MD, slab
  bond carrying, undoable commands, OBJ/glTF and HTML/JSON reports. Phase 15
  is complete: comparative benchmarks and the regression gate, SoA/lazy/
  streaming/parallel paths, split WASM packages, lazy data and cold-start
  measurement. Tasks 9.7.2, 10.5.7, 10.6.4, 10.7.3-10.7.5, 11.4.5, 12.1.4,
  12.3.5, 12.4.1, 14.1.5 and 14.2.6 remain open. The stable API follows the
  SemVer policy and the known subsets below are the only documented reference
  gaps.
- The release workflow is the only GitHub Actions workflow; the CI, rules and
  manual docs workflows were removed by maintainer decision (2026-11-08) to
  avoid runner cost. The release runs on a version tag: three CLI packages,
  PyPI/npm/crates.io from `ubuntu-24.04`, then the docs deployment. The
  `github-pages` environment allows `main` and the `v*` tag pattern.

### Known subsets

The checkboxes record that the documented local subset is implemented, tested
and exposed on the surfaces. The following reference behaviours are subsets or
simplified ports rather than complete reproductions; they are the first
candidates for Phase 8 or a follow-up and are not silent approximations:

- 4.9: the DFTB+ parser is validated against a recorded DFTB+ output
  (`crates/xtalforge-ff/tests/fixtures/dftb/si_scc_detailed.out`), and the xTB
  GFN0/GFN1/GFN2 methods are available in-process through the `sci-form`
  feature. The external launchers still need a real binary, so a complete
  DFTB+/xTB run is not executed in CI; ORCA remains input-only by design.
- 9.6: the Vinet equation of state is not in `EosForm`; the historical
  Birch-Murnaghan and Murnaghan forms are unchanged (Task 9.6.2).
- 5.3-5.4: the porosity grid still has no channel percolation,
  probe-occupiable volume or per-mass area (MOF Tasks 10.3.1/10.3.2); meshes,
  normals and OBJ/glTF are complete.
- 6.1-6.4: SMILES parsing and embedding, sketches, the fragment library,
  substitution and guest placement are Rust-only (no CLI, Python or WASM
  surface), see the open gaps in `docs/api/parity.md`.
- 7.x: indexing-engine agreement and the powder demo remain open; whole-range
  refinement, TCH/Caglioti/FCJ profiles, persistent parameter sets, Pareto
  sweep/knee, named presets and cancellation are complete.
- 15.2.2: streaming readers cover XYZ/LAMMPS/CIF; streaming writers and the
  extxyz/force columns of the streams remain open.
- 15.3.1: the local split WASM packages and budgets are complete; publishing
  the new npm tier names is a release follow-up.
- 15.3.2: the bundles are gzip/lazy in the crates; automatic download-on-first-
  use in WASM remains open.
- 15.2.1/15.2.4: the SoA columns are not fused into the linked-cell search and
  lazy expansion only serves operations that do not need the materialised P1.
- 6.x: the MOF builder uses a curated topology/block subset and a greedy
  orientation instead of PORMAKE's global RMSD fit.
- 9.1: the matcher works on Niggli-reduced input cells and accepts supercell
  relationships; it does not call the standardized primitive-cell path (which
  is a separate API and is being audited), so a conventional/primitive pair is
  matched by expanding the smaller cell.
- 9.3: the Voronoi cutoff must contain every face of the cell (twice the
  longest nearest-neighbour distance is safe); shape measures cover the
  common polyhedra up to 12 vertices; the CrystalNN strategy is a documented
  solid-angle/distance approximation, not pymatgen's full bond-order model.
- 9.4: terminations use midpoint cuts between atomic layers, so cuts inside a
  layer are not generated; adsorption sites are top/bridge/hollow on the
  uppermost layer only.
- 9.5: point defects use Voronoi vertices as interstitial candidates; the
  ordering search caps the decorations (`max_orderings`) and ranks them with
  the in-process Ewald sum; SQS is binary and reports pair/triplet
  correlations, not a full cluster expansion.
- 10.x: the curated six-topology/eight-block builder subset remains until the
  licensed PORMAKE data drop (Task 10.5.7); database readers work on an
  extracted directory, not a zip, and on-demand fetch is Task 10.6.4.
- 10.2: deconstruction uses the all-node rule with metal-adjacent non-carbon
  atoms in the node; MOFid is a deterministic subset of the published syntax
  (unit graphs, topology, interpenetration) and MOFkey is its SHA-256.
- 10.3: PLD/channel dimension are voxelized union-find results from the
  distance grid; they converge with the spacing and are not bit-identical to
  Zeo++.
- 11.1: `--threads` and the Python `threads=`/`set_threads` pool change only
  execution; the parallel implementations keep fixed reductions. WASM has no
  threads without cross-origin isolation, recorded as an exception.
- 11.2-11.3: the facade `write` and the Python file helpers need a filesystem;
  WASM offers the text/JSON forms instead. SMILES embedding builds heavy-atom
  geometry and does not infer hydrogens from aromatic bond orders.
- 11.4: the recipe page is documentation; Task 11.4.5 executes every recipe
  on every surface, and the deprecation aliases are listed in
  `docs/api/naming.md`.
- 10.4: GCMC is single-component with rigid host/guest and an ideal-gas
  chemical potential; guest-guest interactions and Peng-Robinson fugacity are
  left to the RASPA export. Widom uses direct host-guest sums with a 15 A
  cutoff.
- 9.6: the elastic tensor uses the stress finite differences of the native UFF
  calculator (any `Calculator` with stress works), and the EOS/hull fits are
  point estimates with parameter errors, without a global optimizer audit.
- 8.x: every library error exits the CLI with code 3 and Python raises
  `XtalforgeError` (`InputError` for invalid input); the `code` string is the
  discriminator.

## Release milestones

Version rules are defined in `AGENTS.md` section 5. In short:

- **v0.0.1** — Phase 1 complete plus Task 2.1 (moyo integration). The
  "essential core": read a CIF, expand it to P1, detect its space group, build
  supercells and slabs, perceive bonds and rings, write CIF/POSCAR/XYZ, all
  through Rust, CLI, Python and WASM, with docs and demos.
- **v0.0.x** — one PATCH release per completed task group of Phases 2 and 3,
  and of Phase 4 until it is complete.
- **v0.1.0** — Phase 4 complete: force fields and optimizers usable.
- **v0.y.z** — Phases 5 to 7; MINOR bump at each completed phase.
- **v1.0.0** — Phase 8 complete: API frozen, parity matrix closed,
  documentation complete.
- **v1.2.0** — Phase 9 complete: the core features of a standard crystal
  library beyond Crystal-Builder parity (released 2026-11-11).
- **v1.3.0** — Phase 10 complete: the MOF platform (focus of the 2.0 road;
  released 2026-11-12).
- **v1.4.0 … v1.9.0** — Phases 11 to 16, one MINOR release per completed
  phase; every change is backward compatible and anything to be removed is
  deprecated first (Task 11.4.2). `v1.4.0` was released 2026-11-12 with the
  usability audit.
- **v2.0.0** — Phase 17 complete: deprecated surface removed, unified options
  and errors, schema version 2, after a `2.0.0-rc.N` stabilization period.

### Road to 2.0.0 (started 2026-11-11)

The goal is to make xtalforge the default crystal-structure library across
Rust, Python, the CLI and the browser. The plan starts from the core
gap check and usability audit of 2026-11-11 and orders the work by user
impact. **Focus: metal-organic frameworks** (decided 2026-11-11): first the
missing core operations (Phase 9), then the MOF platform (Phase 10), then
remove friction (Phases 11-12), then meet users where they are (Phase 13),
then close the science gaps (Phase 14), then scale (Phases 15-16), and only
then break the API once (Phase 17). No other capability is dropped: the MOF
focus sets priorities inside each phase (table at the end of Phase 10).
Phases 9-16 may overlap; Phase 17 starts when they are
complete or explicitly descoped with a `Note:`.

## Cross-cutting requirements (apply to every phase)

1. Every code file is at most 500 lines (Markdown exempt); data goes to
   `data/`.
2. Language is unrestricted; the former English-only rule is retired
   (2026-11-08).
3. Every user-visible subtask ships four surfaces and docs: Rust, CLI, Python,
   WASM, plus a VitePress API page, and an algorithm page or demo where it
   makes sense.
4. Every ported algorithm is validated against Crystal-Builder or sci-form
   golden data (Task 0.9) with documented tolerances.
5. `CHANGELOG.md`, `docs/api/parity.md` and `AGENTS.md` section 6 are updated
   as part of the work, not afterwards.
6. Each phase ends with a "surfaces and docs" task; it is not optional.
7. Warnings are errors; dead-code, unused and warning-suppression attributes
   are forbidden (AGENTS.md rule 10).
8. Commits are directory-scoped: one directory per commit,
   `feat(<scope>):`/`ref(<scope>):` subject, and the completed tasks listed
   from hardest to easiest in the body (AGENTS.md section 2.1).

## Parity map: Crystal-Builder module to phase

| Crystal-Builder module (`xtal/…`) | Phase |
|---|---|
| `core/` elements, pyykko, lattice, site, structure, molecule, p1, supercell, slab, cellcut, transforms, neighbors, bonding, rings, groups (functional groups), measure, properties, scalars, describe, limits; `io/text.py`, `io/cif_*`, `io/poscar.py`, `io/xyz.py` | 1 |
| `core/` spacegroup, symmetry, subgroups; point-group part of `core/molecule.py`; `analysis/kpath.py` | 2 |
| `io/` all other formats, `io/registry.py`, `io/export.py`, `io/atomic.py`, `io/project.py`, `io/trajectory.py`, `io/pmg.py`; `commands/`; `core/selection.py`, `core/tracking.py`, `core/atom_groups.py`; `workspace.py`, `params.py`, `plugins.py` (registries); `view/`; STL part of `modules/blender.py` | 3 |
| `ff/` (all), `core/prepare.py`, `commands/ff.py`, `commands/prepare.py`, `modules/{forcefield,scan,thermo,qha,dftb,orca,probe,process}.py`, `orca/` | 4 |
| `analysis/` (except kpath), `modules/{pxrd,zeopp,poregrid,report,registry,job,record,stub,net}.py`, `references.py` | 5 |
| `build/`, `carbon/`, `graphene/`, `nanotube/`, `fullerene/`, `polymer/`, `mof/`, `modules/{build,carbon,graphene,nanotube,fullerene,polymer,mof}.py`, `commands/connections.py` | 6 |
| `powder/`, `modules/powder.py` | 7 |
| `agent/`, `cli.py` (completed), ML calculators (`ff/mace`, `ff/orb`, `ff/mattersim`, `ff/ase_engine.py`) | 8 |
| `install.py`, `runtime.py`, Blender rendering in `modules/blender.py`, `xtalapp/` | out of scope (application packaging and GUI) |

## Known risks

- **moyo private APIs.** Niggli, Delaunay and Minkowski reductions exist in
  moyo 0.22 but are crate-private. Plan: propose an upstream PR to export them
  and, until it is released, port the algorithms (MIT OR Apache-2.0 allows it,
  with attribution). Tasks 0.3 and 2.6.
- **Python-only dependencies with no Rust equivalent.** Crystal-Builder uses
  gemmi (CIF, space groups), spglib, scipy (k-d trees, Voronoi, L-BFGS,
  sparse), ASE (band paths, trajectory reading, ML calculators), phonopy
  (phonons, thermodynamics), RDKit (SMILES, sketcher), rietx (indexing,
  Pawley, Rietveld). Each needs a Rust replacement or an explicit decision
  (Task 0.10).
- **Powder refinement.** The numerical core lives in rietx, an external
  Python project with its own kernels (`rietx-kernels`). Phase 7 starts by
  studying it; its license must be verified before porting.
- **Data licensing.** PORMAKE is MIT (verified). The RCSR database snapshots
  (`analysis/data/RCSR*`) and the Zeo++ radii need license checks before
  redistribution (Tasks 5.6.1, 5.2.2 and 6.11.1).
- **500-line limit versus ported code.** Many reference files exceed it
  (Python: `core/prepare.py` 1833, `modules/powder.py` 1893, `mof/build.py`
  1784; sci-form: `space_group_table.rs` 6750, `stereo.rs` 1121, `cell.rs`
  1048). Every port starts with a written split plan.
- **WASM limits.** No threads by default, no filesystem, no subprocesses. The
  parity matrix records each justified exception.
- **Scope size.** Phases 6 and 7 are larger than all the others combined.
  They are staged so each task group ships on its own.
