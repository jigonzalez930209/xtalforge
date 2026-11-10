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

## Current status (2026-11-10)

- Phases 0-8 are implemented, tested and pushed; the Rust API, CLI, Python and
  WASM surfaces build and pass their checks. `v0.0.1`, `v0.0.2` and `v0.0.3`
  are published on crates.io, PyPI (`xtalforge`), npm (`xtalforge-wasm`), the
  GitHub release and GitHub Pages; `1.0.0-rc.1` was published on 2026-11-10
  (13 crates, 7 wheels plus sdist, npm and the hosted CLI assets).
- Phase 8 is complete except the 1.0.0 stabilization: `xtalforge-agent`
  (verbs, session, diagnostics, MCP, skill), plugins, the Python ML bridge,
  distribution, benchmarks, robustness, tutorials and the API freeze.
  `v1.0.0-rc.1` is published and the `v1.0.0` tag follows the fixes-only
  period.
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
- 4.4: held cell lengths/volume (`CellFreedom`) and an external pressure term
  are not implemented; `OptimizeOptions.cell` relaxes the point-group-averaged
  strain only.
- 5.1 and 7.2: the calculated PXRD uses one Gaussian per reflection and the
  `Instrument` has no wavelength or Caglioti `U, V, W`; Rietveld ignores
  `RietveldOptions::kind` (area-normalised Gaussian, K alpha doublet by ratio).
- 5.x/7.x: the Pawley and Le Bail pseudo-Voigt mixing `eta` is a fixed option
  (0.5 by default), not a refined parameter, so a Gaussian pattern fitted with
  the default profile leaves a large Rwp (adding a refinable `eta` would change
  the frozen `PawleyOptions` struct).
- 3.10 and 4.4: `TransformCell`, `MakeSlab`, `AssignWyckoff`, `Invert`,
  `DescendToSubgroup`, `MergeDuplicates`, `InsertMolecules`,
  `ApplyOptimizedGeometry`, `AddHydrogens`, `Prepare` and `SetAtomTypes` are
  not command types, and the Vinet equation of state is missing; the
  functions behind most of them exist and apply through `ReplaceStructure`.
- 5.1-5.4: the PXRD has no anisotropic Debye-Waller factors, no
  pseudo-Voigt/Caglioti width, no monochromator option and a fixed radiation
  set; the porosity grid has no channel percolation, probe-occupiable volume or
  per-mass area; meshes have no normals or OBJ output.
- 5.6-5.10: no structure-to-net contraction, vertex symbols, search by space
  group or transitivity, Blatov interpenetration classes, or HTML/JSON report
  renderers; the module registry has no CLI command.
- 6.1-6.4: SMILES parsing and embedding, sketches, the fragment library,
  substitution and guest placement are Rust-only (no CLI, Python or WASM
  surface), see the open gaps in `docs/api/parity.md`.
- 7.x: no whole-range peak refinement, TCH/Caglioti/FCJ profiles, persistent
  powder parameter sets, Pareto sweep or knee, `run` module steps, named
  instrument presets, indexing engine agreement, cancellation or powder demo.
- 4.4-4.5: no cancellation or progress callback in the optimizers, no
  trajectory output, scan coordinates limited to atom-index distance, angle and
  torsion, and scans that neither write points to disk nor resume.
- 1.8: `Slab.cut` is always zero; severed bonds are not carried through the
  cut.
- 1.11: explicit bonds map to P1 indices and are not propagated to every
  symmetry image; pair-specific bond ranges and a metal-metal switch are not
  exposed.
- 6.x: the MOF builder uses a curated topology/block subset and a greedy
  orientation instead of PORMAKE's global RMSD fit.
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
