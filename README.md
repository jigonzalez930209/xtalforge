# xtalforge

A pure Rust library for building, manipulating, analysing and exporting crystal
structures, with a CLI, Python bindings and WebAssembly bindings developed at
parity with the Rust API.

> Status: Phases 0-7 are complete locally at the documented-subset level
> (`0.0.1`, release prepared, nothing published yet). The core, symmetry, formats and commands,
> force fields, analysis, builders and powder refinement work through the Rust
> API, CLI, Python and WASM packages; the known reference gaps are listed in
> [`docs/roadmap/index.md`](docs/roadmap/index.md). Phase 8 (agent and
> distribution) is next. GitHub Actions remain disabled until the maintainer
> enables them, so the release dry run and the name reservation are still
> pending. The plan is in [`docs/roadmap`](docs/roadmap/index.md); the project
> rules are in [`AGENTS.md`](AGENTS.md).

## What works today

- Crystal structures: elements, lattices, sites, P1 expansion, supercells,
  slabs, transforms, periodic neighbours, bond graphs, rings and measures.
- Symmetry: space-group detection, Wyckoff positions, standard and primitive
  cells, lattice reduction, subgroups, layer and magnetic groups, molecular
  point groups, Brillouin zones and band paths.
- Formats and editing: CIF, POSCAR, XYZ, pymatgen JSON, CSSR, GEN, LAMMPS,
  SHELX and PDB; atomic writes, project archives, undo/redo documents,
  selections, workspace runs and scene data.
- Force fields and thermodynamics: native UFF/UFF4MOF, QEq/EQeq/EEQ charges,
  Ewald, FIRE/L-BFGS, scans, Gamma-point vibrations, IR, thermochemistry and
  quasi-harmonic free energies.
- Analysis and builders: PXRD, porosity, isosurfaces, nets, RCSR matching,
  graphene, nanotubes, fullerenes, SMILES embedding, polymers, MOF assembly
  and powder refinement.
- Four surfaces at parity, with documented exceptions for the filesystem,
  subprocesses and threads in WASM. See `docs/api/parity.md`.

## Install

Nothing is published yet. After the first release the commands will be:

```bash
cargo add xtalforge@0.0.3
pip install xtalforge==0.0.3
npm install xtalforge-wasm@0.0.3
```

## Goals

- Everything the headless core of [Crystal-Builder](https://github.com/JulesOpp/Crystal-Builder)
  does: symmetry-aware structure editing, file formats, force fields and
  optimizers, PXRD, porosity, topology, builders (molecules, carbons,
  polymers, MOFs) and powder refinement.
- Symmetry through [moyo](https://github.com/spglib/moyo).
- Proven code reused from [sci-form](https://github.com/jigonzalez930209/sci-form).
- Four surfaces at parity: Rust (`xtalforge` crate), CLI (`xtalforge`
  binary), Python (`xtalforge` on PyPI), WebAssembly (`xtalforge-wasm` on
  npm).
- Documentation site built with VitePress (`docs/`).

## Working on the project

```bash
scripts/check-rules.sh                 # 500-line limit and warning-suppression checks
cargo clippy --workspace --all-targets -- -D warnings
cd docs && npm install && npm run docs:dev
```

See [`docs/guide/contributing.md`](docs/guide/contributing.md) for the full
checklist and the roadmap workflow.

## License

Licensed under either of

- Apache License, Version 2.0 ([`LICENSE-APACHE`](LICENSE-APACHE))
- MIT license ([`LICENSE-MIT`](LICENSE-MIT))

at your option. Third-party code and data: `THIRD_PARTY.md`.
