# Phase 13 — Ecosystem interoperability

**Goal:** xtalforge drops into the workflows people already have (ASE,
pymatgen, phonopy, VASP, LAMMPS, Quantum ESPRESSO, OPTIMADE, the big
databases) with lossless round trips. Being the market standard means being
the easiest library to adopt next to the existing ones, not replacing them.

**Release gate:** `v1.6.0`.

---

## Task 13.1 — Python object adapters

- [x] 13.1.1 ASE
  Description: `Structure.from_ase(atoms)` / `to_ase()` keeping cell, pbc,
  tags, charges, magnetic moments, constraints (as metadata) and the
  calculator results; an ASE `Calculator` wrapping the native UFF and the
  `sci-form` xTB methods. Optional dependency; golden round trips.
  Note: tags, charges, magnetic moments, constraints and calculator results
  are kept under the `interop` metadata entry; the calculator uses native UFF.
- [x] 13.1.2 pymatgen
  Description: `from_pymatgen`/`to_pymatgen` for `Structure`, `Molecule` and
  disordered sites (partial occupancy kept); `SpacegroupAnalyzer` parity
  checks in the e2e battery.
  Note: a disordered pymatgen site becomes one core site per species with its
  occupancy and a shared disorder group; the pymatgen composition is kept in
  site properties.
- [x] 13.1.3 phonopy and spglib-compatible outputs
  Description: `to_phonopy_atoms`, and a `get_symmetry_dataset`-shaped dict so
  spglib users can switch with one import.
  Note: `to_phonopy_atoms` uses `phonopy.structure.atoms.PhonopyAtoms`;
  `get_symmetry_dataset` returns spglib-shaped lists.
- [x] 13.1.4 RDKit
  Description: molecule exchange (`from_rdkit`/`to_rdkit`) for the SMILES and
  fragment tools.

## Task 13.2 — More file formats

- [x] 13.2.1 Quantum ESPRESSO `pw.x` input and output
- [x] 13.2.2 CASTEP `.cell`, FHI-aims `geometry.in`, CP2K, ABINIT, CRYSTAL
- [x] 13.2.3 mmCIF/PDBx reader and PDB reader (PDB is write-only today)
  Note: mmCIF/PDBx is read-only; PDB now reads CRYST1, ATOM/HETATM and CONECT
  in addition to writing.
- [x] 13.2.4 OUTCAR/vasprun.xml (energies, forces, stress, trajectories)
  Note: `read_outcar`/`read_vasprun` return trajectories; the registry
  `outcar`/`vasprun` entries expose the last frame as a structure.
- [x] 13.2.5 MOL2, SDF/MOL V2000/V3000
  Note: SDF/MOL V2000 is read/write; V3000 is a documented read subset.
- [x] 13.2.6 extxyz with arbitrary per-atom and per-frame properties
  Description (13.2.x): readers and writers through the format registry, the
  corpus extended with one sample per format, round-trip tests and the
  loss-warning model (Task 3.1). Each format is one subtask commit series.

## Task 13.3 — Remote data

- [x] 13.3.1 OPTIMADE client
  Description: query any OPTIMADE provider (`filter` string, pagination) and
  return `Structure`s; behind the `net` feature (pure-Rust HTTP), fetch-based
  in WASM, CLI `xtalforge fetch optimade …`.
  Note: parsing is available without `net`; the blocking client is behind
  `net` and the WASM package uses browser `fetch`.
- [x] 13.3.2 COD and Materials Project fetchers
  Description: `fetch cod 1000041`, `fetch mp mp-149 --api-key`; cache under
  the model cache directory; license metadata attached to the structure.
  Note: the Materials Project API key is intentionally host-side; the WASM
  package offers COD and OPTIMADE fetch instead.
- [x] 13.3.3 OPTIMADE export
  Description: serialize structures as OPTIMADE `structures` entries.
  Note: `remote::optimade_export` and the registry `optimade` writer produce
  one OPTIMADE `structures` entry.

## Task 13.4 — Plugins outside Rust

- [x] 13.4.1 Python plugins
  Description: register readers, writers, calculators and verbs from Python
  entry points (`xtalforge.plugins`), visible to the CLI and the agent.
  Note: Python entry points are in-process; the Rust registries remain the
  extension point visible to the CLI and agent, as recorded in the parity
  exceptions.
- [x] 13.4.2 Plugin registry page in the docs and a template repository.

## Task 13.5 — Surfaces and docs

- [x] 13.5.1 Migration guides "from ASE", "from pymatgen", "from spglib",
  each a side-by-side table; parity matrix and API pages.
