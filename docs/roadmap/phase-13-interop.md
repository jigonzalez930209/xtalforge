# Phase 13 — Ecosystem interoperability

**Goal:** xtalforge drops into the workflows people already have (ASE,
pymatgen, phonopy, VASP, LAMMPS, Quantum ESPRESSO, OPTIMADE, the big
databases) with lossless round trips. Being the market standard means being
the easiest library to adopt next to the existing ones, not replacing them.

**Release gate:** `v1.6.0`.

---

## Task 13.1 — Python object adapters

- [ ] 13.1.1 ASE
  Description: `Structure.from_ase(atoms)` / `to_ase()` keeping cell, pbc,
  tags, charges, magnetic moments, constraints (as metadata) and the
  calculator results; an ASE `Calculator` wrapping the native UFF and the
  `sci-form` xTB methods. Optional dependency; golden round trips.
- [ ] 13.1.2 pymatgen
  Description: `from_pymatgen`/`to_pymatgen` for `Structure`, `Molecule` and
  disordered sites (partial occupancy kept); `SpacegroupAnalyzer` parity
  checks in the e2e battery.
- [ ] 13.1.3 phonopy and spglib-compatible outputs
  Description: `to_phonopy_atoms`, and a `get_symmetry_dataset`-shaped dict so
  spglib users can switch with one import.
- [ ] 13.1.4 RDKit
  Description: molecule exchange (`from_rdkit`/`to_rdkit`) for the SMILES and
  fragment tools.

## Task 13.2 — More file formats

- [ ] 13.2.1 Quantum ESPRESSO `pw.x` input and output
- [ ] 13.2.2 CASTEP `.cell`, FHI-aims `geometry.in`, CP2K, ABINIT, CRYSTAL
- [ ] 13.2.3 mmCIF/PDBx reader and PDB reader (PDB is write-only today)
- [ ] 13.2.4 OUTCAR/vasprun.xml (energies, forces, stress, trajectories)
- [ ] 13.2.5 MOL2, SDF/MOL V2000/V3000
- [ ] 13.2.6 extxyz with arbitrary per-atom and per-frame properties
  Description (13.2.x): readers and writers through the format registry, the
  corpus extended with one sample per format, round-trip tests and the
  loss-warning model (Task 3.1). Each format is one subtask commit series.

## Task 13.3 — Remote data

- [ ] 13.3.1 OPTIMADE client
  Description: query any OPTIMADE provider (`filter` string, pagination) and
  return `Structure`s; behind the `net` feature (pure-Rust HTTP), fetch-based
  in WASM, CLI `xtalforge fetch optimade …`.
- [ ] 13.3.2 COD and Materials Project fetchers
  Description: `fetch cod 1000041`, `fetch mp mp-149 --api-key`; cache under
  the model cache directory; license metadata attached to the structure.
- [ ] 13.3.3 OPTIMADE export
  Description: serialize structures as OPTIMADE `structures` entries.

## Task 13.4 — Plugins outside Rust

- [ ] 13.4.1 Python plugins
  Description: register readers, writers, calculators and verbs from Python
  entry points (`xtalforge.plugins`), visible to the CLI and the agent.
- [ ] 13.4.2 Plugin registry page in the docs and a template repository.

## Task 13.5 — Surfaces and docs

- [ ] 13.5.1 Migration guides "from ASE", "from pymatgen", "from spglib",
  each a side-by-side table; parity matrix and API pages.
