# Phase 18 — Simulation engine bridge

**Goal:** make xtalforge the layer between a structure and a real
calculation: write correct, reproducible inputs for the major simulation
codes, read their outputs back into typed results, analyse electronic
structure, and submit jobs from workflows. The external programs stay
optional (AGENTS.md rule 3); everything that only writes or parses text is
pure Rust and available on every surface.

**Release gate:** `v1.11.0`.

---

## Task 18.1 — Input sets

- [ ] 18.1.1 Input-set framework
  Description: an `InputSet` trait (structure + preset + overrides → named
  files) with presets stored as data files under `data/input-sets/`, each
  versioned and carrying its provenance (Task 16.2). Overrides merge
  deterministically; unknown keys are an error unless explicitly allowed.
  Acceptance: the same inputs give byte-identical files on every surface.
- [ ] 18.1.2 VASP
  Description: INCAR/KPOINTS/POSCAR/POTCAR-spec (POTCAR names only, never the
  licensed files) with relax, static, band, DOS and MOF presets; k-point
  density from a target spacing; MAGMOM and LDAU from site properties.
- [ ] 18.1.3 Quantum ESPRESSO
  Description: `pw.x` input with pseudopotential mapping (SSSP names),
  relax/vc-relax/scf/nscf/bands presets, k-path from Task 2.7.
- [ ] 18.1.4 CP2K
  Description: GPW input with MOLOPT basis/GTH pseudopotential mapping, cell
  and coordinates, DFT-D3; presets for MOF geometry optimisation.
- [ ] 18.1.5 LAMMPS and GULP
  Description: LAMMPS data plus input script with UFF/UFF4MOF parameters from
  the force-field crate; GULP input for lattice relaxations.
- [ ] 18.1.6 RASPA and ORCA refresh
  Description: align the existing RASPA export (Phase 10) and ORCA input
  (Phase 4) with the `InputSet` trait without changing their output.

## Task 18.2 — Output parsers

- [ ] 18.2.1 Common result model
  Description: typed `CalculationResult` (final structure, trajectory,
  energies, forces, stress, magnetic moments, charges, convergence flags,
  timings, warnings) shared by all parsers.
- [ ] 18.2.2 VASP: vasprun.xml (streaming), OUTCAR, OSZICAR, EIGENVAL,
  DOSCAR, PROCAR.
- [ ] 18.2.3 Quantum ESPRESSO: XML output, `pw.x` text output, bands and
  projwfc outputs (extends the Phase 13 QE support).
- [ ] 18.2.4 CP2K, LAMMPS log/dump and GULP outputs.
- [ ] 18.2.5 Fixtures: one recorded output per code and preset under
  `crates/xtalforge-io/tests/fixtures/engines/`, with licenses noted.

## Task 18.3 — Electronic structure analysis

- [ ] 18.3.1 Band structures
  Description: band data on the k-path, band gap (direct/indirect, spin
  resolved), effective masses at the band edges by parabolic fits.
- [ ] 18.3.2 Density of states
  Description: total, element and orbital projected DOS, Gaussian smearing,
  d-band centre; arrays for plotting on every surface and typed buffers in
  WASM (Task 12.3).
- [ ] 18.3.3 Charges and moments: Bader/DDEC file readers and mapping of
  charges and moments back onto the sites.

## Task 18.4 — Job execution

- [ ] 18.4.1 Workflow `compute` step with engines
  Description: a workflow step that writes an input set, runs a configured
  engine and parses the result; failure, restart and resume reuse the runner
  of Task 16.1.2.
- [ ] 18.4.2 Schedulers
  Description: SLURM and PBS submission scripts, submit, poll and collect,
  behind a runtime probe; dry-run mode writes the scripts only.
- [ ] 18.4.3 Error handlers
  Description: known failure patterns (SCF not converged, ZBRENT, walltime)
  detected from the outputs with a documented corrective action per pattern.

## Task 18.5 — Surfaces and docs

- [ ] 18.5.1 Parity matrix (scheduler and process execution are recorded
  WASM exceptions); API page `docs/api/engines.md`; tutorial "MOF from CIF to
  DFT band gap"; algorithm page for band and DOS analysis.
