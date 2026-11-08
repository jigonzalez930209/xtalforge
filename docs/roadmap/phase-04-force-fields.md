# Phase 4 — Force fields, optimization, vibrations and engines

**Goal:** the `Calculator` API with native UFF/UFF4MOF, QEq and EQeq charges,
Ewald sums, symmetry-respecting optimizers (atoms and cell), constraints and
relaxed scans, Gamma-point vibrations and IR, phonons, thermochemistry,
quasi-harmonic free energies, equations of state, hydrogen placement and the
"Prepare for simulation" workflow, external engines (DFTB+, xTB, ORCA input),
and the optional in-process sci-form quantum methods.

**Reference:** `Crystal-Builder/xtal/ff/` (api, registry, record, markers,
uff/{params,typer,terms,calculator,qeq}, charges/{eqeq,sites}, ewald,
optimize 1532 lines, constraints, coordinates, scan 862, vibrations, infrared,
phonons (phonopy), thermo, qha, eos, hydrogens, isotopes, external,
dftb/{calculator,hsd,outputs,params}, xtb/calculator), `xtal/orca/`,
`xtal/core/prepare.py` (1833), `xtal/commands/{ff,prepare}.py`,
`xtal/modules/{forcefield,scan,thermo,qha,dftb,orca,probe,process}.py`;
sci-form `src/forcefield/`, `src/charges_eeq/`,
`src/materials/{geometry_opt,variable_cell}.rs`, `src/xtb/`.

Crate: `xtalforge-ff`. Completing this phase tags **v0.1.0**.

---

## Task 4.1 — Calculator API, registry and run records

- [x] 4.1.1 `Calculator` trait
  Port `ff/api.py`: a calculator is built once for one structure (typing,
  topology, parameters happen at construction) and then answers repeatedly:
  given Cartesian positions and a cell, return energy, forces and, when
  supported, stress. Capability flags: periodic, needs charges, supports
  stress, in-process or external, thread-safe. Units: eV, eV/Angstrom, GPa.
- [x] 4.1.2 Registry and availability
  Port `ff/registry.py`: name to calculator factory, the list read by the CLI
  `--engine`, the Python API and the agent tools. Availability probing returns
  a reason when an engine cannot run (binary missing, feature not compiled,
  model weights absent).
- [x] 4.1.3 Dummy atoms excluded from energies
  Port `ff/markers.py`: dummy atoms (connection points, markers) never enter
  an energy; calculators see a filtered view with index maps back.
- [x] 4.1.4 Run record
  Port `ff/record.py` (339 lines): what a run writes down while it runs
  (engine, options, per-step energies and forces, timings, warnings, final
  result), serializable into the workspace run folder (Task 3.12) and the
  report system (Task 5.10).

## Task 4.2 — UFF and UFF4MOF

- [x] 4.2.1 Parameters
  Port `ff/uff/params.py` (559 lines) including the UFF4MOF extension types.
  Tables to `data/uff.json` with citations (Rappé et al. 1992; Addicoat et
  al. 2014 and Coupry et al. 2016 for UFF4MOF). Compare with sci-form UFF
  parameters and record any difference.
- [x] 4.2.2 Atom typing
  Port `ff/uff/typer.py` (850 lines; split into perception, rules, metals,
  reasons). Typing is where wrong answers come from, so every assignment
  carries a **reason** string (for example "C_R: in an aromatic ring of 6");
  missing hydrogens show up in the reasons. Metal and cluster types for
  UFF4MOF (paddlewheels, Zn4O, Zr6 nodes). User overrides via `SetAtomTypes`.
- [x] 4.2.3 Energy terms and analytic gradients
  Port `ff/uff/terms.py` (783 lines): one type per term (bond, angle with its
  coordination-specific forms, torsion, inversion, van der Waals), each
  holding index arrays and parallel parameter arrays and returning
  `(E, dE/dx)` for the whole crystal at once. Cell gradient (virial) for
  stress. Validate every gradient against finite differences (sci-form
  `fd_gradient.rs`).
- [x] 4.2.4 UFF calculator
  Port `ff/uff/calculator.py` (752 lines): build once (typing, angle and
  torsion enumeration over the periodic bond graph, van der Waals pair lists
  across images with a cutoff and neighbour-list skin), then evaluate cheaply
  as the geometry moves under fixed topology. Optional electrostatics with
  charges from Task 4.3.
- [x] 4.2.5 Validation
  Energies, forces and stress against Crystal-Builder golden data for small
  molecules, molecular crystals and the corpus MOFs; tolerance documented.

## Task 4.3 — Charges and electrostatics

- [x] 4.3.1 QEq
  Port `ff/uff/qeq.py`: Rappé-Goddard charge equilibration with UFF's
  parameter table, for structures without charges.
- [x] 4.3.2 EQeq
  Port `ff/charges/eqeq.py`: extended charge equilibration (Wilmer, Kim and
  Snurr 2012) with measured ionization energies (data file from
  `ff/charges/data`), periodic, validated against `Crystal-Builder/tests/data/eqeq`.
- [ ] 4.3.3 EEQ from sci-form
  Port sci-form EEQ as a third method (`eeq`), with its parameters as a data
  file; compare with EQeq on the corpus and document when each is preferable.
- [x] 4.3.4 Charges per site
  Port `ff/charges/sites.py`: charges computed on the P1 cell become charges
  per asymmetric site (mean over images) so they can be stored
  (`SetCharges`) and written to CIF `_atom_site_charge` or LAMMPS.
- [x] 4.3.5 Ewald summation
  Port `ff/ewald.py`: real-space plus reciprocal-space sums with the
  splitting parameter and cutoffs chosen from a target accuracy; energy,
  forces, stress; deterministic summation order.

## Task 4.4 — Optimizers

- [x] 4.4.1 Symmetry-respecting degrees of freedom
  Port the core idea of `ff/optimize.py`: **the variables are the asymmetric
  unit, not the P1 cell**. The calculator works on P1, forces are projected
  onto the symmetry-allowed displacements of each Wyckoff site, and the space
  group is preserved on every step so the result can be written with its
  symmetry. P1 optimization is the special case.
- [x] 4.4.2 Cell degrees of freedom
  Port `CellFreedom`: which cell quantities are held (any of a, b, c, alpha,
  beta, gamma, or the volume). Nothing held is the ordinary variable cell
  whose projector is the space group's own; holding the volume while the
  shape relaxes is the scan used for flexible frameworks. External pressure.
- [x] 4.4.3 FIRE and L-BFGS
  Port `fire` and `lbfgs` (and compare with sci-form `geometry_opt.rs`,
  `variable_cell.rs`, `minimizer.rs`): convergence on fmax and energy
  change, maximum step, line search for L-BFGS, restarts, cancellation token,
  progress callback with per-step data, trajectory written as extended XYZ.
- [x] 4.4.4 Applying the result
  `ApplyOptimizedGeometry` command: the optimized asymmetric unit and cell
  replace the old ones in one undoable step.
- [x] 4.4.5 Equation of state
  Port `ff/eos.py`: fit Birch-Murnaghan (and Murnaghan, Vinet) to the points of
  a constant-volume scan (each point relaxed with shape and atoms free) to get
  the bulk modulus and its derivative, with fit quality.

## Task 4.5 — Coordinates, constraints and relaxed scans

- [x] 4.5.1 Scan coordinates
  Port `ff/coordinates.py` (571 lines): distance, angle, torsion, angle
  between two planes, and the cell's own parameters and volume. An endpoint is
  a **group of atoms** (a group of one is an atom; a group of many is its
  centroid, recomputed each step), not a dummy atom.
- [x] 4.5.2 Constraints
  Port `ff/constraints.py`: hold a coordinate while everything else relaxes
  by projecting out its gradient rather than freezing the atoms that define it
  (freezing four atoms to hold one dihedral removes twelve degrees of freedom
  to constrain one). Also fixed atoms and fixed Cartesian components.
- [x] 4.5.3 Relaxed scan driver
  Port `ff/scan.py` (862 lines; split) and `modules/scan.py`: for each value of
  an axis set the coordinate, relax the rest, record energy and geometry,
  write each point to disk before starting the next (`scan.csv`, per-point
  CIF), resume after interruption. Rigid scans as an option. Two-dimensional
  scans if Crystal-Builder supports them (verify).

## Task 4.6 — Hydrogens

- [ ] 4.6.1 Hydrogen placement
  Port `ff/hydrogens.py` (513 lines): put back the hydrogens an X-ray
  structure never had, from the bond graph, hybridization and ideal geometry
  (methyl, methylene, aromatic CH, NH, OH with orientation choices), periodic
  aware, as the `AddHydrogens` command. Report each addition.

## Task 4.7 — Prepare for simulation

- [ ] 4.7.1 Workflow steps
  Port `core/prepare.py` (1833 lines; split one file per step). A refined
  structure describes the average over all cells; several things true of an
  average are not true of any one cell. Steps: resolve disorder (choose
  disorder groups or parts, never compute half atoms), resolve partial
  occupancy and mixed sites, remove or keep solvent and counter-ions, add
  missing hydrogens (Task 4.6), check charge balance, assign types and
  charges, check close contacts, produce a P1 or symmetric result as asked.
- [ ] 4.7.2 Report
  `PrepareReport`: every decision, with the atoms involved and the reason, so
  the result is explainable and reproducible; the `Prepare` command applies it
  as one undoable step.
- [ ] 4.7.3 Validation
  Use `Crystal-Builder/resources/samples/prepared/` (inputs and their
  prepared outputs) as golden data.

## Task 4.8 — Vibrations, IR, phonons, thermochemistry

- [ ] 4.8.1 Gamma-point Hessian
  Port `ff/vibrations.py`: each free atom displaced ±delta along each axis
  (6 × free + 1 evaluations), central differences of forces, symmetrized and
  mass-weighted Hessian (isotope masses from Task 1.1.3), projection of
  translations, frequencies and modes. Works with any calculator.
- [ ] 4.8.2 IR spectrum
  Port `ff/infrared.py`: intensities `I_k ∝ |Σ_i Z*_i · e_ik / √m_i|²` with
  Born effective charges approximated from a charge model (the docs page must
  state this approximation as Crystal-Builder does), broadened spectrum.
- [ ] 4.8.3 Phonons (replacing phonopy)
  Crystal-Builder uses phonopy. Implement natively: symmetry-reduced
  displacements in a supercell, force constants from the forces, dynamical
  matrix, frequencies on a q-mesh (Monkhorst-Pack; study sci-form `kmesh.rs`)
  and along the k-path of Task 2.10, phonon DOS. Validate against phonopy
  results recorded in the golden data for simple solids.
- [ ] 4.8.4 Thermochemistry
  Port `ff/thermo.py`: harmonic statistical mechanics over any frequency list
  (Gamma Hessian, q-mesh with weights, or an engine's own modes): zero-point
  energy, H, S, G, Cv versus temperature; handling of imaginary modes.
  `modules/thermo.py` behaviour: relax a copy first when asked.
- [ ] 4.8.5 Quasi-harmonic free energy
  Port `ff/qha.py` and `modules/qha.py`: volume scan, vibrations at each
  converged volume, `F(V, T) = E(V) + F_vib(V, T)` minimized at every
  temperature; thermal expansion and bulk modulus versus temperature; each
  point written to disk before the next starts.

## Task 4.9 — External engines

External programs are optional; nothing in the default build needs them.

- [ ] 4.9.1 Process runner and probes
  Port `ff/external.py`, `modules/process.py` and `modules/probe.py`: write
  inputs into a scratch run folder, launch the binary with a timeout and
  cancellation, capture logs, read results back. A probe checks that a binary
  actually runs (architecture, missing libraries) instead of only checking
  that the file exists. Discovery from `PATH` or configuration. Not available
  in WASM (recorded exception).
- [ ] 4.9.2 DFTB+
  Port `ff/dftb/` (calculator 445 lines, hsd, outputs, params) and
  `modules/dftb.py`: HSD input generation, Slater-Koster parameter set
  location, k-points, parsing `detailed.out` and `results.tag` for energy,
  forces and stress. Tests use the recorded outputs in
  `Crystal-Builder/tests/data/dftb` so CI needs no binary.
- [ ] 4.9.3 xTB
  Port `ff/xtb/calculator.py` (484 lines): GFN1-xTB, GFN2-xTB and GFN-FF
  through the `xtb` binary, periodic where supported.
- [ ] 4.9.4 ORCA input writer
  Port `xtal/orca/` (input 439 lines, catalogue 375 lines plus data) and
  `modules/orca.py`: ORCA is molecular, so the structure is written as a
  cluster with every molecule made whole across cell faces using the bond
  graph; method and basis chosen from a catalogue of ORCA 6.1 functionals
  stored as data. Only writes inputs; never runs ORCA.

## Task 4.10 — In-process quantum methods from sci-form (optional)

- [ ] 4.10.1 Adapter
  Behind feature `sci-form`: implement `Calculator` for sci-form GFN0, GFN1
  and GFN2-xTB, so users get xTB without the external binary (GFN-FF stays
  external). Convert our structures to sci-form inputs; map errors.
- [ ] 4.10.2 Capability limits
  Document which sci-form methods support periodic boundary conditions and
  which are molecular only; refuse invalid use with a clear error;
  cross-check against the xtb binary on molecules.

## Task 4.11 — Machine-learning potentials (decision only)

Implementation is Task 8.4.

- [ ] 4.11.1 Backend decision
  Crystal-Builder runs MACE, ORB-v3 and MatterSim through PyTorch and ASE
  (`ff/ase_engine.py`). Evaluate ONNX Runtime (`ort`, native library) versus
  pure Rust inference (`burn`, `candle`) for these models: operator coverage
  (equivariant layers), accuracy, speed, binary size, WASM feasibility. Write
  `docs/algorithms/ml-backend.md` with the decision.

## Task 4.12 — Surfaces and docs

- [ ] 4.12.1 CLI
  `ff energy`, `ff optimize`, `ff scan`, `ff eos`, `ff vibrations`, `ff ir`,
  `ff phonons`, `ff thermo`, `ff qha`, `ff charges`, `ff types`
  (Crystal-Builder `xtal types`), `prepare`, `hydrogens`, `engines` (lists
  calculators and their availability), `orca input`.
- [ ] 4.12.2 Python and WASM
  Calculator objects; optimizers with progress callbacks (Python) and
  progress events with chunked execution (WASM, to keep the page responsive);
  WASM exposes UFF, QEq, EQeq, Ewald, optimizers, vibrations; external engines
  are recorded exceptions.
- [ ] 4.12.3 Docs and demos
  Algorithm pages: UFF typing with reasons, UFF terms, QEq and EQeq, Ewald
  parameter choice, symmetry-respecting optimization, FIRE and L-BFGS,
  projected constraints, finite-difference Hessian, phonons, harmonic
  thermochemistry, QHA, EOS fits, Prepare steps. Demo: UFF optimization of a
  small MOF in the browser with an energy trace.
- [ ] 4.12.4 Parity and release
  Parity matrix; changelog; tag **v0.1.0** (see `AGENTS.md` section 5).
