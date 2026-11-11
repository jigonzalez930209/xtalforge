# Phase 10 — MOF platform

**Goal:** make xtalforge the reference library for metal-organic frameworks:
read and clean experimental MOF CIFs, understand them (building units,
topology, metal sites, pores), build and modify them, compute the descriptors
used in screening, and run high-throughput studies on databases such as
CoRE MOF and QMOF. Phases 1-8 delivered the MOF **builder** (PORMAKE port,
Tasks 6.11-6.14) and generic porosity/topology (Phase 5); this phase adds
everything around it. The other road-to-2.0 phases stay in scope; MOF use
cases set their priorities (see "MOF priorities" in each).

**Release gate:** `v1.3.0` when Tasks 10.1-10.7 are complete.

**Gap check (2026-11-11):** `crates/*/src` has no open-metal-site detection,
no deconstruction into SBUs and linkers, no MOFid/MOFkey, no adsorption
(Widom, Henry, GCMC) and no native pore-limiting diameter (only through the
Zeo++ bridge, Task 5.5).

**References:** MOFid (Bucior et al. 2019), MOFChecker and MOSAEC (CIF
quality), CoRE MOF 2019/2024, QMOF, Zeo++ (PLD/LCD/ASA), RASPA
(force-field mixing, Widom, GCMC), PORMAKE, ToBaCCo, CrystalNets.jl
(topology). Algorithms are reimplemented from the papers; data sets and
ported code go into `THIRD_PARTY.md` first (rule 8). Golden data from the
reference tools on the CoRE MOF subset already in `tests/external/`.

---

## Task 10.1 — Clean and validate MOF CIFs

- [ ] 10.1.1 MOF checker
  Description: one call reporting the problems of a MOF structure: missing or
  extra hydrogens, overlapping atoms, unbound atoms and fragments, atoms with
  impossible valence, hypervalent carbon, charged frameworks, disordered
  sites, missing metal coordination, isolated solvent; each with severity,
  atoms and a fix verb. Extends `diagnose` (Task 8.1.3). Validated against
  MOFChecker flags on the CoRE MOF subset.
- [ ] 10.1.2 Solvent and counter-ion removal
  Description: detect free and bound solvent molecules (H2O, DMF, EtOH, …,
  from a data file) and remove them, keeping charge-balancing ions unless
  asked; report what was removed. Matches the CoRE MOF "FSR"/"ASR" protocols.
- [ ] 10.1.3 Disorder resolution for MOFs
  Description: pick one ordered configuration of disordered linkers and
  solvent by occupancy groups and clash removal (reuses Task 9.5.2).
- [ ] 10.1.4 Metal oxidation states and charge balance
  Description: per-metal oxidation states from BVS (Task 9.3.3) and the
  linker charges; flags charged frameworks (MOSAEC-style).

## Task 10.2 — Understand the framework

- [ ] 10.2.1 Deconstruction into building units
  Description: split the framework into metal nodes (SBUs), organic linkers
  and their connection points using the "all-node", "single-node" and
  "metal-oxo" algorithms of MOFid; output the units as molecules with
  connection points and the underlying net.
- [ ] 10.2.2 Topology from the structure
  Description: contract the deconstructed units to a periodic net and identify
  it with the RCSR data (Task 5.6), including interpenetration (Task 5.8);
  agreement with CrystalNets.jl/MOFid on the CoRE subset.
- [ ] 10.2.3 MOFid and MOFkey
  Description: canonical SMILES of the units plus topology (`MOFid`) and the
  hashed `MOFkey`; used as the identity of a MOF in datasets (Task 16.3).
- [ ] 10.2.4 Open metal sites
  Description: detect coordinatively unsaturated metals by coordination
  geometry (Task 9.3.2), with type and accessibility to a probe.
- [ ] 10.2.5 Linker and node descriptors
  Description: linker length, functional groups (Task 1.13), metal identity
  and nuclearity, node connectivity, as data for ML.

## Task 10.3 — Porosity descriptors for screening

- [ ] 10.3.1 Native PLD, LCD and accessible surface/volume
  Description: pore-limiting and largest cavity diameters, accessible surface
  area (m²/g, m²/cm³), accessible and probe-occupiable volume, void fraction,
  from the distance grid with channel percolation; Zeo++ results within
  documented tolerances. Closes the porosity item of "Known subsets".
- [ ] 10.3.2 Channel dimensionality and pore-size distribution
- [ ] 10.3.3 Geometric descriptors bundle
  Description: one call returning the standard CoRE MOF descriptor table
  (density, PLD, LCD, ASA, AV, VF, channel dimension) for many structures in
  parallel.

## Task 10.4 — Adsorption basics

- [ ] 10.4.1 Force-field assignment for host-guest
  Description: UFF/DREIDING Lennard-Jones parameters for framework atoms,
  TraPPE guests (CO2, N2, CH4, H2, H2O as data files), Lorentz-Berthelot
  mixing, framework charges from EQeq/EEQ (Task 4.3).
- [ ] 10.4.2 Widom insertion and Henry coefficients
  Description: deterministic (seeded) Widom insertion with Ewald for charged
  guests; Henry coefficient and heat of adsorption at infinite dilution;
  cross-checked with RASPA on a few CoRE MOFs.
- [ ] 10.4.3 GCMC isotherms
  Description: grand-canonical Monte Carlo for rigid frameworks and rigid
  guests, Peng-Robinson fugacity, single-component isotherms, seeded and
  thread-count independent (rule 6). Mixtures are a later subtask.
- [ ] 10.4.4 Energy grids
  Description: precomputed guest energy grids (shared with the porosity grid)
  for fast Widom/GCMC and for visualising binding sites in the playground.
- [ ] 10.4.5 RASPA input export
  Description: write RASPA `simulation.input`, framework and force-field
  files for users who need the full engine (external engine, Task 4.9 style).

## Task 10.5 — Build and modify MOFs

- [ ] 10.5.1 Full PORMAKE and RCSR bundles and global RMSD orientation
  Description: the `mof-data` feature ships the complete topology and block
  libraries; replace the greedy orientation with the global RMSD fit. Closes
  the MOF item of "Known subsets" (moved from Task 14.4.3).
- [ ] 10.5.2 Building units from real MOFs
  Description: blocks extracted by 10.2.1 become builder blocks, so a user can
  rebuild a known MOF on another topology.
- [ ] 10.5.3 Post-synthetic modification
  Description: functionalise linkers (–NH2, –OH, –CH3, –F, …) at chosen or all
  symmetry-equivalent positions, using substitution with clearance (Task 6.3).
- [ ] 10.5.4 Defects in MOFs
  Description: missing linker and missing cluster defects with capping
  (formate, water/hydroxide), at a given concentration with a seed.
- [ ] 10.5.5 Interpenetrated and layered variants
  Description: build n-fold interpenetrated frameworks and 2D MOF stacks
  (AA/AB/offset) with interlayer spacing.
- [ ] 10.5.6 Mixed-metal and mixed-linker frameworks
  Description: ordered and random (seeded) distributions of metals or linkers
  over symmetry-equivalent positions.

## Task 10.6 — MOF databases and screening

- [ ] 10.6.1 CoRE MOF and QMOF readers
  Description: read the published database archives and their property
  tables, with license and citation metadata; fetch on demand (Task 13.3).
- [ ] 10.6.2 Screening pipeline
  Description: a ready workflow (Task 16.1) that cleans, deconstructs and
  computes descriptors and Henry coefficients for a whole database, resumable,
  writing a dataset with MOFkeys.
- [ ] 10.6.3 Similarity and deduplication of MOFs
  Description: by MOFkey, then structure matcher (Task 9.1) for ties.

## Task 10.7 — Surfaces and docs

- [ ] 10.7.1 Four surfaces, parity matrix and API page `docs/api/mof.md`
  extended; algorithm pages for deconstruction, open metal sites, PLD and
  Widom/GCMC.
- [ ] 10.7.2 MOF guide and tutorials: "from a CSD/CoRE CIF to a clean,
  described MOF", "screen CoRE MOF for CO2 Henry coefficients", "build and
  functionalise UiO-66"; playground MOF tab (deconstruction view, open metal
  sites, energy grid).
- [ ] 10.7.3 e2e references: Zeo++, MOFid and RASPA outputs on the CoRE MOF
  subset committed under `tests/e2e/reference/`.

## MOF priorities for the other phases

| Phase | What MOFs need first |
|---|---|
| 9 Core | 9.1 matcher, 9.3 coordination/Voronoi/BVS, 9.5.1-9.5.2 defects and disorder |
| 11 Usability | 11.3 closes SMILES, fragments, substitution and guest placement on all surfaces (linkers and guests) |
| 12 Native bindings | NumPy for grids and isotherms, typed results for descriptors |
| 13 Interop | ASE/pymatgen adapters, CoRE MOF/QMOF/CSD-derived fetchers, RASPA formats |
| 14 Completeness | porosity items not covered here, cell-flexible optimization (breathing) |
| 15 Performance | parallel grids, Widom/GCMC throughput on large cells |
| 16 Workflows | datasets keyed by MOFkey, provenance of every descriptor |
