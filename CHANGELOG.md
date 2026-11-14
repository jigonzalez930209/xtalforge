# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/) (pre-1.0: `0.0.x` until the first
essential release, then `0.x.y`).

## [1.11.0] - 2026-11-14

Phase 18 makes xtalforge the layer between a structure and a real simulation.

### Added

- `xtalforge-engines` crate and facade `engines` feature: versioned input-set
  presets for VASP, Quantum ESPRESSO, CP2K, LAMMPS, GULP, RASPA and ORCA,
  deterministic override merging and byte-identical files.
- Typed `CalculationResult` and parsers for vasprun.xml, OUTCAR, OSZICAR,
  EIGENVAL, DOSCAR, PROCAR, QE text/XML/bands/projwfc, CP2K, LAMMPS log/dump
  and GULP, plus Bader and DDEC charge readers.
- Band structures, direct/indirect/spin gaps, parabolic effective masses,
  total/element/orbital DOS, Gaussian smearing and the d-band centre.
- Workflow compute steps with scratch directories, parser registry, restart
  checkpoints, SLURM/PBS scripts, scheduler probes and the failure catalogue
  (SCF, ZBRENT, walltime, out-of-memory, missing pseudopotential).
- CLI `xtalforge engines presets|inputs|parse|bands|dos|script|failures`,
  Python `engine_*` functions and WASM `engine_*_info` exports.
- The "MOF from CIF to DFT band gap" tutorial, the electronic-structure
  algorithm page, engine fixtures with provenance and the Phase 18 parity rows.

### Notes

- Process execution and scheduler probes are host-only; WASM and Python write
  scripts and parse outputs, recorded as a parity exception.
- The engine fixtures are hand-written representative outputs under the
  repository license, not copied third-party data.
- No licensed POTCAR or pseudopotential content is generated; only names are
  emitted.

## [1.10.0] - 2026-11-13

Phase 17 closes the open subtasks and documented subsets that blocked
production use, before the scope widens in Phases 18-21.

### Added

- Full PORMAKE 0.2.3 topology/block bundle behind `mof-data`: 2,129 quotient
  nets and 867 blocks in lazy gzip files, with the conversion and license in
  `THIRD_PARTY.md`; the curated six-topology/eight-block subset stays embedded.
- On-demand CoRE MOF/QMOF download with a cache directory, SHA-256
  verification and an offline mode (`mof::fetch_core_mof`, CLI
  `mof database fetch`, Python `mof_fetch_database`, WASM
  `fetch_data_bundle`).
- Porosity completeness: probe-occupiable volume in A^3 and cm^3/g,
  accessible surface area in A^2 and m^2/g, percolating directions and a typed
  `PorosityReport`; the CLI/Python/WASM `porosity` output carries the fields.
- Aromatic hydrogen inference for SMILES (`embed_smiles_with_hydrogens`, CLI
  `--hydrogens`, Python `molecule_from_smiles(..., hydrogens=True)`, WASM
  `infer_hydrogens`), including aromatic bond order 1.5.
- Cubic/general indexing agreement (`indexing_agreement` and surfaces).
- YAML workflows (pure-Rust subset), YAML in the agent/MCP verbs, and the
  in-memory WASM runner `run_workflow_memory`; JSON stays canonical.
- Streaming XYZ/extxyz, LAMMPS dump and CIF writers with forces and custom
  per-atom columns; the readers expose the same columns.
- Fused structure-of-arrays neighbor search (no row-array materialisation).
- Intra-layer slab terminations, adsorption sites below the top layer, and
  `TerminationOptions` with defaults that preserve the 1.9 behaviour.
- Standardized-primitive structure matching, so a conventional/primitive pair
  matches directly; canonical orientation keeps the result deterministic.
- CrystalNN-style coordination with distance, solid-angle and
  electronegativity weights, documented against a pymatgen-shaped test set.
- `scripts/gate.sh` (one local gate with a summary table),
  `scripts/publish-npm-packages.sh`, the regenerated/verified golden data, a
  MOF playground tab and the Phase 17 documentation/parity rows.

### Notes

- `mof-data` is off by default and not enabled in the Python/WASM packages.
- The PORMAKE conversion drops self-loop edges and repeated index pairs that
  the quotient builder cannot express; the affected nets are not in the
  bundle.
- The probe-occupiable volume equals the accessible voxel volume; Zeo++'s
  sphere union is not reproduced, as documented on the algorithm page.
- YAML rejects anchors, aliases, tags, block scalars and multi-document files.

## [1.9.0] - 2026-11-13

Phase 16 makes high-throughput runs reproducible, inspectable and shareable.

### Added

- Versioned JSON workflows with `read`/`transform`/`analyse`/`compute`/`write`
  steps, a published JSON Schema, deterministic input expansion, per-item
  error isolation, checkpoints/resume, CLI `xtalforge run`, Python
  `Workflow.run` and agent/MCP `workflow_*` tools.
- `Provenance` records (version, options, input SHA-256, bundle versions, seed)
  with structure-metadata helpers and embedding in CIF, extxyz and project
  files; `xtalforge provenance` and Python `xtalforge.provenance`.
- JSON Lines datasets with typed properties, structure deduplication using the
  matcher plus a fingerprint prefilter, and composition/RDF/Coulomb/SOAP-like
  descriptors; Python `xtalforge.Dataset` with pandas/polars views.
- The high-throughput tutorial over a CoRE MOF subset, the workflow schema
  asset, API pages and parity rows.

### Notes

- Workflows are JSON; `.yaml`/`.yml` files load only when their content is
  JSON. A YAML parser is not a dependency in this release.
- WASM exposes workflow metadata but cannot run workflows that read and write
  files; the filesystem exception is recorded in the parity matrix.
- The SOAP-like descriptor is a documented approximation, not exact SOAP.

## [1.8.0] - 2026-11-13

Phase 15 makes performance and scale a tested contract.

### Added

- Comparative benchmark suite (`benches/compare/`) against gemmi, spglib, ASE
  and pymatgen, a committed baseline, `scripts/check-bench.sh` with a 10 %
  regression gate, and the reproducible performance page.
- Structure-of-arrays `P1Arrays`/`NeighborArrays`, lazy symmetry expansion,
  streaming XYZ/LAMMPS/CIF readers, a parallel UFF energy/force path with a
  fixed-order reduction, and the `huge` limits preset.
- Split WASM packages (`core`, `symmetry`, `analysis`, `build`) with per-package
  budgets, gzip/lazy loading for the UFF, RCSR, ICSD and k-path bundles, and
  `scripts/check-cli-start.sh` (cold `info` median around 5 ms).
- One-million-atom supercell benchmark and a regression test for the corrected
  repeat bound.

### Changed

- The facade now has optional `analysis`, `build`, `ff`, `mof`, `modules` and
  `agent` features; the default `full` keeps the previous surface. Disabling
  default features now selects the base tier (core/io/symmetry/view/commands).
- Several embedded data files are stored as gzip and decompressed on first use.

### Fixed

- `Structure::supercell` no longer truncates repeats above about four
  (50 x 50 x 50 now builds 1,000,000 atoms).

## [1.7.0] - 2026-11-13

Phase 14 closes the scientific subsets that were still simplified in `1.6.0`.

### Added

- Diffraction: Caglioti/pseudo-Voigt/TCH profiles, FCJ asymmetry,
  monochromator polarisation, X-ray/neutron/electron radiation, anisotropic
  `U_ij`, refinable `eta` and profile parameters, single-crystal structure
  factors and HKL export, kinematic electron diffraction, whole-range peak
  refinement, persistent parameter sets and Pareto sweeps.
- Force fields: `CellFreedom`, external pressure, cancellation/progress and
  trajectory output in optimizers, scans and MD; generic scan coordinates
  (cell parameters and collective variables), resumable checkpoints, NVE,
  Langevin, Nose-Hoover and NPT molecular dynamics with explicit seeds.
- Structure operations: slab cuts carry severed bonds and report the real
  `cut` count; explicit bonds propagate through symmetry images; pair-specific
  bond ranges and the metal-metal switch; the remaining undoable command types.
- Meshes: vertex normals, OBJ, glTF 2.0 and GLB export.
- Reports: HTML and JSON renderers, module registry listings and the
  `xtalforge modules` command.
- Data: NIST neutron coherent scattering lengths and named instrument presets.

### Changed

- The calculated PXRD and Rietveld paths use the new profile model; the old
  option structs and numbers remain available.
- The slab `cut` count is no longer a placeholder.
- The WASM package budget moves to 5376 KiB for the Phase 13 additions.

### Notes

- Real DFTB+/xTB execution still needs the corresponding binary (and DFTB+
  Slater-Koster data) locally; the runner reports availability instead of
  pretending the engine ran.
- The indexing-engine agreement and powder demo remain open under 14.1.5.

## [1.6.0] - 2026-11-13

Phase 13 makes xtalforge interoperable with the objects and files users
already work with. The public API is additive.

### Added

- Python object adapters: `Structure.from_ase`/`to_ase`, `from_pymatgen`/
  `to_pymatgen`, `to_phonopy_atoms`, `from_rdkit`/`to_rdkit`, the
  spglib-shaped `get_symmetry_dataset` and the native `AseCalculator`.
  Optional packages are imported only when an adapter is called.
- File formats: Quantum ESPRESSO `pw.x` input/output, CASTEP `.cell`,
  FHI-aims `geometry.in`, CP2K, ABINIT, CRYSTAL, mmCIF/PDBx, a PDB reader,
  OUTCAR, `vasprun.xml`, MOL2, SDF/MOL and extxyz properties.
- Remote data: `net`-gated OPTIMADE, COD and Materials Project fetchers with
  a deterministic cache, `optimade_export`, CLI `fetch`, and browser
  `fetch_optimade`/`fetch_cod` in the WASM package.
- Python entry-point plugins in `xtalforge.plugins`, an installable template
  in `examples/python-plugin`, and the plugin guide.
- Migration guides from ASE, pymatgen and spglib; the interoperability API
  page and 10 parity rows.

### Changed

- The CLI reads and writes through the complete format registry instead of
  the original CIF/POSCAR/XYZ switch.

### Fixed

- PDB is now readable as well as writable, and CIF/mmCIF accepts Cartesian
  `_atom_site` coordinates.

## [1.5.1] - 2026-11-13

### Added

- The xtalforge icon in the README, the docs favicon/logo/Open Graph image,
  the PyPI and npm package READMEs, and every crates.io package page; the
  repository also ships a 1280x640 social-preview image and the AppStream
  metadata lists `1.5.0`.

### Changed

- Refresh the documentation after `1.5.0`: roadmap status and release list,
  API index with the object API page and the split Python/TypeScript stub
  sources, and the open Phase 12 follow-ups in `README.md` and `AGENTS.md`.

## [1.5.0] - 2026-11-13

Phase 12 makes the Python and JavaScript surfaces object-oriented and
array-native while keeping the 1.x text/dict API. The public API is additive.

### Added

- Python `Structure` object API: `to_dict`, `write`, `positions`,
  `frac_positions`, `cell`, `numbers`, `from_arrays`, `set_positions`,
  `ff_energy`, `add_hydrogens` and `_repr_html_`; the common free functions
  accept a `Structure` as their text argument.
- Typed result classes (`SymmetryDataset`, `PxrdPattern`, `PorosityReport`,
  `OptimizeResult`) that remain dict-compatible, and `mypy --strict` stubs
  without `Any` (`scripts/check-python-stubs.sh`).
- Python batch helpers (`xtalforge.batch.read/detect_symmetry/pxrd`) with
  input order preserved, and the `pxrd(..., threads=)` keyword.
- WASM `Structure` handle (parse once, `fromText`/`fromCif`, supercell, slab,
  P1, symmetry, PXRD, energy, hydrogens, positions/cell/numbers arrays) and
  `scene_buffers`, `pxrd_buffers`, `mesh_buffers` typed-array exports.
- `crates/xtalforge-wasm/ts/xtalforge-object.d.ts` typed interfaces checked by
  `scripts/check-ts-any.sh`; `js/structure-registry.js`, `js/worker.js` and
  `js/node.js` for finalization, Web Worker and Node workflows.
- `docs/api/objects.md` and 7 parity rows.

### Fixed

- The `pxrd` Python function accepts a worker count without changing results.

## [1.4.0] - 2026-11-12

Phase 11 closes the 2026-11-11 usability audit. The public API is additive:
the CLI gains accurate help, a real thread option, error categories,
structured warnings, pipelines, `doctor`, molecules, sketches, fragments,
substitution and guests; the facade gains current docs, one-call IO and a
wider prelude; SMILES/sketch/fragment/substitution/guest now exist on all
surfaces.

### Added

- CLI `doctor`, accurate `convert`/`build` help generated from the registered
  lists, `--threads` wired to a scoped Rayon pool, `--spacegroup` overrides,
  `--out -` JSON pipelines and `molecule`, `sketch`, `fragment`, `substitute`
  and `guest` commands.
- `Error::category()` on every crate error and `category` in the CLI/WASM JSON
  error objects; Python `error_category(code)`.
- Structured `StructureWarning` entries (`code`, `message`, `atoms`, `fields`,
  `hint`) in the structure summary, CLI, Python and WASM.
- Facade `read`, `read_str` and `write`, `ReadResult`, and a wider `prelude`.
- Python `read`/`write`, `set_threads`, `deprecation_warning`,
  `molecule_from_smiles`, `sketch_roundtrip`, `fragment_library`,
  `fragment_structure`, `substitute_structure` and `place_guest_structure`;
  the `pxrd` function accepts `threads=`.
- WASM `read`, `molecule_from_smiles`, `sketch_roundtrip`,
  `fragment_library_info`, `fragment_info`, `substitute_info`,
  `place_guest_info` and `deprecation_warning`.
- `docs/api/usability.md`, `docs/api/naming.md`, `docs/guide/deprecation.md`
  and `docs/guide/recipes.md` with 30 task-oriented recipes; 11 new parity
  rows.

### Fixed

- The CLI no longer claims three formats or two builders; the help lists every
  registered format and builder.
- `--out -` always emits the versioned structure document, so the documented
  stdin/stdout pipeline works without temporary files.

## [1.3.0] - 2026-11-12

Phase 10 makes xtalforge a MOF platform: experimental CIFs can be checked and
cleaned, deconstructed into nodes, linkers and a topology, described with the
standard geometric and adsorption descriptors, modified, and screened in
batches. The public API is additive.

### Added

- MOF checking and cleaning: `mof_check` (overlaps, fragments, valence,
  disorder, metal coordination, charge), `mof_solvent_remove` with the
  embedded solvent table and the CoRE FSR/ASR policy, `mof_disorder` and
  `mof_charge_balance`.
- Framework understanding: `mof_deconstruct` (all-node SBU/linker split and
  contracted quotient net), `mof_mofid` (canonical unit SMILES, topology,
  `MOFid` string and SHA-256 `MOFkey`), `mof_open_metal_sites` and
  `mof_descriptors`.
- Native porosity: `mof_pores` with periodic union-find PLD, LCD, accessible
  surface area in m2/g and m2/cm3, accessible volume, void fraction, channel
  dimensionality and the geometric pore-size distribution.
- Adsorption: `mof_guests` (TraPPE CO2, N2, CH4, H2, H2O),
  `mof_forcefield` (UFF/EQeq host, Lorentz-Berthelot mixing),
  `mof_widom` (Henry coefficient and heat of adsorption), `mof_gcmc`
  (single-component isotherm), `mof_energy_grid` and `mof_raspa` input export.
- Modification: `mof_functionalize` (NH2/OH/CH3/F), `mof_defects` (missing
  linker and cluster defects with formate/water/hydroxide caps),
  `mof_interpenetrate`, `mof_stack_layers` and `mof_mixed` (ordered orbit or
  seeded random mixed-metal/linker frameworks).
- Screening: `mof_read_database` for CoRE MOF and QMOF directories,
  `mof_screen` for a resumable descriptor/Henry table with MOFkeys and
  `mof_deduplicate` (MOFkey first, structure matcher for ties).
- CLI `mof` subcommands (`check`, `solvent-remove`, `disorder-resolve`,
  `charge-balance`, `deconstruct`, `mofid`, `open-metal-sites`,
  `descriptors`, `pores`, `adsorption`, `modify`, `database`), Python
  functions and WASM JSON exports. `docs/api/mof.md`,
  `docs/guide/mof-workflows.md` and two new algorithm pages document the
  phase; the parity matrix gains 23 rows.
- `P1Atom` disorder groups, the MOF host-guest data files and the new
  `mof-data` feature flag.

### Fixed

- The MOF checker does not report a charge imbalance for a metal-free
  covalent framework; charge balance is assessed only when metals are
  present.

## [1.2.0] - 2026-11-11

Phase 9 adds the core operations of a standard crystal library beyond the
Crystal-Builder parity of `1.1.0`. Every capability is available in the Rust
API, the CLI, Python and WebAssembly, and the public API is additive.

### Added

- Structure comparison: `match_structures` (Niggli reduction, integer lattice
  mapping, supercell matching, species comparators), `periodic_rmsd`,
  `kabsch_alignment` and `group_structures`. CLI `compare` and `group`.
- Structure edits: `deform` (Voigt strain, deformation gradient, target volume
  and density), seeded `rattle`, `edit_structure` (add, remove, replace, move,
  species maps, partial replacement, sort and wrap) with undoable commands,
  matched `interpolate` and `optimal_supercell` (ASE
  `find_optimal_cell_shape` objective).
- Local environments: periodic `voronoi` with exact half-space clipping and
  solid angles, `coordination` with cutoff/Voronoi/CrystalNN strategies and
  continuous shape measures, `bond_valence_sums`, `oxidation_states` (ICSD
  MAP search with charge balance), `radial_distribution` (total and partial
  RDF/PDF, running coordination), bond-angle distributions and
  `structure_fingerprint`.
- Surfaces: `unique_miller_indices` under the detected point group,
  `slab_terminations` (symmetric slabs, dipole and polar report, minimal
  slab), the `wulff` construction and `adsorption_sites` with
  `place_molecule`.
- Defects and disorder: symmetry-unique `vacancies`, `substitutions` and
  `interstitials` with minimum-image supercells, `orderings` ranked by an
  in-process Ewald sum, `special_quasirandom` with pair/triplet correlation
  reports, `lattice_matches` and `stack_slabs` heterostructures.
- Mechanical and thermodynamic basics: `elastic_tensor` from analytic stress
  finite differences with Voigt-Reuss-Hill moduli and a stability flag,
  `eos_fit` (Birch-Murnaghan, Vinet, Murnaghan with errors) and
  `convex_hull` (formation energies and a two-phase simplex hull distance).
- The CLI gains `compare`, `group`, `deform` and the `features` command group;
  Python gains the corresponding functions; WASM gains JSON-returning
  exports for every operation.
- Documentation: `docs/api/features.md`, algorithm pages for comparison,
  local environments, surfaces, defects and mechanics, and a parity matrix
  with 29 new rows.
- Roadmap to `2.0.0`: Phases 9-17, focused on MOFs (core features beyond
  Crystal-Builder parity, the MOF platform, usability audit and 1.x polish, native bindings, ecosystem interoperability, scientific completeness, performance,
  workflows and the 2.0 redesign), with the 2026-11-11 usability findings.

### Fixed

- `P1Atom` carries the disorder group of its parent site, so orderings and
  SQS generation can find partial occupancies after the P1 expansion.

## [1.1.0] - 2026-11-11

First feature release after `1.0.0`, published on 2026-11-11. It adds the 3D
playground and the end-to-end audit battery, exposes the scene and
text-returning WASM helpers, and folds in the fixes from the functional audit
of every surface (Task 8.9). The public API remains backward compatible: the
new JSON and object fields are additive.

### Added

- The [3D playground](/playground) in the docs: load, transform and build a
  structure with the WASM package and inspect the same Rust scene in an
  orbitable WebGL view.
- `scene_data` on the WASM, Python and CLI (`mesh --scene`) surfaces, and
  `n_atoms` in the WASM structure summary, so a viewer no longer reconstructs
  geometry in JavaScript.
- Text-returning WASM variants of the symmetry transforms and the carbon/MOF
  builders (`to_p1_text`, `primitive_text`, `standardize_text`,
  `asymmetrize_text`, `set_space_group_text`, `change_setting_text`,
  `reduce_text`, `lower_to_subgroup_text`, `build_graphene_text`,
  `build_nanotube_text`, `build_c60_text`, `build_c20_text`, `build_mof_text`,
  `build_ribbon_text`, `build_twist_text`, `build_ztc_text`,
  `build_schwarzite_text`, `build_polymer_text`, `hydrogens_text` and
  `ff_optimize_text`).
- The end-to-end audit battery (`tests/e2e/`): 68 Playwright checks against
  gemmi, spglib, ASE and NumPy references, plus deterministic 3D screenshot
  goldens, and `scripts/e2e.sh` to run it. The playground UI checks cover the
  tabbed panels, atom selection, the Brillouin zone, the powder plot, UFF
  energy, exports, fullscreen, porosity, the MOF builder and all four powder
  refinement methods.
- A guided **Tour** tab in the playground: seventeen steps that load a sample
  or run a builder, set the camera and open the tab where the result is read
  (CIF reading, space group, Brillouin zone, measures, supercell, slab, PXRD,
  nanotube, twisted bilayer, C60, MOF, disordered carbon, polymer, porosity,
  force field, powder refinement and export). The powder plot gains a 2θ axis,
  hkl labels for the strongest lines, I/I0 in the table and recalculates when
  the structure changes; the structure panel adds the disordered-carbon
  builder.
- `build_carbon_text` on the WASM surface returns the disordered carbon as CIF
  for a viewer.

### Fixed

- PXRD intensities used the Z-based Gaussian form factor
  `Z exp(-0.5 Z^(2/3) s^2)` for every element, so relative intensities were off
  by up to 50 % (rutile `(1 1 1)` 28 % of `(1 1 0)` instead of 19 %). The
  neutral-atom International Tables (1992) coefficients are now embedded for H
  to Cf and the Si, NaCl, Fe and rutile patterns match an independent
  calculation to 0.1 %. Rietveld/Le Bail/Pawley stages that start from the
  calculated pattern inherit the change.
- `commensurate_twist` ignored the requested angle: it rounded every twist to
  0 or 60 degrees, and its "60 degree" variant shifted both layers. It now
  builds the commensurate moire cell (`(2, 1)` for 21.79 degrees, 14 atoms per
  layer) and records `twist_angle` and `twist_indices`; `n x m` repeats the
  moire cell.
- `build_polymer` packed the chain into a cube shorter than the chain, so the
  periodic images overlapped (C-C 0.67 Angstrom, coordination 5-7 for
  polyethylene). The repeat units are now joined head to tail along `a` into a
  periodic chain and the cross-section sets the density.
- `embed_smiles` had no bond-angle term, so rings collapsed (benzene with
  1.17 Angstrom 1-3 contacts and non-planar). Angle targets from the
  hybridization (180, 120, 109.47 degrees) are now part of the relaxation and
  the default iteration count is 1500.
- `build_carbon` placed carbons at mesh-patch centres without relaxing them,
  so atoms could sit 0.15 Angstrom apart, and put edge terminations along the
  sheet normal instead of the free in-plane valence. The dual graph is now
  relaxed to 1.42 Angstrom bonds and 120 degree angles, terminations follow
  the free valence, and a termination that would overlap another atom is
  dropped and reported as a bare edge. Termination and graph bonds that
  cross the cell boundary now carry their lattice image; before, 71 of the
  default dia build's bonds were stored without it and drawn 23-33 Angstrom
  long across the cell.
- `P1Cell::expand` and `Lattice::wrap` turned a tiny negative coordinate
  (such as -1e-32 from a Cartesian round trip of 0) into exactly 1.0. The UFF
  view then held a duplicated atom, so `ff vibrations` on silicon reported
  modes of +-14 000 cm^-1 (energy +1028 eV after the Hessian displacements).
  Coordinates now stay in `[0, 1)`; silicon gives three zero acoustic modes and
  degenerate optical modes.
- `le_bail` and `pawley` subtract a SNIP background before extracting
  intensities and `le_bail` uses the profile-weighted partition; on rutile
  with a 50-count baseline R drops from 0.69 to 0.27.
- MOF nodes are oriented with a Kabsch fit over every connection-to-edge
  assignment instead of the greedy mean-direction rule, so planar nodes in
  `srs`/`hcb` no longer force linkers to be stretched 0.55-1.26 times.
- `build_ztc` and `build_schwarzite` return three-coordinate sp2 carbon (the
  relaxed dual of the remeshed gyroid/Schwarz P surface, with explicit bonds)
  instead of mesh vertices with coordination 5-8; their sampled grid now
  includes the closing plane, so `ztc_surface_mesh` and
  `schwarzite_surface_mesh` are closed across the cell faces.
- Wyckoff letters are matched on the input coordinates when the input basis
  is conventional, so they follow the input origin like spglib and the CIF
  (silicon in origin choice 1 is `8a`, not moyo's `8b`).
- `eeq_charges` is now EEQ2019 (multicharge parameters, coordination-dependent
  electronegativity, erf-damped Coulomb) instead of a copy of QEq; isolated
  water gives O -0.598 like the sci-form/dftd4 reference.
- `refine_phase` with `refine_cell` starts from a Pawley cell (new first stage
  `pawley_cell`), so a cell 1 % off no longer traps the refinement.
- `build_mof` refuses two metal-free carboxylate blocks, whose oxygens would
  overlap, instead of returning the overlapping framework.
- `build_ztc` refuses cells below 12 Angstrom and `build_schwarzite` below 8
  Angstrom, where the surfaces cannot hold three-coordinate carbon.
- `eeq_charges` uses an Ewald sum for the periodic Coulomb term (MOF-5 Zn was
  +4.4 with the real-space cutoff).
- Python `ff_charges` and WASM `ff_charges_info` accept `eeq`; the CLI `info`
  JSON and text carry `n_atoms`; WASM errors that were bare strings are now
  `Error` objects with `code = "invalid_input"`; PXRD reflections with equal
  two-theta are ordered identically on native and WebAssembly builds.
- A refinement panel in the playground (Powder tab): load the bundled
  synthetic rutile pattern or a file and run Le Bail, Pawley (with size and
  strain), Rietveld or the automatic run, with an observed/calculated/
  difference plot; the tour has a matching step (seventeen steps).
- The automatic powder result exposed only the quick cubic-index edge as its
  `cell`, so the refinement panel and the CLI table showed the indexing
  candidate as the "best" cell. CLI, Python and WASM now also return
  `best_cell`: the Pawley-refined cell of the leading candidate, or its
  indexed cell when the Pawley stage was skipped (rutile automatics report
  a = 4.5993 Angstrom instead of the indexed 4.5902 Angstrom).
- Running a builder or transform in the playground left the sample selector
  on a stale sample id; it now shows "Built or transformed", and the nanotube
  tour step looks from +x so the tube lies across the view.
- Python `powder_lebail`, WASM `powder_lebail_info` and
  `powder_pawley_width_info`; WASM Pawley results carry `size`, `strain` and
  the fitted `calculated` profile. `scripts/check-surfaces.sh` and the release
  checklist compare 97 fields across CLI, Python and WASM.
- `fit_phase` refines the constant background with the scale, as documented
  (the automatic run's phase Rwp on rutile drops from 1.05 to 0.33).
- The CLI Rietveld uses the library's default iteration count, like Python
  and WASM, and exits quietly (status 141) when its output pipe closes instead
  of panicking.
- The playground radiation selector was empty: seven static `<option>`
  elements were hoisted into one static block that Vue cannot insert into a
  `<select>`. The options are now rendered from a list.
- The playground twisted-bilayer form called the builder with `m = 0`, which
  is rejected; it now defaults to the `(2, 1)` cell at 21.79 degrees.

- `Structure::supercell` copied only the asymmetric sites, so a supercell of a
  symmetric structure lost its symmetry expansion (a `2 x 2 x 2` silicon
  supercell had 8 atoms instead of 64, in the CLI and WASM). It now expands
  the P1 cell first through `Structure::as_p1`; the CLI size estimate and the
  supercell tests use the full-cell count.
- `Structure::is_p1` trusted only the Hall symbol, so a CIF with an explicit
  operation loop and no space-group name was treated as P1 and reported the
  wrong formula, density and charge. It now also requires every operation to
  be the identity.
- The CIF-inspector demo's inline silicon CIF lacked the centring translations
  and reported one atom and 0.29 g/cm3; it now uses the Hall symbol and shows
  `Si8` and 2.329 g/cm3.
- The WASM structure summary now includes the P1 `n_atoms`, matching Python
  and the whole-cell formula work released in `1.0.0`.
- `brillouin_zone_for_structure` detects the primitive cell, so silicon and
  rock salt (fcc), bcc iron and other centred cells no longer return the cube
  of the conventional reciprocal cell. The clipping now starts from a cube
  large enough for anisotropic cells, so graphene returns its hexagonal prism
  instead of a box; the cap ordering also works for non-unit plane normals.

### Changed

- The WASM budget moves from 3904 to 3968 KiB for the audit fixes (IT92
  table, commensurate twist, periodic polymer, embedder angles, carbon
  relaxation, ZTC/Schwarz P dual, Le Bail/Pawley background, MOF quaternion
  fit and `build_carbon_text`; 3943 KiB measured).
- The WASM release budget moves from 3840 to 3904 KiB for `scene_data`,
  `n_atoms` and the text-returning builders (3845 KiB measured).

## [1.0.0] - 2026-11-10

First stable release, published on 2026-11-10 on crates.io (13 crates), PyPI
(7 platform wheels plus an sdist), npm, the GitHub release and GitHub Pages.
The public API is frozen and changes follow the SemVer policy in
`docs/api/stability.md`. This release incorporates the fixes accepted during
the `1.0.0-rc.1` stabilization period.

### Fixed

- `Structure::formula`, `cell_contents`, `cell_mass`, `density`, `z` and
  `total_charge` (and therefore `info`, `structure_info`, the agent `info` verb,
  Python and WASM summaries) counted the asymmetric unit, so every structure
  with symmetry operations reported the wrong formula, density and charge: NaCl
  read from a CIF showed `Na Cl` and 0.55 g/cm3 instead of `Na4 Cl4` and
  2.19 g/cm3. They now use the P1 expansion; `z` is the gcd of the element
  counts.
- The CIF reader records `coincident_atoms_warning` when atoms of different
  sites coincide in the P1 cell with total occupancy above one (a file that
  lists the whole cell under a non-P1 group), and the CLI `info` prints it and
  `symmetry_warning` on standard error.
- Writing to a device such as `/dev/null` works: atomic writes rename the
  temporary file only for regular files, so devices and FIFOs are written in
  place.
- `graphene`: the atoms were stored with Cartesian coordinates as fractional
  ones and the `b` vector used `nx/2` instead of `ny/2`, so the sheet was
  distorted (bonds of 1.07 Angstrom) and not a supercell when `nx != ny`.
- `nanotube`: the strip was not restricted to one circumference, the
  translation vector used `gcd(n, m)` instead of `d_R = gcd(2n + m, 2m + n)`
  and the atoms were not rolled onto a cylinder. (5,5) had 166 atoms and 537
  bonds; it now has the 20 atoms and 30 bonds of the chirality formula. The
  nanotube cell is now `2R + vacuum` on `a` and `c` with the axis along `b`.
- `nanoribbon`: the vacuum stretched the strip instead of adding empty space;
  the bonds are 1.42 Angstrom again.
- `structure_info` (and the CLI `info`, Python and WASM equivalents) rejected
  the versioned `{"schema_version", "structure"}` envelope that `--json`
  commands write, so `xtalforge --json supercell ... | xtalforge info -`
  failed; both forms are accepted now.
- `powder index` (cubic search): a face-centred pattern indexed on a larger
  multiple of the true cell (NaCl returned 7.94 Angstrom); the smallest cell
  among equally good candidates is now chosen.
- `mof build`: the topology cell was fixed, so a long linker was squashed
  (0.45x for stilbene in `pcu`) and thousands of atoms overlapped; the builder
  now fits a uniform cell scale to the linker (`cell_scale` metadata). A node
  with the wrong connection count, a linker with more than two arms, a linker
  scale outside 0.8-1.25 and atom pairs closer than 0.9 Angstrom are reported
  in `warnings`, and the CLI prints them on standard error.
- `powder rietveld --refine` always modelled the Cu K alpha doublet, so a
  single-wavelength pattern drove the ADP multiplier to its bound; the new
  `--no-doublet` flag matches Python's `doublet=False`.
- `powder automatic` / `run_automatic`: the peak search did not subtract the
  background, so on a noisy pattern with a raised baseline every noise maximum
  was a peak, the indexing returned a wrong cell (16.5 Angstrom for NaCl) and
  the run did not finish in five minutes in a release build. It now subtracts
  a SNIP background and applies a noise-based threshold; the same pattern runs
  in 2.4 s and returns Fm-3m first. `PeakOptions::noise_limited` exposes the
  same search and `powder index --background` uses it (a = 5.614 instead of
  16.5).

### Tests

- Phonons of fcc argon with UFF (`phonon_argon.rs`) and the CLI `--no-doublet`
  flag have regression tests.
- Regression tests cover the whole-cell properties and the coincident-atoms
  warning, the versioned `structure_info` envelope, atomic writes to devices,
  the graphene/nanotube/ribbon geometry, cubic indexing, the MOF cell fit,
  background-subtracted peak search and the CLI warning output.
- `xtalforge-wasm` has a browser suite (`tests/browser.rs`) that runs with
  `wasm-pack test --headless --chrome`; the Node suites are unchanged.

### Documentation

- Added 29 diagrams with one shared style under `docs/public/diagrams/`
  (animated SVG; every procedure page: lattice, neighbours, supercell, slab,
  bonding, rings, P1, CIF, formats, prepare, UFF, optimisation, Ewald, PXRD,
  indexing, Le Bail, Pawley, Rietveld, extinction, symmetry, nets, MOF, carbon
  builders, vibrations, limits, validation, agent, measures) and embedded them
  in the matching algorithm pages.
- Audited the algorithm and API pages against the code. Corrected the powder
  instrument description (no wavelength or Caglioti terms), the distance-grid
  sampling and surface area, the optimiser cell freedom, the ellipsoid scale
  factors, the error mapping of the CLI and Python bindings, the quick-start
  version pin and the MOF entry point; recorded the resulting known subsets in
  the roadmap.
- Corrected the quick-start (`structure_info` takes structure JSON; use
  `Structure.from_cif` / `read_cif` for CIF), the WASM size figures, the
  migration table (no `modules` CLI command, no separate powder crate) and the
  stale Phase 0 status and `AGENTS.md` gap list; documented the Pawley
  profile-mismatch caveat.
- Corrected the API examples so every Rust snippet compiles (missing `use`
  lines and the `Calculator` trait, `ParamSet::with_param`, `f64` supercell
  matrices), the `CommandStack::push` signature, the error-code tables
  (`limit_exceeded`, `symmetry`), the isosurface algorithm name (marching
  tetrahedra) and the `external` feature wording.
- Corrected the QEq kernel description (hardness only, truncated 12 Angstrom
  sum) and stated that EEQ equals QEq at the default dielectric, the duplicated
  external-validation bullet and the tutorial commands (`powder lebail/pawley`,
  `build graphene/nanotube` take positional arguments).

## [1.0.0-rc.1] - 2026-11-10

The first release candidate for 1.0.0, published on 2026-11-10 on crates.io
(13 crates), PyPI (7 platform wheels plus sdist), npm, the GitHub release and
GitHub Pages. The public API is frozen; the stabilization period accepts fixes
only, then `1.0.0` follows.

### Added

- Phase 8 Task 8.5: distribution metadata. Every published crate has an
  `include` allow-list (src, its own data bundle, the agent skill), docs.rs
  metadata, and `scripts/check-package.sh` checks the tarballs and the 12 MB
  ceiling. The release workflow builds abi3 wheels for manylinux/musllinux
  (x86_64, aarch64), macOS (x86_64, arm64) and Windows; `xtalforge` gains
  `completions <shell>` and `man`; crates.io publication includes the agent
  crate.
- Phase 8 Task 8.6: regression benchmarks for symmetry, UFF, PXRD and the
  distance grid with `scripts/check-bench-budget.sh`, determinism tests across
  threads, a deterministic mutation harness for every text parser and the
  project archive, and a limits audit at and beyond every budget.
- Phase 8 Task 8.7: tutorials for eight end-to-end workflows, the
  Crystal-Builder migration guide, and consolidated algorithm pages for
  analysis, builders, the MOF builder, the powder pipeline and modules with a
  coverage index.
- Phase 8 Task 8.8: the API stability and SemVer page, a compile-time freeze
  test for representative signatures, and a closed parity matrix with every
  exception recorded.
- Phase 8 Task 8.4: the machine-learning bridge and model cache. Python
  `MlCalculator` wraps any ASE calculator (MACE, ORB-v3, MatterSim) and maps
  energy, forces and stress to eV, eV/Angstrom and GPa; `ml_backends()` probes
  the modules. The CLI `models list|fetch|verify|remove` and
  `xtalforge_ff::models` record URL, license and SHA-256 in a user cache and
  never bundle weights. MACE-MP-0 (small) golden data covers H2, water and
  methane within 1e-6 eV, 1e-5 eV/Angstrom and 1e-6 GPa, with timings.
- Phase 8 Task 8.3: open extension registries and a Python plugin surface.
  `FormatRegistry`, `ScalarRegistry` and `StyleRegistry` join the existing
  calculator and module registries; `examples/plugin` is a buildable crate
  that registers a format, a calculator, a module, a scalar and a style.
  Python `Plugins` registers callables as modules, calculators and formats,
  with JSON in/out and typed errors. The draw-style enum stays closed; the
  open part is the settings preset, recorded in `docs/api/plugins.md`.
- Task 0.9: bundle a 14-structure subset of the QMOF Database (CC BY 4.0,
  Figshare 13147324 v18) and validate the reader against the tabulated atom
  count, volume, density and formula; four P1 files recover their QMOF space
  group at `symprec 0.01 A`, and a CIF round-trip preserves the cell.
- Task 0.9: bundle an 18-complex subset of tmQM 2024 (MIT, one complex per
  metal from Sc to Ir) and validate the XYZ reader against the header
  stoichiometry, the perceived metal coordination and an XYZ round-trip. UFF
  builds, types and evaluates every metal with finite forces.
- Task 0.9: bundle an 18-structure subset of the CoRE MOF 2019 public archive
  (CC BY 4.0, Zenodo 14184621) and validate against it. Every file parses,
  expands to P1 and bonds; four frameworks are checked against their own
  `_chemical_formula_sum`, every declared `_cell_volume` is reproduced within
  0.05 A^3, CIF round-trips preserve sites and volume, and HKUST-1 and ZIF-8
  give bounded PXRD and porosity numbers.
- Phase 8 Task 8.2: `xtalforge mcp` serves the verb catalogue as MCP tools
  over stdio (JSON-RPC 2.0, protocol `2024-11-05`, no external MCP crate),
  with generated input schemas and an `instructions` payload for chat
  clients. Fourteen skill documents ship with the CLI (`xtalforge skill
  path|list|show|install`) and are published as the automation guide. Calls
  run under the standard size limits; cancellation is process-level and the
  window proxy is out of scope, both recorded in
  `docs/algorithms/agent.md` section 11.
- Phase 8 Task 8.1: the `xtalforge-agent` crate: one verb catalogue shared by
  the Rust, CLI, Python and WASM surfaces (22 verbs), a stateful session whose
  mutations use `Origin::Agent`, inspection and diagnostics (overlap,
  occupancy, charge, dangling and under-coordinated atoms, unusual
  coordination, symmetry), command-log change bursts, typed JSON answers and a
  pure-Rust orthographic renderer to SVG and PNG. CLI `capabilities`, `render`
  and `agent`; Python `AgentSession`, `agent_run`, `agent_capabilities`,
  `agent_verbs`; WASM `AgentSession`, `agent_run`, `agent_capabilities`,
  `agent_verbs`. The WASM size budget moves from 3648 to 3840 KiB to carry the
  new exported surface (3775 KiB).
- Phase 0 Task 0.9: external validation suite. EQeq is compared with the
  authors' program on the COD MOF-5 cell (424 charges within 1e-3 e), UFF with
  sci-form on methane and stretched water (within 1e-4 kcal/mol), and the
  DFTB+ parser with a recorded Si SCC output. `docs/algorithms/external-validation.md`
  records the bundled sources, their licences and the data sets still to add.
- Phase 5 Task 5.6: bundle the full RCSR 2019-06-01 snapshot in `data/rcsr/`
  and identify nets by coordination sequence and point symbol against all
  2,931 entries. `RcsrIndex::from_embedded()` returns `RcsrIdentification`
  with the unique name, an ambiguous set, or the nearest candidates; the
  2,909 keyed entries are validated against their canonical keys at index
  build time. The CGD reader now accepts both RCSR edge dialects.
- Phase 6 Task 6.8: the disordered-carbon recipe builder: convolution surface
  around a periodic net, split/collapse/flip remesh with tangential relaxation
  and projection, ribbon paths over the triangle graph, a three-coordinate
  dual and H/F/O edge termination. `CarbonBuild` reports density, ratios,
  ring census, mean C-C, pieces and periodicity; CLI `build carbon`, Python
  `build_carbon`, WASM `build_carbon_recipe_info`.
- Phase 7 Task 7.10.1: every leading-cell row now carries `space_groups`, the
  candidates whose general-position systematic absences every indexed
  reflection obeys. The screen enumerates the cell's crystal system, evaluates
  the structure factor of two generic points per reflection and ranks by
  group order; CLI, Python and WASM expose the candidates.
- Phase 7 Task 7.10.1: the automatic table now ranks extinction classes by a
  Le Bail BIC and Pawley-fits every leading (cell, class), as the reference
  does. A row carries `cell_rank`, `class_rank`, `class_symbol`, `delta_bic`
  and the groups that share the class; the 1 % Rwp tie rule orders the table.
- Phase 4 Task 4.9: validate the DFTB+ energy, force, stress and population
  parser against the recorded `si_scc_detailed.out` fixture. The xTB GFN
  methods remain available in-process through `sci-form` (Task 4.10).
- Phase 1 Task 1.12: `rings_explicit` and `ring_counts_explicit` census a
  builder's stored bonds without distance perception.

### Added

- Phase 4 Task 4.10: link the published MIT crate `sci-form` 2.0.0 behind the
  optional `sci-form` feature and implement `Calculator` for GFN0, GFN1 and
  GFN2-xTB. The methods are molecular only (a space-group structure is
  refused), GFN2 returns central-difference forces, and GFN0/GFN1 refuse
  `energy_forces`. The registry lists `gfn0`, `gfn1` and `gfn2`.

### Fixed

- Phase 5 Task 5.4: order the four crossings of a marching-tetrahedra 2-2 split around their quadrilateral. The fan taken in the generation order crossed itself, leaving non-manifold edges and a surface that was not watertight.
- Phase 7 Task 7.6.3: Le Bail extraction returns intensity esds from the
  profile normal matrix, plus Rwp and the reduced chi-squared; the CLI report
  includes them.
- Phase 7 Task 7.10.1: the automatic run fits a Pawley model of every leading
  indexing candidate and returns a table sorted by Pawley Rwp, with the cell
  esds; CLI, Python and WASM expose the table.
- Phase 7 Task 7.8: run the Rietveld refinement in McCusker order
  (analytic scale/background, zero shift and free cell numbers, width, then
  Wyckoff coordinates, isotropic ADPs, occupancies and March-Dollase texture)
  with a bounded direct search per stage, covariance-based parameter esds and
  the `stages` list. Coordinates move only along the Wyckoff null space
  (`free_directions`), the metric-family cell ties are restored after every
  move, and the copper K alpha doublet can be placed. A rutile golden pattern
  from Crystal-Builder validates the recovered position, cell and zero shift.
- Phase 7 Task 7.6.2: add the Pawley fit with one intensity per reflection,
  exact linear covariance and intensity esds, free cell numbers, zero shift and
  size/strain width model, R-factors and chi-squared; a Python-generated golden
  pattern validates the recovery.
- Phase 7 Task 7.10: the automatic run indexes the five metric families and
  adds a Pawley step of the indexed cell; Python `powder_pawley` and WASM
  `powder_pawley_info` expose the fit.
- Phase 7 Task 7.5.1/7.5.2: add general metric-family indexing (cubic,
  tetragonal, hexagonal, orthorhombic, monoclinic) with cell refinement, de
  Wolff M20, Smith-Snyder F_N and a coverage-ranked table; Python
  `powder_index`, WASM `powder_index_info` and `powder index --general`.
- Phase 7 Task 7.3.1: add the SNIP background envelope, per-line Gaussian plus
  linear background fitting and `fit_peaks`, with `--background` and
  `--refine-peaks` on the CLI.
- Phase 1 Task 1.15.4: map a non-P1 `_geom_bond` row from its site labels,
  operation code and image to P1 contacts (MIL53 matches the 36 golden
  contacts), collapsing symmetry-equivalent images by their minimum-image
  vector. Unmatched records annotate nothing instead of inventing an edge.
- Phase 1 Task 1.15.4: warn when a CIF operation loop's count disagrees with
  the named space group.
- Phase 3 Task 3.6.3: add LAMMPS text dump (element column or type map,
  orthogonal/triclinic, Cartesian/scaled coordinates) and DFTB+ `geo_end`
  readers to the trajectory dispatcher.
- Phase 1 Task 1.15.4: write non-P1 explicit bonds back as a validated
  P1-index CIF loop and compare the full symmetry-operation set with the named
  space group, not only its count.
- Phase 4 Task 4.7: choose overlapping disorder components from the declared
  CIF formula and cell charge, make the planner graph-aware, add hydrogen
  filling, the close-contact list, the prepared-structure undoable command
  path and their tests.
- Phase 4 Task 4.9.2/4.9.3: wire `DftbCalculator` and `XtbCalculator` to the
  runner, parse `detailed.out`/xTB stdout, validate force counts and support
  timeout, cancellation and scratch inspection.
- Phase 4 Task 4.9.1: add the scratch-folder process runner with wall-clock
  timeout, cancellation, capped stdout/stderr logs, relative input files and
  simulated-executable tests.
- Phase 4 Task 4.7: add duplicate merging, disorder-component choice, solvent
  removal, charge-balance reporting and the prepared-sample golden sanity
  pass; `PrepareReport` now returns the prepared structure.
- Phase 4 Task 4.8.3: add supercell force constants, shortest-image folding,
  the acoustic sum rule, a Monkhorst-Pack q-mesh, Gaussian phonon DOS and a
  phonopy spring-model golden comparison. An explicit `q_path` samples a band
  path.

### Changed

- The workspace version is `1.0.0-rc.1`; the facade now forwards `sci-form`
  to `xtalforge-ff` and `xtalforge-agent` instead of declaring a no-op.
- The agent crate joins the published set and the crates.io publication order
  (after the domain crates, before the facade).
- The WASM release budget is 3840 KiB (3767 KiB measured).
- Python wheels are built for manylinux/musllinux (x86_64, aarch64), macOS
  (x86_64, arm64) and Windows (x86_64), plus an sdist.

## [0.0.3] - 2026-11-09

### Changed

- Phase 1 Task 1.11.1/1.11.5: perceive bonds with the reference defaults
  (`scale = 1.15`, metal-metal pairs forbidden) and validate the full golden
  corpus: every available bond set and coordination number matches
  Crystal-Builder after P1 atom mapping.
- Phase 1 Task 1.6/1.15: expand each site on its own duplicate grid, so
  disorder components at the same position both survive, and resolve a CIF
  Hall/H-M symbol into operations when the file omits the operation loop.
- Phase 1 Task 1.12: replace the chordless-cycle scan with the periodic
  Franzblau walk, expose `Ring::path`/`ring_paths` with lattice images on all
  four surfaces, and validate counts and paths against the golden corpus.
- Phase 1 Task 1.15.4: round-trip Uij (matched by label), disorder
  group/assembly, the oxidation-number charge fallback and the explicit-bond
  loop; split the CIF writer and structure extraction under `cif/`.
- Phase 3 Task 3.6.2/3.6.3: add lazy `XyzIndex` frame access and the VASP
  XDATCAR reader with XYZ/XDATCAR detection in the CLI.

## [0.0.2] - 2026-11-09

### Changed

- Phase 1 Task 1.10.1/1.10.2: replace the direct periodic pair scan with a
  linked-cell index in fractional space and exact image enumeration, add the
  10^3/10^4/10^5-atom benchmark target and compare it with Crystal-Builder.
- Phase 1 Task 1.10.2: P1 expansion now deduplicates through a fractional
  grid, removing the quadratic scan that dominated large structures.
- Phase 1 Task 1.7.1: supercell boundary duplicates use the same fractional
  grid, cutting the 2x2x2 MOF-5 build from 37.5 ms to 5.5 ms.

### Fixed

- Rename the Linux desktop entry to the reverse-DNS component ID so
  `appimagetool` finds the AppStream metadata; the AppImage build no longer
  reports the missing-metadata warning, and the `.deb`/`.rpm` payloads carry
  the matching desktop ID.

## [0.0.1] - 2026-11-08

### Changed

- Refresh the documentation against the implemented Phase 1-7 surface: the
  README status, guide, API pages, algorithm pages, validation table, benchmark
  notes, third-party attribution and roadmap current-status/gap list.
- Record the remaining documented subsets (neighbour acceleration, phonon
  q-mesh, preparation steps, general indexing and figures of merit, powder
  profile simplifications) in `docs/roadmap/index.md`.

### Fixed

- Upgrade `zip` to 9.4.2, clearing RUSTSEC-2026-0367, RUSTSEC-2026-0368,
  RUSTSEC-2026-0369 and RUSTSEC-2026-0429.
- Fix the rustdoc intra-doc links that failed `cargo doc -D warnings`.
- Expose C20, twists, ribbons, ZTC/Schwarz P surfaces and polymers in the CLI,
  Python and WASM, closing the Phase 6 surface-parity gaps; add regression
  tests for each surface.

### Added

- Phase 6 Task 6.1/6.5/6.7/6.8/6.9: SMILES 3D embedding, metal coordination,
  the sketch model and SMILES writer, commensurate twists, C20/C60 face
  spirals, gyroid/Schwarz P carbon surfaces and polymer chain packing.
- Phase 7 Task 7.1-7.11: the rietx inventory and acceptance metrics, the
  instrument model, whole-range and structural least-squares refinement,
  free cell parameters and extinction symbols, the combined objective, the
  automatic powder run, and the CLI/Python/WASM surfaces.
- Phase 6 Task 6.11-6.16: the `xtalforge-mof` crate with six topologies,
  eight blocks, framework assembly, rods, layers, merged nets, decoration,
  invented linkers and pore-space partition.
- Phase 6 Task 6.17: `mof` CLI commands, Python/WASM catalogues and builder,
  and the MOF API page.
- Phase 4 Task 4.3.3: native EEQ-style charges.
- Phase 4 Task 4.6: hydrogen placement.
- Phase 4 Task 4.7: the preparation pipeline and report.
- Phase 4 Task 4.9: engine probes plus DFTB+, xTB and ORCA input writers and
  output parsers.
- Phase 4 Task 4.10: the sci-form adapter contract and capability limits.
- Phase 4 Task 4.12: `ff` CLI commands, Python/WASM force-field surfaces and
  the algorithm pages.
- Phase 3 Task 3.7: the `.xtalproj` ZIP container with independent
  `format_version`, readable parts and a newer-version gate.
- Phase 3 Task 3.8.2: `symmetry`, `pxrd` and `uff` parameter presets shared by
  the CLI, Python and WASM schemas.
- Phase 3 Task 3.10.7: `SetCharges`, `Interpenetrate` and
  `MarkConnectionPoints` commands.
- Phase 3 Task 3.12.3: in-memory workspace archives for WASM and on-disk
  workspace binding for Python.
- Phase 3 Task 3.16: `formats`, `params`, `project`, `trajectory` and `script`
  CLI commands; Python `Document`/project/workspace bindings; WASM
  `WasmDocument` plus project and workspace exports; the editor demo.
- Phase 4 Task 4.8: Gamma-point Hessian, vibrational modes, IR spectrum,
  Gamma-point phonons, harmonic/ideal-gas thermochemistry and the
  quasi-harmonic free-energy scan.
- Phase 4 Task 4.3-4.5: QEq/EQeq charges, per-site charges, Ewald energy and
  forces, FIRE/L-BFGS optimization (P1/asymmetric, optional cell), equations of
  state and relaxed scans.
- Phase 4 Task 4.1/4.2: Calculator API, registry, run records and the native
  UFF/UFF4MOF engine with analytic gradients.
- Phase 5 Task 5.1-5.3: calculated PXRD with radiation, reflections,
  structure factors, LP and profiles; periodic distance grids and porosity.
- Phase 6 Task 6.5/6.6: graphene sheets/stacks and nanotube roll-up builders.
- Phase 3: format registry with capabilities/loss warnings, atomic writes,
  pymatgen JSON/CSSR/GEN/LAMMPS/SHELX/PDB, CGD nets, Grid3/cube/XY,
  trajectories, typed parameters, commands/document with undo, selection and
  tracking, workspace/modules, view data, colormaps and ASCII STL.
- Phase 3 Task 3.1-3.3: the format registry with capabilities and loss
  warnings, atomic writes, and pymatgen JSON, CSSR and DFTB+ GEN readers and
  writers.

- Phase 2 Task 2.8: layer-group detection for slabs, experimental collinear
  magnetic groups and molecular point groups.
- Phase 2 Task 2.10.2-2.10.4: Bravais classification, ASE-derived special
  points and band paths in the user's reciprocal basis.
- Phase 2 Task 2.3.5: change of setting by Hall number.

- Project bootstrap: `AGENTS.md` rules and versioning policy, roadmap (phases 0
  to 8, verified against the Crystal-Builder, sci-form and moyo sources),
  rule-check script, and VitePress documentation scaffold.
- Cargo workspace (Task 0.1.1) with `xtalforge-core`, `xtalforge-symmetry`,
  `xtalforge-io`, the facade `xtalforge`, `xtalforge-cli`, `xtalforge-python`
  and `xtalforge-wasm`; shared versions for `nalgebra`, `serde`, `thiserror`,
  `rayon`, `moyo` and `flate2`.
- Toolchain, formatting and lint configuration (Task 0.1.2):
  `rust-toolchain.toml`, `rustfmt.toml` and the workspace lints for unsafe
  code, missing docs, warnings and panic-prone patterns.
- Repository hygiene (Task 0.1.3): `.editorconfig` and `deny.toml` with the
  license allow list, advisory checks, duplicate-version ban and native-build
  bans.
- Facade crate (Task 0.1.4) with stable `core`, `symmetry` and `io` module
  names, a prelude, and the `parallel`, `sci-form`, `ml`, `external` and
  `mof-data` features documented in the crate and in the guide.
- Contributor guide (Task 0.1.5) and the Cargo features guide.
- Rule-check self-test (Task 0.2.1): positive and negative fixtures and a
  runner covering the 500-line boundary and the ban on warning-suppression
  attributes.
- Optional Git hook (Task 0.2.2): `scripts/install-hooks.sh` and a pre-commit
  hook running the rules, formatting and clippy on touched crates.
- CI workflow for the repository rules and their self-test (Task 0.2.4),
  plus blocking formatting and clippy jobs with warnings denied.
- Data bundle loader (Task 0.2.3): `xtalforge_core::data` with embedded gzip
  bundles, SHA-256 verification of the decompressed payload, once-only
  decompression and JSON decoding through `OnceLock`, and lazy external
  bundles for WASM and large databases.
- Data policy page and the `data/` directory convention, including the gzip
  fixture and checksum tests for the loader.
- moyo audit (Task 0.3): requirements matrix, gap decisions, eight-structure
  spglib comparison at four tolerances, wasm/pyo3 build check, and the report
  in `docs/algorithms/symmetry-backend.md`.
- sci-form reuse audit (Task 0.4): verified inventory, port quality gate,
  overlap decisions and the `sci-form` feature policy in
  `docs/algorithms/sci-form-reuse.md`.
- Hello-world slice (Task 0.5): `version()` and `cell_volume()` in the Rust
  API, CLI, Python wheel and WASM package; shared vectors, API pages, parity
  check, binding architecture and error mapping.
- Version and release tooling (Task 0.6): `scripts/sync-version.sh` with a
  self-test, agent-version check, pull request template and release workflow.
- Documentation site (Task 0.7): edit links, search, repository link, API and
  algorithm pages, four-surface code groups, KaTeX, and an offline cell-volume
  WASM demo.
- CI workflows (Task 0.8): Rust matrix, Python, WASM, docs, caching and path
  filters, plus the existing rules job.
- Golden-data oracle (Task 0.9): curated `tests/corpus/` with provenance,
  Phase 1 generator and data, Crystal-Builder test-data inventory,
  `xtalforge-testkit` and the regeneration policy.
- Dependency replacement map (Task 0.10): import survey, replacement
  decisions and license status for every Python dependency.
- Mandatory test layout (AGENTS.md section 2.2): Rust tests live in
  `crates/<crate>/tests/`, independent of the source files.
- Phase 1 Task 1.1: periodic table, deuterium and dummy species, Pyykko radii,
  symbol parsing, `element_info` in Rust, CLI, Python and WASM, the generated
  API table and the shared element vectors.
- Phase 1 Task 1.2: `Lattice` with matrix rows, standard orientation,
  parameters, volume, metric tensor, reciprocal cells, d-spacing and plane
  normals, coordinate conversions, wrapping, exact minimum-image search,
  transforms, strain, volume scaling and validation; `cell --params` in the
  CLI and `lattice_info`/`lattice_minimum_image` in Python and WASM.
- Phase 1 Task 1.3: the limits model with Standard/Generous/WarnOnly
  profiles, soft and hard verdicts, drawn/supercell/grid/phonon/Hessian
  estimates, `Error::LimitExceeded`, and `limits_bounds`/`limits_check` in
  Rust, Python and WASM.
- Phase 1 Task 1.4: `Site` (label, canonical element, occupancy, Uiso/Uij,
  charge, disorder group, props), `Structure`/`StructureBuilder`, molecule in
  a padded box, formula/mass/density/charge properties, the `info` CLI bridge
  and `structure_info` in Python and WASM.
- Phase 1 Task 1.5: `StructureDocument` with `SCHEMA_VERSION`, strict
  top-level fields, version checks and a versioned fixture.
- Phase 1 Task 1.6: `SymOp` parsing/printing/composition/inversion and
  `P1Cell`/`P1Atom` expansion with site/operation/lattice provenance and
  special-position deduplication.
- Phase 1 Task 1.7: diagonal and general supercells as P1, origin shifts,
  metadata/site-property preservation and the supercell limit check.
- Phase 1 Task 1.8: surface bases for arbitrary Miller planes, slab
  construction with layers, vacuum and P1 output, and plane angles.
- Phase 1 Task 1.9: finite-cell images, Cartesian point transforms, plane
  fitting/alignment and displacement ellipsoids.
- Phase 1 Task 1.10: periodic neighbor search with images, distances and
  deterministic ordering; `neighbors` and `neighbors_of`.
- Phase 1 Task 1.11.1/1.11.4: distance-based bond perception with covalent
  radii, dummy exclusion, bond images and the connectivity graph.

### Changed

- License is now `MIT OR Apache-2.0` (confirmed by the maintainer,
  Task 0.1.3), replacing the initial MIT-only `LICENSE` with `LICENSE-MIT` and
  `LICENSE-APACHE`.
- The English-only rule is retired (maintainer decision, 2026-11-08). The
  language checker, allow list, Spanish fixtures and commit-msg hook were
  removed; no language check runs in CI.
- Workspace lints deny warnings and missing docs. `scripts/check-rules.sh`
  rejects `allow`/`expect` attributes for dead code, unused lints and warnings.
- Commit discipline adopted (AGENTS.md section 2.1): one directory per commit,
  `feat(<scope>):` or `ref(<scope>):` subjects, and the completed tasks listed
  from hardest to easiest in the body. The optional pre-commit and commit-msg
  hooks enforce the directory scope and the message format.
- Workspace `unsafe_code` changed from `forbid` to `deny` so the pyo3 binding
  crate can relax it locally with a comment, as the roadmap allows.
- GitHub Actions workflows are committed but disabled at repository level by
  maintainer instruction; they have not run yet.
