# Surface parity matrix

Every public capability must exist in all four surfaces, or have a recorded and
justified exception. Status values: `todo`, `wip`, `done`, `n/a` (with the
reason in the exceptions table). `scripts/check-parity.sh` (Task 0.5.6)
compares this table with the symbols each surface actually exports.

| Capability | Phase | Rust | CLI | Python | WASM | Notes |
|---|---|---|---|---|---|---|
| `version` | 0 | done | done | done | done | Hello-world slice (Task 0.5.7) |
| `cell_volume` | 0 | done | done | done | done | Hello-world slice (Task 0.5.7) |
| `element_info` | 1 | done | done | done | done | Element data (Task 1.1.4) |
| `lattice_info` | 1 | done | done | done | done | Lattice (Task 1.2.5) |
| `lattice_minimum_image` | 1 | done | n/a | done | done | Rust, Python and WASM; the CLI reports geometry through `cell` |
| `limits_bounds` | 1 | done | n/a | done | done | Limits (Task 1.3.3); the CLI has the global `--limits` flag |
| `limits_check` | 1 | done | n/a | done | done | Limits (Task 1.3.3); the CLI has the global `--limits` flag |
| `structure_info` | 1 | done | done | done | done | Structure summary (Task 1.4.5); reads any registered format |
| `read_cif` | 1 | done | done | done | done | CIF reader (Task 1.15); CLI `info` and `convert` |
| `read_structure` | 1 | done | done | done | done | Format dispatch (Task 1.16.3); empty hint auto-detects |
| `supercell` | 1 | done | done | done | done | Supercells (Task 1.7); WASM `supercell_info` |
| `slab` | 1 | done | done | done | done | Slabs (Task 1.8); WASM `slab_info` |
| `bonds` | 1 | done | done | done | done | Bond graph (Task 1.11) |
| `rings` | 1 | done | done | done | done | Primitive rings (Task 1.12) |
| `ring_paths` | 1 | done | done | done | done | Periodic ring paths with lattice images (Task 1.12) |
| `groups` | 1 | done | done | done | done | Functional groups (Task 1.13) |
| `neighbors` | 1 | done | done | done | done | Periodic neighbors (Task 1.10) |
| `distance` | 1 | done | done | done | done | Measures (Task 1.14) |
| `angle` | 1 | done | done | done | done | Measures (Task 1.14) |
| `torsion` | 1 | done | done | done | done | Measures (Task 1.14) |
| `detect_symmetry` | 2 | done | done | done | done | Symmetry search (Task 2.1.4); CLI `symmetry detect` |
| `asymmetrize` | 2 | done | done | done | done | Workflow (Task 2.3.3); WASM `asymmetrize_info` |
| `to_p1` | 2 | done | done | done | done | Workflow (Task 2.3.4); WASM `to_p1_info` |
| `set_space_group` | 2 | done | done | done | done | Workflow (Task 2.3.2); WASM `set_space_group_info` |
| `subgroups` | 2 | done | done | done | done | Enumeration (Task 2.7.1/2.7.4); WASM `subgroups_info` |
| `lower_to_subgroup` | 2 | done | done | done | done | Lowering (Task 2.7.3/2.7.4); WASM `lower_to_subgroup_info` |
| `standardize` | 2 | done | done | done | done | Task 2.5.1 |
| `primitive_cell` | 2 | done | done | done | done | Task 2.5.2 |
| `reduce_structure` | 2 | done | done | done | done | Task 2.6.3 |
| `wyckoff_assignments` | 2 | done | done | done | done | Task 2.4.2 |
| `tolerance_ladder` | 2 | done | done | done | done | Task 2.9.1 |
| `diagnose` | 2 | done | done | done | done | Task 2.9.2 |
| `brillouin_zone` | 2 | done | done | done | done | Task 2.10.1 |
| `detect_layer_symmetry` | 2 | done | done | done | done | Layer groups (Task 2.8.1); WASM `layer_group_info` |
| `detect_magnetic_symmetry` | 2 | done | done | done | done | Magnetic groups (Task 2.8.2, experimental); WASM `magnetic_group_info` |
| `change_setting` | 2 | done | done | done | done | Setting change (Task 2.3.5); WASM `change_setting_info` |
| `band_path` | 2 | done | done | done | done | Bravais/k-path (Task 2.10.2-2.10.4); WASM `kpath_info` |
| `molecular_point_group` | 2 | done | done | done | done | Point groups (Task 2.8.3); WASM `pointgroup_info` |
| `scene_summary` | 3 | done | done | done | done | View data (Task 3.14.4); CLI `mesh` |
| `scene_data` | 3 | done | done | done | done | Full scene geometry for a viewer (Task 3.14.4); CLI `mesh --scene` |
| `mesh_stl` | 3 | done | done | done | done | ASCII STL (Task 3.15.1); CLI `mesh` |
| `pxrd` | 5 | done | done | done | done | Calculated PXRD (Task 5.1); Rust `calculate_pxrd`, WASM `pxrd_info` |
| `param_schema` | 3 | done | done | done | done | Parameter presets (Task 3.8.2); Rust `preset`, CLI `params` |
| `project_pack` | 3 | done | done | done | done | Project archives (Task 3.7); Rust `write_project`, CLI `project` |
| `workspace_archive` | 3 | done | n/a | done | done | In-memory workspace (Task 3.12.3); the CLI reads `run.json` directly |
| `prepare` | 4 | done | done | done | done | Preparation pipeline (Task 4.7); Python `ff_prepare`, WASM `ff_prepare_info` |
| `add_hydrogens` | 4 | done | done | done | done | Hydrogen placement (Task 4.6) |
| `porosity` | 5 | done | done | done | done | Distance grid (Task 5.2/5.3); Rust `distance_grid`; 17.1.2 adds POV, ASA per mass and percolating directions |
| `net_invariants` | 5 | done | done | done | done | Quotient graphs and the full RCSR fingerprint identification (Task 5.6/5.7); Rust `quotient_graph` and `RcsrIndex` |
| `zeo_parse` | 5 | done | done | done | done | Zeo++ output parser (Task 5.5); Rust `parse_zeo_channels` |
| `build_graphene` | 6 | done | done | done | done | Carbon builders (Task 6.5); Rust `graphene`, WASM `build_graphene_info` |
| `build_nanotube` | 6 | done | done | done | done | Carbon builders (Task 6.6); Rust `nanotube`, WASM `build_nanotube_info` |
| `build_c60` | 6 | done | done | done | done | C60 (Task 6.7); CLI `build c60`, WASM `build_c60_info` |
| `build_c20` | 6 | done | done | done | done | C20 (Task 6.7); Rust `face_spiral(20)`, CLI `build c20` |
| `build_mof` | 6 | done | done | done | done | MOF assembly (Task 6.12-6.16); WASM `build_mof_info` |
| `twist` | 6 | done | done | done | done | Commensurate twist (Task 6.5.3); CLI `build twist` |
| `ribbon` | 6 | done | done | done | done | Nanoribbon (Task 6.8.3); CLI `build ribbon` |
| `build_ztc` | 6 | done | done | done | done | Gyroid carbon model (Task 6.8); CLI `build ztc` |
| `build_schwarzite` | 6 | done | done | done | done | Schwarz P carbon model (Task 6.8); CLI `build schwarzite` |
| `build_carbon` | 6 | done | done | done | done | Disordered carbon recipe with convolution surface, remesh, ribbons and termination (Task 6.8); CLI `build carbon`, Python `build_carbon`, WASM `build_carbon_recipe_info` |
| `build_polymer` | 6 | done | done | done | done | Polymer chain (Task 6.9); CLI `build polymer` |
| `run_automatic` | 7 | done | done | done | done | Automatic run with the Le Bail-BIC extinction-class ranking and the Pawley-ranked (cell, class) table (Task 7.10); Python `powder_automatic` |
| `refine_phase` | 7 | done | done | done | done | McCusker-order structural refinement with esds and stages (Task 7.8); Python `powder_refine` |
| `powder_index` | 7 | done | done | done | done | General metric-family indexing with M20/F_N (Task 7.5); Python `powder_index`, WASM `powder_index_info`; 17.1.7 adds the cubic/general agreement verdict |
| `powder_lebail` | 7 | done | done | done | done | Le Bail intensity extraction with a SNIP background (Task 7.6.1; Python/WASM added by the 2026-11-11 audit); Rust `le_bail`, CLI `powder lebail`, WASM `powder_lebail_info` |
| `powder_pawley` | 7 | done | done | done | done | Pawley fit with esds and free cell/instrument parameters (Task 7.6.2); Python `powder_pawley`, WASM `powder_pawley_info` and `powder_pawley_width_info` (size/strain refinement, added by the 2026-11-11 audit) |
| `agent_capabilities` | 8 | done | done | done | done | Runtime feature, engine, bundle and verb report (Task 8.1.1); Rust `capability_report`, CLI `capabilities` |
| `agent_verbs` | 8 | done | done | done | done | The shared verb catalogue (Task 8.1.1); Rust `catalogue`, CLI `capabilities` |
| `agent_run` | 8 | done | done | done | done | Batch session over a JSON script (Task 8.1.2); Rust `run_verbs`, CLI `agent` |
| `ml_backends` | 8 | n/a | done | done | n/a | The ML bridge is Python-only by the Task 4.11 decision; the CLI `models` command lists the catalogue (Task 8.4) |
| `models` | 8 | done | done | n/a | n/a | Model cache with SHA-256 and license (Task 8.4.2); Rust `fetch`, CLI `models` |
| `match_structures` | 9 | done | done | done | done | Structure matcher (Task 9.1.1); Rust `match_structures`, CLI `compare`; 17.3.3 standardizes the primitive cell by default |
| `periodic_rmsd` | 9 | done | done | done | done | Periodic RMSD under a mapping (Task 9.1.2) |
| `kabsch_alignment` | 9 | done | done | done | done | Kabsch alignment (Task 9.1.2) |
| `group_structures` | 9 | done | done | done | done | Equality classes (Task 9.1.3); CLI `group` |
| `deform` | 9 | done | done | done | done | Strain, gradient, volume and density (Task 9.2.1); CLI `deform` |
| `rattle` | 9 | done | done | done | done | Seeded displacements (Task 9.2.2) |
| `edit_structure` | 9 | done | done | done | done | Site edits and species maps (Task 9.2.3); Rust `replace_species` and the site helpers, Python `edit_structure` |
| `interpolate` | 9 | done | done | done | done | Matched interpolation (Task 9.2.4); Python `interpolate_structures` |
| `optimal_supercell` | 9 | done | done | done | done | Cube-like supercell (Task 9.2.5) |
| `voronoi` | 9 | done | done | done | done | Periodic Voronoi cells (Task 9.3.1) |
| `coordination` | 9 | done | done | done | done | Coordination shells and shape measures (Task 9.3.2) |
| `bond_valence_sums` | 9 | done | done | done | done | Bond-valence sums (Task 9.3.3) |
| `oxidation_states` | 9 | done | done | done | done | Oxidation-state estimation (Task 9.3.3) |
| `radial_distribution` | 9 | done | done | done | done | Total and partial RDF (Task 9.3.4) |
| `structure_fingerprint` | 9 | done | done | done | done | Composition, cell and RDF fingerprint (Task 9.3.4); Rust `fingerprint`, Python `structure_fingerprint` |
| `unique_miller_indices` | 9 | done | done | done | done | Symmetry-unique Miller indices (Task 9.4.1) |
| `slab_terminations` | 9 | done | done | done | done | All slab terminations (Task 9.4.2); 17.3.2 adds intra-layer cuts and lower-layer adsorption sites |
| `wulff` | 9 | done | done | done | done | Wulff construction (Task 9.4.3) |
| `adsorption_sites` | 9 | done | done | done | done | Top, bridge and hollow sites (Task 9.4.4) |
| `vacancies` | 9 | done | done | done | done | Symmetry-unique vacancies (Task 9.5.1) |
| `substitutions` | 9 | done | done | done | done | Symmetry-unique substitutions (Task 9.5.1) |
| `interstitials` | 9 | done | done | done | done | Interstitial candidates (Task 9.5.1) |
| `orderings` | 9 | done | done | done | done | Ordered decorations ranked by Ewald energy (Task 9.5.2) |
| `special_quasirandom` | 9 | done | done | done | done | SQS generation with pair and triplet correlations (Task 9.5.3) |
| `lattice_matches` | 9 | done | done | done | done | Zur-McGill 2D lattice matching (Task 9.5.4) |
| `stack_slabs` | 9 | done | done | done | done | Heterostructure stacking (Task 9.5.4) |
| `elastic_tensor` | 9 | done | done | done | done | Stress finite differences, Voigt-Reuss-Hill moduli (Task 9.6.1) |
| `eos_fit` | 9 | done | done | done | done | Birch-Murnaghan, Vinet and Murnaghan fits (Task 9.6.2); Rust `fit_eos`, CLI `features eos` |
| `convex_hull` | 9 | done | done | done | done | Formation energies and hull distances (Task 9.6.3) |
| `mof_check` | 10 | done | done | done | done | MOF quality checker (Task 10.1.1); CLI `mof check` |
| `mof_solvent_remove` | 10 | done | done | done | done | Solvent/counter-ion removal (Task 10.1.2); CLI `mof solvent-remove` |
| `mof_disorder` | 10 | done | done | done | done | Disorder resolution (Task 10.1.3); CLI `mof disorder-resolve` |
| `mof_charge_balance` | 10 | done | done | done | done | Metal states and charge balance (Task 10.1.4); CLI `mof charge-balance` |
| `mof_deconstruct` | 10 | done | done | done | done | Building-unit deconstruction and contracted net (Task 10.2.1-10.2.2); CLI `mof deconstruct` |
| `mof_mofid` | 10 | done | done | done | done | MOFid and MOFkey (Task 10.2.3); CLI `mof mofid` |
| `mof_open_metal_sites` | 10 | done | done | done | done | Open metal sites (Task 10.2.4); CLI `mof open-metal-sites` |
| `mof_descriptors` | 10 | done | done | done | done | Linker/node/framework descriptors (Task 10.2.5); CLI `mof descriptors` |
| `mof_pores` | 10 | done | done | done | done | Native PLD, LCD, ASA, AV and channel dimension (Task 10.3.1-10.3.3); CLI `mof pores` |
| `mof_guests` | 10 | done | done | done | done | TraPPE guest table (Task 10.4.1); CLI `mof adsorption guests` |
| `mof_forcefield` | 10 | done | done | done | done | Host UFF/EQeq and guest mixing (Task 10.4.1); CLI `mof adsorption force-field` |
| `mof_widom` | 10 | done | done | done | done | Widom insertion and Henry coefficient (Task 10.4.2); CLI `mof adsorption widom` |
| `mof_gcmc` | 10 | done | done | done | done | Single-component GCMC isotherm (Task 10.4.3); CLI `mof adsorption gcmc` |
| `mof_energy_grid` | 10 | done | done | done | done | Guest energy grid (Task 10.4.4); CLI `mof adsorption grid` |
| `mof_raspa` | 10 | done | done | done | done | RASPA input export (Task 10.4.5); CLI `mof adsorption raspa` |
| `mof_functionalize` | 10 | done | done | done | done | Post-synthetic functionalisation (Task 10.5.3); CLI `mof modify functionalize` |
| `mof_defects` | 10 | done | done | done | done | Missing linker/cluster defects (Task 10.5.4); CLI `mof modify defects` |
| `mof_interpenetrate` | 10 | done | done | done | done | Interpenetrated variants (Task 10.5.5); CLI `mof modify interpenetrate` |
| `mof_stack_layers` | 10 | done | done | done | done | Layered variants (Task 10.5.5); CLI `mof modify stack-layers` |
| `mof_mixed` | 10 | done | done | done | done | Mixed-metal/linker frameworks (Task 10.5.6); CLI `mof modify mixed` |
| `mof_read_database` | 10 | done | done | done | n/a | CoRE MOF/QMOF directory readers (Task 10.6.1); no filesystem in WASM; 17.1.4 adds `fetch_core_mof`/`fetch_qmof` with cache, checksum and offline mode |
| `mof_screen` | 10 | done | done | done | n/a | Screening pipeline (Task 10.6.2); no filesystem in WASM |
| `mof_deduplicate` | 10 | done | done | done | n/a | MOFkey + matcher deduplication (Task 10.6.3); no filesystem in the browser |
| `read` | 11 | done | done | done | n/a | One-call IO (Task 11.2.3); WASM has no filesystem, use `read(text, hint)` |
| `write` | 11 | done | done | done | n/a | Atomic writer (Task 11.2.3); no filesystem in WASM |
| `set_threads` | 11 | done | done | done | n/a | Process-wide worker pool (Task 11.1.2); WASM threads need cross-origin isolation |
| `error_category` | 11 | done | done | done | done | Error categories in JSON objects (Task 11.1.3) |
| `deprecation_warning` | 11 | done | n/a | done | done | Deprecation mechanism (Task 11.4.2); the CLI uses `warn_deprecated` |
| `molecule_from_smiles` | 11 | done | done | done | done | SMILES parse and embed (Task 11.3.1); CLI `molecule` |
| `sketch_roundtrip` | 11 | done | done | done | done | Sketch round trip (Task 11.3.2); CLI `sketch` |
| `fragment_library` | 11 | done | done | done | done | Fragment lookup (Task 11.3.3); CLI `fragment`, WASM `fragment_library_info` |
| `fragment_structure` | 11 | done | done | done | done | Embed a fragment (Task 11.3.3); CLI `fragment NAME` |
| `substitute_structure` | 11 | done | done | done | done | Substitution with clearance (Task 11.3.4); CLI `substitute` |
| `place_guest_structure` | 11 | done | done | done | done | Guest placement with clearance (Task 11.3.5); CLI `guest` |
| `structure_handle` | 12 | done | n/a | done | done | Parsed structure handle (Task 12.1.1/12.3.1); CLI has no long-lived handle |
| `structure_from_arrays` | 12 | n/a | n/a | done | done | NumPy/typed-array input (Task 12.2.2); WASM `setPositions` |
| `typed_results` | 12 | n/a | n/a | done | done | Typed result classes/interfaces (Task 12.1.2/12.3.2); legacy text API keeps loose types |
| `typed_arrays` | 12 | n/a | n/a | done | done | Positions, cell, scene, PXRD and mesh buffers (Task 12.2.1/12.3.3); Python lists/NumPy |
| `batch_api` | 12 | n/a | n/a | done | n/a | Python batch helpers (Task 12.2.3); the browser has no filesystem |
| `wasm_worker` | 12 | n/a | n/a | n/a | done | Optional Web Worker wrapper (Task 12.3.4) |
| `node_entry` | 12 | n/a | n/a | n/a | done | Node file helpers (Task 12.3.4) |
| `ase_adapter` | 13 | n/a | n/a | done | n/a | ASE Atoms round trip and UFF calculator (Task 13.1.1); Python-only optional package |
| `pymatgen_adapter` | 13 | n/a | n/a | done | n/a | pymatgen Structure/Molecule round trip (Task 13.1.2); Python-only optional package |
| `phonopy_adapter` | 13 | n/a | n/a | done | n/a | PhonopyAtoms export (Task 13.1.3); Python-only optional package |
| `rdkit_adapter` | 13 | n/a | n/a | done | n/a | RDKit molecule exchange (Task 13.1.4); Python-only optional package |
| `get_symmetry_dataset` | 13 | n/a | n/a | done | n/a | spglib-shaped symmetry dict (Task 13.1.3); Python-only adapter |
| `interop_formats` | 13 | done | done | done | done | QE, CASTEP, FHI-aims, CP2K, ABINIT, CRYSTAL, mmCIF/PDBx, PDB, OUTCAR, vasprun.xml, MOL2, SDF/MOL and extxyz (Task 13.2) |
| `optimade_fetch` | 13 | done | done | n/a | done | OPTIMADE query, parse and export (Tasks 13.3.1/13.3.3); CLI `fetch optimade`, WASM `fetch_optimade` |
| `cod_fetch` | 13 | done | done | n/a | done | COD CIF fetch (Task 13.3.2); CLI `fetch cod`, WASM `fetch_cod` |
| `materials_project_fetch` | 13 | done | done | n/a | n/a | Materials Project fetch (Task 13.3.2); API-key handling is host-side |
| `python_plugins` | 13 | n/a | n/a | done | n/a | `xtalforge.plugins` entry points (Task 13.4.1); Rust registries cover the CLI and agent |
| `diffraction_profiles` | 14 | done | done | n/a | n/a | Caglioti/pseudo-Voigt/FCJ profiles and single-crystal factors (Task 14.1); native option structs and file writers |
| `md` | 14 | done | n/a | n/a | n/a | NVE/Langevin/Nose-Hoover/NPT molecular dynamics (Task 14.2.5); host-side trajectory and seed control |
| `net_contraction` | 14 | done | n/a | n/a | n/a | Structure-to-net contraction, vertex symbols and RCSR search (Task 14.4.2); native graph API |
| `mesh_export` | 14 | done | done | done | done | OBJ and glTF/GLB with normals (Task 14.4.4) |
| `report_renderers` | 14 | done | n/a | n/a | n/a | HTML and JSON report renderers (Task 14.5.1); module listing is the CLI surface |
| `workflow` | 16 | done | done | done | done | Versioned workflow runner and `xtalforge run` (Task 16.1); 17.2.1 adds YAML and 17.2.5 the in-memory WASM runner |
| `provenance` | 16 | done | done | done | n/a | Provenance record, embedding and `provenance` command (Task 16.2); WASM has no file embedding |
| `dataset` | 16 | done | n/a | done | n/a | JSON Lines dataset and pandas/polars view (Task 16.3.1); CLI stays file-oriented |
| `dataset_dedup` | 16 | done | n/a | n/a | n/a | Matcher + fingerprint prefilter (Task 16.3.2); Rust analysis API |
| `ml_descriptors` | 16 | done | n/a | n/a | n/a | Composition, RDF, Coulomb and SOAP-like arrays (Task 16.3.3); Rust analysis API |
| `porosity_report` | 17 | done | done | done | done | POV, ASA per mass, percolating directions and PLD from one grid (Task 17.1.2); CLI/Python/WASM extend the existing `porosity` output |
| `indexing_agreement` | 17 | done | done | done | done | Cubic/general engine verdict (Task 17.1.7); CLI `powder index --general`, Python `powder_index_agreement`, WASM `powder_index_agreement_info` |
| `database_fetch` | 17 | done | done | done | done | On-demand CoRE MOF/QMOF download with cache, SHA-256 and offline mode (Task 17.1.4); WASM `fetch_data_bundle` leaves storage to the page |
| `mof_data` | 17 | done | done | n/a | n/a | Full PORMAKE 0.2.3 bundle behind `mof-data` (Task 17.1.3); Python/WASM leave it off for package size |
| `engine_inputs` | 18 | done | done | done | done | 20 versioned input-set presets (Task 18.1); CLI `engines inputs`, Python `engine_input_set`, WASM `engine_input_set_info` |
| `engine_outputs` | 18 | done | done | done | done | Typed result model and VASP/QE/CP2K/LAMMPS/GULP parsers (Task 18.2); Bader/DDEC readers |
| `engine_electronic` | 18 | done | done | done | done | Bands, gaps, effective masses and DOS (Task 18.3); CLI `engines bands/dos` |
| `engine_jobs` | 18 | done | done | done | done | Scheduler scripts and failure patterns (Task 18.4); process execution is a host-side exception |
| `magnetic_moments` | 19 | done | done | done | done | Collinear and vector moments on sites, preserved by P1/supercell/slab (Task 19.1.1); mCIF and the JSON schema |
| `mcif` | 19 | done | done | done | done | magCIF reader/writer with magnetic operations and moments (Task 19.1.3) |
| `magnetic_orderings` | 19 | done | done | done | done | Ferro/antiferro/ferri enumeration ranked by magnetic symmetry (Task 19.1.4) |
| `disorder` | 19 | done | done | done | done | Partial-occupancy groups, composition, density and scattering weights (Task 19.2.1) |
| `derivative_structures` | 19 | done | done | done | done | HNF supercells and symmetry-distinct decorations (Task 19.2.2) |
| `neutron_powder` | 19 | done | done | done | done | Nuclear neutron pattern, constant wavelength and TOF (Task 19.3.1) |
| `magnetic_powder` | 19 | done | done | done | done | Magnetic neutron structure factors from site moments (Task 19.3.2) |
| `electron_diffraction` | 19 | done | done | done | done | Kinematic SAED with an image buffer (Task 19.3.3) |
| `raman` | 19 | done | done | done | done | Powder-averaged Raman from engine tensors (Task 19.4.2) |
| `phonon_dos_neutron` | 19 | done | done | done | done | Phonon DOS and inelastic-neutron weights (Task 19.4.3) |

## Recorded exceptions

| Capability | Surface | Reason | Alternative offered |
|---|---|---|---|
| Reading and writing paths, atomic writes | WASM | No filesystem in the browser | String and `Uint8Array` APIs |
| On-disk workspace of runs | WASM | No filesystem | In-memory workspace exported as a zip |
| External engines (DFTB+, xtb, Zeo++) and process probes | WASM | No subprocesses in the browser | Native engines (UFF, QEq, EQeq, grid porosity) |
| In-process GFN0/GFN1/GFN2 (`sci-form` feature) | CLI, Python, WASM | Optional build-time feature whose methods are molecular only; the Rust `Calculator` and the registry listing are the supported surface, and no binding negotiation for the feature exists yet | Rust `quantum::adapter`; CLI `ff engines` lists the registry entries |
| Multi-threading | WASM | Threads need cross-origin isolation | Optional threaded build |
| Process execution and scheduler probes | WASM, Python | No subprocesses in the browser, and the Python package does not drive a scheduler; script writing and parsing stay available | CLI `engines script`; WASM/Python script text |
| Python progress callbacks | Python | The profile refinement entry points that accept `RunControl` are not bound as separate Python functions yet; the Python surface exposes the results and cancellation is available in Rust | Rust `RunControl`/`Cancellation`; Python `Structure` context manager |
| MCP server over stdio | Python, WASM | It is a process, served by the CLI | `xtalforge mcp` |
| Machine-learning potentials (`MlCalculator`, `ml_backends`) | Rust, CLI, WASM | The Task 4.11 decision is the Python/ASE bridge; torch is not a WASM target and the Rust `ml` feature is a declared no-op | Python `MlCalculator`; the CLI `models` command manages weights |
| Plugins (`Plugins`, runtime registries) | CLI, WASM | Python callables and Rust registries are host-side extension points; a browser has no dynamic library or Python runtime | Rust `FormatRegistry`, `ScalarRegistry`, `StyleRegistry`, `CalculatorRegistry`, `ModuleRegistry`; Python `Plugins` |
| Model cache (`models`) | Python, WASM | It is a filesystem command with `curl`; the browser and Python package leave the download to the user | CLI `models`; `xtalforge_ff::models` |
| Full MOF topology bundle | WASM | Package size | The builder always carries the curated six-topology/eight-block subset and the RCSR fingerprint index; the full PORMAKE bundle is a Rust/CLI feature (`mof-data`) and stays out of the browser/Python packages for size |
| lattice_minimum_image, limits_bounds, limits_check | CLI | The CLI exposes geometry through `cell` and the global `--limits` flag | Use `cell` and `--limits` |
| workspace_archive | CLI | The CLI has no `modules`/`run` command; the workspace is a library API | Rust `Workspace` and Python `workspace_*`; `project` for archives |
| ASE, pymatgen, phonopy and RDKit adapters | Rust, CLI, WASM | They wrap optional Python objects and optional third-party packages | Rust readers/writers and the typed structure model |
| Materials Project fetch | WASM | The API key must not be embedded in a browser bundle | COD and OPTIMADE fetch from the browser; `xtalforge fetch mp` from the CLI |
| Python entry-point plugins | Rust, CLI, WASM | Entry points are discovered by the Python interpreter; a browser or the Rust binary cannot import them | Rust plugin registries are visible to the CLI and the agent; Python plugins stay in-process |
| Diffraction profiles, MD, NET contraction and report renderers | Python, WASM | The option structs carry host callbacks and trajectory/file controls; 1.7 exposes the stable serializable parts and keeps the full model native | CLI for profiles; the public Rust API for the rest |
| Workflow runner | WASM | Workflows read filesystem inputs and write checkpoint/output files; the browser has no filesystem | `workflow_schema`/`workflow_validate` style metadata can be added later |
| Provenance embedding | WASM | CIF/extxyz/project files are written by the host | The text/JSON provenance record is available in Rust and Python |
| Dataset files | CLI, WASM | Dataset files are a host filesystem container; CLI analysis commands stay per-structure | Rust `Dataset` and the Python pandas/polars view |

## Open gaps: Rust-only capabilities

None. Every capability that had a Rust-only row in 1.9 now has a CLI, Python and
WASM surface or a recorded exception in the table above. The former
SMILES/sketch/fragment/substitution/guest gaps were closed in Phase 11 and the
aromatic hydrogen inference is the opt-in `hydrogens`/`infer_hydrogens` flag
(Task 17.1.1); the RCSR and PORMAKE bundles are data, not an operation.
