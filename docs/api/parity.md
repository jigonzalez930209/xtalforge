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
| `porosity` | 5 | done | done | done | done | Distance grid (Task 5.2/5.3); Rust `distance_grid` |
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
| `powder_index` | 7 | done | done | done | done | General metric-family indexing with M20/F_N (Task 7.5); Python `powder_index`, WASM `powder_index_info` |
| `powder_pawley` | 7 | done | done | done | done | Pawley fit with esds and free cell/instrument parameters (Task 7.6.2); Python `powder_pawley`, WASM `powder_pawley_info` |
| `agent_capabilities` | 8 | done | done | done | done | Runtime feature, engine, bundle and verb report (Task 8.1.1); Rust `capability_report`, CLI `capabilities` |
| `agent_verbs` | 8 | done | done | done | done | The shared verb catalogue (Task 8.1.1); Rust `catalogue`, CLI `capabilities` |
| `agent_run` | 8 | done | done | done | done | Batch session over a JSON script (Task 8.1.2); Rust `run_verbs`, CLI `agent` |
| `ml_backends` | 8 | n/a | done | done | n/a | The ML bridge is Python-only by the Task 4.11 decision; the CLI `models` command lists the catalogue (Task 8.4) |
| `models` | 8 | done | done | n/a | n/a | Model cache with SHA-256 and license (Task 8.4.2); Rust `fetch`, CLI `models` |

## Recorded exceptions

| Capability | Surface | Reason | Alternative offered |
|---|---|---|---|
| Reading and writing paths, atomic writes | WASM | No filesystem in the browser | String and `Uint8Array` APIs |
| On-disk workspace of runs | WASM | No filesystem | In-memory workspace exported as a zip |
| External engines (DFTB+, xtb, Zeo++) and process probes | WASM | No subprocesses in the browser | Native engines (UFF, QEq, EQeq, grid porosity) |
| In-process GFN0/GFN1/GFN2 (`sci-form` feature) | CLI, Python, WASM | Optional build-time feature whose methods are molecular only; the Rust `Calculator` and the registry listing are the supported surface, and no binding negotiation for the feature exists yet | Rust `quantum::adapter`; CLI `ff engines` lists the registry entries |
| Multi-threading | WASM | Threads need cross-origin isolation | Optional threaded build |
| MCP server over stdio | Python, WASM | It is a process, served by the CLI | `xtalforge mcp` |
| Machine-learning potentials (`MlCalculator`, `ml_backends`) | Rust, CLI, WASM | The Task 4.11 decision is the Python/ASE bridge; torch is not a WASM target and the Rust `ml` feature is a declared no-op | Python `MlCalculator`; the CLI `models` command manages weights |
| Plugins (`Plugins`, runtime registries) | CLI, WASM | Python callables and Rust registries are host-side extension points; a browser has no dynamic library or Python runtime | Rust `FormatRegistry`, `ScalarRegistry`, `StyleRegistry`, `CalculatorRegistry`, `ModuleRegistry`; Python `Plugins` |
| Model cache (`models`) | Python, WASM | It is a filesystem command with `curl`; the browser and Python package leave the download to the user | CLI `models`; `xtalforge_ff::models` |
| Full MOF topology bundle | WASM | Package size | The builder carries the curated six-topology/eight-block subset (Tasks 6.11-6.16); the RCSR fingerprint index (35 KB) is embedded on every surface |
| lattice_minimum_image, limits_bounds, limits_check | CLI | The CLI exposes geometry through `cell` and the global `--limits` flag | Use `cell` and `--limits` |
| workspace_archive | CLI | The CLI has no `modules`/`run` command; the workspace is a library API | Rust `Workspace` and Python `workspace_*`; `project` for archives |

## Open gaps: Rust-only capabilities

These capabilities exist in the Rust API but have **no** CLI, Python or WASM
surface and no justification in the table above, so under rule 4 of
`AGENTS.md` they are not finished. They are recorded here instead of being
hidden; closing them is post-1.0 work.

| Capability | Rust | CLI | Python | WASM |
|---|---|---|---|---|
| SMILES parser (`parse_smiles`) and 3D embedding (`embed_smiles`), Tasks 6.1.1-6.1.3 | done | open | open | open |
| Sketch data model (`Sketch`, `from_smiles`, `to_smiles`), Task 6.1.6 | done | open | open | open |
| Fragment library (`find_fragment`), Task 6.2 | done | open | open | open |
| Substitution and clearance (`substitute`), Task 6.3 | done | open | open | open |
| Guest placement (`guest` module), Task 6.4 | done | open | open | open |
