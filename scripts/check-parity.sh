#!/usr/bin/env bash
# Checks docs/api/parity.md against the symbols the four surfaces export
# (Task 0.5.6).
#
# The Python and WASM checks use the generated stub/declaration when it is
# present and fall back to the source that produces it. The CLI check runs
# `cargo run --help`, so it needs a built workspace.
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT" || exit 1

MATRIX="docs/api/parity.md"
PY_STUB="crates/xtalforge-python/python/xtalforge"
WASM_TS="crates/xtalforge-wasm/pkg/xtalforge_wasm.d.ts"
WASM_SRC="crates/xtalforge-wasm/src/lib.rs"
status=0

trim() {
  printf '%s' "$1" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//'
}

check_rust() {
  local name=$1
  case "$name" in
    bonds) name=bond_graph ;;
    pxrd) name=calculate_pxrd ;;
    build_graphene) name=graphene ;;
    build_nanotube) name=nanotube ;;
    build_c60) name=fullerene ;;
    build_mof) name=build_mof ;;
    build_c20) name=face_spiral ;;
    twist) name=commensurate_twist ;;
    ribbon) name=nanoribbon ;;
    build_ztc) name=build_ztc ;;
    build_schwarzite) name=build_schwarzite ;;
    build_polymer) name=build_polymer ;;
    run_automatic) name=run_automatic ;;
    refine_phase) name=refine_phase ;;
    powder_index) name=index_pattern ;;
    powder_pawley) name=pawley_refine ;;
    powder_lebail) name=le_bail ;;
    param_schema) name=preset ;;
    project_pack) name=write_project ;;
    workspace_archive) name=archive_records ;;
    prepare) name=prepare ;;
    add_hydrogens) name=add_hydrogens ;;
    porosity) name=distance_grid ;;
    net_invariants) name=quotient_graph ;;
    zeo_parse) name=parse_zeo_channels ;;
    groups) name=find_groups ;;
    ring_paths) name=ring_paths ;;
    subgroups) name=translationengleiche_subgroups ;;
    slab) name=make_slab ;;
    agent_capabilities) name=capability_report ;;
    agent_verbs) name=catalogue ;;
    agent_run) name=run_verbs ;;
    models) name=fetch ;;
    edit_structure) name=replace_species ;;
    structure_fingerprint) name=fingerprint ;;
    eos_fit) name=fit_eos ;;
    mof_check) name=check_mof ;;
    mof_solvent_remove) name=remove_solvent ;;
    mof_disorder) name=resolve_disorder ;;
    mof_charge_balance) name=charge_balance ;;
    mof_deconstruct) name=deconstruct ;;
    mof_mofid) name=mofid ;;
    mof_open_metal_sites) name=open_metal_sites ;;
    mof_descriptors) name=framework_descriptors ;;
    mof_pores) name=pore_descriptors ;;
    mof_guests) name=guests ;;
    mof_forcefield) name=force_field ;;
    mof_widom) name=widom_insertion ;;
    mof_gcmc) name=gcmc_isotherm ;;
    mof_energy_grid) name=energy_grid ;;
    mof_raspa) name=raspa_input ;;
    mof_functionalize) name=functionalize ;;
    mof_defects) name=missing_defects ;;
    mof_interpenetrate) name=interpenetrate ;;
    mof_stack_layers) name=stack_layers ;;
    mof_mixed) name=mixed_framework ;;
    mof_read_database) name=read_core_mof ;;
    mof_screen) name=screen_database ;;
    mof_deduplicate) name=deduplicate ;;
    read) name=read ;;
    write) name=write ;;
    set_threads) name=set_global_threads ;;
    error_category) name=category ;;
    molecule_from_smiles) name=parse_smiles ;;
    sketch_roundtrip) name=from_smiles ;;
    fragment_library) name=fragments ;;
    fragment_structure) name=find_fragment ;;
    substitute_structure) name=substitute ;;
    place_guest_structure) name=place_guest ;;
    structure_handle) name=structure_info ;;
    structure_from_arrays) name=add_sites ;;
    typed_results) name=structure_info ;;
    typed_arrays) name=radial_distribution ;;
    batch_api) name=read ;;
    wasm_worker) name=build_scene ;;
    optimade_fetch) name=optimade_query ;;
    cod_fetch) name=cod_fetch ;;
    materials_project_fetch) name=materials_project_fetch ;;
    interop_formats) grep -Rq 'name: "castep"' crates/xtalforge-io/src; return ;;
  esac
  grep -Rq "pub fn $name(" crates/*/src
}

check_python() {
  local name=$1
  case "$name" in
    structure_handle) grep -Rq --include='*.pyi' "class Structure" "$PY_STUB"; return ;;
    structure_from_arrays) grep -Rq --include='*.pyi' "def from_arrays" "$PY_STUB"; return ;;
    typed_results) grep -Rq --include='*.pyi' "class SymmetryDataset" "$PY_STUB"; return ;;
    typed_arrays) grep -Rq --include='*.pyi' "positions: list" "$PY_STUB"; return ;;
    batch_api) grep -q "def read(" "$ROOT/crates/xtalforge-python/python/xtalforge/batch.pyi"; return ;;
    ase_adapter) grep -Rq --include='*.pyi' "def from_ase" "$PY_STUB"; return ;;
    pymatgen_adapter) grep -Rq --include='*.pyi' "def from_pymatgen" "$PY_STUB"; return ;;
    phonopy_adapter) grep -Rq --include='*.pyi' "def to_phonopy_atoms" "$PY_STUB"; return ;;
    rdkit_adapter) grep -Rq --include='*.pyi' "def from_rdkit" "$PY_STUB"; return ;;
    get_symmetry_dataset) grep -Rq --include='*.pyi' "def get_symmetry_dataset" "$PY_STUB"; return ;;
    python_plugins) test -f "$PY_STUB/plugins.py"; return ;;
    interop_formats) grep -q "def read(" "$PY_STUB/__init__.pyi"; return ;;
  esac
  case "$name" in
    prepare) name=ff_prepare ;;
    ring_paths) name=ring_paths ;;
    run_automatic) name=powder_automatic ;;
    powder_index) name=powder_index ;;
    powder_pawley) name=powder_pawley ;;
    refine_phase) name=powder_refine ;;
  esac
  grep -Rq --include='*.pyi' "def $name(" "$PY_STUB"
}

check_wasm() {
  local name=$1
  case "$name" in
    structure_handle) grep -q "class Structure" "$WASM_TS"; return ;;
    structure_from_arrays) grep -q "set_positions(" "$WASM_TS"; return ;;
    typed_results) grep -q "toSummary(" "$WASM_TS"; return ;;
    typed_arrays) grep -q "function scene_buffers(" "$WASM_TS"; return ;;
    wasm_worker) test -f "$ROOT/crates/xtalforge-wasm/js/worker.js"; return ;;
    node_entry) test -f "$ROOT/crates/xtalforge-wasm/js/node.js"; return ;;
    optimade_fetch) grep -q "pub async fn fetch_optimade" "$ROOT/crates/xtalforge-wasm/src/remote_export.rs"; return ;;
    cod_fetch) grep -q "pub async fn fetch_cod" "$ROOT/crates/xtalforge-wasm/src/remote_export.rs"; return ;;
    interop_formats) grep -q "function read(" "$WASM_TS"; return ;;
  esac
  case "$name" in
    supercell) name=supercell_info ;;
    slab) name=slab_info ;;
    asymmetrize) name=asymmetrize_info ;;
    to_p1) name=to_p1_info ;;
    set_space_group) name=set_space_group_info ;;
    subgroups) name=subgroups_info ;;
    lower_to_subgroup) name=lower_to_subgroup_info ;;
    standardize) name=standardize_info ;;
    primitive_cell) name=primitive_info ;;
    reduce_structure) name=reduce_info ;;
    wyckoff_assignments) name=wyckoff_info ;;
    tolerance_ladder) name=tolerance_ladder_info ;;
    diagnose) name=diagnose_info ;;
    detect_layer_symmetry) name=layer_group_info ;;
    detect_magnetic_symmetry) name=magnetic_group_info ;;
    change_setting) name=change_setting_info ;;
    band_path) name=kpath_info ;;
    molecular_point_group) name=pointgroup_info ;;
    pxrd) name=pxrd_info ;;
    build_graphene) name=build_graphene_info ;;
    build_nanotube) name=build_nanotube_info ;;
    build_c60) name=build_c60_info ;;
    build_mof) name=build_mof_info ;;
    build_c20) name=build_c20_info ;;
    twist) name=twist_info ;;
    ribbon) name=ribbon_info ;;
    build_ztc) name=ztc_info ;;
    build_schwarzite) name=schwarzite_info ;;
    build_carbon) name=build_carbon_recipe_info ;;
    build_polymer) name=polymer_info ;;
    run_automatic) name=powder_automatic_info ;;
    powder_index) name=powder_index_info ;;
    powder_pawley) name=powder_pawley_info ;;
    powder_lebail) name=powder_lebail_info ;;
    refine_phase) name=powder_refine_info ;;
    param_schema) name=param_schema ;;
    project_pack) name=project_pack ;;
    workspace_archive) name=workspace_archive ;;
    prepare) name=ff_prepare_info ;;
    ring_paths) name=rings ;;
    add_hydrogens) name=hydrogens_info ;;
    porosity) name=porosity_info ;;
    net_invariants) name=net_invariants_info ;;
    zeo_parse) name=zeo_parse_info ;;
    fragment_library) name=fragment_library_info ;;
    fragment_structure) name=fragment_info ;;
    substitute_structure) name=substitute_info ;;
    place_guest_structure) name=place_guest_info ;;
  esac
  if [ -f "$WASM_TS" ]; then
    grep -q "function $name(" "$WASM_TS"
  else
    grep -q "pub fn $name(" "$WASM_SRC"
  fi
}

check_docs() {
  grep -Rq --exclude=parity.md -- "$1" docs/api
}

cli_help=$(cargo run --quiet -p xtalforge-cli -- --help 2>/dev/null) || {
  echo "FAIL: could not run the CLI to read --help" >&2
  exit 1
}

check_cli() {
  local command=${1//_/-}
  case "$1" in
    element_info) command=element ;;
    lattice_info) command=cell ;;
    lattice_minimum_image) command=cell ;;
    structure_info) command=info ;;
    read_cif|read_structure) command=info ;;
    detect_symmetry|asymmetrize|to_p1|set_space_group|subgroups|lower_to_subgroup|standardize|primitive_cell|reduce_structure|wyckoff_assignments|tolerance_ladder|diagnose|brillouin_zone|detect_layer_symmetry|detect_magnetic_symmetry|change_setting|band_path|molecular_point_group) command=symmetry ;;
    scene_summary|mesh_stl|scene_data) command=mesh ;;
    pxrd) command=pxrd ;;
    build_graphene|build_nanotube|build_c20|build_c60|twist|ribbon|build_ztc|build_schwarzite|build_carbon|build_polymer) command=build ;;
    build_mof) command=mof ;;
    run_automatic|refine_phase|powder_index|powder_pawley|powder_lebail) command=powder ;;
    param_schema) command=params ;;
    project_pack) command=project ;;
    prepare) command=ff ;;
    add_hydrogens) command=hydrogens ;;
    porosity|net_invariants|zeo_parse|ring_paths) command=analysis ;;
    agent_capabilities|agent_verbs) command=capabilities ;;
    agent_run) command=agent ;;
    ml_backends|models) command=models ;;
    mof_check|mof_solvent_remove|mof_disorder|mof_charge_balance|mof_deconstruct|mof_mofid|mof_open_metal_sites|mof_descriptors|mof_pores|mof_guests|mof_forcefield|mof_widom|mof_gcmc|mof_energy_grid|mof_raspa|mof_functionalize|mof_defects|mof_interpenetrate|mof_stack_layers|mof_mixed|mof_read_database|mof_screen|mof_deduplicate) command=mof ;;
    read) command=info ;;
    write) command=convert ;;
    set_threads|error_category) command=doctor ;;
    molecule_from_smiles) command=molecule ;;
    sketch_roundtrip) command=sketch ;;
    fragment_library|fragment_structure) command=fragment ;;
    substitute_structure) command=substitute ;;
    place_guest_structure) command=guest ;;
    structure_handle|structure_from_arrays|typed_results|typed_arrays|batch_api|wasm_worker|node_entry) command=version ;;
    optimade_fetch|cod_fetch|materials_project_fetch) command=fetch ;;
    interop_formats) command=convert ;;
    match_structures|periodic_rmsd|kabsch_alignment|group_structures) command=compare ;;
    deform) command=deform ;;
    rattle|edit_structure|interpolate|optimal_supercell|voronoi|coordination|bond_valence_sums|oxidation_states|radial_distribution|structure_fingerprint|unique_miller_indices|slab_terminations|wulff|adsorption_sites|vacancies|substitutions|interstitials|orderings|special_quasirandom|lattice_matches|stack_slabs|elastic_tensor|eos_fit|convex_hull) command=features ;;
  esac
  printf '%s\n' "$cli_help" | grep -q -- "$command"
}

rows=$(grep -E '^\| `[^`]+` \|' "$MATRIX" || true)
if [ -z "$rows" ]; then
  echo "FAIL: no capability rows found in $MATRIX" >&2
  exit 1
fi

while IFS='|' read -r _ capability _ rust cli python wasm _; do
  capability=$(trim "$capability")
  capability=${capability//\`/}
  rust=$(trim "$rust")
  cli=$(trim "$cli")
  python=$(trim "$python")
  wasm=$(trim "$wasm")

  for value in "$rust" "$cli" "$python" "$wasm"; do
    case "$value" in
      todo|wip|done|n/a) ;;
      *)
        echo "FAIL: $capability has invalid status '$value'" >&2
        status=1
        ;;
    esac
  done

  [ "$rust" = "done" ] && ! check_rust "$capability" && {
    echo "FAIL: Rust does not export '$capability'" >&2
    status=1
  }
  [ "$cli" = "done" ] && ! check_cli "$capability" && {
    echo "FAIL: CLI does not expose '$capability'" >&2
    status=1
  }
  [ "$python" = "done" ] && ! check_python "$capability" && {
    echo "FAIL: Python stub does not declare '$capability'" >&2
    status=1
  }
  [ "$wasm" = "done" ] && ! check_wasm "$capability" && {
    echo "FAIL: WASM does not export '$capability'" >&2
    status=1
  }
  [ "$rust" = "done" ] && ! check_docs "$capability" && {
    echo "FAIL: no API page documents '$capability'" >&2
    status=1
  }
done <<< "$rows"

if [ "$status" -eq 0 ]; then
  echo "OK: parity matrix matches the exported symbols"
fi
exit "$status"
