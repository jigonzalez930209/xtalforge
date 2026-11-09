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
PY_STUB="crates/xtalforge-python/python/xtalforge/__init__.pyi"
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
  esac
  grep -Rq "pub fn $name(" crates/*/src
}

check_python() {
  local name=$1
  case "$name" in
    prepare) name=ff_prepare ;;
    ring_paths) name=ring_paths ;;
    run_automatic) name=powder_automatic ;;
    refine_phase) name=powder_refine ;;
  esac
  grep -q "def $name(" "$PY_STUB"
}

check_wasm() {
  local name=$1
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
    build_polymer) name=polymer_info ;;
    run_automatic) name=powder_automatic_info ;;
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
    scene_summary|mesh_stl) command=mesh ;;
    pxrd) command=pxrd ;;
    build_graphene|build_nanotube|build_c20|build_c60|twist|ribbon|build_ztc|build_schwarzite|build_polymer) command=build ;;
    build_mof) command=mof ;;
    run_automatic|refine_phase) command=powder ;;
    param_schema) command=params ;;
    project_pack) command=project ;;
    prepare) command=ff ;;
    add_hydrogens) command=hydrogens ;;
    porosity|net_invariants|zeo_parse|ring_paths) command=analysis ;;
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
