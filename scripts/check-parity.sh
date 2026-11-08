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
    groups) name=find_groups ;;
    slab) name=make_slab ;;
  esac
  grep -Rq "pub fn $name(" crates/*/src
}

check_python() {
  grep -q "def $1(" "$PY_STUB"
}

check_wasm() {
  local name=$1
  case "$name" in
    supercell) name=supercell_info ;;
    slab) name=slab_info ;;
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
