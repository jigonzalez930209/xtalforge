#!/usr/bin/env bash
# Builds the split WASM packages (Task 15.3.1).
#
# `full` is the npm package `xtalforge-wasm` and keeps the 1.7.0 surface in
# `crates/xtalforge-wasm/pkg`. The other four packages live in subdirectories
# and are named `xtalforge-wasm-<tier>`, so a viewer can install only the tier
# it uses:
#
#   core      structure reading, conversion, supercells/slabs, scene and mesh
#             buffers, documents, the object handle
#   symmetry  core plus space groups, Wyckoff, standardization, reduction,
#             subgroups, the Brillouin zone and the matcher
#   analysis  core plus PXRD, porosity, nets, defects, orderings, mechanics
#             and powder refinement
#   build     core plus graphene, nanotubes, fullerenes, polymers, SMILES,
#             fragments and hydrogen filling
#   full      every export, including force fields, MOF and the agent surface
#
# Usage: scripts/build-wasm-packages.sh [core] [symmetry] [analysis] [build] [full]
# With no arguments every package is built (full last).
#
# Each package is loaded and exercised by
# crates/xtalforge-wasm/tests/packages.mjs; `scripts/check-wasm-size.sh`
# enforces the per-package budget.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
WASM_DIR="$ROOT/crates/xtalforge-wasm"

if ! command -v wasm-pack >/dev/null 2>&1; then
  echo "build-wasm-packages: wasm-pack is not installed" >&2
  exit 2
fi
if ! command -v node >/dev/null 2>&1; then
  echo "build-wasm-packages: node is not installed" >&2
  exit 2
fi

# Rewrites the generated package.json name and description for a split tier.
rename_package() {
  local manifest=$1 name=$2
  node -e '
    const fs = require("node:fs");
    const [manifest, name] = process.argv.slice(1);
    const doc = JSON.parse(fs.readFileSync(manifest, "utf8"));
    const tier = name.split("-").pop();
    doc.name = name;
    doc.description = `WebAssembly bindings for xtalforge (${tier} tier).`;
    fs.writeFileSync(manifest, `${JSON.stringify(doc, null, 2)}\n`);
  ' "$manifest" "$name"
}

build_split() {
  local tier=$1
  local out="$WASM_DIR/pkg/$tier"
  local out_name="xtalforge_wasm_$tier"
  echo "== building $tier"
  (
    cd "$WASM_DIR"
    wasm-pack build --release --target web --out-dir "pkg/$tier" --out-name "$out_name" \
      -- --no-default-features --features "$tier" >/dev/null
  )
  rename_package "$out/package.json" "xtalforge-wasm-$tier"
  node "$WASM_DIR/tests/packages.mjs" "$out" "$tier"
  local bytes
  bytes=$(wc -c < "$out/${out_name}_bg.wasm")
  printf 'OK: %-8s %5s KiB\n' "$tier" "$((bytes / 1024))"
}

build_full() {
  echo "== building full"
  (cd "$WASM_DIR" && wasm-pack build --release --target web >/dev/null)
  node "$WASM_DIR/tests/packages.mjs" "$WASM_DIR/pkg" "full"
  local bytes
  bytes=$(wc -c < "$WASM_DIR/pkg/xtalforge_wasm_bg.wasm")
  printf 'OK: %-8s %5s KiB\n' "full" "$((bytes / 1024))"
}

tiers=("$@")
if [ "${#tiers[@]}" -eq 0 ]; then
  tiers=(core symmetry analysis build full)
fi

for tier in "${tiers[@]}"; do
  case "$tier" in
    core|symmetry|analysis|build) build_split "$tier" ;;
    full) build_full ;;
    *)
      echo "build-wasm-packages: unknown tier $tier" >&2
      exit 2
      ;;
  esac
done
