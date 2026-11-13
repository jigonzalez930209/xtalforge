#!/usr/bin/env bash
# Checks the typed object interface has no `any` and the generated package
# exports the Structure handle and the typed array helpers (Task 12.3.2).
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TS="$ROOT/crates/xtalforge-wasm/ts/xtalforge-object.d.ts"
DTS="$ROOT/crates/xtalforge-wasm/pkg/xtalforge_wasm.d.ts"
status=0

if grep -nE '\bany\b' "$TS" >/dev/null; then
  echo "FAIL: $TS contains 'any'" >&2
  status=1
fi
if [ -f "$DTS" ]; then
  for symbol in "class Structure" "scene_buffers" "pxrd_buffers" "mesh_buffers" "fetch_optimade" "fetch_cod"; do
    grep -q "$symbol" "$DTS" || { echo "FAIL: $DTS is missing $symbol" >&2; status=1; }
  done
fi
[ "$status" -eq 0 ] && echo "OK: object API interface is typed and exported"
exit "$status"
