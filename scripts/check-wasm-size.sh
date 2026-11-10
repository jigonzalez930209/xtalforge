#!/usr/bin/env bash
# Checks the release WASM size against the budget (Task 0.5.5).
#
# 3840 KiB since Phase 8: the agent session, its diagnostics and the PNG
# renderer add functions every surface exports, and the 3648 KiB budget from
# Phase 7 no longer fits (3775 KiB before this bump).
#
# Build first:
#   (cd crates/xtalforge-wasm && wasm-pack build --release --target web)
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
WASM="$ROOT/crates/xtalforge-wasm/pkg/xtalforge_wasm_bg.wasm"
BUDGET_KIB=3840

if [ ! -f "$WASM" ]; then
  echo "check-wasm-size: build the package first: $WASM" >&2
  exit 2
fi

size=$(wc -c < "$WASM")
budget=$((BUDGET_KIB * 1024))
human=$((size / 1024))
if [ "$size" -gt "$budget" ]; then
  echo "FAIL: ${human} KiB exceeds the ${BUDGET_KIB} KiB WASM budget" >&2
  exit 1
fi
echo "OK: ${human} KiB (budget ${BUDGET_KIB} KiB)"
