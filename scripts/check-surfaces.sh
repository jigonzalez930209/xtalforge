#!/usr/bin/env bash
# Runs the same operations through the CLI, the Python bindings and the WASM
# package on the eight docs samples and compares every field (functional
# audit 2026-11-11). Needs a release CLI, the wheel installed in the Python
# environment and a built `crates/xtalforge-wasm/pkg` (wasm-pack, target web).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PYTHON="${PYTHON:-$ROOT/.venv/bin/python}"
[ -x "$PYTHON" ] || PYTHON="$(command -v python3)"
OUT="$(mktemp -d)"
trap 'rm -rf "$OUT"' EXIT
cd "$ROOT"
cargo build --release -p xtalforge-cli >/dev/null
"$PYTHON" tests/surfaces/cli_surface.py "$OUT/cli.json"
"$PYTHON" tests/surfaces/python_surface.py "$OUT/python.json"
node tests/surfaces/wasm_surface.mjs "$OUT/wasm.json"
"$PYTHON" tests/surfaces/compare.py "$OUT/cli.json" "$OUT/python.json" "$OUT/wasm.json"
