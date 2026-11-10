#!/usr/bin/env bash
# Builds the WASM package and the docs, regenerates the independent reference
# data and runs the Playwright end-to-end battery (Task 8.9).
#
# Usage: scripts/e2e.sh [--skip-build] [--update-goldens]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PYTHON="${PYTHON:-$ROOT/.venv/bin/python}"
if [ ! -x "$PYTHON" ]; then
  PYTHON="$(command -v python3)"
fi
SKIP_BUILD=0
UPDATE_GOLDENS=0
for argument in "$@"; do
  case "$argument" in
    --skip-build) SKIP_BUILD=1 ;;
    --update-goldens) UPDATE_GOLDENS=1 ;;
    *)
      echo "unknown argument: $argument" >&2
      exit 2
      ;;
  esac
done

if [ "$SKIP_BUILD" -eq 0 ]; then
  command -v wasm-pack >/dev/null || {
    echo "wasm-pack is required: cargo install wasm-pack --version 0.13.1 --locked" >&2
    exit 1
  }
  (cd "$ROOT/crates/xtalforge-wasm" && wasm-pack build --release --target web)
  if [ ! -d "$ROOT/docs/node_modules" ]; then
    (cd "$ROOT/docs" && npm install --no-audit --no-fund)
  fi
  (cd "$ROOT/docs" && npm run docs:build)
fi

"$PYTHON" - <<'PY'
import importlib
for name in ("ase", "gemmi", "numpy", "PIL", "playwright", "pytest", "scipy", "spglib"):
    try:
        importlib.import_module(name)
    except ImportError as error:
        raise SystemExit(
            f"missing {name}: {error}\n"
            "install the battery with: "
            ".venv/bin/python -m pip install -r tests/e2e/requirements.txt"
        )
PY

"$PYTHON" "$ROOT/tests/e2e/reference/generate.py" --check-samples
if [ "$UPDATE_GOLDENS" -eq 1 ]; then
  XF_UPDATE_GOLDENS=1 "$PYTHON" -m pytest "$ROOT/tests/e2e" -q
else
  "$PYTHON" -m pytest "$ROOT/tests/e2e" -q
fi
