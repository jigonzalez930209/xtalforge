#!/usr/bin/env bash
# Checks the release WASM size against the budget (Task 0.5.5).
#
# 3904 KiB since the 3D playground: scene_data, n_atoms, the text-returning
# builders and the Brillouin-zone fix add exported surface (3845 KiB measured);
# the Phase 8 budget of 3840 KiB no longer fits.
#
# 3968 KiB since the functional audit of 2026-11-11: the International Tables
# form factors, the commensurate twist, the periodic polymer, the embedder
# angle term, the carbon relaxation, the remeshed ZTC/Schwarz P dual, the
# Le Bail/Pawley background and the MOF quaternion fit (3943 KiB measured).
#
# 4736 KiB since 1.2.0: the Phase 9 matcher, Voronoi/shape kernels, bond
# valence data, surface and defect generators, ordering/SQS search and
# elasticity/EOS/hull fits (4603 KiB measured).
#
# Build first:
#   (cd crates/xtalforge-wasm && wasm-pack build --release --target web)
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
WASM="$ROOT/crates/xtalforge-wasm/pkg/xtalforge_wasm_bg.wasm"
BUDGET_KIB=4736

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
