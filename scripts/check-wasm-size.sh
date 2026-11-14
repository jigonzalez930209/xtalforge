#!/usr/bin/env bash
# Checks the release WASM size against the budgets (Tasks 0.5.5 and 15.3.1).
#
# The full package is the npm `xtalforge-wasm` build in pkg/. Task 15.3.1 adds
# four split packages, each with its own budget:
#
#   core      viewer surface (reading, conversion, scene and mesh buffers)
#   symmetry  core plus space groups, Wyckoff, standardization and the matcher
#   analysis  core plus PXRD, porosity, nets, defects and powder refinement
#   build     core plus carbon, polymer, SMILES, fragment and hydrogen builders
#   full      every export, including force fields, MOF and the agent surface
#
# Measured on 2026-11-13 (release build): core 1957, symmetry 2807,
# analysis 3057, build 2367, full 5346 KiB. The budgets leave about 5 % of
# headroom over those numbers; the full budget keeps the historical 1.7.0
# ceiling.
# History of the full-package budget:
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
# 5120 KiB since 1.3.0: the MOF platform checker, deconstruction, MOFid,
# native percolation, TraPPE host-guest fields, Widom/GCMC/grids and
# modification (4964 KiB measured).
#
# 5376 KiB since 1.6.0: the Phase 13 formats, PDB/mmCIF readers, VASP
# trajectories, extxyz properties and browser remote fetchers (5301 KiB
# measured).
#
# 5504 KiB since 1.7.0: the diffraction profile model, MD/scan control, mesh
# export and net contraction (5398 KiB measured).
#
# 5632 KiB since 1.9.0: workflows, provenance and dataset descriptors
# 5888 KiB since 1.10.0: Phase 17 parity surfaces, YAML, MOF data loader and
#   the in-memory workflow runner
# (5606 KiB measured).
#
# Build first:
#   scripts/build-wasm-packages.sh
# or, for the full package only:
#   (cd crates/xtalforge-wasm && wasm-pack build --release --target web)
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PKG="$ROOT/crates/xtalforge-wasm/pkg"
status=0

# Prints the package.json name and size of one build.
check_package() {
  local label=$1 wasm=$2 budget_kib=$3
  if [ ! -f "$wasm" ]; then
    echo "FAIL: $label is missing $wasm" >&2
    status=1
    return
  fi
  local size budget human
  size=$(wc -c < "$wasm")
  budget=$((budget_kib * 1024))
  human=$((size / 1024))
  if [ "$size" -gt "$budget" ]; then
    echo "FAIL: ${label} ${human} KiB exceeds the ${budget_kib} KiB budget" >&2
    status=1
  else
    echo "OK: ${label} ${human} KiB (budget ${budget_kib} KiB)"
  fi
}

check_package "full" "$PKG/xtalforge_wasm_bg.wasm" 5888

split=(core symmetry analysis build)
present=0
for tier in "${split[@]}"; do
  [ -d "$PKG/$tier" ] && present=$((present + 1))
done

if [ "$present" -eq 0 ]; then
  echo "NOTE: split packages not built; run scripts/build-wasm-packages.sh to check them"
elif [ "$present" -ne "${#split[@]}" ]; then
  echo "FAIL: only ${present} of ${#split[@]} split packages exist under $PKG" >&2
  status=1
else
  check_package "core" "$PKG/core/xtalforge_wasm_core_bg.wasm" 2048
  check_package "symmetry" "$PKG/symmetry/xtalforge_wasm_symmetry_bg.wasm" 2944
  check_package "analysis" "$PKG/analysis/xtalforge_wasm_analysis_bg.wasm" 3264
  check_package "build" "$PKG/build/xtalforge_wasm_build_bg.wasm" 2496
fi

exit "$status"
