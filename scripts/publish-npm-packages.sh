#!/usr/bin/env bash
# Publishes the split WASM packages to npm (Task 17.2.4 follow-up).
#
# `scripts/build-wasm-packages.sh` writes the five local packages:
#   xtalforge-wasm, xtalforge-wasm-core, xtalforge-wasm-symmetry,
#   xtalforge-wasm-analysis and xtalforge-wasm-build.
#
# This script checks that every manifest is present and runs `npm publish` in
# the requested order: the tier packages first, the default `full` package
# last, so a partial upload never leaves the default pointing at a missing
# tier. The default is a dry run; pass --publish to upload for real.
#
# Usage: scripts/publish-npm-packages.sh [--publish] [--tag next]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
WASM_DIR="$ROOT/crates/xtalforge-wasm"
PUBLISH=0
TAG=latest

while [ "$#" -gt 0 ]; do
  case "$1" in
    --publish) PUBLISH=1 ;;
    --tag)
      shift
      TAG="${1:-latest}"
      ;;
    -h|--help)
      sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *)
      echo "publish-npm-packages: unknown argument: $1" >&2
      exit 2
      ;;
  esac
  shift
done

if ! command -v npm >/dev/null 2>&1; then
  echo "publish-npm-packages: npm is not installed" >&2
  exit 2
fi

entries=(
  "core:$WASM_DIR/pkg/core"
  "symmetry:$WASM_DIR/pkg/symmetry"
  "analysis:$WASM_DIR/pkg/analysis"
  "build:$WASM_DIR/pkg/build"
  "full:$WASM_DIR/pkg"
)

for entry in "${entries[@]}"; do
  name=${entry%%:*}
  directory=${entry#*:}
  manifest="$directory/package.json"
  if [ ! -f "$manifest" ]; then
    echo "publish-npm-packages: missing $manifest; run scripts/build-wasm-packages.sh first" >&2
    exit 2
  fi
  package=$(node -p "require('$manifest').name")
  version=$(node -p "require('$manifest').version")
  if [ "$PUBLISH" -eq 1 ]; then
    echo "== publishing $package@$version (tag $TAG)"
    (cd "$directory" && npm publish --access public --tag "$TAG")
  else
    echo "would publish $package@$version from $directory (tag $TAG)"
  fi
done

if [ "$PUBLISH" -eq 0 ]; then
  echo "dry run only; pass --publish to upload"
fi
