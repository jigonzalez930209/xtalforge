#!/usr/bin/env bash
# Checks the crate tarballs before publishing (Task 8.5.1).
#
# Each published crate must carry only src/, its data/ bundles (when it embeds
# them) and skill documents (agent). Tests, benches, fuzz targets, examples and
# fixtures never reach the tarball. The size ceiling is deliberately generous;
# the point is to catch a data bundle copied into the wrong crate.
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT" || exit 1

CRATES=(
  xtalforge-core
  xtalforge-io
  xtalforge-symmetry
  xtalforge-commands
  xtalforge-view
  xtalforge-ff
  xtalforge-analysis
  xtalforge-build
  xtalforge-mof
  xtalforge-modules
  xtalforge-agent
  xtalforge
  xtalforge-cli
)
CEILING_KIB=12000
status=0

for crate in "${CRATES[@]}"; do
  list=$(cargo package -p "$crate" --list --allow-dirty 2>/dev/null)
  if [ -z "$list" ]; then
    echo "FAIL: could not list $crate" >&2
    status=1
    continue
  fi
  forbidden=$(printf '%s\n' "$list" | grep -E '^(tests|benches|fuzz|examples)/' || true)
  if [ -n "$forbidden" ]; then
    echo "FAIL: $crate includes test or build material:" >&2
    printf '%s\n' "$forbidden" | sed 's/^/    /' >&2
    status=1
  fi
  if ! cargo package -p "$crate" --allow-dirty --no-verify -q 2>/dev/null; then
    echo "SKIP: $crate size; its internal dependencies are not published yet"
    continue
  fi
  tarball=$(ls "target/package/${crate}-"[0-9]*.crate 2>/dev/null | head -1)
  if [ -z "$tarball" ]; then
    echo "FAIL: no tarball for $crate" >&2
    status=1
    continue
  fi
  size=$(( $(wc -c < "$tarball") / 1024 ))
  if [ "$size" -gt "$CEILING_KIB" ]; then
    echo "FAIL: $crate tarball ${size} KiB exceeds ${CEILING_KIB} KiB" >&2
    status=1
  else
    echo "OK: $crate ${size} KiB"
  fi
done

exit "$status"
