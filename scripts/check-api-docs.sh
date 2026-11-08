#!/usr/bin/env bash
# Checks that every public Rust item of the core crate has an API page
# (Task 0.7.2). The scan is intentionally simple: public functions, structs,
# enums, constants and type aliases; modules and re-exports are covered by the
# facade page.
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT" || exit 1

status=0
items=$(grep -RhoE '^pub (fn|struct|enum|const|type) [A-Za-z_][A-Za-z0-9_]*' \
  crates/xtalforge-core/src | awk '{print $3}' | sort -u)

if [ -z "$items" ]; then
  echo "FAIL: no public items found in the core crate" >&2
  exit 1
fi

for item in $items; do
  if ! grep -Rq --exclude=parity.md -- "$item" docs/api; then
    echo "FAIL: no API page mentions '$item'" >&2
    status=1
  fi
done

if [ "$status" -eq 0 ]; then
  echo "OK: every public core item has an API page"
fi
exit "$status"
