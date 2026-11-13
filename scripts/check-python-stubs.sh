#!/usr/bin/env bash
# Checks that the public Python stub has no `Any` in a signature (Task 12.1.3).
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
STUB="$ROOT/crates/xtalforge-python/python/xtalforge"
status=0
if grep -RIn --include='*.pyi' '\bAny\b' "$STUB"; then
  echo "FAIL: public Python stubs contain Any" >&2
  status=1
fi
[ "$status" -eq 0 ] && echo "OK: Python stubs have no Any"
exit "$status"
