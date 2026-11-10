#!/usr/bin/env bash
# Runs the process-level benchmarks and fails on a large regression
# (Task 8.6.1). Ceilings are deliberately loose: they catch a 10x regression,
# not normal machine noise.
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT" || exit 1

status=0

check() {
    local label=$1 value=$2 ceiling=$3
    if [ -z "$value" ]; then
        echo "FAIL: no timing for $label" >&2
        status=1
        return
    fi
    if awk -v v="$value" -v c="$ceiling" 'BEGIN { exit !(v > c) }'; then
        echo "FAIL: $label ${value} ms exceeds ${ceiling} ms" >&2
        status=1
    else
        echo "OK: $label ${value} ms (ceiling ${ceiling} ms)"
    fi
}

symmetry=$(cargo bench -q -p xtalforge-symmetry --bench detection 2>/dev/null)
value=$(printf '%s\n' "$symmetry" | awk '/ 64 atoms/ {print $(NF-1)}')
check "symmetry 64 atoms" "$value" 100

uff=$(cargo bench -q -p xtalforge-ff --bench uff 2>/dev/null)
value=$(printf '%s\n' "$uff" | awk '/methane/ {print $(NF-1)}')
check "UFF methane" "$value" 5
value=$(printf '%s\n' "$uff" | awk '/sites/ {print $(NF-3)}')
check "UFF MOF-5" "$value" 500

analysis=$(cargo bench -q -p xtalforge-analysis --bench pxrd 2>/dev/null)
value=$(printf '%s\n' "$analysis" | awk '/PXRD/ {print $(NF-1)}')
check "PXRD MOF-5" "$value" 2000
value=$(printf '%s\n' "$analysis" | awk '/porosity/ {print $(NF-1)}')
check "porosity MOF-5" "$value" 15000

exit "$status"
