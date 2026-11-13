#!/usr/bin/env bash
# Regression gate for the comparative benchmark suite (Task 15.1.2).
#
# Runs the xtalforge rows of `benches/compare/run.py`, compares each median with
# `benches/compare/baseline.json`, and fails when a tracked case regresses by
# more than 10 %. The baseline is machine-specific: regenerate it on the same
# host with `--update` after an intentional change.
#
# Usage:
#   scripts/check-bench.sh                 # gate
#   scripts/check-bench.sh --update        # rewrite the baseline
#   scripts/check-bench.sh --self-test     # compare logic only, no timings
#   scripts/check-bench.sh --reps 1        # faster local smoke run
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PYTHON="${PYTHON:-$ROOT/.venv/bin/python}"
[ -x "$PYTHON" ] || PYTHON="$(command -v python3)"
BASELINE="$ROOT/benches/compare/baseline.json"
COMPARE="$ROOT/benches/compare/compare.py"
RUN="$ROOT/benches/compare/run.py"

mode=check
reps=()
while [ $# -gt 0 ]; do
    case "$1" in
        --update) mode=update ;;
        --self-test) mode=self-test ;;
        --reps) reps=(--reps "$2"); shift ;;
        --threshold) threshold=$2; shift ;;
        -h|--help)
            sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *)
            echo "unknown option: $1" >&2
            exit 2
            ;;
    esac
    shift
done

if [ "$mode" = "self-test" ]; then
    exec "$PYTHON" "$COMPARE" --self-test
fi

"$PYTHON" -c "import xtalforge" >/dev/null 2>&1 || {
    echo "xtalforge is not importable by $PYTHON" >&2
    echo "build and install the wheel first (maturin build --release -m crates/xtalforge-python/Cargo.toml)" >&2
    exit 2
}

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "== benchmark gate: xtalforge rows =="
"$PYTHON" "$RUN" --only xtalforge --quiet --output "$TMP/current.json" "${reps[@]}"

if [ "$mode" = "update" ]; then
    "$PYTHON" "$COMPARE" --baseline "$BASELINE" --current "$TMP/current.json" --update
else
    args=(--baseline "$BASELINE" --current "$TMP/current.json")
    if [ -n "${threshold:-}" ]; then
        args+=(--threshold "$threshold")
    fi
    if ! "$PYTHON" "$COMPARE" "${args[@]}"; then
        echo "retrying once to separate load noise from a regression" >&2
        "$PYTHON" "$RUN" --only xtalforge --quiet \
            --output "$TMP/current-2.json" "${reps[@]}"
        retry=(--baseline "$BASELINE"
               --current "$TMP/current.json"
               --current "$TMP/current-2.json")
        if [ -n "${threshold:-}" ]; then
            retry+=(--threshold "$threshold")
        fi
        "$PYTHON" "$COMPARE" "${retry[@]}"
    fi
fi
