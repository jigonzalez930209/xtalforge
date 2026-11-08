#!/usr/bin/env bash
# Enforces the repository rules from AGENTS.md:
#   1. No code file longer than 500 lines (Markdown is exempt).
#   2. No dead-code, unused or warnings suppression attributes in Rust code.
#
# Usage: scripts/check-rules.sh [ROOT]
#   ROOT defaults to the repository root (the parent of this script).
#
# The language rule was retired by the maintainer on 2026-11-08; there is no
# language check in this script or in CI.
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="${1:-$(cd "$SCRIPT_DIR/.." && pwd)}"
MAX_LINES=500
status=0

prune=( -path '*/target' -o -path '*/node_modules' -o -path '*/.git'
        -o -path '*/.vitepress/cache' -o -path '*/.vitepress/dist'
        -o -path '*/pkg' -o -path '*/.venv' -o -path '*/__pycache__' )

list_files() {
  # $@: -name patterns joined with -o
  find "$ROOT" \( "${prune[@]}" \) -prune -o -type f \( "$@" \) -print
}

echo "== Rule 1: code files must not exceed ${MAX_LINES} lines =="
while IFS= read -r f; do
  n=$(wc -l < "$f")
  if [ "$n" -gt "$MAX_LINES" ]; then
    echo "TOO LONG ($n lines): ${f#"$ROOT"/}"
    status=1
  fi
done < <(list_files -name '*.rs' -o -name '*.py' -o -name '*.pyi' \
  -o -name '*.ts' -o -name '*.mts' -o -name '*.js' -o -name '*.mjs' \
  -o -name '*.vue' -o -name '*.sh' -o -name '*.toml' -o -name '*.yml' \
  -o -name '*.yaml')

echo "== Rule 2: no warning-suppression attributes in Rust code =="
# Warnings are errors in the workspace. Dead-code and unused warnings must be
# fixed, not hidden behind allow/expect attributes. The only exception is the
# test-only clippy configuration at the crate root.
suppression='#!?\[[[:space:]]*(allow|expect)\([[:space:]]*(dead_code|unused[a-z_]*|warnings)\b'
while IFS= read -r f; do
  if grep -n -E -m5 "$suppression" "$f" > /dev/null; then
    echo "SUPPRESSED WARNING in ${f#"$ROOT"/}:"
    grep -n -E -m5 "$suppression" "$f" | sed 's/^/    /'
    status=1
  fi
done < <(list_files -name '*.rs')

if [ "$status" -eq 0 ]; then echo "OK: all rules satisfied"; fi
exit "$status"
