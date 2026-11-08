#!/usr/bin/env bash
# Checks that AGENTS.md section 6 carries the workspace version (Task 0.6.4).
#
# Intended for tagged commits in CI. It reads only the first `- Version:` line.
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="${1:-$(cd "$SCRIPT_DIR/.." && pwd)}"
CARGO="$ROOT/Cargo.toml"
AGENTS="$ROOT/AGENTS.md"

version=$(sed -n 's/^version = "\([^"]*\)"/\1/p' "$CARGO" | head -n 1)
line=$(grep -m 1 '^- Version:' "$AGENTS" || true)
found=$(printf '%s\n' "$line" | sed -n 's/.*`\([^`]*\)`.*/\1/p')

if [ -z "$version" ] || [ -z "$found" ]; then
  echo "FAIL: could not read the workspace or AGENTS.md version" >&2
  exit 2
fi

if [ "$version" != "$found" ]; then
  echo "FAIL: AGENTS.md says $found but the workspace is $version" >&2
  exit 1
fi

echo "OK: AGENTS.md and the workspace both say $version"
