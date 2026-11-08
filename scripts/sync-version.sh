#!/usr/bin/env bash
# Keeps every derived version location in lockstep with
# `[workspace.package].version` in the root Cargo.toml (Task 0.6.1).
#
# Usage: scripts/sync-version.sh [--check] [ROOT]
#   --check   report drift and exit non-zero; never write
#   ROOT      repository root for tests; defaults to the parent of this script
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CHECK=0

for arg in "$@"; do
  case "$arg" in
    --check) CHECK=1 ;;
    *) ROOT="$arg" ;;
  esac
done

CARGO="$ROOT/Cargo.toml"
DOCS_PKG="$ROOT/docs/package.json"
DOCS_LOCK="$ROOT/docs/package-lock.json"
README="$ROOT/README.md"
PYPROJECT="$ROOT/crates/xtalforge-python/pyproject.toml"
WASM_PKG="$ROOT/crates/xtalforge-wasm/pkg/package.json"
status=0

if [ ! -f "$CARGO" ]; then
  echo "sync-version: no Cargo.toml in $ROOT" >&2
  exit 2
fi

version=$(sed -n 's/^version = "\([^"]*\)"/\1/p' "$CARGO" | head -n 1)
if [ -z "$version" ]; then
  echo "sync-version: no [workspace.package] version in $CARGO" >&2
  exit 2
fi
case "$version" in
  *[!0-9A-Za-z.+-]*)
    echo "sync-version: '$version' is not a valid version" >&2
    exit 2
    ;;
esac

json_version() {
  sed -n 's/^[[:space:]]*"version": "\([^"]*\)".*/\1/p' "$1" | head -n 1
}

report() {
  echo "sync-version: $1" >&2
  status=1
}

check_pins() {
  local pins
  pins=$(grep -E '^xtalforge(-[a-z]+)* = \{ version = "=' "$CARGO" \
    | sed -E 's/.*version = "=([^"]+)".*/\1/')
  if [ -z "$pins" ]; then
    report "no internal crate pins found in $CARGO"
    return
  fi
  for pin in $pins; do
    [ "$pin" = "$version" ] || report "internal pin =$pin should be =$version"
  done
}

check_json() {
  local file="$1" label="$2"
  [ -f "$file" ] || return 0
  local found
  found=$(json_version "$file")
  [ "$found" = "$version" ] || report "$label version $found should be $version"
}

check_lock() {
  [ -f "$DOCS_LOCK" ] || return 0
  local values
  values=$(sed -n 's/^[[:space:]]*"version": "\([^"]*\)".*/\1/p' "$DOCS_LOCK" | head -n 3)
  for found in $values; do
    [ "$found" = "$version" ] || report "docs/package-lock.json version $found should be $version"
  done
}

check_readme() {
  [ -f "$README" ] || return 0
  local found
  found=$(grep -oE 'xtalforge@[0-9]+\.[0-9]+\.[0-9]+' "$README" | head -n 1 | cut -d@ -f2)
  [ -z "$found" ] || [ "$found" = "$version" ] || report "README cargo version $found should be $version"
  found=$(grep -oE 'xtalforge==[0-9]+\.[0-9]+\.[0-9]+' "$README" | head -n 1 | sed 's/.*==//')
  [ -z "$found" ] || [ "$found" = "$version" ] || report "README pip version $found should be $version"
  found=$(grep -oE 'xtalforge-wasm@[0-9]+\.[0-9]+\.[0-9]+' "$README" | head -n 1 | cut -d@ -f2)
  [ -z "$found" ] || [ "$found" = "$version" ] || report "README npm version $found should be $version"
}

check_pyproject() {
  [ -f "$PYPROJECT" ] || return 0
  grep -q 'dynamic = \["version"\]' "$PYPROJECT" \
    || report "crates/xtalforge-python/pyproject.toml must keep a dynamic version"
}

write_pins() {
  sed -i -E "s#^(xtalforge(-[a-z]+)* = \{ version = \")=[^\"]+\"#\1=$version\"#" "$CARGO"
}

write_json_versions() {
  local file="$1" count="$2"
  [ -f "$file" ] || return 0
  python3 - "$file" "$version" "$count" <<'PY'
import re
import sys

path, version, count = sys.argv[1], sys.argv[2], int(sys.argv[3])
with open(path, encoding="utf-8") as stream:
    text = stream.read()
text, _ = re.subn(
    r'("version":\s*")[^"]*(")',
    lambda match: match.group(1) + version + match.group(2),
    text,
    count=count,
)
with open(path, "w", encoding="utf-8") as stream:
    stream.write(text)
PY
}

write_readme() {
  [ -f "$README" ] || return 0
  sed -i -E "s/(xtalforge@)[0-9]+\.[0-9]+\.[0-9]+/\1$version/g" "$README"
  sed -i -E "s/(xtalforge==)[0-9]+\.[0-9]+\.[0-9]+/\1$version/g" "$README"
  sed -i -E "s/(xtalforge-wasm@)[0-9]+\.[0-9]+\.[0-9]+/\1$version/g" "$README"
}

if [ "$CHECK" -eq 1 ]; then
  check_pins
  check_json "$DOCS_PKG" "docs/package.json"
  check_lock
  check_readme
  check_pyproject
  check_json "$WASM_PKG" "crates/xtalforge-wasm/pkg/package.json"
  if [ "$status" -eq 0 ]; then
    echo "OK: all versions are $version"
  fi
  exit "$status"
fi

write_pins
write_json_versions "$DOCS_PKG" 1
write_json_versions "$DOCS_LOCK" 3
write_readme
write_json_versions "$WASM_PKG" 1
check_pyproject

echo "sync-version: wrote version $version to the derived locations"
exit "$status"
