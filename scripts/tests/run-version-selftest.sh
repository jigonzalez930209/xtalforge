#!/usr/bin/env bash
# Self-test for scripts/sync-version.sh (Task 0.6.1).
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SYNC="$SCRIPT_DIR/../sync-version.sh"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
failures=0

mkdir -p "$tmp/docs" "$tmp/crates/xtalforge-python" "$tmp/crates/xtalforge-wasm/pkg"

cat > "$tmp/Cargo.toml" <<'EOF'
[workspace.package]
version = "1.2.3"

[workspace.dependencies]
xtalforge = { version = "=0.0.0", path = "crates/xtalforge" }
xtalforge-core = { version = "=0.0.0", path = "crates/xtalforge-core" }
xtalforge-io = { version = "=0.0.0", path = "crates/xtalforge-io" }
xtalforge-symmetry = { version = "=0.0.0", path = "crates/xtalforge-symmetry" }
xtalforge-ff = { version = "=0.0.0", path = "crates/xtalforge-ff" }
xtalforge-mof = { version = "=0.0.0", path = "crates/xtalforge-mof" }
EOF

cat > "$tmp/docs/package.json" <<'EOF'
{
  "name": "xtalforge-docs",
  "version": "0.0.0"
}
EOF

cat > "$tmp/docs/package-lock.json" <<'EOF'
{
  "name": "xtalforge-docs",
  "version": "0.0.0",
  "lockfileVersion": 3,
  "packages": {
    "": {
      "name": "xtalforge-docs",
      "version": "0.0.0"
    },
    "../crates/xtalforge-wasm/pkg": {
      "name": "xtalforge-wasm",
      "version": "0.0.0"
    }
  }
}
EOF

cat > "$tmp/README.md" <<'EOF'
cargo add xtalforge@0.0.0
pip install xtalforge==0.0.0
npm install xtalforge-wasm@0.0.0
EOF

cat > "$tmp/crates/xtalforge-python/pyproject.toml" <<'EOF'
[project]
name = "xtalforge"
dynamic = ["version"]
EOF

cat > "$tmp/crates/xtalforge-wasm/pkg/package.json" <<'EOF'
{
  "name": "xtalforge-wasm",
  "version": "0.0.0"
}
EOF

pass() { echo "ok - $1"; }
fail() {
  echo "FAIL - $1"
  failures=$((failures + 1))
}

if "$SYNC" "$tmp" > /dev/null; then pass "write mode runs"; else fail "write mode runs"; fi

grep -q 'version = "=1.2.3"' "$tmp/Cargo.toml" && pass "internal pins updated" \
  || fail "internal pins updated"
grep -q 'xtalforge-ff = { version = "=1.2.3"' "$tmp/Cargo.toml" && pass "all internal pins updated" \
  || fail "all internal pins updated"
grep -q '"version": "1.2.3"' "$tmp/docs/package.json" && pass "docs package updated" \
  || fail "docs package updated"
grep -q '"version": "1.2.3"' "$tmp/crates/xtalforge-wasm/pkg/package.json" \
  && pass "wasm package updated" || fail "wasm package updated"
grep -q 'xtalforge@1.2.3' "$tmp/README.md" && pass "README cargo snippet updated" \
  || fail "README cargo snippet updated"
grep -q 'xtalforge==1.2.3' "$tmp/README.md" && pass "README pip snippet updated" \
  || fail "README pip snippet updated"
grep -q 'xtalforge-wasm@1.2.3' "$tmp/README.md" && pass "README npm snippet updated" \
  || fail "README npm snippet updated"

if "$SYNC" --check "$tmp" > /dev/null; then pass "check mode accepts synced files"; \
  else fail "check mode accepts synced files"; fi

sed -i 's/"version": "1.2.3"/"version": "0.0.0"/' "$tmp/docs/package.json"
if "$SYNC" --check "$tmp" > /dev/null 2>&1; then
  fail "check mode rejects drift"
else
  pass "check mode rejects drift"
fi

if [ "$failures" -ne 0 ]; then
  echo "FAILED: $failures sync-version case(s)"
  exit 1
fi
echo "OK: sync-version self-test passed"
