#!/usr/bin/env bash
# Self-test for scripts/check-rules.sh (Tasks 0.2.1 and the warning policy).
#
# Each case builds a small tree in a temporary directory, runs the rule check
# against it and asserts both the exit status and the message.
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RULES="$SCRIPT_DIR/../check-rules.sh"
FIXTURES="$SCRIPT_DIR/rules"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
failures=0

make_lines() {
  local path="$1" count="$2" i=1
  : > "$path"
  while [ "$i" -le "$count" ]; do
    printf '// fixture line %d\n' "$i" >> "$path"
    i=$((i + 1))
  done
}

pass_case() {
  local label="$1" root="$2" out="$tmp/pass.log"
  if "$RULES" "$root" > "$out" 2>&1; then
    echo "ok - $label"
  else
    echo "FAIL - $label (expected success)"
    sed 's/^/    /' "$out"
    failures=$((failures + 1))
  fi
}

fail_case() {
  local label="$1" root="$2" needle="$3" out="$tmp/fail.log"
  if "$RULES" "$root" > "$out" 2>&1; then
    echo "FAIL - $label (expected failure)"
    failures=$((failures + 1))
  elif ! grep -q -F -- "$needle" "$out"; then
    echo "FAIL - $label (missing message: $needle)"
    sed 's/^/    /' "$out"
    failures=$((failures + 1))
  else
    echo "ok - $label"
  fi
}

# 1. Ordinary docs and an exactly 500-line code file pass.
mkdir -p "$tmp/pass"
cp "$FIXTURES/pass_docs.md.fixture" "$tmp/pass/conventions.md"
make_lines "$tmp/pass/exactly_500.rs" 500
pass_case "ordinary docs and exactly 500 lines" "$tmp/pass"

# 2. A 501-line code file fails with the line count in the message.
mkdir -p "$tmp/long"
make_lines "$tmp/long/too_long.rs" 501
fail_case "501-line code file is rejected" "$tmp/long" "TOO LONG (501 lines)"

# 3. A dead-code allow attribute fails.
mkdir -p "$tmp/dead_code"
cp "$FIXTURES/fail_allow_dead_code.rs.fixture" "$tmp/dead_code/suppressed.rs"
fail_case "allow(dead_code) is rejected" "$tmp/dead_code" \
  "SUPPRESSED WARNING in suppressed.rs"

# 4. An unused-variable expect attribute fails.
mkdir -p "$tmp/unused"
cp "$FIXTURES/fail_expect_unused.rs.fixture" "$tmp/unused/suppressed.rs"
fail_case "expect(unused_variables) is rejected" "$tmp/unused" \
  "SUPPRESSED WARNING in suppressed.rs"

# 5. A file with exactly 500 lines is the boundary: one line fewer passes.
mkdir -p "$tmp/boundary"
make_lines "$tmp/boundary/exactly_500.rs" 500
pass_case "the 500-line boundary is inclusive" "$tmp/boundary"

if [ "$failures" -ne 0 ]; then
  echo "FAILED: $failures rule self-test case(s)"
  exit 1
fi
echo "OK: rules self-test passed"
