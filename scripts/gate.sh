#!/usr/bin/env bash
# Local quality gate (Task 17.4.1).
#
# One command that runs the checks the release checklist requires (AGENTS.md
# section 5.4) and prints a summary table. It replaces the removed CI workflows
# at no runner cost; the maintainer runs it on the exact commit that will be
# tagged.
#
# Usage:
#   scripts/gate.sh            full gate (default)
#   scripts/gate.sh --quick    skip the WASM build/tests, Python wheel, docs
#                              build and benchmark regression
#   scripts/gate.sh --list     print the commands without running them
#
# Every step runs even after a failure, so the table shows the complete state.
# The exit status is 1 when any step failed and 0 only when all of them passed.
set -u

root=$(git rev-parse --show-toplevel 2>/dev/null || pwd)
cd "$root" || exit 2

quick=0
list_only=0
for argument in "$@"; do
  case "$argument" in
    --quick) quick=1 ;;
    --list) list_only=1 ;;
    -h|--help)
      sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *)
      echo "gate.sh: unknown argument: $argument" >&2
      exit 2
      ;;
  esac
done

python_command() {
  if [ -x "$root/.venv/bin/python" ] && "$root/.venv/bin/python" -c 'import mypy' >/dev/null 2>&1; then
    printf '%s' "$root/.venv/bin/python"
  elif command -v python3 >/dev/null 2>&1; then
    printf 'python3'
  else
    printf 'python'
  fi
}

names=()
commands=()
add() {
  names+=("$1")
  commands+=("$2")
}

add "rules" "scripts/check-rules.sh"
add "version" "scripts/sync-version.sh --check"
add "agent-version" "scripts/check-agent-version.sh"
add "parity" "scripts/check-parity.sh"
add "api-docs" "scripts/check-api-docs.sh"
add "format" "cargo fmt --check"
add "clippy" "cargo clippy --workspace --all-targets --all-features -- -D warnings"
add "rustdoc" "RUSTDOCFLAGS='-D warnings' cargo doc --workspace --no-deps"
add "rust-tests" "cargo test --workspace --all-features"
add "deny" "cargo deny check"
add "wasm-check" "cargo check -p xtalforge-wasm --target wasm32-unknown-unknown"
add "cli-start" "scripts/check-cli-start.sh"
add "surfaces" "scripts/check-surfaces.sh"
if [ "$quick" -eq 0 ]; then
  add "wasm-packages" "scripts/build-wasm-packages.sh"
  add "wasm-size" "scripts/check-wasm-size.sh"
  add "wasm-tests" "(cd crates/xtalforge-wasm && wasm-pack test --node --release)"
  add "python-wheel" "(cd crates/xtalforge-python && maturin build --release)"
  add "python-types" "$(python_command) -m mypy --strict crates/xtalforge-python/python/xtalforge"
  add "docs-build" "(cd docs && npm run docs:build)"
  add "bench" "scripts/check-bench.sh"
fi

if [ "$list_only" -eq 1 ]; then
  for index in "${!names[@]}"; do
    printf '%-16s %s\n' "${names[$index]}" "${commands[$index]}"
  done
  exit 0
fi

results=()
durations=()
failures=0
started_all=$(date +%s)

printf 'xtalforge local gate\n'
printf 'commit: %s\n' "$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
printf 'mode:   %s\n\n' "$([ "$quick" -eq 1 ] && echo quick || echo full)"

for index in "${!names[@]}"; do
  name=${names[$index]}
  command=${commands[$index]}
  printf '[%d/%d] %s\n' "$((index + 1))" "${#names[@]}" "$name"
  started=$(date +%s)
  if bash -c "$command" >"/tmp/xtalforge-gate-${name}.log" 2>&1; then
    results+=("PASS")
  else
    results+=("FAIL")
    failures=$((failures + 1))
  fi
  finished=$(date +%s)
  durations+=("$((finished - started))")
  printf '      %s in %ss (log: %s)\n' "${results[$index]}" "${durations[$index]}" "/tmp/xtalforge-gate-${name}.log"
done

total=$(( $(date +%s) - started_all ))
printf '\nsummary\n'
printf '%-16s %-6s %s\n' "step" "status" "seconds"
for index in "${!names[@]}"; do
  printf '%-16s %-6s %s\n' "${names[$index]}" "${results[$index]}" "${durations[$index]}"
done
printf '%-16s %-6s %s\n' "total" "$([ "$failures" -eq 0 ] && echo PASS || echo FAIL)" "$total"

if [ "$failures" -ne 0 ]; then
  printf '\n%d step(s) failed; inspect the logs above.\n' "$failures" >&2
  exit 1
fi
printf '\nAll gate steps passed.\n'
