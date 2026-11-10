#!/usr/bin/env bash
# Publishes every publishable workspace crate to crates.io in dependency order.
#
# crates.io has no atomic batch endpoint: every crate is uploaded separately.
# New crates have a burst of five, then one new crate every ten minutes, so
# this script waits for the retry time in the 429 response instead of hammering
# the index.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VERSION=$(sed -n 's/^version = "\([^"]*\)"/\1/p' "$ROOT/Cargo.toml" | head -n 1)
USER_AGENT="xtalforge-release/${VERSION} (https://github.com/jigonzalez930209/xtalforge)"

CRATES=(
  xtalforge-core
  xtalforge-symmetry
  xtalforge-io
  xtalforge-view
  xtalforge-commands
  xtalforge-ff
  xtalforge-analysis
  xtalforge-build
  xtalforge-mof
  xtalforge-modules
  xtalforge-agent
  xtalforge
  xtalforge-cli
)

published() {
  curl --fail --silent --show-error --user-agent "$USER_AGENT" \
    "https://crates.io/api/v1/crates/$1/$VERSION" >/dev/null 2>&1
}

retry_after() {
  local output=$1 retry target now delay
  retry=$(printf '%s\n' "$output" |
    sed -n 's/.*Please try again after \(.*\) GMT.*/\1/p' | head -n 1)
  if [ -z "$retry" ]; then
    echo 600
    return
  fi
  target=$(date -d "$retry GMT" +%s 2>/dev/null || echo 0)
  now=$(date +%s)
  delay=$((target - now + 30))
  if [ "$delay" -lt 30 ]; then
    delay=600
  fi
  echo "$delay"
}

for crate in "${CRATES[@]}"; do
  if published "$crate"; then
    echo "${crate} ${VERSION} is already published; skipping"
    continue
  fi

  echo "::group::${crate} ${VERSION}"
  while true; do
    set +e
    output=$(cargo publish --locked --no-verify -p "$crate" 2>&1)
    status=$?
    set -e
    printf '%s\n' "$output"

    if [ "$status" -eq 0 ]; then
      break
    fi

    if printf '%s\n' "$output" | grep -qE 'already (exists|uploaded)|crate version .* is already'; then
      echo "${crate} ${VERSION} is already published; skipping"
      break
    fi

    if printf '%s\n' "$output" | grep -qiE '429|too many|rate limit|try again after'; then
      delay=$(retry_after "$output")
      echo "crates.io rate limit: sleeping ${delay}s before retrying ${crate}"
      sleep "$delay"
      continue
    fi

    echo "could not publish ${crate}" >&2
    exit "$status"
  done
  echo "::endgroup::"

  # Let the sparse index expose the new version to the next dependency.
  sleep 30
done
