#!/usr/bin/env bash
# Regenerates the embedded RCSR fingerprint index (Task 5.6.4).
#
# Usage: scripts/golden/rcsr_index.sh
# Requires the bundled snapshots in data/rcsr/ and a release build of the
# analysis crate. The generator rebuilds every keyed entry from its canonical
# key and fails on any fingerprint mismatch.
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
json=$(mktemp "${TMPDIR:-/tmp}/rcsr-XXXXXX.json")
trap 'rm -f "$json"' EXIT

gzip -dc "$root/data/rcsr/rcsr-2019-06-01.json.gz" > "$json"
cargo run --release --manifest-path "$root/Cargo.toml" \
  -p xtalforge-analysis --example build_rcsr_index -- \
  "$json" "$root/crates/xtalforge-analysis/data/rcsr_index.bin"
