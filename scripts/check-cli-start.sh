#!/usr/bin/env bash
# Measures the cold start of `xtalforge info` (Task 15.3.3).
#
# The roadmap budget is 20 ms for the whole process: `fork`/`exec`, the dynamic
# loader, clap parsing, reading and parsing one small CIF, and exit. The script
# reports the first invocation separately (it pays the page-cache cost of a
# fresh process image) and then the median, minimum, mean and p95 of RUNS
# repetitions. Only the median is compared with the ceiling; the p95 is
# reported so a noisy machine is visible.
#
# Reproducibility:
#   * the release binary is built from the current tree before measuring;
#   * the fixture is the tracked docs sample docs/public/samples/Si.cif;
#   * timings use the bash `EPOCHREALTIME` builtin, so no helper process runs
#     inside the measured interval;
#   * RUNS (default 30) and CEILING_MS (default 20) are environment variables.
#
# A true cold page cache would need root (`drop_caches`); set COLD=1 to try,
# and the script falls back to a warm cache with a note when it cannot.
#
# Usage: scripts/check-cli-start.sh [CEILING_MS]
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT" || exit 1

BIN="$ROOT/target/release/xtalforge"
SAMPLE="$ROOT/docs/public/samples/Si.cif"
RUNS="${RUNS:-30}"
CEILING_MS="${1:-${CEILING_MS:-20}}"

case "$RUNS" in
  ''|*[!0-9]*) echo "FAIL: RUNS must be a positive integer" >&2; exit 2 ;;
esac
if [ "$RUNS" -lt 5 ]; then RUNS=5; fi
if [ ! -f "$SAMPLE" ]; then
  echo "FAIL: fixture $SAMPLE is missing" >&2
  exit 2
fi

echo "== building the release CLI"
if ! cargo build --release -q -p xtalforge-cli; then
  echo "FAIL: could not build xtalforge-cli" >&2
  exit 2
fi

if [ "${COLD:-0}" = "1" ]; then
  if [ -w /proc/sys/vm/drop_caches ] && sync && echo 3 > /proc/sys/vm/drop_caches 2>/dev/null; then
    echo "== page cache dropped"
  else
    echo "NOTE: cannot drop the page cache without root; measuring with a warm cache"
  fi
fi

# Milliseconds between two EPOCHREALTIME values.
elapsed_ms() {
  awk -v start="$1" -v end="$2" 'BEGIN { printf "%.3f", (end - start) * 1000 }'
}

first_start=$EPOCHREALTIME
"$BIN" info "$SAMPLE" >/dev/null 2>&1
first_status=$?
first_end=$EPOCHREALTIME
if [ "$first_status" -ne 0 ]; then
  echo "FAIL: xtalforge info exited with status $first_status" >&2
  exit 1
fi
first=$(elapsed_ms "$first_start" "$first_end")

values=()
for _ in $(seq 1 "$RUNS"); do
  start=$EPOCHREALTIME
  "$BIN" info "$SAMPLE" >/dev/null 2>&1
  status=$?
  end=$EPOCHREALTIME
  if [ "$status" -ne 0 ]; then
    echo "FAIL: xtalforge info exited with status $status" >&2
    exit 1
  fi
  values+=("$(elapsed_ms "$start" "$end")")
done

summary=$(printf '%s\n' "${values[@]}" | sort -n | awk '
  { samples[NR] = $1; total += $1 }
  END {
    n = NR;
    median = (n % 2) ? samples[(n + 1) / 2] : (samples[n / 2] + samples[n / 2 + 1]) / 2;
    p95_index = int((95 * n + 99) / 100);
    if (p95_index < 1) p95_index = 1;
    if (p95_index > n) p95_index = n;
    printf "%.3f %.3f %.3f %.3f", median, samples[1], total / n, samples[p95_index];
  }
')
read -r median minimum mean p95 <<<"$summary"

echo "first run: ${first} ms"
printf 'xtalforge info: median %s ms (min %s, mean %s, p95 %s; %s runs)\n' \
  "$median" "$minimum" "$mean" "$p95" "$RUNS"

if awk -v value="$median" -v ceiling="$CEILING_MS" 'BEGIN { exit !(value > ceiling) }'; then
  echo "FAIL: median ${median} ms exceeds the ${CEILING_MS} ms cold-start budget" >&2
  exit 1
fi
echo "OK: median ${median} ms (budget ${CEILING_MS} ms)"
