"""Compares a benchmark result document with the stored baseline.

The baseline tracks only xtalforge rows, so the file stays small and a
regression in a reference library cannot fail the local gate. A row regresses
when its fastest sample (minimum) is more than `threshold_percent` above the
baseline minimum; missing rows also fail. The minimum is used instead of the
median because CPU-frequency changes move the median by 20-30 % on a laptop,
while the fastest sample stays within a few percent between runs.

Examples:
    python compare.py --baseline baseline.json --current latest.json
    python compare.py --baseline baseline.json --current run1.json --current run2.json
    python compare.py --baseline baseline.json --current latest.json --update
    python compare.py --self-test

With two `--current` documents the comparison keeps the fastest sample for each
row, which lets the local gate retry once when a shared machine is busy.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SCHEMA_VERSION = 1
THRESHOLD_PERCENT = 10.0


def _key(record: dict) -> str:
    return f"{record['task']}/{record['library']}"


def _baseline_entries(document: dict) -> dict[str, dict]:
    entries: dict[str, dict] = {}
    for record in document.get("results", []):
        if record.get("library") != "xtalforge":
            continue
        if "median_ms" not in record:
            continue
        entries[_key(record)] = {
            "median_ms": record["median_ms"],
            "min_ms": record.get("min_ms", record["median_ms"]),
            "max_ms": record.get("max_ms", record["median_ms"]),
            "repetitions": record.get("repetitions"),
            "output": record.get("output", {}),
        }
    return dict(sorted(entries.items()))


def _baseline_document(document: dict) -> dict:
    return {
        "schema_version": SCHEMA_VERSION,
        "suite": "xtalforge-compare",
        "release": document.get("release"),
        "threshold_percent": THRESHOLD_PERCENT,
        "environment": {
            "cpu_model": document.get("environment", {}).get("cpu_model"),
            "cpu_threads": document.get("environment", {}).get("cpu_threads"),
            "memory_gib": document.get("environment", {}).get("memory_gib"),
            "os": document.get("environment", {}).get("os"),
            "versions": document.get("environment", {}).get("versions", {}),
        },
        "benchmarks": _baseline_entries(document),
    }


def _current_entries(document: dict) -> dict[str, dict]:
    entries: dict[str, dict] = {}
    for record in document.get("results", []):
        if "median_ms" in record:
            entries[_key(record)] = record
    return entries


def _merge_documents(documents: list[dict]) -> dict:
    """Merges repeated runs, keeping the fastest sample for each row."""
    merged = dict(documents[0])
    best: dict[str, dict] = {}
    for document in documents:
        for record in document.get("results", []):
            if "median_ms" not in record:
                continue
            key = _key(record)
            previous = best.get(key)
            if previous is None or float(record.get("min_ms", record["median_ms"])) < float(
                previous.get("min_ms", previous["median_ms"])
            ):
                best[key] = record
    merged["results"] = list(best.values())
    return merged


def compare(
    baseline: dict, current: dict, threshold: float
) -> tuple[bool, list[str], list[str]]:
    """Returns (ok, lines, failures) for a baseline/current pair."""
    recorded = _current_entries(current)
    lines: list[str] = []
    failures: list[str] = []
    for key, expected in sorted(baseline.get("benchmarks", {}).items()):
        actual = recorded.get(key)
        if actual is None:
            message = f"MISSING {key}: no current benchmark row"
            lines.append(message)
            failures.append(message)
            continue
        base = float(expected.get("min_ms", expected["median_ms"]))
        value = float(actual.get("min_ms", actual["median_ms"]))
        if base <= 0:
            message = f"INVALID {key}: baseline minimum must be positive"
            lines.append(message)
            failures.append(message)
            continue
        change = (value / base - 1.0) * 100.0
        status = "FAIL" if change > threshold else "OK"
        lines.append(
            f"{status:<4} {key:<20} {base:9.3f} -> {value:9.3f} ms "
            f"({change:+.1f}%)"
        )
        if change > threshold:
            failures.append(
                f"{key} regressed {change:.1f}% (limit {threshold:.1f}%)"
            )
    return not failures, lines, failures


def _self_test() -> int:
    baseline = {
        "benchmarks": {
            "a/xtalforge": {"median_ms": 10.0, "min_ms": 10.0},
            "b/xtalforge": {"median_ms": 20.0, "min_ms": 20.0},
            "c/xtalforge": {"median_ms": 30.0, "min_ms": 30.0},
        }
    }
    current = {
        "results": [
            {"task": "a", "library": "xtalforge", "median_ms": 10.5, "min_ms": 10.2},
            {"task": "b", "library": "xtalforge", "median_ms": 20.1, "min_ms": 19.8},
            {"task": "c", "library": "xtalforge", "median_ms": 33.1, "min_ms": 34.0},
        ]
    }
    ok, lines, failures = compare(baseline, current, 10.0)
    assert not ok
    assert len(failures) == 1, failures
    assert "c/xtalforge" in failures[0]
    assert any("MISSING" in line for line in lines) is False
    assert len(lines) == 3

    improved = {
        "results": [
            {"task": "a", "library": "xtalforge", "median_ms": 5.0},
            {"task": "b", "library": "xtalforge", "median_ms": 20.0},
            {"task": "c", "library": "xtalforge", "median_ms": 30.0},
        ]
    }
    ok, _, failures = compare(baseline, improved, 10.0)
    assert ok and not failures

    missing = {"results": []}
    ok, _, failures = compare(baseline, missing, 10.0)
    assert not ok and len(failures) == 3

    document = {
        "release": "1.8.0",
        "environment": {"cpu_model": "test", "versions": {"xtalforge": "1.7.0"}},
        "results": [
            {"task": "a", "library": "xtalforge", "median_ms": 10.0},
            {"task": "b", "library": "gemmi", "median_ms": 1.0},
        ],
    }
    stored = _baseline_document(document)
    assert list(stored["benchmarks"]) == ["a/xtalforge"]
    assert stored["threshold_percent"] == 10.0

    slower = {
        "results": [
            {"task": "a", "library": "xtalforge", "median_ms": 12.0, "min_ms": 12.0},
        ]
    }
    faster = {
        "results": [
            {"task": "a", "library": "xtalforge", "median_ms": 9.0, "min_ms": 9.0},
        ]
    }
    merged = _merge_documents([slower, faster])
    assert merged["results"][0]["min_ms"] == 9.0
    print("self-test: OK")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=Path("baseline.json"))
    parser.add_argument("--current", type=Path, action="append", required=False)
    parser.add_argument("--update", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--threshold", type=float, default=None)
    args = parser.parse_args(argv)

    if args.self_test:
        return _self_test()
    if args.current is None:
        parser.error("--current is required unless --self-test is used")

    current = _merge_documents(
        [json.loads(path.read_text(encoding="utf-8")) for path in args.current]
    )
    if args.update:
        args.baseline.write_text(
            json.dumps(_baseline_document(current), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"updated baseline: {args.baseline}")
        return 0

    if not args.baseline.is_file():
        print(
            f"missing baseline {args.baseline}; run scripts/check-bench.sh --update",
            file=sys.stderr,
        )
        return 2
    baseline = json.loads(args.baseline.read_text(encoding="utf-8"))
    threshold = (
        args.threshold
        if args.threshold is not None
        else float(baseline.get("threshold_percent", THRESHOLD_PERCENT))
    )
    hardware = baseline.get("environment", {})
    observed = current.get("environment", {})
    if hardware.get("cpu_model") and hardware["cpu_model"] != observed.get("cpu_model"):
        print(
            "warning: baseline CPU differs from the current CPU; "
            "timing comparison is advisory",
            file=sys.stderr,
        )
    ok, lines, failures = compare(baseline, current, threshold)
    for line in lines:
        print(line)
    if not ok:
        print("regression gate: FAILED", file=sys.stderr)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr)
        return 1
    print(f"regression gate: OK (limit {threshold:.1f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
