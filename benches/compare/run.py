"""Runs the comparative benchmark suite and writes deterministic JSON.

Usage:
    .venv/bin/python benches/compare/run.py --output results/v1.8.0.json
    .venv/bin/python benches/compare/run.py --only xtalforge --reps 5

The JSON contains the environment, versions, methodology, repetitions and one
record per task/library pair. Timing values are medians of the recorded samples;
everything else is stable across runs so diffs are readable.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

from suite import IMPORT_ERRORS, REPETITIONS, TASK_INFO, cases

ROOT = Path(__file__).resolve().parents[2]
SUITE_VERSION = 1
RELEASE = "1.8.0"
LIBRARIES = ["xtalforge", "gemmi", "spglib", "ase", "pymatgen"]


def _cpu_model() -> str:
    try:
        for line in Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
            if line.lower().startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def _memory_gib() -> float | None:
    try:
        for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("MemTotal:"):
                kib = float(line.split()[1])
                return round(kib / 1024 / 1024, 2)
    except (OSError, IndexError, ValueError):
        pass
    return None


def _package_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def _rustc_version() -> str | None:
    try:
        completed = subprocess.run(
            ["rustc", "--version"],
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError:
        return None
    if completed.returncode != 0:
        return None
    return completed.stdout.strip()


def _cli_version() -> str | None:
    cli = ROOT / "target" / "release" / "xtalforge"
    if not cli.is_file():
        return None
    completed = subprocess.run(
        [str(cli), "version"], text=True, capture_output=True, check=False
    )
    if completed.returncode != 0:
        return None
    return completed.stdout.strip().removeprefix("xtalforge ").strip()


def _environment() -> dict:
    import xtalforge

    versions = {
        "python": platform.python_version(),
        "numpy": _package_version("numpy"),
        "xtalforge": xtalforge.version(),
        "xtalforge_cli": _cli_version(),
        "gemmi": _package_version("gemmi"),
        "spglib": _package_version("spglib"),
        "ase": _package_version("ase"),
        "pymatgen": _package_version("pymatgen"),
        "rustc": _rustc_version(),
    }
    return {
        "cpu_model": _cpu_model(),
        "cpu_threads": os.cpu_count(),
        "memory_gib": _memory_gib(),
        "os": platform.platform(),
        "machine": platform.machine(),
        "versions": {key: value for key, value in versions.items() if value},
        "import_errors": IMPORT_ERRORS,
    }


def _selected(args: argparse.Namespace) -> list:
    from dataclasses import replace

    wanted_tasks = set(args.tasks.split(",")) if args.tasks else None
    wanted_libraries = set(args.only.split(",")) if args.only else None
    selected = []
    for case in cases():
        if wanted_tasks is not None and case.task not in wanted_tasks:
            continue
        if wanted_libraries is not None and case.library not in wanted_libraries:
            continue
        repetitions = args.reps if args.reps else case.repetitions
        selected.append(replace(case, repetitions=max(1, repetitions)))
    order = {task: index for index, task in enumerate(TASK_INFO)}
    return sorted(
        selected,
        key=lambda case: (order.get(case.task, 99), LIBRARIES.index(case.library)),
    )


def _record(case, measured: dict) -> dict:
    entry = {
        "task": case.task,
        "library": case.library,
        "median_ms": measured["median_ms"],
        "min_ms": measured["min_ms"],
        "max_ms": measured["max_ms"],
        "samples_ms": measured["samples_ms"],
        "repetitions": case.repetitions,
        "warmups": case.warmups,
        "output": measured["output"],
    }
    if case.external:
        entry["source"] = "rust-bench"
    return entry


def _summary(records: list[dict]) -> str:
    lines = [f"{'task':<12} {'library':<10} {'median ms':>10} {'reps':>5}  output"]
    for record in records:
        output = json.dumps(record["output"], sort_keys=True)
        lines.append(
            f"{record['task']:<12} {record['library']:<10} "
            f"{record['median_ms']:>10.3f} {record['repetitions']:>5}  {output}"
        )
    return "\n".join(lines)


def _document(args: argparse.Namespace, records: list[dict]) -> dict:
    unavailable = {
        task: info["unavailable"]
        for task, info in TASK_INFO.items()
        if info["unavailable"]
    }
    return {
        "schema_version": SUITE_VERSION,
        "suite": "xtalforge-compare",
        "release": RELEASE,
        "environment": _environment(),
        "methodology": {
            "statistic": "median",
            "timer": "time.perf_counter_ns in the same process",
            "warmups": 1,
            "repetitions": {
                task: (
                    args.reps if args.reps else REPETITIONS[task]
                )
                for task in TASK_INFO
            },
            "input": "the same in-memory CIF text is parsed by every library inside the timed region",
            "supercell": "in-memory build; the 10^6-atom xtalforge case invokes the release core benchmark",
            "notes": [
                "The reference numbers are end-to-end from CIF text and include that library's parser.",
                "xtalforge returns each periodic pair once; ASE and pymatgen return both directions.",
                "pymatgen's PXRD calculator returns merged peaks, xtalforge returns reflections.",
                "No installed reference package provides grid-based porosity; that row is xtalforge only.",
            ],
        },
        "tasks": TASK_INFO,
        "unavailable_references": unavailable,
        "results": records,
        "skipped": [
            {"task": task, "library": library, "reason": reason}
            for task, entries in unavailable.items()
            for library, reason in entries.items()
        ],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write JSON to this path")
    parser.add_argument("--only", help="comma-separated libraries")
    parser.add_argument("--tasks", help="comma-separated task ids")
    parser.add_argument("--reps", type=int, help="override repetitions for every case")
    parser.add_argument("--list", action="store_true", help="list tasks and exit")
    parser.add_argument("--quiet", action="store_true", help="no progress or table")
    args = parser.parse_args(argv)

    if args.list:
        for task, info in TASK_INFO.items():
            print(f"{task:<12} {info['title']}")
        return 0

    selected = _selected(args)
    if not selected:
        print("no benchmark cases selected", file=sys.stderr)
        return 2

    from suite import _measure

    if not args.quiet:
        print(_summary_header_for(selected), file=sys.stderr)
    records = []
    for index, case in enumerate(selected, start=1):
        measured = _measure(case)
        records.append(_record(case, measured))
        if not args.quiet:
            print(
                f"[{index}/{len(selected)}] {case.task}/{case.library} "
                f"{records[-1]['median_ms']:.3f} ms",
                file=sys.stderr,
            )

    document = _document(args, records)
    text = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
        if not args.quiet:
            print(_summary(records))
    else:
        print(text, end="")
    return 0


def _summary_header_for(selected: list) -> str:
    tasks = sorted({case.task for case in selected})
    return f"running {len(selected)} cases: {', '.join(tasks)}"


if __name__ == "__main__":
    raise SystemExit(main())
