#!/usr/bin/env python3
"""Golden extinction classes for the rutile pattern (Task 7.10.1).

Runs Crystal-Builder's indexing and rietx's extinction-class ranking on the
same rutile XY fixture the Rietveld golden uses. The class list is a property
of the tetragonal catalogue and the observed lines; the Rust port compares the
group partition and that P4_2/mnm is the leading class, not the reference's
absolute delta-BIC values.

    PYTHONPATH=/path/to/Crystal-Builder python scripts/golden/extinction.py
"""

import json
from pathlib import Path

from xtal.powder.data import PowderData, Radiation
from xtal.powder.index import IndexOptions, index
from xtal.powder.peaks import PeakOptions, fit_peaks

ROOT = Path(__file__).resolve().parents[2]
PATTERN = ROOT / "tests/golden/phase-7/rietveld/rutile.xy"
TARGET = ROOT / "tests/golden/phase-7/extinction/rutile-classes.json"


def main() -> None:
    data = PowderData.from_xy(str(PATTERN))
    peaks = fit_peaks(data, Radiation("cu"), PeakOptions())
    options = IndexOptions(
        bravais=frozenset({"tP"}),
        longest_axis=6.0,
        budget=20.0,
        zero_error=0.0,
        rank_groups=1,
    )
    result = index(peaks.for_indexing(), data, Radiation("cu"), options)
    row = result.rows[0]
    payload = {
        "source": "Crystal-Builder index + rietx determine_extinction_symbol",
        "cell": list(row.cell),
        "classes": [
            {
                "symbol": group.symbol,
                "representative": group.representative,
                "groups": list(group.space_groups),
                "delta_bic": group.delta_bic,
                "refuted": group.refuted,
                "conditions": list(group.conditions),
            }
            for group in row.classes
        ],
    }
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    TARGET.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"wrote {TARGET}: {len(payload['classes'])} classes")


if __name__ == "__main__":
    main()
