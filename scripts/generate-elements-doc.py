#!/usr/bin/env python3
"""Generate the element API page from the embedded reference table."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "crates" / "xtalforge-core" / "data" / "elements.json.gz"
OUT = ROOT / "docs" / "api" / "elements.md"


def fmt(value: object, digits: int = 3) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def main() -> None:
    data = json.loads(gzip.decompress(DATA.read_bytes()))
    elements = data["elements"]
    lines = [
        "# Elements",
        "",
        "Reference data for the periodic table, generated from",
        "`crates/xtalforge-core/data/elements.json.gz`. The bundle records its",
        "Crystal-Builder source and a SHA-256 checksum; `all()` verifies it on",
        "first use. `D` is deuterium and `X` is the dummy connection point.",
        "",
        "## API",
        "",
        "::: code-group",
        "",
        "```rust",
        "let iron = xtalforge::core::element_info(\"Fe\")?;",
        "let z = xtalforge::core::elements::atomic_number(\"Fe\")?;",
        "```",
        "",
        "```bash",
        "xtalforge element Fe",
        "xtalforge --json element 26",
        "```",
        "",
        "```python",
        "info = xtalforge.element_info(\"Fe\")",
        "info = xtalforge.element_info(26)",
        "```",
        "",
        "```js",
        "import { element_info } from \"xtalforge-wasm\";",
        "const info = element_info(\"Fe\");",
        "```",
        "",
        ":::",
        "",
        "The Rust module also exposes `element_by_z`, `element_by_symbol`,",
        "`parse_symbol`, `canonical_symbol`, `is_symbol`, `symbol_from_z`,",
        "`atomic_number`, `isotopes`, `isotope`, `mass_of`, `schema_version`",
        "and `data_source`, plus the `Element` and `Isotope` records.",
        "",
        "## Table",
        "",
        "| Z | Symbol | Name | Mass | Covalent | vdW | Pyykko single | Double | Aromatic | Valence | EN | Color |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for element in elements:
        color = "#{:02X}{:02X}{:02X}".format(*element["color"])
        lines.append(
            "| "
            + " | ".join(
                [
                    str(element["z"]),
                    element["symbol"],
                    element["name"],
                    fmt(element["mass"], 4),
                    fmt(element["covalent_radius"]),
                    fmt(element["vdw_radius"]),
                    fmt(element["pyykko_single"]),
                    fmt(element["pyykko_double"]),
                    fmt(element["pyykko_aromatic"]),
                    fmt(element["valence"]),
                    fmt(element["electronegativity"]),
                    color,
                ]
            )
            + " |"
        )
    lines.append("")
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)} ({len(elements)} rows)")


if __name__ == "__main__":
    main()
