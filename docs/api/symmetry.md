# Symmetry

`xtalforge-symmetry` wraps `moyo` behind our own types. A search returns a
`SymmetryDataset` with the space-group number, Hall number, H-M symbol, Pearson
symbol, operations in the input cell, orbits, Wyckoff letters, site-symmetry
symbols and the standardized cell.

## API

::: code-group

```rust
use xtalforge::symmetry::{SymmetryOptions, detect_symmetry};

let dataset = detect_symmetry(&structure, SymmetryOptions::default())?;
println!("{} ({})", dataset.hm_symbol, dataset.number);
println!("{:?}", dataset.wyckoffs);
```

```bash
xtalforge symmetry detect structure.cif --symprec 0.001
xtalforge --json symmetry detect structure.cif
```

```python
dataset = xtalforge.detect_symmetry(text, symprec=0.001)
print(dataset["number"], dataset["wyckoffs"])
```

```js
import { detect_symmetry } from "xtalforge-wasm";
const dataset = detect_symmetry(cifText, 0.001);
console.log(dataset.number, dataset.wyckoffs);
```

:::

## Options

`SymmetryOptions { symprec, angle_tolerance, setting }`:

- `symprec` is a position tolerance in Angstrom; the default is `1e-5` and the
  CLI/Python default is `1e-3` to match Crystal-Builder's corpus options.
- `angle_tolerance` is in radians; `None` uses moyo's symprec-driven criterion.
- `setting` is `Spglib` (default), `Standard` or `HallNumber(n)` for a specific
  Hall setting.

The structure is expanded to P1 first; the input is never modified. Atom types
are `(element, occupancy)` pairs rounded to three decimals, so a mixed site is
not equivalent to a pure one and dummy `X` has its own type.

## Space groups

`SpaceGroup::from_number(225)`, `SpaceGroup::from_hall_number(525)` and
`SpaceGroup::parse("Fm-3m")` accept numbers, H-M symbols with optional setting
(`Fd-3m:1`), Hall symbols and plain numbers. The value exposes `number`,
`hall_number`, `hm_short`, `hm_full`, `hall_symbol`, `setting`, `centering`,
`crystal_system`, `lattice_system`, `operations`, `order`,
`is_centrosymmetric` and `is_chiral`. `all_settings()` lists the 530 Hall
settings.

`WyckoffPosition::list(hall_number)` returns multiplicity, letter,
site-symmetry symbol and coordinate template for every position of a setting.

## Errors

A failed search returns `Error::Symmetry` with code `symmetry` and a message
that includes the tolerance. Unknown symbols return the same code.

## Corpus note

The small-corpus test compares the detected number with the Crystal-Builder
golden header for the `simple/` structures. The full corpus sweep is an
ignored test because the large MOFs take minutes; six files are recorded as
exceptions on the algorithm page.
