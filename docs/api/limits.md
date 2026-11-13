# Limits

Size checks are arithmetic, made before anything is allocated. Two levels:
over the **soft** limit the caller is warned; over the **hard** limit the
operation is refused and the message says what would fit.

| Profile | Soft | Hard |
|---|---|---|
| `standard` | table value | table value |
| `generous` | 2x | 2x |
| `huge` | million-atom work | 64 GB class, per-operation budgets |
| `warn_only` | table value | never refuses |

The table (Standard) is: drawn 100,000 / 400,000; supercell 50,000 / 200,000;
grid 1,500,000 / 4,000,000; carbon 10,000 / 30,000; nanotube and graphene
50,000 / 200,000; phonons 2,000 / 8,000; analytic Hessian 100 / 300.

## API

::: code-group

```rust
use xtalforge::core::{Limits, Profile, limits_bounds, limits_check};

let limits = Limits::new(Profile::Standard);
let verdict = xtalforge::core::limits::check_supercell(
    &limits, 10_000, [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
)?;
let (soft, hard) = limits_bounds("standard", "supercell")?;
let report = limits_check("standard", "supercell", 80_000)?;
```

```bash
xtalforge --limits standard cell --params 3 4 5 60 70 80
```

```python
soft, hard = xtalforge.limits_bounds("standard", "supercell")
report = xtalforge.limits_check("standard", "supercell", 80_000)
```

```js
import { limits_bounds, limits_check } from "xtalforge-wasm";
const bounds = limits_bounds("standard", "supercell");
const report = limits_check("standard", "supercell", 80000);
```

:::

## Estimates

- `count(value)` — formats a count as `52,000` or `1.4 M`.
- `drawn_atoms(n_atoms, ranges)` — cell atoms times the range volume.
- `multiplier(p)` and `supercell_atoms(n_atoms, p)` — `|det P|`.
- `grid_shape(lattice, spacing)` and `grid_points(lattice, spacing)` — from
  the cell lengths, minimum 8 points per axis.
- `check_drawn`, `check_supercell`, `check_grid`, `check_phonons` and
  `check_hessian` return a `Verdict { level, estimate, sentence }`.
- `Limits::bounds(budget)` returns `(soft, hard)`; `hard` is `None` under
  `WarnOnly`.
- `Verdict::refused` and `Verdict::warned` are the predicates; a refused
  operation returns `Error::LimitExceeded(verdict)` with code
  `limit_exceeded`.
- `Budget` names the estimate; `LimitLevel` is `Ok`, `Warn` or `Refuse`.
  `UNDO_ATOMS`, `UNDO_KEEP` and `HESSIAN_MIB_PER_ATOM` are the model constants.

The carbon, nanotube, graphene, phonon and Hessian budgets are part of the
table. The builders allocate through the generic supercell/grid paths; a
dedicated per-builder check is not exposed, so callers use
`check_supercell`, `check_grid`, `check_phonons` or `check_hessian` with the
matching budget.

## Surfaces

The CLI has the global `--limits {standard,generous,huge,warn-only}` flag. Python
takes `limits_check(profile, budget, estimate)`; WASM takes the same and
returns an object. WASM defaults to `standard` because browser memory is small.
