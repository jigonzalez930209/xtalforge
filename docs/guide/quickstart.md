# Quick start

## Rust

```rust
use xtalforge::io::read_cif;

let text = std::fs::read_to_string("structure.cif")?;
let structure = xtalforge::io::read_cif(&text)?;
println!("{}", structure.summary()?.formula);
let p1 = structure.p1()?;
let bonds = structure.bond_graph()?;
let rings = xtalforge::core::rings(&structure, 8)?;
```

Add the facade crate to `Cargo.toml`:

```toml
[dependencies]
xtalforge = "0.0.0"
```

## CLI

```bash
xtalforge info structure.cif
xtalforge --json info structure.cif
xtalforge cell --params 3 4 5 60 70 80
xtalforge element Fe
xtalforge convert structure.cif --format xyz
```

The binary reads a path or `-` for standard input and writes stdout, or a path
with `--out`. `--json` switches every command to a JSON object.

## Python

```python
import xtalforge

summary = xtalforge.structure_info(open("structure.cif").read())
info = xtalforge.element_info("Fe")
volume = xtalforge.cell_volume(3, 4, 5, 60, 70, 80)
```

Install the wheel with `pip install xtalforge` after the first release; the
local development wheel is built with `maturin build --release`.

## WebAssembly

```js
import init, { structure_info, version } from "xtalforge-wasm";
await init();
const summary = structure_info(cifText);
```

The browser demo in [Demos](/demos/) runs the same package offline. The
package is built with `wasm-pack build --release --target web`.

## Next

- [Concepts](./concepts)
- [File formats](./formats)
- [Data bundles](./data)
- [API reference](/api/)
