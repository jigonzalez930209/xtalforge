# Hello-world API

`version()` and `cell_volume()` are implemented in all four surfaces and are
tested against the same vectors in `tests/vectors/hello_world.json`.

## Version

::: code-group

```rust
let version: &str = xtalforge::core::version();
```

```bash
xtalforge version
xtalforge --json version
```

```python
import xtalforge
version: str = xtalforge.version()
```

```js
import { version } from "xtalforge-wasm";
const value = version();
```

:::

## Cell volume

`cell_volume(a, b, c, alpha, beta, gamma)` takes lengths in Angstrom and
angles in degrees and returns cubic Angstrom.

::: code-group

```rust
let volume = xtalforge::core::cell_volume(3.0, 4.0, 5.0, 60.0, 70.0, 80.0)?;
```

```bash
xtalforge cell-volume 3 4 5 60 70 80
xtalforge --json cell-volume 3 4 5 60 70 80
```

```python
volume = xtalforge.cell_volume(3.0, 4.0, 5.0, 60.0, 70.0, 80.0)
```

```js
import { cell_volume } from "xtalforge-wasm";
const volume = cell_volume(3, 4, 5, 60, 70, 80);
```

:::

## Errors

Every surface reports the stable `invalid_input` code for a non-positive
length, an angle outside `(0, 180)` degrees, a non-finite value, or a
degenerate cell.

- Rust: `Err(xtalforge::core::Error::InvalidInput(..))`, `error.code()`.
- CLI: exit code `3` and, with `--json`, an object under `error`.
- Python: `xtalforge.InputError`, a subclass of `xtalforge.XtalforgeError`;
  `args[0]` is the code and `args[1]` the message.
- WASM: a JavaScript `Error` whose `code` property is `invalid_input`.
