# Data bundles

Large tables, parameter sets copied from papers and databases never become
Rust literals. They live in data files, are compressed when large, and are
loaded through `xtalforge_core::data` (Task 0.2.3).

## When a table becomes a data file

Use a data file when any of these is true:

- the table has more than about 100 rows;
- the values come from a publication or a third-party database;
- the table is shared by more than one algorithm or surface;
- embedding it as Rust code would push a file toward the 500-line limit.

Small hand-written constants (basis names, unit strings, a handful of
thresholds) stay in code.

## Formats and location

| Format | Use |
|---|---|
| `.json` | small or medium tables that must stay readable in review |
| `.json.gz` | the default for large tables; reproducible compression |
| `.bin` | compact numeric arrays with a documented layout and endianness |

The canonical location is `data/` at the workspace root. A bundle that is
embedded with `include_bytes!` lives inside the crate that embeds it
(`crates/<crate>/data/`) so `cargo package` carries it; the root `data/`
directory holds canonical files, external bundles and data used by tests,
documentation and tooling. A copy that exists in both places must stay
byte-identical; the release tooling verifies it before publishing.

Every data file records its origin and license in `THIRD_PARTY.md`, and the
columns, units and units conversions at the top of the Rust module that reads
it. Regenerating a file is a reviewed change with the reason in the commit
message (Task 0.9.5).

## Embedded bundles

`EmbeddedBundle` takes the gzip bytes with `include_bytes!` and the SHA-256
digest of the **decompressed** payload. The first access decompresses and
verifies the payload; the result is cached in a `OnceLock`, so later accesses
cost one branch. `JsonBundle<T>` adds one-time JSON decoding into a typed
value.

```rust
use xtalforge_core::data::{EmbeddedBundle, JsonBundle};

static ELEMENTS: JsonBundle<Vec<ElementRecord>> = JsonBundle::new(EmbeddedBundle::new(
    "elements",
    include_bytes!("../data/elements.json.gz"),
    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
));
```

Compute the digest after compressing, from the uncompressed payload:

```bash
gzip -9 -n -c data/elements.json > data/elements.json.gz
gzip -dc data/elements.json.gz | sha256sum
```

The `-n` flag removes the original file name and timestamp, so the compressed
file is reproducible. The digest constant is part of the bundle declaration,
next to `include_bytes!`, and the test for the owning module calls `bytes()` or
`get()` so a corrupted or edited bundle fails the suite.

## External bundles

`load_external_gzip` and `load_external_json` accept bytes supplied by the
caller. Native code reads a file from `data/`; WebAssembly fetches the bytes
from the host or a URL. External bundles are not cached by the loader: the
caller keeps the decoded value in its own `OnceLock` or handle, which avoids a
hidden global cache in short-lived WASM calls.

Large optional databases, such as the MOF bundle (Task 6.11.3), are external
for every surface. The default build never downloads anything.

## Checksum scope

The checksum covers the decompressed payload, not the gzip stream. The gzip
container already has a CRC-32 for accidental corruption; SHA-256 catches an
edited or substituted payload and is stable across compression settings. It is
an integrity check, not a signature: bundles are reviewed in the pull request
that adds them.
