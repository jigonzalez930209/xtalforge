# Core API

The `xtalforge::core` module is the source of truth for the other surfaces.

## `version() -> &'static str`

Returns the workspace version. The Python package exposes the same value as
`xtalforge.__version__` and the WASM package as `version()`.

## `cell_volume(a, b, c, alpha, beta, gamma) -> Result<f64, Error>`

Lengths in Angstrom, angles in degrees, result in cubic Angstrom. See the
[hello-world page](./hello-world) for the four surfaces and the shared
vectors.

## `Error`

One enum for the crate. `code()` returns a stable `snake_case` string and
`message()` the human-readable form.

| Variant | Code |
|---|---|
| `InvalidInput` | `invalid_input` |
| `DataCompression` | `data_compression` |
| `DataChecksum` | `data_checksum` |
| `DataDigest` | `data_digest` |
| `DataJson` | `data_json` |
| `LimitExceeded` | `limit_exceeded` |
| `Symmetry` | `symmetry` |

## Data bundles

`xtalforge_core::data` loads compressed tables:

- `EmbeddedBundle::new(name, gzip, sha256)` and `bytes()`: decompress once,
  verify the SHA-256 digest of the payload, cache the result.
- `JsonBundle<T>::new(bundle)` and `get()`: decode the payload once into `T`.
- `load_external_gzip(name, gzip, sha256)`: runtime bytes for WASM and large
  databases.
- `load_external_json<T>(name, gzip, sha256)`: the same with JSON decoding.

The policy and the generator command are in the [data guide](../guide/data).
