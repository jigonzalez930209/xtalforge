# Plugins and extension points

xtalforge exposes open registries instead of a global plugin loader. A plugin
is a crate (or a Python object) that registers itself with the registry the
caller owns. There is no hidden mutable state: a server, a test and a front end
each build their own registries.

## Rust: formats, calculators, modules, scalars and styles

| Registry | Crate | Registers |
|---|---|---|
| `FormatRegistry` | `xtalforge-io` | a reader and/or writer with a name and extensions |
| `CalculatorRegistry` | `xtalforge-ff` | an energy engine implementing `Calculator` |
| `ModuleRegistry` | `xtalforge-modules` | a `Module` with JSON parameters and result |
| `ScalarRegistry` | `xtalforge-core` | a named `ScalarFn` and its `ScalarInfo` |
| `StyleRegistry` | `xtalforge-view` | a named `ViewSettings` preset |

Every registry starts from the built-ins (`new`) or empty (`empty`), rejects a
duplicate name with a typed error, and exposes its entries for a front end to
list.

The `FormatRegistry` mirrors the static `FORMATS` table:

```rust
use xtalforge::io::registry::{Capabilities, FormatRegistry, FormatSpec};

let mut formats = FormatRegistry::new();
formats.register(FormatSpec {
    name: "two-line",
    extensions: &[".2l"],
    read: Some(my_reader),
    write: Some(my_writer),
    capabilities: Capabilities { cell: true, ..Capabilities::default() },
})?;
# Ok::<(), xtalforge::core::Error>(())
```

The scalar registry carries the four built-ins (`bond_lengths`,
`coordination_numbers`, `atom_charges`, `atom_masses`); `ScalarRegistry::register`
adds a plugin function and `compute` runs it.

The draw style itself is a closed enum: a renderer must agree on the geometry
it draws. What is open is the settings preset (radii, colours, cell,
hydrogens), which is what `StyleRegistry` stores.

A complete, buildable example is `examples/plugin`:

```bash
cargo test -p xtalforge-plugin-example
```

It registers a two-line format, a zero-energy calculator, a site-count module,
a scalar and a large-atom style, and checks every registry.

The full API is in the rustdoc of `xtalforge_io::registry`,
`xtalforge_ff::CalculatorRegistry`, `xtalforge_modules::ModuleRegistry`,
`xtalforge_core::ScalarRegistry` and `xtalforge_view::StyleRegistry`.

## Python: callables

```python
import json
import xtalforge

plugins = xtalforge.Plugins()
plugins.register_module("echo", lambda parameters: json.dumps({"echo": json.loads(parameters)}))
plugins.register_calculator(
    "zero",
    lambda structure: json.dumps({"energy": 0.0,
                                  "n_atoms": len(json.loads(structure)["structure"]["sites"])}),
)
plugins.register_format("two-line", [".2l"], reader, writer)

print(plugins.names_json())
print(plugins.run_module("echo", '{"n": 2}'))
```

The contract is JSON in, JSON out:

| Kind | Callable |
|---|---|
| module | `parameters_json -> result_json` |
| calculator | `structure_json -> {"energy": ..., "forces": [...]}` |
| format reader | `text -> structure_json` |
| format writer | `structure_json -> text` |

The GIL is released around the Rust work and reacquired to call the plugin.
A Python exception propagates unchanged; a non-string return is a `TypeError`;
a duplicate or empty name is an `InputError` with code `invalid_input`.

A Python plugin is visible to the Python code that owns the `Plugins` object.
It is not inserted into the Rust registries, so the CLI, MCP and the optimizer
do not see it. A plugin that must be visible there is a Rust crate; this is a
deliberate boundary, not a missing feature.

## Python entry points

`xtalforge.plugins` discovers the `xtalforge.plugins` entry-point group. Each
entry point returns a `Plugins` registry, a callable that fills one, or a
mapping with `modules`, `calculators` and `formats`:

```python
import xtalforge.plugins

registry = xtalforge.plugins.default()
print(registry.names_json())
```

The repository includes the installable template in `examples/python-plugin`;
see [Python plugins](../guide/plugins.md) for the package layout. The entry
points remain in-process: the CLI and agent use the Rust registries above.
