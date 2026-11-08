# Typed parameters

`ParamSpec` describes a value a module, engine or builder needs before it
runs: name, `ParamType`, default, help text, required flag and optional group.
`ParamSet` owns the descriptions and the current values.

```rust
use xtalforge::core::{ParamSet, ParamSpec, ParamType};
use serde_json::json;

let set = ParamSet::new().add(ParamSpec {
    name: "symprec".to_owned(),
    kind: ParamType::Float { min: Some(0.0), max: Some(1.0), unit: Some("Angstrom".to_owned()) },
    default: json!(0.001),
    help: "Position tolerance".to_owned(),
    required: false,
    group: Some("Symmetry".to_owned()),
})?;
```

`ParamType` is `Integer`, `Float`, `Boolean`, `Choice`, `Text`, `Vector3` or
`Path`. `validate` checks one value, `set` validates before storing,
`check_required` reports missing required parameters and `schema` returns a
JSON-schema-like object with the defaults and help text.

`ParamSet::to_json`/`from_json` persist and revalidate a parameter set.

The CLI, Python and WASM flag generation from the same description is Task
3.8.2; the schema is already the input for that generation.
