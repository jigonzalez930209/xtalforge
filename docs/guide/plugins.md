# Python plugin template

A Python package can extend the in-process registry without changing
xtalforge. Declare an entry point in the `xtalforge.plugins` group:

```toml
[project.entry-points."xtalforge.plugins"]
demo = "my_package:register"
```

```python
def register(registry):
    registry.register_module("demo", lambda parameters: '{"ok": true}')
```

Then:

```python
import xtalforge.plugins

registry = xtalforge.plugins.default()
print(registry.names_json())
```

The repository ships a working template in `examples/python-plugin`: one
module and one `.demo` format. Python entry points are in-process; the Rust
plugin registries are the extension point for the CLI and the agent.
