# xtalforge Python plugin template

Install this package next to `xtalforge`, then load the registry:

```bash
pip install -e examples/python-plugin
python -c "import xtalforge.plugins; print(xtalforge.plugins.default().names_json())"
```

The entry point lives in the `xtalforge.plugins` group and returns a callable
that fills a `Plugins` registry. See `docs/api/interop.md` for the contract.
