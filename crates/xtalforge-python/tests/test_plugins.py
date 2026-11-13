"""Python-side plugins (Task 8.3.2)."""

import json

import pytest
import xtalforge


def test_module_and_calculator_plugins() -> None:
    plugins = xtalforge.Plugins()
    plugins.register_module("echo", lambda parameters: json.dumps({"echo": json.loads(parameters)}))
    plugins.register_calculator(
        "zero",
        lambda structure: json.dumps({"energy": 0.0, "n_atoms": len(json.loads(structure)["structure"]["sites"])}),
    )
    names = json.loads(plugins.names_json())
    assert names["modules"] == ["echo"]
    assert names["calculators"] == ["zero"]
    assert json.loads(plugins.run_module("echo", '{"n": 2}')) == {"echo": {"n": 2}}
    structure = xtalforge.Structure("2\nNaCl\nNa 0 0 0\nCl 0.5 0.5 0.5\n", "xyz").to_json()
    assert json.loads(plugins.calculator_energy("zero", structure))["energy"] == 0.0


def test_format_plugins() -> None:
    plugins = xtalforge.Plugins()

    def reader(text: str) -> str:
        assert text.startswith("cell ")
        return xtalforge.Structure("1\nNa\nNa 0 0 0\n", "xyz").to_json()

    def writer(structure: str) -> str:
        assert "structure" in json.loads(structure)
        return "cell 4.0\nNa 0 0 0\n"

    plugins.register_format("two-line", [".2l"], reader, writer)
    names = json.loads(plugins.names_json())
    assert names["formats"][0]["name"] == "two-line"
    assert names["formats"][0]["write"] is True
    structure = plugins.read_format("two-line", "cell 4.0\nNa 0 0 0\n")
    text = plugins.write_format("two-line", structure)
    assert text.startswith("cell 4.0")


def test_plugin_errors_are_mapped() -> None:
    plugins = xtalforge.Plugins()
    with pytest.raises(xtalforge.InputError) as error:
        plugins.run_module("missing", "{}")
    assert error.value.args[0] == "unknown_module"
    plugins.register_module("bad", lambda parameters: 42)
    with pytest.raises(TypeError):
        plugins.run_module("bad", "{}")
    plugins.register_module("dup", lambda parameters: "{}")
    with pytest.raises(xtalforge.InputError):
        plugins.register_module("dup", lambda parameters: "{}")


def test_entry_point_discovery(monkeypatch: pytest.MonkeyPatch) -> None:
    import xtalforge.plugins as plugins_module

    class Point:
        name = "demo"

        def load(self) -> object:
            def register(registry: xtalforge.Plugins) -> None:
                registry.register_module("demo", lambda parameters: '{"ok": true}')

            return register

    monkeypatch.setattr(plugins_module, "_entry_points", lambda group: [Point()])
    plugins_module.reset()
    registry = plugins_module.default()
    assert json.loads(registry.names_json())["modules"] == ["demo"]
    assert json.loads(registry.run_module("demo", "{}")) == {"ok": True}


def test_entry_point_failure_is_wrapped(monkeypatch: pytest.MonkeyPatch) -> None:
    import xtalforge.plugins as plugins_module

    class Point:
        name = "broken"

        def load(self) -> object:
            raise RuntimeError("boom")

    monkeypatch.setattr(plugins_module, "_entry_points", lambda group: [Point()])
    with pytest.raises(plugins_module.PluginError) as error:
        plugins_module.discover()
    assert "broken" in str(error.value)
