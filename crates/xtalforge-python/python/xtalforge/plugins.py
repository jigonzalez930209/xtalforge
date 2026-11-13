"""Python entry-point plugins (Task 13.4).

Plugins register through the ``xtalforge.plugins`` entry-point group. Each
entry point returns one of:

- a :class:`xtalforge.Plugins` instance;
- a callable that receives a fresh registry and registers into it;
- a mapping with ``modules``, ``calculators`` and ``formats`` callables.

The registry is in-process, like the Rust extension's ``Plugins`` class. The
CLI and the agent use the Rust plugin registry; a Python plugin can expose its
operations to those surfaces through a small wrapper script that prints
structure JSON.
"""

from __future__ import annotations

from importlib import metadata
from typing import Callable, Iterable, cast

from . import Plugins

PluginFactory = Callable[[Plugins], None]
ENTRY_POINT_GROUP = "xtalforge.plugins"


class PluginError(RuntimeError):
    """A plugin entry point could not be loaded or did not match the contract."""


def _entry_points(group: str) -> Iterable[metadata.EntryPoint]:
    points = metadata.entry_points()
    if hasattr(points, "select"):
        return points.select(group=group)
    return cast(
        Iterable[metadata.EntryPoint],
        points.get(group, ()),  # type: ignore[attr-defined]
    )


def _module_runner(source: Plugins, name: str) -> Callable[[str], str]:
    def run(parameters: str) -> str:
        return source.run_module(name, parameters)

    return run


def _calculator_runner(source: Plugins, name: str) -> Callable[[str], str]:
    def run(structure: str) -> str:
        return source.calculator_energy(name, structure)

    return run


def _reader_runner(source: Plugins, name: str) -> Callable[[str], str]:
    def run(text: str) -> str:
        return source.read_format(name, text)

    return run


def _writer_runner(source: Plugins, name: str) -> Callable[[str], str]:
    def run(structure: str) -> str:
        return source.write_format(name, structure)

    return run


def _apply(registry: Plugins, plugin: object, name: str) -> None:
    if isinstance(plugin, Plugins):
        _merge(registry, plugin)
        return
    if callable(plugin):
        plugin(registry)
        return
    if isinstance(plugin, dict):
        for key, callable_plugin in plugin.get("modules", {}).items():
            registry.register_module(str(key), callable_plugin)
        for key, callable_plugin in plugin.get("calculators", {}).items():
            registry.register_calculator(str(key), callable_plugin)
        for key, (extensions, reader, writer) in plugin.get("formats", {}).items():
            registry.register_format(str(key), list(extensions), reader, writer)
        return
    raise PluginError(
        f"plugin {name!r} returned {type(plugin).__name__}; "
        "expected Plugins, a callable or a mapping"
    )


def _merge(target: Plugins, source: Plugins) -> None:
    import json

    parsed = json.loads(source.names_json())
    for module in parsed.get("modules", []):
        target.register_module(module, _module_runner(source, module))
    for calculator in parsed.get("calculators", []):
        target.register_calculator(
            calculator, _calculator_runner(source, calculator)
        )
    for fmt in parsed.get("formats", []):
        name = fmt["name"]
        target.register_format(
            name,
            list(fmt["extensions"]),
            _reader_runner(source, name),
            _writer_runner(source, name) if fmt["write"] else None,
        )


def discover(group: str = ENTRY_POINT_GROUP) -> Plugins:
    """Loads every entry point in *group* into a new registry."""
    registry = Plugins()
    for entry in _entry_points(group):
        try:
            plugin = entry.load()
            _apply(registry, plugin, entry.name)
        except PluginError:
            raise
        except Exception as error:  # noqa: BLE001 - entry points are user code
            raise PluginError(f"could not load plugin {entry.name!r}: {error}") from error
    return registry


_DEFAULT: Plugins | None = None


def default() -> Plugins:
    """Returns the process-wide entry-point registry, loading it once."""
    global _DEFAULT
    if _DEFAULT is None:
        _DEFAULT = discover()
    return _DEFAULT


def reset() -> None:
    """Forgets the cached registry; the next :func:`default` reloads it."""
    global _DEFAULT
    _DEFAULT = None


__all__ = [
    "ENTRY_POINT_GROUP",
    "PluginError",
    "default",
    "discover",
    "reset",
]
