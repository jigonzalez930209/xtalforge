"""Template plugin: register one module and one format."""

import json

import xtalforge


def register(registry: xtalforge.Plugins) -> None:
    """Registers a demo module and a minimal two-line format."""
    registry.register_module(
        "demo",
        lambda parameters: json.dumps({"plugin": "demo", "parameters": json.loads(parameters)}),
    )

    def read(text: str) -> str:
        lines = [line for line in text.splitlines() if line.strip()]
        if not lines or lines[0] != "demo":
            raise ValueError("expected a 'demo' header")
        element, x, y, z = lines[1].split()
        return xtalforge.Structure(
            f"1\n{element}\n{element} {x} {y} {z}\n", "xyz"
        ).to_json()

    def write(structure: str) -> str:
        document = json.loads(structure)["structure"]
        site = document["sites"][0]
        return f"demo\n{site['element']} {site['frac'][0]} {site['frac'][1]} {site['frac'][2]}\n"

    registry.register_format("demo", [".demo"], read, write)
