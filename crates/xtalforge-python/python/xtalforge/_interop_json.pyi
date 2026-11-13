"""Shared JSON helpers for the optional adapters (Task 13.1)."""

from ._object_api import Structure

INTEROP: str
MOLECULE_PADDING: float

def jsonable(value: object) -> object: ...
def document(structure: Structure) -> dict[str, object]: ...
def interop(document_data: dict[str, object]) -> dict[str, object]: ...
def filled(values: list[float | None]) -> list[float] | None: ...
def cart(frac: list[float], lattice: list[list[float]]) -> list[float]: ...
def build_document(
    name: str,
    entries: list[dict[str, object]],
    periodic: bool,
    lattice: list[list[float]] | None,
    interop_data: dict[str, object],
    cartesian: bool = False,
) -> dict[str, object]: ...
