from typing import Final

class XtalforgeError(Exception):
    """Base class for every xtalforge error.

    ``args`` carries the stable error code first and the human-readable
    message second, for example ``("invalid_input", "invalid input: ...")``.
    """

class InputError(XtalforgeError):
    """A geometry or value does not satisfy the API preconditions."""

__version__: Final[str]

def version() -> str: ...
def element_info(input: str | int) -> dict[str, object]: ...
def lattice_info(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> dict[str, object]: ...
def lattice_minimum_image(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
    delta: list[float],
) -> dict[str, object]: ...
def limits_bounds(profile: str, budget: str) -> tuple[int, int | None]: ...
def limits_check(profile: str, budget: str, estimate: int) -> dict[str, object]: ...
def structure_info(text: str) -> dict[str, object]: ...
def cell_volume(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> float: ...
