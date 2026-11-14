from typing import Any

def magnetic_orderings(
    text: str,
    supercell: int = 1,
    max_orderings: int = 8,
    format_hint: str | None = None,
) -> list[dict[str, object]]: ...
def disorder_groups(
    text: str, format_hint: str | None = None
) -> dict[str, object]: ...
def derivative_structures(
    text: str,
    max_index: int = 2,
    max_structures: int = 8,
    format_hint: str | None = None,
) -> list[dict[str, object]]: ...
def neutron_powder(
    text: str,
    wavelength: float = 1.54,
    step: float = 0.05,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def magnetic_powder(
    text: str,
    wavelength: float = 1.54,
    step: float = 0.05,
    allow_missing_moments: bool = False,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def electron_diffraction(
    text: str,
    zone: tuple[int, int, int] = (0, 0, 1),
    voltage: float = 200.0,
    size: int = 256,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def raman_powder_average(
    tensors: list[list[list[float]]],
) -> list[float]: ...
def born_diagonal_average(
    tensors: list[list[list[float]]],
) -> list[float]: ...
