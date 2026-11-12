"""Phase 9 comparison, edit, analysis and mechanics stubs."""

from typing import Final

StructureDict = dict[str, object]

def match_structures(
    first: str,
    second: str,
    format_hint: str | None = None,
    ltol: float = 0.2,
    stol: float = 0.3,
    angle_tol: float = 5.0,
    allow_supercell: bool = True,
) -> dict[str, object]: ...
def periodic_rmsd(
    first: str, second: str, mapping: list[int], format_hint: str | None = None
) -> float: ...
def kabsch_alignment(
    points: list[list[float]], reference: list[list[float]]
) -> dict[str, object]: ...
def group_structures(
    structures: list[str], ltol: float = 0.2, stol: float = 0.3
) -> dict[str, object]: ...
def deform(
    text: str,
    format_hint: str | None = None,
    voigt: list[float] | None = None,
    gradient: list[list[float]] | None = None,
    volume: float | None = None,
    density: float | None = None,
) -> StructureDict: ...
def rattle(
    text: str,
    format_hint: str | None = None,
    amplitude: float = 0.05,
    seed: int = 0,
    bounded: bool = False,
) -> StructureDict: ...
def edit_structure(
    text: str,
    format_hint: str | None = None,
    add: list[tuple[str, list[float]]] | None = None,
    remove: list[int] | None = None,
    replace: list[tuple[int, str]] | None = None,
    move_sites: list[tuple[int, list[float]]] | None = None,
    species: list[tuple[str, str]] | None = None,
    fraction: float = 1.0,
    seed: int = 0,
    sort: str | None = None,
    wrap: bool = False,
) -> StructureDict: ...
def interpolate(
    first: str, second: str, images: int = 5, format_hint: str | None = None
) -> list[StructureDict]: ...
def interpolate_structures(
    first: str, second: str, images: int = 5, format_hint: str | None = None
) -> list[StructureDict]: ...
def optimal_supercell(
    text: str, cells: int, format_hint: str | None = None
) -> StructureDict: ...
def voronoi(
    text: str, cutoff: float = 6.0, format_hint: str | None = None
) -> dict[str, object]: ...
def coordination(
    text: str,
    index: int,
    strategy: str = "voronoi",
    cutoff: float = 5.0,
    solid_angle: float = 0.1,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def bond_valence_sums(
    text: str, cutoff: float = 5.0, format_hint: str | None = None
) -> dict[str, object]: ...
def oxidation_states(
    text: str, cutoff: float = 5.0, format_hint: str | None = None
) -> dict[str, object]: ...
def radial_distribution(
    text: str,
    rmax: float = 10.0,
    bin: float = 0.05,
    smoothing: int = 0,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def structure_fingerprint(
    text: str, rmax: float = 10.0, bin: float = 0.1, format_hint: str | None = None
) -> dict[str, object]: ...
def fingerprint_similarity(
    first: str, second: str, rmax: float = 10.0, bin: float = 0.1
) -> float: ...
def unique_miller_indices(
    text: str,
    max_index: int = 2,
    symprec: float = 1e-5,
    format_hint: str | None = None,
) -> list[dict[str, object]]: ...
def slab_terminations(
    text: str,
    h: int,
    k: int,
    l: int,
    layers: int = 2,
    vacuum: float = 10.0,
    symmetric: bool = False,
    format_hint: str | None = None,
) -> list[StructureDict]: ...
def wulff(
    text: str, surfaces: list[list[float]], format_hint: str | None = None
) -> dict[str, object]: ...
def adsorption_sites(
    text: str, tolerance: float = 0.05, format_hint: str | None = None
) -> list[dict[str, object]]: ...
def place_molecule(
    slab: str,
    molecule: str,
    site_index: int = 0,
    offset: float = 2.0,
    format_hint: str | None = None,
) -> StructureDict: ...
def vacancies(
    text: str, min_image: float = 10.0, format_hint: str | None = None
) -> list[StructureDict]: ...
def substitutions(
    text: str, species: str, min_image: float = 10.0, format_hint: str | None = None
) -> list[StructureDict]: ...
def interstitials(
    text: str,
    species: str,
    min_image: float = 10.0,
    cutoff: float = 5.0,
    clearance: float = 1.0,
    format_hint: str | None = None,
) -> list[StructureDict]: ...
def orderings(
    text: str,
    supercell: tuple[int, int, int] = (1, 1, 1),
    max_orderings: int = 64,
    format_hint: str | None = None,
) -> list[StructureDict]: ...
def special_quasirandom(
    text: str,
    supercell: tuple[int, int, int] = (2, 2, 2),
    seed: int = 0,
    cutoff: float = 6.0,
    iterations: int = 5000,
    format_hint: str | None = None,
) -> StructureDict: ...
def lattice_matches(
    film: str,
    substrate: str,
    max_index: int = 3,
    max_area: float = 200.0,
    tolerance: float = 0.05,
    format_hint: str | None = None,
) -> list[dict[str, object]]: ...
def stack_slabs(
    film: str,
    substrate: str,
    gap: float = 2.0,
    vacuum: float = 15.0,
    shift: tuple[float, float] = (0.0, 0.0),
    format_hint: str | None = None,
) -> StructureDict: ...
def elastic_tensor(
    text: str, strain: float = 0.005, format_hint: str | None = None
) -> dict[str, object]: ...
def eos_fit(
    points: list[tuple[float, float]], model: str = "birch"
) -> dict[str, object]: ...
def convex_hull(
    entries: list[tuple[str, dict[str, float], float]],
) -> dict[str, object]: ...
