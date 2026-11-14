def generate_crystal(
    composition: dict[str, int],
    space_group: int = 225,
    seed: int = 0,
    volume_per_atom: float = 20.0,
    min_distance: float = 1.0,
    attempts: int = 64,
) -> dict[str, object]: ...
def generate_crystal_range(
    composition: dict[str, int],
    groups: list[int],
    seed: int = 0,
    volume_per_atom: float = 20.0,
    min_distance: float = 1.0,
    attempts: int = 64,
) -> list[dict[str, object]]: ...
def generate_molecule_crystal(
    text: str,
    space_group: int = 1,
    copies: int = 1,
    seed: int = 0,
    volume_per_atom: float = 20.0,
    min_distance: float = 1.0,
    attempts: int = 64,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def generate_layer_crystal(
    composition: dict[str, int],
    layer_group: int = 1,
    seed: int = 0,
    area_per_atom: float = 10.0,
    min_distance: float = 1.0,
    attempts: int = 64,
) -> dict[str, object]: ...
def generate_rod_crystal(
    composition: dict[str, int],
    rod_group: int = 1,
    seed: int = 0,
    repeat: float = 5.0,
    min_distance: float = 1.0,
    attempts: int = 64,
) -> dict[str, object]: ...
def layer_groups() -> list[dict[str, object]]: ...
def prototype_library() -> list[dict[str, object]]: ...
def build_prototype(label: str, scale: float = 1.0) -> dict[str, object]: ...
def identify_prototype(
    text: str, stol: float = 0.3, format_hint: str | None = None
) -> dict[str, object]: ...
def substitution_candidates(
    text: str,
    radius_tolerance: float = 0.15,
    limit: int = 32,
    same_charge: bool = True,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def random_search(
    composition: dict[str, int],
    space_groups: list[int],
    structures: int = 8,
    seed: int = 0,
    relax: bool = True,
    cell: bool = False,
    steps: int = 200,
) -> dict[str, object]: ...
def evolutionary_search(
    composition: dict[str, int],
    space_groups: list[int],
    population: int = 8,
    generations: int = 5,
    seed: int = 0,
    cell: bool = False,
    steps: int = 150,
) -> dict[str, object]: ...
def search_workflow(
    kind: str,
    seed_path: str,
    options: object,
    checkpoint: str | None = None,
    output: str = "search-results.json",
    name: str = "search",
) -> dict[str, object]: ...
def mof_dataset(
    topologies: list[str] = [],
    nodes: list[str] = [],
    linkers: list[str] = [],
    limit: int = 32,
    relax: bool = False,
    descriptors: bool = True,
    seed: int = 0,
) -> dict[str, object]: ...
def functionalize_sweep(
    text: str,
    groups: list[str] = ["amino"],
    fractions: list[float] = [1.0],
    seed: int = 0,
    descriptors: bool = True,
    format_hint: str | None = None,
) -> dict[str, object]: ...
