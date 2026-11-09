from typing import Any, Final

class XtalforgeError(Exception):
    """Base class for every xtalforge error.

    ``args`` carries the stable error code first and the human-readable
    message second, for example ``("invalid_input", "invalid input: ...")``.
    """

class InputError(XtalforgeError):
    """A geometry or value does not satisfy the API preconditions."""

StructureDict = dict[str, object]
BondDict = dict[str, object]
GroupDict = dict[str, object]
NeighborDict = dict[str, object]

__version__: Final[str]

def version() -> str: ...

class Lattice:
    def __init__(self, a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> None: ...
    @staticmethod
    def from_matrix(matrix: list[list[float]]) -> Lattice: ...
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float
    volume: float
    parameters: list[float]
    matrix: Any
    def to_cart(self, point: list[float]) -> Any: ...
    def to_frac(self, point: list[float]) -> Any: ...
    def d_spacing(self, hkl: list[float]) -> float: ...
    def distance(self, first: list[float], second: list[float]) -> float: ...
    def minimum_image(self, delta: list[float]) -> tuple[Any, float]: ...
    def scaled_to_volume(self, volume: float) -> Lattice: ...

class Site:
    def __init__(self, label: str, element: str, x: float, y: float, z: float) -> None: ...
    label: str
    element: str
    frac: Any
    occupancy: float
    u_iso: float | None
    charge: float | None
    def with_occupancy(self, occupancy: float) -> Site: ...
    def with_charge(self, charge: float) -> Site: ...
    def to_dict(self) -> dict[str, object]: ...

class Structure:
    def __init__(self, text: str, format_hint: str = "") -> None: ...
    @staticmethod
    def from_cif(text: str) -> Structure: ...
    name: str
    formula: str
    hall: str
    lattice: Lattice
    sites: list[Site]
    n_sites: int
    volume: float
    density: float
    is_p1: bool
    def p1(self) -> P1Cell: ...
    def supercell(self, na: int, nb: int, nc: int) -> Structure: ...
    def slab(self, h: int, k: int, l: int, layers: int, vacuum: float) -> Structure: ...
    def bonds(self) -> list[dict[str, object]]: ...
    def rings(self, max_size: int) -> list[list[int]]: ...
    def ring_paths(
        self, max_size: int
    ) -> list[list[tuple[int, tuple[int, int, int]]]]: ...
    def summary(self) -> dict[str, object]: ...
    def to_format(self, format: str) -> str: ...
    def to_json(self) -> str: ...

class P1Cell:
    n_atoms: int
    lattice: Lattice
    elements: list[str]
    frac: Any
    occupancy: list[float]
    atoms: list[dict[str, object]]

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
def read_structure(text: str, format_hint: str) -> StructureDict: ...
def read_cif(text: str) -> StructureDict: ...
def supercell(
    text: str,
    na: int,
    nb: int,
    nc: int,
    format_hint: str | None = None,
) -> StructureDict: ...
def slab(
    text: str,
    h: int,
    k: int,
    l: int,
    layers: int,
    vacuum: float,
    format_hint: str | None = None,
) -> StructureDict: ...
def bonds(text: str, format_hint: str | None = None) -> list[BondDict]: ...
def rings(
    text: str,
    max_size: int,
    format_hint: str | None = None,
) -> list[list[int]]: ...
def ring_paths(
    text: str,
    max_size: int,
    format_hint: str | None = None,
) -> list[list[tuple[int, tuple[int, int, int]]]]: ...
def groups(text: str, format_hint: str | None = None) -> list[GroupDict]: ...
def neighbors(
    text: str,
    cutoff: float,
    format_hint: str | None = None,
) -> list[NeighborDict]: ...
def distance(
    text: str,
    i: int,
    j: int,
    format_hint: str | None = None,
) -> float: ...
def angle(
    text: str,
    i: int,
    j: int,
    k: int,
    format_hint: str | None = None,
) -> float: ...
def torsion(
    text: str,
    i: int,
    j: int,
    k: int,
    l: int,
    format_hint: str | None = None,
) -> float: ...

def detect_symmetry(
    text: str,
    format_hint: str | None = None,
    symprec: float = 1e-3,
    angle_tolerance: float | None = None,
    hall_number: int | None = None,
) -> dict[str, object]: ...

def asymmetrize(
    text: str,
    format_hint: str | None = None,
    symprec: float = 1e-3,
) -> StructureDict: ...
def to_p1(text: str, format_hint: str | None = None) -> StructureDict: ...
def set_space_group(
    text: str,
    symbol: str,
    mode: str = "reinterpret",
    format_hint: str | None = None,
) -> StructureDict: ...

def subgroups(hall_number: int, epsilon: float = 1e-5) -> list[dict[str, object]]: ...
def lower_to_subgroup(
    text: str,
    subgroup: int,
    tolerance: float = 0.5,
    format_hint: str | None = None,
) -> StructureDict: ...

def standardize(
    text: str,
    symprec: float = 1e-3,
    format_hint: str | None = None,
) -> StructureDict: ...
def primitive_cell(
    text: str,
    symprec: float = 1e-3,
    format_hint: str | None = None,
) -> StructureDict: ...
def reduce_structure(
    text: str,
    algorithm: str = "niggli",
    format_hint: str | None = None,
) -> StructureDict: ...
def wyckoff_assignments(
    text: str,
    symprec: float = 1e-3,
    format_hint: str | None = None,
) -> list[dict[str, object]]: ...
def brillouin_zone(text: str, format_hint: str | None = None) -> dict[str, object]: ...
def tolerance_ladder(
    text: str,
    format_hint: str | None = None,
) -> list[dict[str, object]]: ...
def diagnose(
    text: str,
    symbol: str,
    format_hint: str | None = None,
) -> dict[str, object]: ...

def detect_layer_symmetry(
    text: str,
    symprec: float = 1e-3,
    format_hint: str | None = None,
) -> dict[str, object]: ...

def detect_magnetic_symmetry(
    text: str,
    symprec: float = 1e-3,
    format_hint: str | None = None,
) -> dict[str, object]: ...

def change_setting(
    text: str,
    hall_number: int,
    format_hint: str | None = None,
) -> StructureDict: ...

def band_path(text: str, format_hint: str | None = None) -> dict[str, object]: ...

def molecular_point_group(
    text: str,
    tolerance: float = 1e-3,
    format_hint: str | None = None,
) -> dict[str, object]: ...

def scene_summary(
    text: str,
    style: str = "ball-stick",
    scale: float = 1.0,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def mesh_stl(
    text: str,
    style: str = "ball-stick",
    scale: float = 1.0,
    format_hint: str | None = None,
) -> str: ...

def pxrd(
    text: str,
    radiation: str = "cu-ka",
    max_two_theta: float = 80.0,
    step: float = 0.02,
    format_hint: str | None = None,
) -> dict[str, object]: ...

def build_graphene(
    nx: int,
    ny: int,
    layers: int = 1,
    interlayer: float = 3.35,
    vacuum: float = 15.0,
) -> StructureDict: ...
def build_nanotube(
    n: int,
    m: int,
    length: int = 1,
    vacuum: float = 10.0,
) -> StructureDict: ...
def build_c60() -> StructureDict: ...


class Document:
    def __init__(self, text: str, format_hint: str | None = None) -> None: ...
    def move_sites(self, sites: list[int], delta: list[float]) -> None: ...
    def set_element(self, site: int, element: str) -> None: ...
    def set_charges(self, charges: list[float]) -> None: ...
    def interpenetrate(self, shift: list[float]) -> None: ...
    def undo(self) -> bool: ...
    def redo(self) -> bool: ...
    can_undo: bool
    can_redo: bool
    dirty: bool
    history: list[str]
    def to_cif(self) -> str: ...
    def to_json(self) -> str: ...
    def summary(self) -> StructureDict: ...

def param_schema(name: str) -> dict[str, object]: ...
def project_pack(text: str, format_hint: str | None = None) -> bytes: ...
def project_read(bytes: bytes) -> StructureDict: ...
def workspace_open(path: str) -> str: ...
def workspace_write_run(path: str, record_json: str) -> str: ...
def workspace_runs(path: str) -> list[str]: ...
def workspace_read_run(path: str, id: str) -> str: ...
def workspace_archive(records_json: str) -> bytes: ...
def workspace_read_archive(bytes: bytes) -> list[str]: ...


def ff_energy(text: str, format_hint: str | None = None) -> float: ...
def ff_optimize(
    text: str,
    minimizer: str = "fire",
    max_steps: int = 500,
    force_tolerance: float = 0.05,
    cell: bool = False,
    format_hint: str | None = None,
) -> StructureDict: ...
def ff_charges(
    text: str,
    method: str = "qeq",
    total_charge: float = 0.0,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def ff_vibrations(
    text: str,
    step: float = 0.005,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def ff_ir(
    text: str,
    fwhm: float = 10.0,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def ff_thermo(
    text: str,
    temperature: float = 298.15,
    pressure: float = 1.0,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def ff_types(text: str, format_hint: str | None = None) -> list[dict[str, object]]: ...


def porosity(
    text: str,
    spacing: float = 0.35,
    probe: float = 1.2,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def zeo_parse(text: str) -> dict[str, object]: ...
def zeo_run(
    text: str,
    probe: float = 1.2,
    samples: int = 1000,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def net_invariants(text: str, tolerance: float = 0.01) -> list[dict[str, object]]: ...
def net_overlay(text: str, tolerance: float = 0.01) -> dict[str, object]: ...
def rcsr_names() -> list[str]: ...
def add_hydrogens(
    text: str,
    distance: float = 1.09,
    format_hint: str | None = None,
) -> StructureDict: ...

def ff_prepare(
    text: str,
    minimize: bool = True,
    charges: bool = True,
    ewald: bool = False,
    vibrations: bool = False,
    format_hint: str | None = None,
    merge_duplicates: bool = False,
    order_disorder: bool = False,
    remove_solvent: bool = False,
    check_charge_balance: bool = False,
    hydrogens: bool = False,
    check_close_contacts: bool = False,
) -> dict[str, object]: ...

def mof_topologies() -> list[dict[str, object]]: ...
def mof_blocks() -> list[dict[str, object]]: ...
def build_mof(topology: str, node: str, linker: str) -> StructureDict: ...

def powder_automatic(
    text: str,
    pattern_text: str,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def powder_index(
    pattern_text: str,
    wavelength: float = 1.54184,
    max_index: int = 3,
    tolerance: float = 0.05,
    background: bool = False,
    refine: bool = False,
) -> dict[str, object]: ...
def powder_refine(
    text: str,
    pattern_text: str,
    refine_positions: bool = True,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def powder_cells(number: int) -> dict[str, object]: ...
def powder_supported_groups() -> list[int]: ...
def rietx() -> dict[str, object]: ...

def build_c20() -> StructureDict: ...
def twist(n: int, m: int, theta: float, layers: int = 2) -> StructureDict: ...
def ribbon(width: int, length: int, vacuum: float = 12.0) -> StructureDict: ...
def build_ztc(
    samples: int = 16, cell: float = 12.0, isovalue: float = 0.0
) -> StructureDict: ...
def build_schwarzite(
    samples: int = 16, cell: float = 12.0, isovalue: float = 0.0
) -> StructureDict: ...
def build_polymer(
    smiles: str, monomers: int, density: float = 1.0
) -> StructureDict: ...
