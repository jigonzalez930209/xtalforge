"""Structure builder stubs (Phase 6)."""

StructureDict = dict[str, object]

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

def build_c20() -> StructureDict: ...
def twist(n: int, m: int, theta: float, layers: int = 2) -> StructureDict: ...
def ribbon(width: int, length: int, vacuum: float = 12.0) -> StructureDict: ...
def build_ztc(
    samples: int = 16, cell: float = 12.0, isovalue: float = 0.0
) -> StructureDict: ...
def build_schwarzite(
    samples: int = 16, cell: float = 12.0, isovalue: float = 0.0
) -> StructureDict: ...
def build_carbon(
    net: str = "dia",
    repeat: tuple[int, int, int] = (1, 1, 1),
    density: float = 0.42,
    radius_ratio: float = 0.25,
    coverage: float = 0.38,
    seed: int = 0,
    hydrogen: float = 0.07,
    fluorine: float = 0.29,
    oxygen: float = 0.044,
) -> StructureDict: ...
def build_polymer(
    smiles: str, monomers: int, density: float = 1.0
) -> StructureDict: ...
