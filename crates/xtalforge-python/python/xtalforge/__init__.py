"""Python bindings for xtalforge."""

import json

from ._xtalforge import (
    AgentSession,
    Document,
    InputError,
    Lattice,
    MlCalculator,
    P1Cell,
    Plugins,
    Site,
    Structure,
    XtalforgeError,
    __version__,
    add_hydrogens,
    agent_capabilities,
    agent_run,
    agent_verbs,
    angle,
    asymmetrize,
    band_path,
    bonds,
    brillouin_zone,
    build_c60,
    build_graphene,
    build_c20,
    build_carbon,
    build_mof,
    build_polymer,
    build_schwarzite,
    build_ztc,
    build_nanotube,
    cell_volume,
    change_setting,
    detect_layer_symmetry,
    detect_magnetic_symmetry,
    detect_symmetry,
    diagnose,
    distance,
    element_info,
    ff_charges,
    ff_energy,
    ff_ir,
    ff_optimize,
    ff_prepare,
    ff_thermo,
    ff_types,
    ff_vibrations,
    groups,
    lattice_info,
    lattice_minimum_image,
    limits_bounds,
    limits_check,
    ml_backends,
    lower_to_subgroup,
    mesh_gltf,
    mesh_obj,
    mesh_stl,
    mof_blocks,
    mof_topologies,
    molecular_point_group,
    net_invariants,
    net_overlay,
    neighbors,
    porosity,
    rcsr_names,
    param_schema,
    primitive_cell,
    project_pack,
    project_read,
    powder_automatic,
    powder_index,
    powder_index_agreement,
    powder_lebail,
    powder_pawley,
    powder_cells,
    powder_refine,
    powder_supported_groups,
    pxrd,
    rietx,
    read_cif,
    read_structure,
    reduce_structure,
    ribbon,
    ring_paths,
    rings,
    scene_data,
    scene_summary,
    set_space_group,
    slab,
    standardize,
    structure_info,
    subgroups,
    supercell,
    to_p1,
    tolerance_ladder,
    torsion,
    twist,
    version,
    set_threads,
    error_category,
    deprecation_warning,
    molecule_from_smiles,
    sketch_roundtrip,
    fragment_library,
    fragment_names,
    fragment_structure,
    substitute_structure,
    place_guest_structure,
    read,
    write,
    workflow_run,
    workflow_schema,
    workflow_validate,
    provenance_json,
    workspace_archive,
    workspace_open,
    workspace_read_archive,
    workspace_read_run,
    workspace_runs,
    workspace_write_run,
    wyckoff_assignments,
    zeo_parse,
    mof_check,
    mof_solvent_remove,
    mof_disorder,
    mof_charge_balance,
    mof_deconstruct,
    mof_mofid,
    mof_open_metal_sites,
    mof_descriptors,
    mof_pores,
    mof_guests,
    mof_forcefield,
    mof_widom,
    mof_gcmc,
    mof_energy_grid,
    mof_raspa,
    mof_functionalize,
    mof_defects,
    mof_interpenetrate,
    mof_stack_layers,
    mof_mixed,
    mof_fetch_database,
    mof_read_database,
    mof_screen,
    mof_deduplicate,
    adsorption_sites,
    bond_valence_sums,
    convex_hull,
    coordination,
    deform,
    edit_structure,
    elastic_tensor,
    eos_fit,
    fingerprint_similarity,
    group_structures,
    interstitials,
    interpolate,
    interpolate_structures,
    kabsch_alignment,
    lattice_matches,
    match_structures,
    optimal_supercell,
    orderings,
    oxidation_states,
    periodic_rmsd,
    place_molecule,
    radial_distribution,
    rattle,
    slab_terminations,
    special_quasirandom,
    structure_fingerprint,
    substitutions,
    unique_miller_indices,
    vacancies,
    voronoi,
    wulff,
    zeo_run,
)

def _text_input(value):
    """Returns structure JSON for a Structure, otherwise the value unchanged."""
    return value.to_json() if isinstance(value, Structure) else value


def _accepts_structure(function):
    def wrapped(first, *args, **kwargs):
        return function(_text_input(first), *args, **kwargs)

    wrapped.__name__ = getattr(function, "__name__", "wrapped")
    wrapped.__doc__ = getattr(function, "__doc__", None)
    return wrapped


class _TypedResult(dict):
    """Read-only mapping with the legacy dict shape."""

    def to_dict(self) -> dict:
        return dict(self)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({dict.__repr__(self)})"

    def _repr_html_(self) -> str:
        rows = "".join(
            f"<tr><th>{key}</th><td>{value}</td></tr>" for key, value in self.items()
        )
        return f"<table>{rows}</table>"


class SymmetryDataset(_TypedResult):
    """Typed symmetry dataset; also behaves like the old dict."""


class PxrdPattern(_TypedResult):
    """Typed PXRD pattern; also behaves like the old dict."""


class PorosityReport(_TypedResult):
    """Typed porosity report; also behaves like the old dict."""


class OptimizeResult(_TypedResult):
    """Typed force-field optimization result; also behaves like the old dict."""


_raw_detect_symmetry = detect_symmetry
_raw_pxrd = pxrd
_raw_porosity = porosity
_raw_ff_optimize = ff_optimize


def detect_symmetry(*args, **kwargs):
    return SymmetryDataset(_raw_detect_symmetry(_text_input(args[0]), *args[1:], **kwargs))


def pxrd(*args, **kwargs):
    return PxrdPattern(_raw_pxrd(_text_input(args[0]), *args[1:], **kwargs))


def porosity(*args, **kwargs):
    return PorosityReport(_raw_porosity(_text_input(args[0]), *args[1:], **kwargs))


def ff_optimize(*args, **kwargs):
    return OptimizeResult(_raw_ff_optimize(_text_input(args[0]), *args[1:], **kwargs))


for _name in [
    "add_hydrogens",
    "asymmetrize",
    "bonds",
    "ff_charges",
    "ff_energy",
    "ff_vibrations",
    "groups",
    "mof_check",
    "mof_deconstruct",
    "mof_mofid",
    "neighbors",
    "primitive_cell",
    "reduce_structure",
    "rings",
    "slab",
    "standardize",
    "supercell",
    "to_p1",
]:
    if _name in globals() and callable(globals()[_name]):
        globals()[_name] = _accepts_structure(globals()[_name])

from . import plugins as plugins
from .workflow import Workflow as Workflow
from .dataset import Dataset as Dataset


def provenance(text: str) -> dict[str, object] | None:
    """Returns the embedded provenance record, when present."""
    value = provenance_json(text)
    return None if value is None else json.loads(value)

__all__ = [
    "AgentSession",
    "AseCalculator",
    "Document",
    "InputError",
    "Lattice",
    "MlCalculator",
    "P1Cell",
    "Plugins",
    "Site",
    "Structure",
    "SymmetryDataset",
    "PxrdPattern",
    "PorosityReport",
    "OptimizeResult",
    "XtalforgeError",
    "__version__",
    "add_hydrogens",
    "agent_capabilities",
    "agent_run",
    "agent_verbs",
    "angle",
    "asymmetrize",
    "band_path",
    "bonds",
    "brillouin_zone",
    "build_c60",
    "build_graphene",
    "build_c20",
    "build_carbon",
    "build_mof",
    "build_polymer",
    "build_schwarzite",
    "build_ztc",
    "build_nanotube",
    "cell_volume",
    "change_setting",
    "detect_layer_symmetry",
    "detect_magnetic_symmetry",
    "detect_symmetry",
    "diagnose",
    "distance",
    "element_info",
    "ff_charges",
    "ff_energy",
    "ff_ir",
    "ff_optimize",
    "ff_prepare",
    "ff_thermo",
    "ff_types",
    "ff_vibrations",
    "get_symmetry_dataset",
    "groups",
    "lattice_info",
    "lattice_minimum_image",
    "limits_bounds",
    "limits_check",
    "lower_to_subgroup",
    "ml_backends",
    "mesh_gltf",
    "mesh_obj",
    "mesh_stl",
    "mof_blocks",
    "mof_topologies",
    "molecular_point_group",
    "net_invariants",
    "net_overlay",
    "neighbors",
    "porosity",
    "plugins",
    "rcsr_names",
    "param_schema",
    "primitive_cell",
    "project_pack",
    "project_read",
    "powder_automatic",
    "powder_index",
    "powder_index_agreement",
    "powder_lebail",
    "powder_pawley",
    "powder_cells",
    "powder_refine",
    "powder_supported_groups",
    "pxrd",
    "rietx",
    "read_cif",
    "read_structure",
    "reduce_structure",
    "ribbon",
    "ring_paths",
    "rings",
    "scene_data",
    "scene_summary",
    "set_space_group",
    "slab",
    "standardize",
    "structure_info",
    "subgroups",
    "supercell",
    "to_p1",
    "tolerance_ladder",
    "torsion",
    "twist",
    "version",
    "set_threads",
    "error_category",
    "deprecation_warning",
    "molecule_from_smiles",
    "sketch_roundtrip",
    "fragment_library",
    "fragment_names",
    "fragment_structure",
    "substitute_structure",
    "place_guest_structure",
    "read",
    "write",
    "Workflow",
    "Dataset",
    "provenance",
    "workspace_archive",
    "workspace_open",
    "workspace_read_archive",
    "workspace_read_run",
    "workspace_runs",
    "workspace_write_run",
    "wyckoff_assignments",
    "zeo_parse",
    "mof_check",
    "mof_solvent_remove",
    "mof_disorder",
    "mof_charge_balance",
    "mof_deconstruct",
    "mof_mofid",
    "mof_open_metal_sites",
    "mof_descriptors",
    "mof_pores",
    "mof_guests",
    "mof_forcefield",
    "mof_widom",
    "mof_gcmc",
    "mof_energy_grid",
    "mof_raspa",
    "mof_functionalize",
    "mof_defects",
    "mof_interpenetrate",
    "mof_stack_layers",
    "mof_mixed",
    "mof_fetch_database",
    "mof_read_database",
    "mof_screen",
    "mof_deduplicate",
    "adsorption_sites",
    "bond_valence_sums",
    "convex_hull",
    "coordination",
    "deform",
    "edit_structure",
    "elastic_tensor",
    "eos_fit",
    "fingerprint_similarity",
    "group_structures",
    "interstitials",
    "interpolate",
    "interpolate_structures",
    "kabsch_alignment",
    "lattice_matches",
    "match_structures",
    "optimal_supercell",
    "orderings",
    "oxidation_states",
    "periodic_rmsd",
    "place_molecule",
    "radial_distribution",
    "rattle",
    "slab_terminations",
    "special_quasirandom",
    "structure_fingerprint",
    "substitutions",
    "unique_miller_indices",
    "vacancies",
    "voronoi",
    "wulff",
    "zeo_run",
]

from .adapters import get_symmetry_dataset as get_symmetry_dataset


def __getattr__(name):
    """Resolves the optional ASE calculator without importing ASE eagerly."""
    if name == "AseCalculator":
        from .adapters import AseCalculator

        return AseCalculator
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
