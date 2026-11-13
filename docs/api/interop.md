# Interoperability (Phase 13)

xtalforge is additive: it reads and writes the formats and objects the wider
ecosystem already uses. Optional Python packages are imported only when an
adapter is called, so the base install stays dependency-free.

## Python object adapters

| Object | Import | Export | Keeps |
|---|---|---|---|
| ASE `Atoms` | `Structure.from_ase(atoms)` | `structure.to_ase()` | cell, pbc, tags, charges, magnetic moments, fixed-atom constraints |
| pymatgen `Structure` / `Molecule` | `Structure.from_pymatgen(obj)` | `structure.to_pymatgen()` | lattice, species, Cartesian/fractional positions, occupancy |
| phonopy | `structure.to_phonopy_atoms()` | — | cell and scaled positions |
| spglib | `xtalforge.get_symmetry_dataset(structure)` | — | spglib-shaped dict |
| RDKit `Mol` | `Structure.from_rdkit(mol)` | `structure.to_rdkit()` | atoms, conformer coordinates, bonds |

```python
from ase.io import read
import xtalforge

structure = xtalforge.Structure.from_ase(read("POSCAR"))
dataset = xtalforge.get_symmetry_dataset(structure)
print(dataset["number"], dataset["international"])
```

The ASE calculator wrapper uses the native UFF implementation:

```python
from ase import Atoms
from xtalforge.adapters import UffAseCalculator

atoms = Atoms("Si2", positions=[[0, 0, 0], [1.3575, 1.3575, 1.3575]], cell=[5.43] * 3, pbc=True)
atoms.calc = UffAseCalculator()
print(atoms.get_potential_energy())
```

## File formats

The registry now includes the crystallographic and electronic-structure
formats used by the common engines. Every entry declares its extensions and
the data it preserves; conversion warns about lost fields.

| Family | Formats |
|---|---|
| Electronic structure | Quantum ESPRESSO `pw.x` input/output, CASTEP `.cell`, FHI-aims `geometry.in`, CP2K, ABINIT, CRYSTAL |
| Trajectories | OUTCAR, `vasprun.xml`, extended XYZ |
| Crystallography | CIF, mmCIF/PDBx, PDB, CSSR, SHELX |
| Molecules | MOL2, SDF/MOL V2000, XYZ |

```bash
xtalforge convert pw.in output.cif
xtalforge trajectory OUTCAR --frame 0 --out final.cif
xtalforge convert structure.cif --format optimade
```

## Remote data

The `net` feature adds blocking HTTP; the WASM package uses browser `fetch`.
Parsing is available without the feature for offline workflows.

```bash
xtalforge fetch optimade --base https://optimade.materialsproject.org \
  --filter 'elements HAS "Si"' --limit 5 --out si.json
xtalforge fetch cod 1000041 --out cod.cif
xtalforge fetch mp mp-149 --api-key "$MP_API_KEY" --out mp.cif
```

The CLI prints a JSON object with `id`, `source`, `license` and the structure
document. Downloads are cached under `~/.cache/xtalforge/remote`.
`xtalforge::remote::optimade_export` writes a structure as one OPTIMADE
`structures` entry for sharing with the same providers.

## Python plugins

Install a package that declares an entry point in the `xtalforge.plugins`
group. The entry point returns a `Plugins` registry, a callable that fills one,
or a mapping with `modules`, `calculators` and `formats`.

```toml
[project.entry-points."xtalforge.plugins"]
my-plugins = "my_package.plugin:register"
```

```python
import xtalforge.plugins

registry = xtalforge.plugins.default()
print(registry.names_json())
```

Rust plugins, which are visible to the CLI and the agent, are described in
[Plugins and extension points](./plugins.md).

## Parity keys

The Phase 13 rows are `ase_adapter`, `pymatgen_adapter`, `phonopy_adapter`,
`rdkit_adapter`, `get_symmetry_dataset`, `interop_formats`, `optimade_fetch`,
`cod_fetch`, `materials_project_fetch` and `python_plugins`. The matrix
records which surfaces are intentionally `n/a` and why.
