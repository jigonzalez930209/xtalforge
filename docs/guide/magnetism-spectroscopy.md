# Magnetism and spectroscopy from a CIF

This tutorial walks from a magnetic CIF to an ordering, a neutron pattern and
an IR/Raman spectrum. The magnetic and diffraction parts are pure Rust on all
four surfaces; the spectra need modes and tensors from an engine or the native
force field.

## 1. Read and inspect the moments

```bash
xtalforge info magnetic.cif --json
```

`xtalforge_io::read_mcif` maps `_atom_site_moment_*` onto
`Site::magnetic_moment`. A collinear moment is an axis plus a signed magnitude;
a non-collinear one is a Cartesian vector. The JSON structure document carries
the same field.

## 2. Detect and symmetrise

```bash
xtalforge symmetry magnetic magnetic.cif --json
```

The result is the UNI number, BNS/OG symbols, the magnetic operations and the
stabilizer count. Python uses `detect_magnetic_symmetry` and
`symmetrize_moments`; the WASM export returns the same JSON.

## 3. Enumerate orderings

```bash
xtalforge features magnetic-orderings magnetic.cif --supercell 2 --max 8 --json
```

The list is ranked by magnetic symmetry and carries a ready `MAGMOM` string.
Use it in the Phase 18 VASP input:

```bash
xtalforge engines inputs vasp/static magnetic.cif --out work/afm \
  --overrides '{"incar":{"MAGMOM":"1 -1 1 -1"}}'
```

## 4. Neutron and electron patterns

```bash
xtalforge analysis neutron magnetic.cif --geometry tof --path-length 15 --json
xtalforge analysis electron magnetic.cif --zone 0,0,1 --json
```

The neutron pattern uses the embedded coherent scattering lengths; the
magnetic pattern adds the moment contribution. The electron pattern is
kinematic SAED with an optional image buffer for the playground.

## 5. IR and Raman

```bash
xtalforge ff vibrations magnetic.cif --json
```

`ir_spectrum_with_born` takes Born tensors from a file or the `qeq`/`eqeq`
diagonal. Raman needs tensors from an engine output:

```python
import xtalforge

spectrum = xtalforge.raman_spectrum(
    modes=modes,
    tensors=xtalforge.RamanTensorSet.from_text(open("raman.txt").read()),
)
```

Phonon DOS and inelastic-neutron weights are available from the same phonon
result. See [Spectroscopy](../algorithms/spectroscopy.md) and
[Neutron/electron diffraction](../algorithms/neutron-electron-diffraction.md)
for the formulas and limits.
