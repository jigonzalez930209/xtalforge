# MOF workflows

Phase 10 turns a CIF into a checked, described framework, and then into
adsorption descriptors. The same functions are available in Rust, the CLI,
Python and WASM; the examples below use the CLI because it shows every
intermediate file.

## From a CIF to a clean, described MOF

```bash
xtalforge mof check framework.cif
xtalforge mof solvent-remove framework.cif --out clean.cif
xtalforge mof disorder-resolve clean.cif --out ordered.cif
xtalforge mof charge-balance ordered.cif
xtalforge mof deconstruct ordered.cif
xtalforge mof mofid ordered.cif
xtalforge mof descriptors ordered.cif
xtalforge mof open-metal-sites ordered.cif
```

`mof check` reports overlaps, fragments, valence problems, disorder, missing
metal coordination and the charge verdict. `mof solvent-remove` removes free
water, DMF, ethanol and the other formulas in the embedded table while
keeping bound solvent and, by default, counter-ions. The deconstruction lists
nodes, linkers, connection points and the contracted net; `mof mofid` prints
the canonical identity string and its SHA-256 MOFkey.

## Screen a database for CO2 Henry coefficients

```bash
xtalforge mof database read core-mof/ --kind core
xtalforge mof database screen core-mof/ --guest CO2 --limit 100 --samples 2000
```

The screening record contains the MOFkey, descriptors, PLD/LCD/ASA/AV/VF,
channel dimension, the Widom Henry coefficient and the checker counts. Core
MOF and QMOF directories are read with their license and citation metadata;
the reader accepts the published directory tree (extract a zip first).

## Build and functionalise UiO-66

```bash
xtalforge mof build pcu zr6 benzene_di --out uio.cif
xtalforge mof modify functionalize uio.cif --group amino --fraction 0.5 --seed 7 --out uio-nh2.cif
xtalforge mof modify defects uio.cif --linker-fraction 0.1 --cap formate --seed 7 --out uio-defect.cif
xtalforge mof modify interpenetrate uio.cif --folds 2 --out uio-2fold.cif
xtalforge mof modify stack-layers uio.cif --layers 3 --offset 0.5,0.5 --out uio-stack.cif
```

## Adsorption

```bash
xtalforge mof adsorption guests
xtalforge mof adsorption force-field clean.cif
xtalforge mof adsorption widom clean.cif --guest CO2 --temperature 298 --samples 5000
xtalforge mof adsorption gcmc clean.cif --guest CO2 --pressures 1,5,10 --steps 2000
xtalforge mof adsorption grid clean.cif --guest CO2 --spacing 1.0
xtalforge mof adsorption raspa clean.cif --guest CO2 --dir raspa-run/
```

The Widom result is a dimensionless Henry coefficient (the average Boltzmann
factor) and a heat of adsorption. GCMC is a single-component, rigid-host,
rigid-guest engine with the ideal-gas chemical potential; use the RASPA export
when the full engine, mixtures or framework flexibility are required.
