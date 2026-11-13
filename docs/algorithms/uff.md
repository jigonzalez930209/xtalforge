# UFF and UFF4MOF

`xtalforge-ff` implements the Universal Force Field natively. The parameter
table is `data/uff.json.gz`: 218 rows ported from Crystal-Builder
`xtal/ff/uff/params.py` (Rappé 1992, Addicoat 2014, Coupry 2016), with the
reference units (Angstrom, degrees, kcal/mol). The public boundary converts to
eV and eV/Angstrom.

## Typing

An atom type is `(element, geometry, oxidation state)`. The typer combines
the element table, the perceived bond graph and the local geometry:

- explicit rules for H, C, N, O, F, Cl, Br, I, S, P and the common metals;
- aromatic ring bonds use the resonant carbon/nitrogen types;
- terminal atoms refine their geometry from the number of neighbours;
- framework metals prefer the UFF4MOF rows when the set is `Uff4Mof`;
- a dummy `X` atom has no parameters and never contributes.

`Typing` records the chosen `AtomType` and a `Confidence`, so a caller can see
why an atom got a type. `UffModel` exposes the prepared topology.

## Energy terms

![UFF energy terms](/diagrams/uff-terms.svg)


- Bond: `0.5 k (r - r0)^2`, `k = 664.12 z1 z2 / r0^3`.
- Angle: the reference cosine expansion whose curvature at `theta0` is the UFF
  `K_IJK`.
- Torsion and inversion: the reference Fourier terms.
- van der Waals: Lennard-Jones with the UFF `x` and `D` parameters, 1-2 and
  1-3 exclusions, a cutoff and a shifted potential.

Every term has an analytic gradient. The finite-difference test moves each
coordinate by a small Cartesian step and compares forces within `1e-4 eV/A`.

## Validation

- Ethylene energy `0.1113230794 eV`, matching the reference term-by-term sum.
- Diamond periodic energy `1.7200457169 eV`.
- Argon dimer shifted well `-0.0080043731 eV`.
- Methane at its UFF bond length, water at its angle, a benzene resonant ring
  and a Zn4O UFF4MOF node.

The Crystal-Builder golden comparison over the full corpus is not run in CI;
the values above are the regression anchors.
