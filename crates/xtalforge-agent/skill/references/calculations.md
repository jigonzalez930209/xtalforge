# Reference: calculations

## Symmetry

`symmetry` runs the moyo backend on the P1 expansion. `symprec` is a Cartesian
tolerance in Angstrom (default 1e-5). The answer has the ITA number, the
Hermann-Mauguin and Pearson symbols, the operation count, the Wyckoff letter
and site-symmetry symbol of every input atom. A larger `symprec` merges sites
that a deposited file left slightly off; report the value you used.

## Bonds and rings

`bonds` perceives bonds with the covalent-radius rule (scale 1.15, minimum
0.4 A) unless the structure stores explicit bonds. Images are lattice shifts
of the second atom. `rings` counts primitive rings up to `max_size` (default
12) with the periodic shortest-path search. A ring count is a graph property,
not a chemical aromaticity assignment.

## PXRD

`pxrd` sums Gaussian peaks over the allowed reflections of the space group.
Defaults: Cu Ka, 5-80 degrees, 0.02 degree step, 0.1 degree FWHM. Intensities
are normalized to 100 and the first 240 reflections are returned with the
decimated profile. The pattern is calculated from the model; it has no
instrument broadening beyond the Gaussian FWHM.

## UFF

`uff_energy` builds the UFF/UFF4MOF model over the P1 cell, including periodic
van der Waals images inside the default cutoff. The answer is the total energy
in eV, the largest force in eV/Angstrom and the stress in GPa. `optimize`
relaxes the geometry with FIRE (default) or L-BFGS and replaces the document.
The run reports the final energy, steps, convergence flag and largest force.

Caveats: UFF is a classical force field; do not use it for bond breaking or
electronic properties. An element without UFF parameters is an error, not a
silent zero.
