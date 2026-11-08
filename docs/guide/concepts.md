# Concepts

## Asymmetric unit and P1

A `Structure` is a lattice, a space-group identity (the Hall symbol) and the
asymmetric-unit sites. The P1 cell is derived: every atom carries its parent
site, the operation that generated it and the lattice shift. An edit to a P1
atom therefore maps back to the site that owns the orbit.

## Hall symbol

The Hall symbol is the identity of a space group. A number alone does not fix
the setting: 227 has two origin choices, 146 has hexagonal and rhombohedral
axes. Until moyo supplies the operations, a structure carries the explicit
operation list read from the file or supplied by a builder.

## Coordinates

- Fractional coordinates are `[f64; 3]`.
- The lattice matrix rows are the vectors a, b, c.
- Fractional to Cartesian is `frac @ M`.
- The standard orientation puts a along +x and b in the xy plane.
- The reciprocal lattice uses the crystallographic convention: `a . a* = 1`
  and `d(hkl) = 1 / |hkl @ reciprocal|`.

## Units

Angstrom and degrees at every API boundary; radians internally. Energies in
eV, forces in eV/Angstrom and stress in GPa in Phase 4.

## Limits

Every large allocation is estimated first (drawn atoms, supercell atoms, grid
points). The soft limit warns; the hard limit refuses with a sentence that
says what would fit. The profile is Standard, Generous or WarnOnly.

## Bonds and rings

Bonds are perceived from covalent-radius sums on the P1 cell, then explicit
additions/removals and topology edges are applied. Rings are chordless
primitive cycles on the bond graph. Dummy `X` atoms never bond.
