# Magnetic symmetry and orderings

## Magnetic moments

A collinear moment is an axial scalar with a sign; a non-collinear moment is a
Cartesian vector. The site stores the moment next to occupancy and ADP data, so
the P1 expansion carries it unchanged. The mCIF writer emits crystal-axis
components because that is the conventional magnetic CIF frame.

## Magnetic space groups

`detect_magnetic_symmetry_with_options` passes the P1 cell and the per-atom
collinear moments to moyo's magnetic cell. The returned `MagneticDataset`
carries:

- the UNI number in `1..=1651`;
- the BNS and OG symbols when the backend provides them;
- the magnetic operations, including the time-reversal part;
- the action on the moments.

A site without a moment is treated as non-magnetic (`0`). The search tolerance
is the ordinary `symprec` in Angstrom; the magnetic part adds no separate
length tolerance.

## Symmetrisation

`symmetrize_moments` applies every magnetic operation to the moment array and
keeps the value when the operation maps the site onto itself or onto a
symmetry-equivalent site. The result reports how many operations leave the
array invariant. A structure whose moments are already symmetric has the full
operation count; a lower count means the input breaks part of the group. The
function refuses a non-collinear array when the detection was collinear and
vice versa.

## Orderings

The enumeration starts from a supercell of the requested size. For every
symmetry-distinct assignment of the available magnetic species it builds the
moment array, detects the magnetic group and ranks by:

1. number of magnetic operations (higher is more symmetric);
2. the kind order ferro, ferri, antiferro (a deterministic tie-break, not a
   physical energy);
3. the assignment bit pattern.

The search is bounded; `max_orderings` is honoured and the result is a
deterministic prefix, never a random sample. The `MAGMOM` helper writes the
collinear values in POSCAR order, which is what the Phase 18 VASP preset reads.

## Limits

- No magnetic exchange constants or spin dynamics; ordering energies are not
  evaluated here.
- mCIF propagation-vector symmetry is read from the operations but no
  propagation-vector refinement is performed.
- The rank is symmetry-first, not energy-first.
