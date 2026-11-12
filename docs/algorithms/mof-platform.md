# MOF deconstruction, open metal sites and native porosity

**Tasks 10.1-10.3.** References: MOFid (Bucior et al., *Mol. Syst. Des. Eng.*
2019), MOFChecker/MOSAEC quality flags, CrystalNets.jl (topology), Zeo++
(PLD/LCD/ASA), CoRE MOF 2019/2024 and QMOF.

## Cleaning

The checker starts from the P1 expansion and the perceived bond graph. It
reports overlapping atoms (distance below 0.4 Angstrom), unbound atoms and
fragments (connected components other than the largest), impossible valences
from the sum of bond orders, hypervalent carbon, missing hydrogens on
under-coordinated carbon, partial occupancies, metals below four bonds and,
for metal-containing frameworks only, the charge balance of Task 9.3.3.
Every finding has a severity (`info`, `warning`, `error`), the involved atoms
and a suggested fix verb.

Solvent removal takes the largest component as the framework. Each other
component is matched against the embedded formula table (water, alcohols,
DMF, DMF analogues, DMSO, acetonitrile, acetone, THF, ethers, chlorinated
solvents, benzene, toluene, hexane, cyclohexane); a small organic fragment
without a metal is treated as solvent even when it is not in the table. Free
counter-ions are kept unless `remove_ions` is set. Bound fragments are never
removed because they belong to the framework component.

## Deconstruction

The "all-node" rule assigns every metal atom and every non-carbon atom
directly bonded to a metal (plus its hydrogens) to a node. Metals that share a
non-carbon bridge are joined into one node (union-find). The remaining atoms
of the main component form linkers. A linker with two connections becomes an
edge of the contracted net; a node or a higher-connected linker becomes a net
vertex. Edges carry the lattice image that makes the two node centroids
nearest neighbours. The quotient graph is identified against the embedded
RCSR index by its coordination sequence and point symbol; ties are reported
as a name set. The number of bond-graph components that contain a metal is
the interpenetration count.

MOFid is a deterministic subset of the published syntax:
`MOFid:{topology}:{nodes}|{linkers}|{n}`. Each unit is turned into a graph
with external bonds cut, atom colours are refined five times from
`(element, degree, neighbour colours)`, the atoms are sorted by colour and a
depth-first SMILES writer emits the string. MOFkey is the SHA-256 of MOFid.

## Open metal sites

The coordination shell comes from the Voronoi tessellation with a 3.5
Angstrom cutoff. A metal is open when its coordination number is below the
typical value of its family (6 for first-row transition metals, 7 for
alkaline earths, 8 for lanthanides) and the vector sum of its neighbour
directions leaves a non-zero free direction. The clearance along that
direction is read from the distance grid; a probe of 1.2 Angstrom fits when
the clearance reaches that radius.

## Native PLD, LCD and ASA

The distance grid stores `distance to the nearest atom surface`. LCD is
`2 * max(field)`. PLD is the largest probe radius for which a connected
component of samples spans the cell in at least one axis. Samples are added in
decreasing field order; a union-find joins active periodic neighbours and
records which of the six faces the component touches. The first threshold at
which a component touches opposite faces is the PLD radius, so the PLD is
`2 * radius`. Channel dimensionality repeats the union-find at a fixed probe
radius and counts the axes with a percolating component.

Accessible volume is the number of samples with field at least the probe
radius times the sample volume. Void fraction is that volume over the cell
volume. Accessible surface area counts the grid faces whose two sides straddle
the probe threshold and multiplies by the face area; the area is reported per
mass (`A / (m * 1.66054) * 1e4`) and per volume (`A / V * 1e4`). These are
voxelized estimates that converge from below as the spacing shrinks, exactly
like the generic distance grid of Task 5.2.
