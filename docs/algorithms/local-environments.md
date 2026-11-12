# Local environments and chemistry

**Tasks 9.3.1-9.3.4.** References: pymatgen `analysis.local_env`,
`analysis.bond_valence`; O'Keefe and Brese, *JACS* 113 (1991) 3226; ICSD
bond-valence data mining as shipped by pymatgen (MIT).

## Periodic Voronoi

The Voronoi cell of an atom is the intersection of the bisector half-spaces

```
2 (r_j - r_i) . x <= |r_j - r_i|^2
```

for all neighbors within the cutoff, in coordinates relative to the central
atom. The half-spaces are clipped from a large cube face by face
(`xtalforge_core::polyhedron`); every new cap face is ordered around the plane
normal. The volume follows from the signed tetrahedra of the outward faces,
and the solid angle of a face is `area * cos(theta) / d^2`. Face planes are
matched back to the neighbor that produced them by the normal cosine, so each
neighbor carries its area and solid angle. The cutoff must include every face
of the cell; at least twice the longest nearest-neighbor distance is safe.

## Coordination and shape measures

Three strategies are available:

- `cutoff`: every neighbor within a distance.
- `voronoi`: faces with a solid angle above a threshold.
- `crystalnn`: a CrystalNN-like weight `exp(-((d - d_min)/0.4)^2) * sqrt(Omega/pi)`.

Continuous shape measures compare the neighbor vectors with ideal polyhedra
of the same size. Both sets are centred and normalized to unit RMS radius; the
best rotation comes from Kabsch alignment and the vertex correspondence from
the Hungarian algorithm on the squared-distance cost matrix, iterated four
times. The reported value is

```
CShM = 100 * sum |q - p|^2 / sum |p|^2
```

which is zero for the ideal shape. The ideal sets are linear, trigonal planar,
tetrahedral, square planar, trigonal bipyramid, square pyramid, octahedral,
pentagonal bipyramid, cube, square antiprism and icosahedral.

## Bond valence and oxidation states

The bond valence of a pair is

```
v_ij = exp((R_ij - d_ij) / 0.31)
R_ij = r_i + r_j - r_i r_j (sqrt(c_i) - sqrt(c_j))^2 / (c_i r_i + c_j r_j)
```

with the 1991 `r` and `c` parameters. Only pairs where at least one element is
in the O'Keefe electronegative set contribute, and the sign follows the
electronegativity order. The site sum is signed.

Oxidation states use a maximum-a-posteriori search: for every element the
mean bond-valence sum is compared with the ICSD Gaussian for each candidate
state, multiplied by its occurrence prior. A depth-first search over the
element states keeps the combination with the smallest negative log posterior
plus `10 * |total charge|`. The state is then assigned to every atom of that
element. Cells that are not neutral are reported, never hidden.

## RDF, PDF and fingerprints

The RDF is a pair histogram with the normalization
`g(r) = n(r) / (N rho 4 pi r^2 dr)`; partial functions use the concentration
of the second element. A moving average can smooth every curve. The running
coordination number integrates `rho g(r)` from zero. Bond angles use all
neighbor pairs of every atom. The fingerprint stores the element fractions,
the six cell parameters and the smoothed total RDF; its similarity is
`0.4 * composition + 0.2 * cell + 0.4 * RDF` cosine terms.
