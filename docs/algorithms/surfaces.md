# Crystal surfaces

**Tasks 9.4.1-9.4.4.** References: pymatgen `core.surface`,
`analysis.wulff`, `analysis.adsorption`; Wulff (1901); Zur and McGill (1984).

## Symmetry-unique Miller indices

The point group comes from a symmetry search of the structure. Every reduced
`(h k l)` inside the index bound is transformed by each rotation and by its
negative (the Laue class always contains inversion); the family is the set of
distinct transformed triples. The representative is the family member with the
largest positive component count and the largest absolute components, which
picks `(1 0 0)`, `(1 1 0)` and `(1 1 1)` for a cubic crystal. Multiplicity is
the size of the family and the d-spacing follows the reciprocal metric
tensor.

## Slab terminations

The atomic layers are the fractional heights in the surface basis
(`surface_basis`, Task 1.8), taken modulo one. Every gap between consecutive
layers defines a cut shift at the midpoint; `make_slab` builds the slab for
each shift. A signature over `(element, z, x, y)` removes cuts that produce the
same decoration. A symmetric slab mirrors every atom across the middle of the
film. The net dipole is `sum q_i z_i` along the surface normal, with the
oxidation states of Task 9.3.3 as charges when available, and a dipole above
`1e-4 e Angstrom` marks a polar surface. The minimal slab for a thickness
uses `ceil(thickness / d_hkl)` layers.

## Wulff construction

For a surface energy `gamma(h k l)`, the facet plane lies at distance
`gamma / |g|` from the origin, where `g` is the reciprocal-lattice vector of
`(h k l)`. The Wulff polyhedron is the intersection of the half-spaces
`n . x <= gamma / |g|` for every direction; a finite set with positive
energies is bounded and contains the origin. The half-space clipper of the
Voronoi module reports the vertices, faces, volume and area. Each facet is
matched to the direction with the largest normal cosine, and its fraction is
the facet area over the total surface area. Lower energies therefore produce
larger facets, as the construction requires.

## Adsorption sites

The top layer is the set of atoms within a fractional tolerance of the maximum
`z`. Top sites sit on each surface atom; bridge sites at the midpoint of two
surface atoms within `1.25 * (r_i + r_j)` covalent radii; hollow sites at the
centroid of three pairwise-bonded surface atoms. Duplicate sites are removed
by position. `place_molecule` translates the molecule centroid to the site
plus an offset along the surface normal and appends its atoms in the slab
cell.
