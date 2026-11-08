# Slab construction

`surface_basis(hkl)` performs Euclid's algorithm on the three Miller indices as
column operations on the identity. The result is an integer, right-handed,
unimodular basis whose first two rows lie in `(hkl)` and whose third row
satisfies `h . v3 = 1`.

The slab lattice is `[a, b, normal * sign(rise) * height]`, where `height =
layers * |rise| + vacuum`. Atoms are placed in surface coordinates: one layer
is one unit along the third basis vector, the in-plane component of `c` is
removed by solving it in the `a,b` plane, and the vacuum becomes the final
fractional range. `SLAB_FLOOR` keeps the slab off the cell face so atoms on
the face are not drawn twice.

The cut is P1. Bond carrying through the cut needs the bond graph of Task 1.11
and is not in this pass; `Slab.cut` is zero until then.
