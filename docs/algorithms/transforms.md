# Cartesian transforms and ellipsoids

Point transforms are free functions over Cartesian points: translation,
rotation about an axis through a centre (Rodrigues), mirroring through a
plane, inversion through a centre, alignment of one vector to another, and
plane fitting from at least three non-collinear points.

The finite cell duplicates every atom whose fractional coordinate is within
the tolerance of a face, edge or corner, using a coordinate-wise duplicate
test (not periodic distance, which would merge 0 and 1).

Displacement ellipsoids diagonalize the symmetric Uij tensor with Jacobi
rotations, sort the eigenvalues, take square roots for the semi-axis lengths
and scale by `sqrt(chi2_3(p))`, the square root of the chi-square quantile for
three degrees of freedom: 1.538, 2.500, 2.796 and 3.368 at `p` = 0.5, 0.9, 0.95
and 0.99, linearly interpolated between them and clamped to the end values
outside `[0.5, 0.99]`. A non-positive-semidefinite tensor is
rejected.
