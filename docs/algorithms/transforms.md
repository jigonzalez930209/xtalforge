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
and scale by the chi-square quantile for three degrees of freedom (0.5, 0.9,
0.95, 0.99 with linear interpolation). A non-positive-semidefinite tensor is
rejected.
