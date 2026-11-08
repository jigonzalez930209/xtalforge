# Measures

Distance, angle and torsion use the minimum image of the P1 cell, so a
measurement across a boundary equals the same local geometry inside the cell.

For a pivot `j`, the arms are `minimum_image(x_i - x_j)` and
`minimum_image(x_k - x_j)`, converted to Cartesian. The angle is
`acos(u·v / |u||v|)`. The torsion uses the two plane normals
`n1 = b1 × b2` and `n2 = b2 × b3` and the signed angle around `b2`.

Degenerate arms, collinear segments and zero normals return
`Error::InvalidInput` rather than a silent NaN. Centroids and interplanar
angles are plain Cartesian operations.
