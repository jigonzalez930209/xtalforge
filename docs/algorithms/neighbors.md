# Periodic neighbor search

For each pair `(i, j)` the search enumerates integer images `n` and computes
`|(x_j + n - x_i) M|`. The image range comes from the smallest interplanar
spacing `d_min = min(1/|a*|, 1/|b*|, 1/|c*|)`:
`range = ceil(cutoff / d_min) + 1`, which is conservative for skewed cells.
Self-images with a nonzero image are included; the zero image is not.

Results are sorted by `(i, j, image, distance)` for determinism. The direct
search is exact but O(N^2 * range^3); the linked-cell acceleration of Task
1.10.2 must preserve this output byte for byte. A rayon build will use a fixed
reduction order.
