# Periodic neighbor search

The search combines a linked-cell index with exact image enumeration. It
returns every contact within the cutoff once per lattice image, including the
self-images of a small cell.

## Linked-cell index

![Linked-cell neighbour search](/diagrams/neighbors-linked-cell.svg)


1. The structure is expanded to P1.
2. The perpendicular width of each axis comes from the reciprocal cell:
   `d_k = 1 / |k*|`. Edge lengths are not enough for a skewed cell.
3. Each axis is divided into `n_k = max(1, floor(d_k / cutoff))` bins in
   fractional coordinates. A pair within `cutoff` satisfies
   `|delta x_k| <= cutoff / d_k <= 1 / n_k`, so it always falls in the same or
   an adjacent bin, wraps included.
4. Each atom is inserted into one bin. Candidate pairs come from the 27
   neighbouring bins (deduplicated when an axis has one or two bins).
5. For every candidate pair the integer images that can reach the cutoff are
   enumerated per axis, and the exact Cartesian distance filters the result.

Self-images (`i == j`, nonzero image) are included; the zero image is not.
When the cutoff is at least one perpendicular width the index degenerates to
one bin and the search is quadratic, which is the only case where this happens.

## Determinism

Results are sorted by `(i, j, image, distance)`; the order does not depend on
hash-map iteration or platform. The same structure and cutoff return the same
list on every surface.

## Benchmarks

The table compares the Rust search with the Crystal-Builder Python
implementation (`xtal.core.neighbors.neighbor_pairs`, scipy cKDTree) on the
same coordinates: a cubic cell at 0.05 atoms/Angstrom^3, cutoff 3.0 Angstrom,
median of five Rust runs and three Python runs. Rust is a release build; Python
is 3.14.4 with numpy 2.4.4 and scipy 1.17.1 on the same Linux x86-64 machine.

| Atoms | Rust | Python | Neighbors |
|---|---|---|---|
| 1,000 | 2.9 ms | 10.2 ms | 2,838 |
| 10,000 | 44.1 ms | 133.1 ms | 28,409 |
| 100,000 | 845.4 ms | 3,563.1 ms | 282,450 |

The Rust side is reproducible with:

```bash
cargo bench -p xtalforge-core --bench neighbors
```

The P1 expansion that every query performs also uses a fractional grid for
duplicate detection, so its cost is linear in the number of generated atoms
instead of the previous quadratic scan.
