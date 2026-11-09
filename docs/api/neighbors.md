# Neighbors

`neighbors(structure, cutoff)` returns every periodic pair within `cutoff`
Angstrom as `Neighbor { i, j, image, distance }`, where `image` is the lattice
translation of atom `j` relative to atom `i`. `neighbors_of` filters the list
for one atom. The structure is expanded to P1 first.

```rust
let found = neighbors(&structure, 3.2)?;
for neighbor in &found {
    println!("{} - {} at {:.3} A", neighbor.i, neighbor.j, neighbor.distance);
}
```

The search uses a linked-cell index in fractional space. Bins are sized from
the perpendicular widths of the cell, so triclinic cells are exact; each
candidate pair then enumerates the integer images that can reach the cutoff.
Small cells and skewed cells include the repeated images inside the cutoff.
Results are sorted by `(i, j, image)`, which makes them deterministic.
Self-images (`i == j`, nonzero image) are included; the zero image is not.

The benchmark target covers 10^3, 10^4 and 10^5 atoms:

```bash
cargo bench -p xtalforge-core --bench neighbors
```

See [Periodic neighbor search](../algorithms/neighbors.md) for the index
construction, the complexity and the comparison with the Crystal-Builder
implementation.
