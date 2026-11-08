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

The integer image range comes from the smallest interplanar spacing, so small
cells and skewed cells include the repeated images inside the cutoff. Results
are sorted by `(i, j, image)`, which makes them deterministic. Self-images
(`i == j`, nonzero image) are included; the zero image is not.

## Status

The direct periodic search is complete and exact. The linked-cell
acceleration and criterion benchmarks from Task 1.10.2 are pending; they must
preserve this output exactly. A rayon build will use a fixed reduction order.
