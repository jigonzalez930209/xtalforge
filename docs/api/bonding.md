# Bonding

`Structure::bond_graph()` perceives distance-based bonds on the P1 cell with
`BondRules::default()`: a bond exists when the distance is at least 0.4 A and
at most 1.2 times the sum of the two covalent radii. Dummy `X` atoms never
bond. Each `Bond` carries `i`, `j`, the lattice `image`, `distance`, `order`
and the `explicit` flag.

```rust
let graph = structure.bond_graph()?;
for bond in &graph.bonds {
    println!("{} - {} at {:.3} A", bond.i, bond.j, bond.distance);
}
let components = graph.components();
```

`BondGraph` exposes `bonds`, `n_atoms`, `degree`, `bonds_of` and
`components`. The graph is built over the expanded P1 cell, so periodic bonds
keep their image. `perceive_bonds(structure, rules)` is the configurable
entry point.

## Status

Distance perception and the connectivity graph are complete. Explicit bonds
with symmetry propagation (1.11.2), topology bonds (1.11.3) and the full
corpus comparison (1.11.5) follow with the CIF reader and the command layer.
