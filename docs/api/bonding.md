# Bonding

`Structure::bond_graph()` perceives distance-based bonds on the P1 cell with
`BondRules::default()`: a bond exists when the distance is at least 0.4 A and
at most 1.15 times the sum of the two covalent radii. Dummy `X` atoms and
metal-metal pairs never bond. Each `Bond` carries `i`, `j`, the lattice `image`, `distance`, `order`
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

## Explicit and topology bonds

`ExplicitBond { i, j, image, order, remove }` is applied after perception: a
removal deletes the matching perceived bond, and an addition inserts it (or
updates its order) and marks it `explicit`. `TopologyBond { i, j, image }` is
stored separately in `Structure::topology_bonds()`; it is a network edge, not
chemistry, and never enters the chemical bond graph. Both types are
serialized with the structure.

## Status

Distance perception, explicit additions/removals, topology edges and the
connectivity graph are complete. The full golden corpus comparison of bond
sets and coordination numbers is in `crates/xtalforge-io/tests/golden_bonds.rs`
(Task 1.11.5). Bond symmetry propagation to every image (1.11.2) remains a
later refinement.
