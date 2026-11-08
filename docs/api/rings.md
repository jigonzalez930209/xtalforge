# Rings

`rings(structure, max_size)` finds chordless simple cycles on the perceived
periodic bond graph. A cycle with a chord is a sum of smaller rings, so it is
not primitive; this is the Franzblau criterion applied to the molecular graph.
`ring_counts` aggregates by size.

```rust
for ring in rings(&structure, 8)? {
    println!("{}: {:?}", ring.size(), ring.atoms);
}
```

The search is periodic-aware because the graph is built on the P1 cell. It
returns sorted, duplicated-free `Ring { atoms }` values. Aromaticity flags and
the full framework ring statistics follow with the molecule and framework
tests.
