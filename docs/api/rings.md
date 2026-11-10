# Rings

`rings(structure, max_size)` finds the primitive (shortest-path) rings of
Franzblau on the periodic bond graph. A cycle with a shortcut through the graph
is a sum of smaller rings, so it is not primitive. `ring_counts` aggregates by
size, `ring_membership` counts rings per atom and `is_aromatic` is a labelled
heuristic (five- or six-membered, all C/N/O) until bond orders and planarity
are available. `rings_explicit` and `ring_counts_explicit` use only the stored
explicit bonds, which is the census a builder states: distance perception is
not allowed to add a close contact across a pore.

```rust
use xtalforge::core::rings;

for ring in rings(&structure, 8)? {
    println!("{} atoms:", ring.size());
    for step in &ring.path {
        println!("  atom {} image {:?}", step.atom, step.image);
    }
}
```

`Ring { atoms, path }` carries both forms: `atoms` is the sorted set of unique
P1 indices for membership tests, and `path` is the ring in walking order as
`RingStep { atom, image }`. Add `image` to the atom's fractional position to
lay the path out in Cartesian space; a small cell can therefore show the same
atom more than once, and `size()` is the path length. `ring_paths` returns only
the paths.

The search is periodic-aware because the graph is built on the P1 cell and the
cycle walk carries lattice images. Results are sorted deterministically, and
the same structure and limit return the same list on every surface.

See [Primitive rings](../algorithms/rings.md) for the walk, the Franzblau
filter, the path budget and the golden-corpus validation.
