# Functional groups

`find_groups(structure)` matches a small catalogue on the perceived bond
graph and returns `GroupMatch { kind, atoms }`; `group_counts` aggregates by
`GroupKind`. The current catalogue is hydroxyl (O with an H neighbour),
carboxylate (C with two O neighbours), amine (N with one to three neighbours)
and fluoride (C-F).

```rust
use xtalforge::core::find_groups;

for group in find_groups(&structure)? {
    println!("{}: {:?}", group.kind.name(), group.atoms);
}
```

The matcher is periodic-aware because it runs on the P1 graph. The catalogue
is a fixed enum with these four patterns; the remaining patterns from
Crystal-Builder `core/groups.py` and the framework corpus tests are not
implemented. Each current pattern is tested on a small molecule.
