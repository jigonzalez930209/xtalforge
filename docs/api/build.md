# Carbon builders

`xtalforge-build` starts the carbon family:

```rust
use xtalforge_build::{graphene, nanotube};

let sheet = graphene(4, 4, 1, 3.35, 15.0)?;
let tube = nanotube(5, 5, 2, 10.0)?;
```

`graphene(nx, ny, layers, interlayer, vacuum)` builds an AA-stacked honeycomb
sheet in a hexagonal cell (C-C = 1.42 Angstrom). `nanotube(n, m, length,
vacuum)` rolls the honeycomb along the chiral vector; the current builder
uses a direct lattice scan and buckets atoms into the tube cell, which is
correct for the common armchair/zigzag cases and is validated by atom counts
and cell dimensions.

`fullerene(60)` builds C60 from the embedded ASE-derived coordinates; the
result is a padded molecule with 60 atoms and 90 perceived bonds.

The remaining builders (SMILES molecules, fragments, fullerenes, amorphous
carbon, polymers) are later Phase 6 tasks.
