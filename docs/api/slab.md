# Slabs

`make_slab(structure, hkl, layers, vacuum, shift)` cuts `layers` spacings of
the plane `(hkl)` and adds `vacuum` Angstrom above. The surface cell is the
plane's own lattice: `surface_basis(hkl)` returns an integer, right-handed,
unimodular basis with the first two rows in the plane and `h . v3 = 1`.

```rust
use xtalforge::core::{make_slab, surface_basis};

let basis = surface_basis([1, 1, 1])?;
let slab = make_slab(&structure, [1, 1, 1], 2, 15.0, 0.0)?;
assert!(slab.structure.is_p1());
```

The returned `Slab` carries `structure`, the reduced `hkl`, `layers`,
`thickness`, `vacuum` and `cut`. The structure is P1. `SLAB_FLOOR` keeps the
slab a little above the bottom of the cell so atoms on the face are not drawn
twice. `plane_angle(first, second)` returns the angle between two Miller
planes.

## Terminations and adsorption depth

`slab_cut_shifts(structure, hkl, &options)` returns the sorted cut shifts.
`slab_terminations(structure, hkl, layers, vacuum)` returns the midpoint cuts
of Task 9.4 as a `SlabTerminationReport { hkl, terminations }`, where each
`SlabTermination { shift, slab }` pairs a cut with its slab.
`slab_terminations_with_options(..., &TerminationOptions)` adds
`cuts` (extra shifts in interplanar-spacing units), `include_intra_layer`
(quarter and three-quarter cuts inside each gap) and `layer_tolerance`.
`adsorption_sites_with_options(slab, tolerance, &options)` reads
`TerminationOptions::site_layers` to enumerate sites down to that many atomic
layers below the top; the default of one reproduces `adsorption_sites`. The
CLI exposes `--intra-layer`, `--cut` and `--site-layers`.

## Bond carrying

The slab builder maps every perceived or explicit bond onto its layer copies
using the surface basis. A bond whose partner lies outside the cut is carried
on the surface copies and counted in `cut`; a bond crossing the periodic
boundary counts at both surfaces. Explicit bond orders and suppressions
survive. No atom capping is added: the caller decides which chemistry to apply
to the dangling valences. The CIF inspector demo downloads the same slabs
through the WASM bindings.
