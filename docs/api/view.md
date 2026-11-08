# View data

`xtalforge-view` has no windowing dependency: a viewer consumes `SceneData`.

```rust
use xtalforge::view::{DrawStyle, ViewSettings, build_scene, build_mesh, write_stl};

let settings = ViewSettings { style: DrawStyle::SpaceFill, ..Default::default() };
let scene = build_scene(&structure, &settings)?;
let stl = write_stl(&build_mesh(&scene), "structure");
```

## Settings

`ViewSettings { style, atom_scale, bond_radius, show_cell, show_axes,
background, colormap, hide_hydrogens }` with `sanitize`. `DrawStyle` is
`BallStick`, `SpaceFill`, `Wireframe` or `Sticks`.

## Colormaps

`Colormap::{Viridis, Plasma, Inferno, Magma, CoolWarm, Grayscale}` samples a
0-1 value into 0-255 RGB.

## Scene

`scene_summary` returns atom, bond and cell-edge counts plus the element
list and radii for the CLI, Python and WASM. `mesh_stl` returns the full ASCII
STL string.

`SceneData { atoms, bonds, cell_edges, settings }`: atoms carry element,
Cartesian position, display radius, colour, index and a hidden flag; bonds
carry endpoints, radius, colour and order; cell edges are 12 Cartesian
segments.

## STL

`build_mesh` turns visible atom octahedra and bond prisms into triangles and
`write_stl` writes ASCII STL.

## CLI

```bash
xtalforge mesh structure.cif --style space-fill --out structure.stl
```

`--style` accepts `ball-stick`, `space-fill`, `wireframe` or `sticks`;
`--scale` multiplies the atom radii. Without `--out` the ASCII STL goes to
standard output.
