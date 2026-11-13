# Meshes, normals and export

The view crate builds triangle meshes for isosurfaces, molecular surfaces and
the scene. Phase 14 adds vertex normals and the two exchange formats.

## Normals

Each triangle contributes its area-weighted normal to its three vertices; the
accumulated vertex normals are normalised after the pass. Degenerate triangles
(zero area) contribute nothing and keep the vertex normal at zero, which is
reported rather than hidden. The result is deterministic: triangle order,
summation order and the final normalisation are fixed.

```text
n_t = (p1 - p0) x (p2 - p0)
n_v = normalize(sum_{t incident to v} n_t)
```

## OBJ

`write_obj(mesh)` emits `v`, `vn` and `f` records with one-based indices and
`f v//vn` corners. Vertex order is preserved. The output carries no material
or texture references, so it imports into Blender, MeshLab and PyVista without
extra files.

## glTF

`write_gltf(mesh)` emits glTF 2.0 JSON with one mesh, one primitive
(`TRIANGLES`), `POSITION` and `NORMAL` accessors and a `UNSIGNED_INT` index
buffer. `write_glb(mesh)` embeds the same JSON and binary buffer in the GLB
container. Both use the normalised normals above; the bounding box is computed
from the positions.

## Validation

- `sum |n_v| = 1` for every non-degenerate vertex.
- OBJ and glTF round trip the position, normal and index arrays.
- A GLB file has the `glTF` magic, version 2 and a valid JSON chunk length.
- Degenerate triangles are excluded and the mesh still parses.
