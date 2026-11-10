# MOF builder

Reference: PORMAKE 0.2.3 as vendored in Crystal-Builder (MIT). The crate keeps
the algorithm shape but ships a curated topology/block subset; it is not the
full PORMAKE database.

## Data model

- **Topology**: a finite periodic quotient. Nodes are fractional coordinates in
  a cell `[a, b, c, alpha, beta, gamma]`; edges are undirected index pairs.
- **Block**: a small molecule in Cartesian Angstrom with an origin and outward
  connection-point vectors. A block is a node or a linker; its `cell` matrix is
  a display bounding box, not a placement cell.
- **Build**: copy a node to every topology vertex, orient it, place a linker on
  every edge, then scale the result.

## Node orientation

![MOF assembly](/diagrams/mof-assembly.svg)


PORMAKE fits the connection directions to a slot by RMSD over all rotations.
This port uses a deterministic subset:

1. Local reference = normalised mean of the node connection vectors. If that
   mean is zero (octahedral or trigonal), use the lowest-index connection.
2. Target = normalised mean of the incident edge directions; if zero, use the
   lowest-index edge.
3. Rotation = minimal proper rotation between those directions. The
   antiparallel case uses a fixed perpendicular axis, so it is deterministic.
4. Greedy matching: the connection/edge pair with the largest cosine is
   assigned first; ties use the index order.

The difference is recorded because an RMSD fit picks the global best over all
rotations, while the greedy fit is cheaper and reproducible. For the curated
symmetric topologies the two agree on the resulting connectivity; the geometry
can differ by a rotation of a node.

## Linker placement and scale

For an edge `i-j`, the node connection points `p_i` and `p_j` are known. The
linker connection pair `(a, b)` with the largest separation is aligned to
`p_j - p_i`; the linker is translated so `a` lands on `p_i`. The edge scale is
the ratio of the node-to-node distance to the linker span, clamped to the
documented range. A scale outside the range is reported, not silently applied.

## Assembly

`build_mof` returns a P1 structure with the node and linker atoms,
explicit topology bonds and metadata naming the topology and blocks. It does
not relax the geometry. Interpenetration, if requested, is a separate operation
over the assembled document.

The advanced subset adds rods, layers, merged nets, decoration and invented
linkers. Each is a geometric construction with its own documented limitations;
the pore partition reports the accessible volume per framework rather than
claiming a topology assignment.

## Validation

- `crates/xtalforge-mof/tests/mof.rs` checks every curated block's connection
  count and span, every topology's edge/node consistency, a build for each
  topology, determinism, and that a bad block name or a scale out of range is a
  typed error.
- `crates/xtalforge-mof/tests/advanced.rs` covers rods, layers, merged nets and
  decoration.
- The builder is not validated against PORMAKE's full output; the subset is
  illustrative and the test records that.

## Differences from PORMAKE

- Curated topology and block subset, not the full database.
- Greedy deterministic orientation instead of a global RMSD fit.
- The result is P1 and unrelaxed; there is no automatic symmetry assignment.
- The builder refuses a block with the wrong connection count instead of
  repairing it.
