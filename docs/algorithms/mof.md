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

## Cell scale and linker placement

The catalogue cell of a topology fixes the node-to-node distances, but the
linker decides how long an edge has to be. Before placing anything the builder
therefore **scales the whole cell uniformly** by the factor `s` that minimises
`sum_edges (gap_e(s) / span - 1)^2`, where `gap_e(s) = |s d_e + (R_j p_j - R_i p_i)|`
is the distance the linker has to bridge between the two oriented node
connection points (`d_e` is the unscaled node-to-node vector) and `span` is
the largest distance between two connection points of the linker. The search
is a golden-section on `[0.1, 20]` and is deterministic; the factor is stored
in the structure metadata as `cell_scale`. A longer linker (stilbene instead of
benzene in `pcu`) therefore gives a larger cell instead of a squashed linker.

For each edge `i-j` the linker pair `(a, b)` with the largest separation is then
aligned to `p_j - p_i`, translated so `a` lands on `p_i`, and scaled by
`gap_e / span`. After the fit that residual scale is close to one when every
edge of the topology has the same length. When it is not (for example `srs`
with a paddlewheel), the builder records one warning if the scales fall outside
`[0.8, 1.25]` (`LINKER_SCALE_MIN`, `LINKER_SCALE_MAX`); the scale is still
applied, so the warning is the signal that the linker is distorted.

## Assembly

`build_mof` returns a P1 structure with the node and linker atoms,
explicit topology bonds and metadata naming the topology, the blocks and the
`cell_scale`, plus the list of warnings. It does
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
- A node whose connection count differs from the vertex degree is **not**
  refused: the edges without a connection point are skipped and one warning
  names the number of vertices affected. A linker with more than two
  connection points spans each edge with its widest pair and is reported.
- A final check lists atom pairs closer than 0.9 Angstrom (`CLOSE_CONTACT`) in
  one warning; every combination of the curated blocks builds, and the
  combinations that do not fit chemically carry at least one of these warnings.
