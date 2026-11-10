# Primitive rings

## Problem

A crystal's ring statistics must count the rings a chemist sees, not every
simple cycle. Graphene's 576 atoms have 1728 ten-cycles that are sums of
hexagons; those are not primitive rings. A ring may also close through a cell
face, and in a two-atom cell it can contain several copies of the same atom.

## Formulation

![Primitive ring chord test](/diagrams/rings.svg)


Build the periodic bond graph on the P1 cell, with both directions of every
bond carrying its lattice image. For each start atom, walk simple paths of
`(atom, image)` steps up to `max_size`; a step that returns to the start with
net translation zero is a candidate cycle. The walk only visits atoms whose
index is at least the start, which makes the start the smallest atom in the
ring. Each candidate is translated home so its least step has zero image and
then compared as a set, so the two directions and every lattice translate name
the same ring once.

A candidate is kept only if it is primitive: for every pair of ring atoms that
are not neighbours in the path, the graph distance between them (computed by a
breadth-first search up to `max_size / 2` bonds, carrying images) must not be
shorter than their distance round the ring. A chord is the common shortcut and
splits the cycle into two smaller rings, so the cycle is not primitive.

The path budget is `1_000_000 + 500 * n_atoms`; a graph that exhausts it is
refused with an input error instead of hanging. The reference uses the same
floor, per-atom share and depth.

## Complexity

Exponential in the worst case, bounded by `max_size`, the graph degree and the
path budget. With `max_size = 8`, the golden corpus (including 16,000-atom
MIL-101 and the 13,552-atom MIL-100) runs in well under a second per structure
in a release build. The result is deterministic because the bond graph, the
adjacency order and the path sort are deterministic.

## Caveats

- A small cell can repeat an atom in the path. `Ring::atoms` is the sorted set
  of unique indices; `Ring::size()` and `Ring::path` describe the periodic
  path.
- `is_aromatic` remains a heuristic (five/six-membered, C/N/O) until bond
  orders and planarity are available.

## Validation

`crates/xtalforge-io/tests/golden_rings.rs` compares every golden structure
that has ring data against Crystal-Builder at `max_size = 8`. Counts must match
exactly; when the golden list is complete, every periodic path must match after
P1 atom mapping, translation home and forgetting direction. Structures where
the reference refused a too-dense graph must also be refused here. The six
prepared structures above the golden generator's atom limit have no ring list
and are skipped.
