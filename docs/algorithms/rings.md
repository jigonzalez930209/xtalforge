# Primitive rings

## Problem

A crystal's ring statistics must count the rings a chemist sees, not every
simple cycle. Graphene's 576 atoms have 1728 ten-cycles that are sums of
hexagons; those are not primitive rings.

## Formulation

Build the periodic bond graph on the P1 cell. For each start atom, depth-first
search simple paths up to `max_size`; a path that closes to the start is a
candidate cycle. Keep only the start-minimum representative and drop any cycle
with a chord: a chord splits the cycle into two smaller cycles, so the cycle is
not primitive (Franzblau).

## Complexity

Exponential in the worst case; bounded by `max_size` and the graph degree.
Typical frameworks with `max_size = 8` finish in milliseconds. The search is
deterministic because the bond graph and the start order are deterministic.

## Caveats

- Ring membership and aromaticity are counted on atom indices; a ring crossing
  a cell boundary is still one ring.
- `is_aromatic` is a heuristic (five/six-membered, C/N/O) until bond orders
  and planarity are available.
