# Functional-group matching

The matcher runs on the perceived periodic bond graph. For every atom it reads
the element and its bonded neighbours, then applies the catalogue:

| Group | Pattern |
|---|---|
| hydroxyl | O with an H neighbour |
| carboxylate | C with at least two O neighbours |
| amine | N with one to three neighbours |
| fluoride | F with one neighbour |

A match is `GroupMatch { kind, atoms }`, sorted and deduplicated.
`group_counts` aggregates by kind. The catalogue is a fixed enum now; it
becomes data-driven with the remaining Crystal-Builder patterns and their
framework tests.

The matcher is deterministic because the bond graph and atom order are
deterministic. It runs on the P1 cell, so a group crossing a cell boundary is
found once.
