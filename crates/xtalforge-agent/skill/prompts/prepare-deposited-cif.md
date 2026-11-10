# Prompt: prepare a deposited CIF

Goal: turn a raw deposited CIF into a model that is safe to analyse.

1. `load` the CIF. Keep the deposited file as the source; never overwrite it.
2. `info`: read the formula, the space group symbol, the occupancy list and
   the diagnostics.
3. `diagnostics` decides the next step:
   - `partial_occupancy`: list the affected sites and ask which component is
     wanted. Do not silently delete atoms.
   - `missing_hydrogens`: run `hydrogens` only for an X-ray model that should
     have them; record the X-H distance used (default 1.09 A).
   - `no_symmetry`: the file may have lost its space group. Run `symmetry`
     with a larger `symprec` (for example 1e-3) and compare.
   - `overlapping_atoms`: the file is probably disordered or in the wrong
     setting; report it rather than editing.
4. `wrap` after any edit that moves sites.
5. `symmetry` to record the space group and the standardized cell. Report the
   Hall number when two settings share a symbol.
6. `save` the prepared CIF under a new name and report every step from
   `changes`.

The full preparation pipeline (tasks, order and reports) is in
`references/prepare.md`. The agent verbs above are the interactive subset.
