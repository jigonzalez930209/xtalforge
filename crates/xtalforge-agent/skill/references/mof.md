# Reference: MOF builder

The MOF builder assembles a framework from a topology and building blocks.
The full builder lives in `xtalforge-mof` and the CLI `mof` command; the agent
loads and analyses the result.

## Curated subset

The published crate carries an illustrative subset of topologies and blocks,
not the full PORMAKE database. `xtalforge mof topologies` lists the names and
`xtalforge mof blocks` the blocks, each with its connection points. A name
outside the list is an error.

## Workflow

1. Read `xtalforge mof topologies` and choose by the coordination of the node
   and the linker.
2. Build with `xtalforge mof build --topology <name> --node <name> --linker
   <name> --out framework.cif`.
3. `load` the CIF in the agent, then `info` and `diagnostics`.
4. `symmetry` to see whether the assembly landed in the intended group;
   `rings` to check the smallest ring; `porosity` for the pore numbers.
5. `render` and look at the picture.

## Caveats

- The builder places blocks on the topology edges; it does not relax the
  geometry. Run `optimize` only on a small cell, or use the preparation and
  force-field surfaces deliberately.
- The topology assignment of a loaded framework is a net-identification
  question. Use the RCSR index and report ambiguity rather than guessing from
  the coordination alone.
- Interpenetration is a separate operation; do not infer it from a density.
