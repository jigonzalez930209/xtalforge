# Prompt: build a MOF

Goal: assemble a framework from a topology and a set of building blocks, then
check it.

1. `capabilities`: confirm `mof` verbs and the RCSR bundle.
2. `load` the metal cluster or node CIF if you have one, and the linker CIF.
   Name them `node` and `linker`.
3. If the blocks are connection points rather than full molecules, use the
   MOF builder with the topology name and the block names it accepts. The
   curated subset and the topology list are in `references/mof.md`.
4. `info` on the result. Check `diagnostics`; `dangling_atoms` often means a
   connection point was not consumed, and `overlapping_atoms` means two blocks
   were placed on top of each other.
5. `symmetry` to confirm the expected space group and `rings` to see the
   smallest rings.
6. `porosity` with `spacing: 0.5` for a first look, then `0.3` when the answer
   is needed for a report. Record the probe radius (default 1.2 A).
7. `render` a PNG and inspect it visually. A folded linker or a cluster in the
   pore is easier to see than to read from the site table.
8. `save` the CIF. `changes` lists every edit the assistant made; `undo`
   reverts the last burst.

Do not claim a topology assignment from the builder alone: run the net
identification and report the RCSR name and its ambiguity, if any.
