# Reference: polymers

`build polymer` packs chains from a SMILES repeat unit into a periodic cell.

Options: `smiles`, `monomers`, `density`. The builder embeds one monomer,
polymerizes it along the connection points and packs chains until the target
density is met. The result is a P1 structure with explicit bonds.

Use it when the user needs an amorphous or semi-crystalline polymer model, not
a crystal. For a known crystal, `load` the deposited CIF instead.

Checks after building:

- `info`: formula, density and cell. Compare the density with the target.
- `diagnostics`: `overlapping_atoms` means the packing failed;
  `dangling_atoms` means a chain end was not terminated.
- `rings`: the repeat unit's rings should appear.
- `uff_energy` and a short `optimize` only for a small cell.

The build is deterministic for the same SMILES, monomer count and density.
Report those three values with any result.
