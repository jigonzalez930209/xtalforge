# Prompt: porosity report

Goal: a defensible porosity number for one structure.

1. `load` the CIF and run `info`. If `diagnostics` reports
   `overlapping_atoms`, `empty_structure` or a charge problem, stop and report
   the model problem before measuring.
2. Remove guests or solvent only if the user asked; record the edit with
   `changes` so the report can say what was removed.
3. Run `porosity` twice:
   - `spacing: 0.6` for a quick value;
   - `spacing: 0.3` for the reported value.
   If the two differ by more than a few percent, report both and say the grid
   is not converged. The grid spacing is part of the result.
4. Report the probe radius (default 1.2 A), the accessible volume, the
   accessible fraction, the surface area and the largest included sphere
   diameter. Do not convert the sphere diameter into a pore size without
   saying it is the largest included sphere, not a channel diameter.
5. `render` a picture along a channel axis if the user wants one.
6. `save` nothing unless asked; a measurement does not modify the document.

The reference behaviour is Crystal-Builder's `zeo` bridge and the distance
grid in `docs/algorithms/agent.md`; the numbers are geometric, not
force-field derived.
