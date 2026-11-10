# Reference: porosity

`porosity` builds a periodic distance grid and reads four numbers from it.

- `spacing`: grid step in Angstrom (default 0.4). Cost scales with
  `1/spacing^3`; 0.3-0.4 is a report value, 0.6-1.0 a quick look.
- `probe_radius`: sphere rolled through the pore (default 1.2 A, the usual
  nitrogen probe).

Answer fields:

| Field | Meaning |
|---|---|
| `cell_volume` | volume of the cell in Angstrom^3 |
| `accessible_volume` | volume the probe centre can reach |
| `accessible_fraction` | accessible volume over cell volume |
| `surface_area` | accessible surface area in Angstrom^2 |
| `largest_included_diameter` | largest sphere that fits, in Angstrom |

The grid uses the element van der Waals radii; it is a geometric measurement,
not a force-field or a zeolite-specific probe. Report the spacing and the probe
radius with the numbers, and run two spacings when the value matters: if they
differ by more than a few percent, the grid is not converged.

The limits check refuses a grid above the standard budget before allocation;
increase `spacing` or change the `--limits` profile deliberately.
