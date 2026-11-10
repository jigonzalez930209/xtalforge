# Reference: diagnostics

`diagnose` returns a report with `n_atoms`, a list of findings and counts per
severity. Findings are sorted by severity then by code.

| Code | Severity | Rule | Usual response |
|---|---|---|---|
| `empty_structure` | error | no atoms in P1 | the file or the read failed |
| `overlapping_atoms` | error | periodic distance < 0.4 A | disorder, wrong setting or wrong conversion |
| `partial_occupancy` | info | occupancy < 0.999 | choose the component before measuring |
| `charge_imbalance` | warning | abs(total charge) > 1e-6 e | check whether the model is ionic |
| `dangling_atoms` | warning | non-metal atom with no bond | bond rules or a missing connection point |
| `missing_hydrogens` | info | C/N/O below valence | `hydrogens` for an X-ray model |
| `unusual_coordination` | warning | more than 12 bonds | bond rules too generous |
| `no_symmetry` | info | P1 with more than three atoms | run `symmetry` with a larger tolerance |
| `symmetry_failed` | info | search did not finish | the cell is too distorted for the tolerance |

Atom indices refer to the P1 cell. Empty atom lists mean the finding is about
the structure as a whole. `suggestion` names the verb that would address it
when one exists.

Severity is a hint, not a verdict: an ionic framework legitimately has a
non-zero stored charge, and a P1 model is a valid choice. Report the finding
and the model context together.
