# Hand-made edge cases

These CIFs were written for the test corpus (CC0, no third-party data) to cover
the cases the sample structures do not:

| File | Case |
|---|---|
| `triclinic.cif` | triclinic cell, P1 |
| `hexagonal.cif` | hexagonal metric, P1 |
| `partial_occupancy.cif` | two sites at occupancy 0.5 |
| `disorder.cif` | two alternative carbon positions at occupancy 0.5 |
| `special_position.cif` | sites at the origin and body centre of P-1 |
| `malformed.cif` | truncated cell item and atom loop; readers must return an error |

They contain one or two atoms on purpose: they test parsing, symmetry and
occupancy handling, not chemistry.
