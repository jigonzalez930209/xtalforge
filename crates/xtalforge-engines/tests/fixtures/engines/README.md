# Engine output fixtures

Every file in this directory is a **hand-written representative output**
authored for the xtalforge test suite on 2026-11-13. None of them is a copy of
a third-party output file, a published data set or a program distribution.
They contain only small, physically plausible values (a two-atom MgO or Si2
cell) chosen to exercise each parser field; they are **not** results of a real
calculation.

Provenance and licensing:

- Author: the xtalforge maintainers.
- License: the repository license, `MIT OR Apache-2.0`.
- Format documentation used while writing them: the VASP Wiki pages for
  `vasprun.xml`, `OUTCAR`, `OSZICAR`, `EIGENVAL`, `DOSCAR` and `PROCAR`; the
  Quantum ESPRESSO `pw.x`, `bands.x`, `projwfc.x` and `data-file-schema.xml`
  documentation; the CP2K output manual; the LAMMPS `thermo`/dump
  documentation; the GULP output manual; the Henkelman group `ACF.dat` layout
  and the Chargemol DDEC6 net-atomic-charge XYZ layout.

Because the files are representative rather than recorded, tests assert the
values written here, not values from an external run. A parser must stay
tolerant of the header variations described in its module documentation, and
must return a typed `Error::Parse` when a file like these is truncated.

## Layout

| Path | Engine | Parser |
|---|---|---|
| `vasp/vasprun.xml` | VASP 6 | `parse::parse_vasp_vasprun` |
| `vasp/OUTCAR` | VASP 6 | `parse::parse_vasp_outcar` |
| `vasp/OSZICAR` | VASP 6 | `parse::parse_vasp_oszicar` |
| `vasp/EIGENVAL` | VASP 6 | `parse::parse_vasp_eigenval` |
| `vasp/DOSCAR` | VASP 6 | `parse::parse_vasp_doscar` |
| `vasp/PROCAR` | VASP 6 | `parse::parse_vasp_procar` |
| `qe/pw.out` | Quantum ESPRESSO 7 | `parse::parse_qe` |
| `qe/data-file-schema.xml` | Quantum ESPRESSO 7 | `parse::parse_qe_xml` |
| `qe/bands.dat` | Quantum ESPRESSO 7 | `parse::parse_qe_bands` |
| `qe/projwfc.out` | Quantum ESPRESSO 7 | `parse::parse_qe_projwfc` |
| `cp2k/geo_opt.out` | CP2K 2024.1 | `parse::parse_cp2k` |
| `lammps/minimize.log` | LAMMPS 2023 | `parse::parse_lammps` (log) |
| `lammps/final.dump` | LAMMPS 2023 | `parse::parse_lammps_dump` |
| `gulp/optimise.out` | GULP 6 | `parse::parse_gulp` |
| `bader/ACF.dat` | Bader | `result::bader_charges` |
| `bader/ddec6_charges.xyz` | Chargemol DDEC6 | `result::ddec_charges` |

Malformed and truncated cases are inline strings in `tests/parsers.rs`; they
are deliberately not stored here so the fixture directory stays a catalogue of
valid representative files.
