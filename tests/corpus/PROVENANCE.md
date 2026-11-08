# Test corpus provenance

Copied from Crystal-Builder commit `78d5476fd9e851fda00cde57b7ab8da3195dcf54`,
`resources/samples/`. The full per-file provenance and citations are in
[`PROVENANCE-Crystal-Builder.md`](./PROVENANCE-Crystal-Builder.md); this page
records what was copied, what was skipped and why.

## Copied

| Path | Files | License | Notes |
|---|---|---|---|
| `simple/` | 8 CIF | CC0 (COD) | Textbook solids, fetched from the COD. |
| `cod/` | 21 CIF | CC0 (COD) | MOF frameworks listed in the source provenance. |
| `polymer/` | 3 CIF | COD CC0 / project-generated | The two written files carry their construction recipe in the header. |
| `prepared/` | 21 CIF | CC0-derived project models | Made from the CC0 COD files by Crystal-Builder's `prepare` pipeline. |
| `CFA1.cif` | 1 CIF | CC0 (COD 7021651) | Replaced the CCDC copy in Crystal-Builder. |
| `MIL53.cif` | 1 CIF | project-generated | PORMAKE build recorded in the source provenance. |
| `NiHITP.cif` | 1 CIF | project-generated | Pawley-refined model recorded in the source provenance. |
| `edge/` | hand-made | CC0 | Partial occupancy, disorder, special position, triclinic, hexagonal and malformed CIF cases. |

The COD data are CC0: they may be shipped, changed and passed on freely. The
`prepared/` files are models derived from CC0 data and carry no additional
license.

## Skipped

| File | Reason |
|---|---|
| `MFU4l.cif` | CCDC download with the CCDC access policy; redistribution is not clearly allowed. |
| `Ni2Cl2BTDD.cif` | CCDC export with no license statement. |
| `MOF-5.cif`, `ZIF-8.cif`, `zn_oac.cif` | Model exports of unrecorded origin. |
| `HKUST1.cif`, `UIO66.cif` | RASPA exports with no recorded redistribution license. |

## Reproducing the copy

```bash
cp -r Crystal-Builder/resources/samples/{simple,cod,polymer,prepared} tests/corpus/
cp Crystal-Builder/resources/samples/{CFA1,MIL53,NiHITP}.cif tests/corpus/
```

The copy is deliberate and reviewed; CI never regenerates it (Task 0.9.5).
