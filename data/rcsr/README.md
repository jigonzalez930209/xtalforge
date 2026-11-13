# RCSR snapshot (2019-06-01)

Bundled with the maintainer's authorisation (2026-11-09). The RCSR
(Reticular Chemistry Structure Resource, <https://rcsr.ou.edu>) is freely
available for use by the general public; cite

> M. O'Keeffe, M. A. Peskov, S. J. Ramsden and O. M. Yaghi,
> *Accounts of Chemical Research* **41** (2008) 1782-1789.

## Files

| File | SHA-256 | Use |
|---|---|---|
| `RCSRnets-2019-06-01.cgd.gz` | `15503a7cb3cdedcdeeeab9036382ea6428446a4ce88f1e5d0053da7728f696f9` | Source nets (Systre CGD), provenance |
| `rcsr-2019-06-01.json.gz` | `731c91f5a045ee572309e9cc0d3e9a59584804f2ec73948ed40cd8d6e74f2fa1` | Fingerprints and canonical keys, index source |
| `rcsr-transitivity-2026-09-24.json.gz` | `2368c316596acd4f181f276a916e634f103fa538000c96a6bcc06f3f5a895d6f9` | Vertex/edge/face/tile transitivity, not used at runtime |

The files are the snapshots shipped by Crystal-Builder
(`xtal/analysis/data/`), copied unchanged.

## Regenerating the embedded index

The crate embeds `crates/xtalforge-analysis/data/rcsr_index.bin.gz`, a compact
gzip hash index derived from the JSON fingerprints. The generator prints the
SHA-256 of the decompressed payload, which is the manifest constant in
`crates/xtalforge-analysis/src/rcsr.rs` (Task 15.3.2):

```bash
scripts/golden/rcsr_index.sh
```

The generator rebuilds every keyed entry from its canonical key, compares the
result with the stored coordination sequences and point symbols, and fails on
any mismatch.
