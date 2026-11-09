# Nets and RCSR identification

## Quotient graph

A periodic net is stored as a quotient graph: `nodes` are fractional
vertices, `edges` carry the lattice shift of the second endpoint, and the
labels are only decoration. `components` finds disconnected pieces and
`interpenetrate` duplicates a graph at a shift. The CGD reader accepts both
dialects of the RCSR file: an edge as a pair of points, or as a node label
followed by one neighbour position. The file gives the asymmetric unit; a
CGD entry read here must already be the quotient cell (or a `P1` entry), since
space-group expansion is a separate step not performed by `quotient_graph`.

## Fingerprint

`NetFingerprint` is:

- the periodicity: the integer rank of the lattice the net's own cycles
  generate (spanning-forest voltages, not the width of the cell),
- the distinct **coordination sequences** to ten shells, computed over
  `(vertex, lattice shift)` pairs so a six-neighbour vertex in a one-vertex
  cell counts six different atoms,
- the distinct **point symbols**: for every pair of edges at a vertex, the
  smallest ring through that angle, with its multiplicity (`4^12.6^3` for
  `pcu`). A ring longer than twelve is written `*`.

The token is `periodicity|sequences|symbols`, the form the RCSR index is
built on. `QuotientGraph::from_canonical_key` parses the RCSR canonical key
(`d:n:a,b,sx,sy,sz;...`) back into one minimal cell, which is how the index
is validated.

## RCSR catalogue

`crates/xtalforge-analysis/data/rcsr_index.bin` holds a 35 KB FNV-1a hash per
catalogue entry: the periodicity, the distinct coordination numbers and the
fingerprint hash, for all 2,931 nets of the 2019-06-01 snapshot. The full snapshots and their checksums live in
`data/rcsr/` at the repository root.

`RcsrIndex::from_embedded()` reads the index and `identify` computes the input
fingerprint:

- one name with the same hash is `named`,
- several names are `ambiguous` (they share both invariants),
- no name gives the `nearest` candidates with the same periodicity and
  coordination numbers.

The reference names 2,665 of the 2,679 three-periodic nets uniquely this way.
The index contains 35 shared fingerprint hashes (including 2D and polyhedral
entries); the corresponding names are all returned rather than one being
chosen silently.

## Validation

`examples/build_rcsr_index.rs` rebuilds every keyed entry from its canonical
key, computes the fingerprint and compares it with the stored coordination
sequences and point symbols. The committed index reports 2,909 keyed entries,
22 unkeyed entries and no mismatch; the unkeyed entries are taken from the
stored fingerprint because their canonical key is absent.

`tests/rcsr.rs` identifies `tbo`, `pcu`, `fcu`, `sod`, `msf`, `srs`, `dia`
and `sql` from their canonical keys and checks an unknown graph. `tests/cgd.rs`
reads the bundled snapshot (2,931 entries, one known incomplete block,
`moo-a`) and both edge dialects.

## Documented deviations

- Canonical keys are consumed from the RCSR JSON when the index is built, but
  an input graph's key is not recomputed; nets that share both invariants are
  reported as an ambiguous set.
- Space-group expansion of a CGD entry is not implemented; the reader accepts
  the expanded quotient cell and the `P1` form.
- The 22 unkeyed catalogue entries and the two-periodic/zero-periodic entries
  are matched by the stored fingerprint, not by a rebuilt graph.
