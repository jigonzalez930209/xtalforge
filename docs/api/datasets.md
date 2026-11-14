# Datasets and descriptors

`xtalforge_analysis::Dataset` is a JSON Lines container:
`DatasetRecord { id, structure, properties }`, stable `DATASET_VERSION`, and
`to_json_lines`/`from_json_lines` plus `save`/`load`. `column` and
`property_names` expose tabular properties; iteration preserves file order.

`deduplicate(dataset, &DedupOptions)` uses the structure matcher with a
composition/supercell prefilter and returns a greedy `DedupReport` with the
kept and removed groups. The report's `apply` returns a duplicate-free dataset.

`composition_features`, `rdf_descriptor`, `coulomb_matrix` and
`soap_like_descriptor` return numeric descriptors for ML. The SOAP-like
descriptor is a Gaussian radial basis times Legendre angular moments and is
documented as an approximation, not exact SOAP.

Python exposes the same container as `xtalforge.Dataset` with JSONL I/O and
optional pandas/polars views.

Parity keys: `dataset`, `dataset_dedup`, `ml_descriptors`.

## Phase 21: splits, imports and Parquet

`split_dataset(dataset, SplitOptions { train, validation, test, seed, strategy })`
returns one `SplitLabel` per record and the group key. Strategies are
`random`, `composition` (whole composition groups), `cluster` (whole matcher
clusters) and `prototype` (space group plus Wyckoff signature). Fractions are
largest-remainder rounded and groups are shuffled with the explicit seed.

`import_dataset(text, ImportOptions { format, name, id_prefix })` reads
extxyz (including `Lattice`, `Properties` and `pbc`), JSON Lines or a JSON
array, and OPTIMADE JSON. `import_cif_directory` reads a directory of CIFs for
CoRE MOF or QMOF. The report carries the dataset, the imported count and the
skipped count.

With the `ml-parquet` feature, `write_parquet` and `read_parquet` store and
load the same dataset as a columnar file. The schema is `id`, `structure`
(versioned structure document JSON) and `properties` (JSON object); the
conversion is lossless and JSON Lines stays the default.

Parity keys: `crystal_graph`, `soap_descriptor`, `acsf_descriptor`,
`mbtr_descriptor`, `mof_ml_features`, `ml_split`, `ml_import`, `ml_predict`,
`ml_select`, `ml_active_workflow`, `dataset_parquet`, `onnx_inference`.
