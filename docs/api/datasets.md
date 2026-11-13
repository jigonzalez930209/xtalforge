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
