"""Phase 21 runtime re-exports and NumPy helpers (representations and ML)."""

from ._xtalforge import (
    acsf_descriptor,
    crystal_graph_arrays,
    mbtr_descriptor,
    ml_active_workflow,
    ml_import,
    ml_predict,
    ml_predict_dataset,
    ml_select,
    ml_split,
    mof_ml_features,
    soap_descriptor,
)


def crystal_graph_numpy(
    text: str,
    cutoff: float = 5.0,
    max_neighbors: int | None = None,
    format_hint: str | None = None,
) -> dict[str, object]:
    """Returns the graph arrays as NumPy arrays."""
    import numpy

    arrays = crystal_graph_arrays(text, cutoff, max_neighbors, format_hint)
    return {key: numpy.asarray(value) for key, value in arrays.items()}


__all__ = [
    "acsf_descriptor",
    "crystal_graph_arrays",
    "crystal_graph_numpy",
    "mbtr_descriptor",
    "ml_active_workflow",
    "ml_import",
    "ml_predict",
    "ml_predict_dataset",
    "ml_select",
    "ml_split",
    "mof_ml_features",
    "soap_descriptor",
]
