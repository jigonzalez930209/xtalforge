def crystal_graph_arrays(
    text: str,
    cutoff: float = 5.0,
    max_neighbors: int | None = None,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def crystal_graph_numpy(
    text: str,
    cutoff: float = 5.0,
    max_neighbors: int | None = None,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def soap_descriptor(
    text: str,
    cutoff: float = 5.0,
    n_max: int = 4,
    l_max: int = 4,
    sigma: float = 0.5,
    average: bool = False,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def acsf_descriptor(
    text: str, cutoff: float = 5.0, average: bool = False, format_hint: str | None = None
) -> dict[str, object]: ...
def mbtr_descriptor(
    text: str,
    cutoff: float = 5.0,
    grid: int = 64,
    sigma: float = 0.1,
    angular: bool = True,
    format_hint: str | None = None,
) -> dict[str, object]: ...
def mof_ml_features(
    text: str, max_path: int = 3, pores: bool = True, format_hint: str | None = None
) -> dict[str, object]: ...
def ml_split(
    text: str,
    strategy: str = "random",
    train: float = 0.6,
    validation: float = 0.2,
    test: float = 0.2,
    seed: int = 0,
) -> dict[str, object]: ...
def ml_import(
    text: str, format: str = "extxyz", name: str = "imported"
) -> dict[str, object]: ...
def ml_predict(model_json: str, features: list[float]) -> dict[str, object]: ...
def ml_predict_dataset(
    dataset_jsonl: str,
    model_json: str,
    target: str = "predicted",
    feature: str = "soap",
) -> dict[str, object]: ...
def ml_select(
    features: list[list[float]],
    uncertainties: list[float],
    batch: int = 8,
    strategy: str = "uncertainty",
    diversity_radius: float = 0.0,
) -> list[int]: ...
def ml_active_workflow(
    name: str = "active-learning",
    dataset: str = "dataset.jsonl",
    model: str = "model.json",
    batch: int = 8,
    rounds: int = 4,
    checkpoint: str | None = None,
) -> dict[str, object]: ...
def ml_write_parquet(dataset_jsonl: str, path: str) -> int: ...
def ml_read_parquet(path: str) -> str: ...
def onnx_predict(
    model_path: str, features: list[float], n_features: int
) -> list[float]: ...
def onnx_energy(
    text: str,
    model_path: str,
    n_features: int,
    displacement: float = 0.001,
    format_hint: str | None = None,
) -> float: ...
