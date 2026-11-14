"""Phase 21 Python surfaces: ML representations, splits and inference."""

from pathlib import Path

import pytest

xtalforge = pytest.importorskip("xtalforge")

ROOT = Path(__file__).resolve().parents[3]
NACL = (ROOT / "tests" / "corpus" / "simple" / "NaCl.cif").read_text()


def test_graph_and_numpy_helper() -> None:
    graph = xtalforge.crystal_graph_arrays(NACL, cutoff=6.0, max_neighbors=4)
    assert graph["edge_index"]
    assert len(graph["atomic_numbers"]) == 8
    arrays = xtalforge.crystal_graph_numpy(NACL, cutoff=6.0, max_neighbors=4)
    assert arrays["edge_index"].shape[1] == 2


def test_soap_acsf_mbtr() -> None:
    soap = xtalforge.soap_descriptor(NACL, n_max=2, l_max=1, average=True)
    assert soap["n_features"] == 6
    acsf = xtalforge.acsf_descriptor(NACL, average=True)
    assert acsf["n_features"] > 0
    mbtr = xtalforge.mbtr_descriptor(NACL, grid=8)
    assert len(mbtr["features"]) > 0


def test_mof_features() -> None:
    features = xtalforge.mof_ml_features(NACL, pores=False)
    assert len(features["features"]) == 12


def test_split_and_import() -> None:
    records = []
    for index in range(4):
        document = xtalforge.read_structure(NACL, "cif")
        records.append({"id": f"r{index}", "structure": document, "properties": {"e": index}})
    dataset = xtalforge.Dataset(records).to_jsonl()
    split = xtalforge.ml_split(dataset, strategy="random", seed=1)
    assert len(split["labels"]) == 4

    extxyz = (
        '1\nLattice="5 0 0 0 5 0 0 0 5" Properties=species:S:1:pos:R:3 pbc="T T T"\n'
        "C 0 0 0\n"
    )
    report = xtalforge.ml_import(extxyz, format="extxyz", name="x")
    assert report["n_records"] == 1


def test_model_prediction_and_selection() -> None:
    model = '{"name":"m","weights":[1.0,2.0],"bias":0.5,"mean":[],"scale":[]}'
    value = xtalforge.ml_predict(model, [1.0, 2.0])
    assert value["value"] == pytest.approx(5.5)
    selected = xtalforge.ml_select(
        [[0.0], [1.0], [5.0]], [0.1, 0.9, 0.2], batch=2, strategy="uncertainty"
    )
    assert len(selected) == 2
    workflow = xtalforge.ml_active_workflow("al", "d.jsonl", "m.json", 4, 2)
    assert workflow["steps"][2]["verb"] == "ml_predict"


def test_widget_optional() -> None:
    try:
        widget = xtalforge.StructureWidget(NACL)
    except ImportError:
        pytest.skip("anywidget is not installed")
    else:
        assert widget.frames == 1


def test_parquet_round_trip(tmp_path) -> None:
    document = xtalforge.read_structure(NACL, "cif")
    dataset = xtalforge.Dataset(
        [{"id": "a", "structure": document, "properties": {"e": -1.5}}]
    ).to_jsonl()
    path = tmp_path / "data.parquet"
    count = xtalforge.ml_write_parquet(dataset, str(path))
    assert count == 1
    restored = xtalforge.ml_read_parquet(str(path))
    assert '"id":"a"' in restored.replace(" ", "")


def test_onnx_fixture() -> None:
    model = ROOT / "crates" / "xtalforge-ml" / "data" / "onnx_linear.onnx"
    values = xtalforge.onnx_predict(str(model), [0.0, 0.0, 0.0, 2.0, 0.0, 0.0], 6)
    assert abs(values[0] - 0.2) < 1e-5
    poscar = (
        "C2\n1.0\n20 0 0\n0 20 0\n0 0 20\nC\n2\nDirect\n0 0 0\n0.1 0 0\n"
    )
    energy = xtalforge.onnx_energy(poscar, str(model), 6)
    assert abs(energy - 0.2) < 1e-5
