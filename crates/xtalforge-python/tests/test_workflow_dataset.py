"""Workflow and dataset surfaces (Phase 16)."""

import json

import xtalforge

CIF = """data_Si
_cell_length_a 5.43
_cell_length_b 5.43
_cell_length_c 5.43
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0 0 0
"""


def test_workflow_runs_and_validates(tmp_path) -> None:
    source = tmp_path / "Si.cif"
    source.write_text(CIF, encoding="utf-8")
    document = json.dumps(
        {
            "workflow_version": 1,
            "name": "inspect",
            "steps": [
                {"kind": "read", "path": "{{input}}"},
                {"kind": "analyse", "name": "info", "verb": "info"},
            ],
        }
    )
    workflow = xtalforge.Workflow(document)
    summary = workflow.validate()
    assert summary["name"] == "inspect"
    report = workflow.run([str(source)])
    assert report["succeeded"] == 1
    assert report["failed"] == 0
    assert report["items"][0]["status"] == "ok"
    assert "schema" in xtalforge.Workflow.schema()


def test_dataset_jsonl_and_pandas(tmp_path) -> None:
    dataset = xtalforge.Dataset(
        [{"id": "a", "formula": "Si"}, {"id": "b", "formula": "NaCl"}]
    )
    text = dataset.to_jsonl()
    back = xtalforge.Dataset.from_jsonl(text)
    assert [record["id"] for record in back] == ["a", "b"]
    path = tmp_path / "data.jsonl"
    back.save(path)
    assert xtalforge.Dataset.load(path).records == back.records
    frame = back.to_pandas()
    assert list(frame["id"]) == ["a", "b"]
