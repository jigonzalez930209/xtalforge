# Train and run a MOF property model

This tutorial builds a property model from MOF features, validates it with a
leakage-safe split and runs it on a new framework. It uses the native
`xtalforge-ml` surface; the Python ML bridge from Phase 8 is still the route
for PyTorch message-passing models.

## 1. Describe a framework

```bash
xtalforge ml mof-features framework.cif --out features.json
xtalforge ml soap framework.cif --average --out soap.json
xtalforge ml graph framework.cif --k 12 --out graph.json
```

`mof-features` concatenates the RACs and the pore vector. `soap` is the exact
power spectrum; `graph` is the periodic graph with image offsets.

## 2. Build the dataset

```bash
xtalforge ml import core-mof.extxyz --format extxyz --name core --out core.jsonl
xtalforge ml split core.jsonl --strategy composition --seed 7 --out split.json
```

Composition, cluster and prototype splits keep whole groups on one side, so a
model cannot memorise the same composition in training and test.

## 3. Fit a model

The native `PropertyModel` is linear. Fit the weights with any tool and save
the JSON:

```json
{"name": "co2-uptake", "weights": [0.1, -0.2, 0.05], "bias": 1.5,
 "mean": [0.0, 0.0, 0.0], "scale": [1.0, 1.0, 1.0]}
```

```bash
xtalforge ml predict core.jsonl --model model.json --features soap --target uptake --out predicted.jsonl
```

`predict_dataset` and `append_predictions` do the same in Rust; Python exposes
`ml_predict` and `ml_predict_dataset`.

## 4. Select the next experiments

```bash
xtalforge ml select features.json --batch 8 --strategy uncertainty --out selected.json
```

The feature table carries `features` and optional `uncertainties`. Strategies
are `uncertainty`, `farthest_point` and `diversity`; the output is a sorted
index list.

## 5. Run the active loop

```bash
xtalforge ml active-workflow \
  --name al --dataset core.jsonl --model model.json \
  --batch 8 --rounds 4 --checkpoint al.checkpoint.json \
  --out al.workflow.json
```

The document has read, search, predict, select, append and write steps. The
`search` step calls the Phase 20 checkpointed search verb; the ML steps call
the `ml_*` verbs. Running it with `xtalforge run al.workflow.json` keeps the
provenance.

## 6. From Python

```python
import xtalforge

graph = xtalforge.crystal_graph_numpy(open("framework.cif").read(), cutoff=5.0, max_neighbors=12)
soap = xtalforge.soap_descriptor(open("framework.cif").read(), average=True)
model = '{"name":"m","weights":[0.1,0.2],"bias":0.0}'
value = xtalforge.ml_predict(model, soap["rows"][0][:2])
print(graph["edge_index"].shape, value)
```

## Limits

- The native model is linear; fit it with scikit-learn or another tool and
  export the weights.
- ONNX inference is included in the default wheels, CLI and WASM package and
  needs a single-input graph; message-passing models use the Python bridge.
- Parquet is included in the default wheels, CLI and WASM package. WASM uses
  byte buffers because it has no filesystem.
- SOAP, ACSF and MBTR match the recorded DScribe fixtures within `5e-6`
  (`5e-4` for SOAP); the fast analysis approximation stays available for large
  screens.
- WASM has no filesystem: use JSONL text and the in-memory byte APIs.
