//! ML verb definitions.

use crate::verbs::ParameterKind::{Array, Integer, Number, String as Text};
use crate::verbs::VerbInfo;
use crate::verbs::param;

/// ML verbs appended to the catalogue.
pub(crate) const ML_VERBS: [VerbInfo; 3] = [
    VerbInfo {
        name: "ml_predict",
        summary: "Predict a scalar property for the active document with a linear model",
        parameters: &[
            param(
                "document",
                Text,
                false,
                "Document name; defaults to the active document",
                Some("\"active\""),
            ),
            param("model", Text, true, "Model JSON text or path", None),
            param(
                "features",
                Text,
                false,
                "`soap` or `acsf`",
                Some("\"soap\""),
            ),
        ],
        result: "ml_prediction",
        example: r#"{"model": "model.json"}"#,
        mutation: false,
    },
    VerbInfo {
        name: "ml_select",
        summary: "Select a batch from a feature table with uncertainty, farthest-point or diversity",
        parameters: &[
            param("features", Array, true, "Feature rows", None),
            param(
                "uncertainties",
                Array,
                false,
                "One uncertainty per row",
                None,
            ),
            param("batch", Integer, false, "Batch size", Some("8")),
            param(
                "strategy",
                Text,
                false,
                "`uncertainty`, `farthest_point` or `diversity`",
                Some("\"uncertainty\""),
            ),
        ],
        result: "ml_selection",
        example: r#"{"features": [[0.0], [1.0]], "uncertainties": [0.1, 0.9], "batch": 1}"#,
        mutation: false,
    },
    VerbInfo {
        name: "ml_append",
        summary: "Append a numeric metadata property to the active document",
        parameters: &[
            param(
                "document",
                Text,
                false,
                "Document name; defaults to the active document",
                Some("\"active\""),
            ),
            param("key", Text, true, "Metadata key", None),
            param("value", Number, true, "Numeric value", None),
        ],
        result: "ml_append",
        example: r#"{"key": "predicted", "value": 1.5}"#,
        mutation: true,
    },
];
