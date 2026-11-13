//! The verb catalogue (Task 8.1.1).
//!
//! One table lists every operation an assistant can call. The CLI help, the
//! Python and WebAssembly bindings and the MCP server read this table, so a
//! capability is declared once. `capabilities` reports what the current
//! installation can actually run; this table declares what the layer offers.

mod types;

pub use types::{ParameterInfo, ParameterKind, VerbInfo};

use types::param;

use ParameterKind::{Array, Boolean, Integer, Object, String as Text};

const DOCUMENT: ParameterInfo = param(
    "document",
    ParameterKind::String,
    false,
    "Document name; defaults to the active document",
    Some("\"active\""),
);

const PARAM_WORKFLOW: ParameterInfo = param("workflow", Object, false, "Document object", None);
const PARAM_TEXT: ParameterInfo = param("text", Text, false, "Document text", None);
const PARAM_PATH: ParameterInfo = param("path", Text, false, "Workflow file", None);
const PARAM_INPUTS: ParameterInfo = param("inputs", Array, true, "Inputs to process", None);
const PARAM_JOBS: ParameterInfo = param("jobs", Integer, false, "Items at once", Some("1"));
const PARAM_CHECKPOINT: ParameterInfo = param("checkpoint", Text, false, "Checkpoint file", None);
const PARAM_RESUME: ParameterInfo = param("resume", Boolean, false, "Resume the checkpoint", None);

/// Source arguments shared by `workflow_validate` and `workflow_run`.
const WORKFLOW_SOURCE: &[ParameterInfo] = &[PARAM_WORKFLOW, PARAM_TEXT, PARAM_PATH];

/// Arguments of `workflow_run`.
const WORKFLOW_RUN_PARAMETERS: &[ParameterInfo] = &[
    PARAM_WORKFLOW,
    PARAM_TEXT,
    PARAM_PATH,
    PARAM_INPUTS,
    PARAM_JOBS,
    PARAM_CHECKPOINT,
    PARAM_RESUME,
];

/// Every verb, in the order the skill document presents them.
pub const VERBS: &[VerbInfo] = &[
    VerbInfo {
        name: "capabilities",
        summary: "Report what this installation can do: features, engines, data bundles and verbs",
        parameters: &[],
        result: "capabilities",
        example: "{}",
        mutation: false,
    },
    VerbInfo {
        name: "list_documents",
        summary: "List the documents open in the session",
        parameters: &[],
        result: "documents",
        example: "{}",
        mutation: false,
    },
    VerbInfo {
        name: "load",
        summary: "Read a structure from a file or from text and open it as a document",
        parameters: &[
            param(
                "path",
                ParameterKind::String,
                false,
                "File to read; omit when `text` is given",
                None,
            ),
            param(
                "text",
                ParameterKind::String,
                false,
                "Structure content",
                None,
            ),
            param(
                "format",
                ParameterKind::String,
                false,
                "Format name or extension; sniffed when omitted",
                None,
            ),
            param(
                "name",
                ParameterKind::String,
                false,
                "Document name; defaults to the file stem",
                None,
            ),
        ],
        result: "loaded",
        example: r#"{"path": "MOF-5.cif", "name": "mof5"}"#,
        mutation: false,
    },
    VerbInfo {
        name: "save",
        summary: "Write a document to a file in a chosen format",
        parameters: &[
            DOCUMENT,
            param("path", ParameterKind::String, true, "Output file", None),
            param(
                "format",
                ParameterKind::String,
                false,
                "Format name or extension",
                None,
            ),
        ],
        result: "saved",
        example: r#"{"path": "out.cif"}"#,
        mutation: false,
    },
    VerbInfo {
        name: "info",
        summary: "Name, formula, cell, density, charge and symmetry of a document",
        parameters: &[DOCUMENT],
        result: "info",
        example: "{}",
        mutation: false,
    },
    VerbInfo {
        name: "diagnostics",
        summary: "Problems in a structure with severity, atoms and a suggested verb",
        parameters: &[DOCUMENT],
        result: "diagnostics",
        example: "{}",
        mutation: false,
    },
    VerbInfo {
        name: "supercell",
        summary: "Replace a document with a supercell",
        parameters: &[
            DOCUMENT,
            param(
                "repeats",
                ParameterKind::Array,
                true,
                "Three positive integers",
                None,
            ),
        ],
        result: "structure",
        example: r#"{"repeats": [2, 2, 2]}"#,
        mutation: true,
    },
    VerbInfo {
        name: "wrap",
        summary: "Wrap every site into the unit cell",
        parameters: &[DOCUMENT],
        result: "structure",
        example: "{}",
        mutation: true,
    },
    VerbInfo {
        name: "hydrogens",
        summary: "Add the hydrogens the perceived coordination is missing",
        parameters: &[
            DOCUMENT,
            param(
                "distance",
                ParameterKind::Number,
                false,
                "X-H distance in Angstrom",
                Some("1.09"),
            ),
        ],
        result: "structure",
        example: r#"{"distance": 1.09}"#,
        mutation: true,
    },
    VerbInfo {
        name: "symmetry",
        summary: "Detect the space group, Wyckoff letters and standardized cell",
        parameters: &[
            DOCUMENT,
            param(
                "symprec",
                ParameterKind::Number,
                false,
                "Position tolerance in Angstrom",
                Some("1e-5"),
            ),
        ],
        result: "symmetry",
        example: r#"{"symprec": 1e-3}"#,
        mutation: false,
    },
    VerbInfo {
        name: "bonds",
        summary: "Perceive bonds and report the graph, including periodic images",
        parameters: &[DOCUMENT],
        result: "bonds",
        example: "{}",
        mutation: false,
    },
    VerbInfo {
        name: "rings",
        summary: "Count primitive rings up to a size limit",
        parameters: &[
            DOCUMENT,
            param(
                "max_size",
                ParameterKind::Integer,
                false,
                "Largest ring",
                Some("12"),
            ),
        ],
        result: "rings",
        example: r#"{"max_size": 8}"#,
        mutation: false,
    },
    VerbInfo {
        name: "pxrd",
        summary: "Calculate a powder pattern",
        parameters: &[
            DOCUMENT,
            param(
                "radiation",
                ParameterKind::String,
                false,
                "cu-ka, mo-ka, ...",
                Some("\"cu-ka\""),
            ),
            param(
                "two_theta_min",
                ParameterKind::Number,
                false,
                "Lowest 2 theta",
                Some("5.0"),
            ),
            param(
                "two_theta_max",
                ParameterKind::Number,
                false,
                "Highest 2 theta",
                Some("80.0"),
            ),
            param(
                "fwhm",
                ParameterKind::Number,
                false,
                "Gaussian FWHM in degrees",
                Some("0.1"),
            ),
        ],
        result: "pxrd",
        example: r#"{"radiation": "cu-ka", "two_theta_max": 60.0}"#,
        mutation: false,
    },
    VerbInfo {
        name: "porosity",
        summary: "Surface area, pore volume and largest included sphere from a distance grid",
        parameters: &[
            DOCUMENT,
            param(
                "spacing",
                ParameterKind::Number,
                false,
                "Grid spacing in Angstrom",
                Some("0.4"),
            ),
            param(
                "probe_radius",
                ParameterKind::Number,
                false,
                "Probe radius",
                Some("1.2"),
            ),
        ],
        result: "porosity",
        example: r#"{"spacing": 0.5}"#,
        mutation: false,
    },
    VerbInfo {
        name: "uff_energy",
        summary: "UFF energy, forces and stress of a document",
        parameters: &[DOCUMENT],
        result: "energy",
        example: "{}",
        mutation: false,
    },
    VerbInfo {
        name: "optimize",
        summary: "Relax the geometry with UFF and replace the document",
        parameters: &[
            DOCUMENT,
            param(
                "max_steps",
                ParameterKind::Integer,
                false,
                "Step limit",
                Some("200"),
            ),
            param(
                "algorithm",
                ParameterKind::String,
                false,
                "fire or lbfgs",
                Some("\"fire\""),
            ),
        ],
        result: "optimization",
        example: r#"{"max_steps": 50}"#,
        mutation: true,
    },
    VerbInfo {
        name: "render",
        summary: "Draw the scene to SVG or PNG for an assistant to look at",
        parameters: &[
            DOCUMENT,
            param(
                "format",
                ParameterKind::String,
                false,
                "svg or png",
                Some("\"svg\""),
            ),
            param(
                "width",
                ParameterKind::Integer,
                false,
                "Image width in pixels",
                Some("800"),
            ),
            param(
                "height",
                ParameterKind::Integer,
                false,
                "Image height in pixels",
                Some("600"),
            ),
            param(
                "path",
                ParameterKind::String,
                false,
                "Write the image here",
                None,
            ),
        ],
        result: "picture",
        example: r#"{"format": "png", "path": "mof5.png"}"#,
        mutation: false,
    },
    VerbInfo {
        name: "undo",
        summary: "Undo the most recent edit",
        parameters: &[DOCUMENT],
        result: "history",
        example: "{}",
        mutation: true,
    },
    VerbInfo {
        name: "redo",
        summary: "Redo the most recently undone edit",
        parameters: &[DOCUMENT],
        result: "history",
        example: "{}",
        mutation: true,
    },
    VerbInfo {
        name: "history",
        summary: "The document's command log with origins",
        parameters: &[DOCUMENT],
        result: "history",
        example: "{}",
        mutation: false,
    },
    VerbInfo {
        name: "changes",
        summary: "Bursts of edits made through this agent session",
        parameters: &[DOCUMENT],
        result: "changes",
        example: "{}",
        mutation: false,
    },
    VerbInfo {
        name: "close",
        summary: "Close a document and drop it from the session",
        parameters: &[DOCUMENT],
        result: "documents",
        example: "{}",
        mutation: false,
    },
    VerbInfo {
        name: "workflow_schema",
        summary: "Return the JSON Schema of a workflow document",
        parameters: &[],
        result: "workflow_schema",
        example: "{}",
        mutation: false,
    },
    VerbInfo {
        name: "workflow_validate",
        summary: "Parse and validate a workflow document without running it",
        parameters: WORKFLOW_SOURCE,
        result: "workflow_validation",
        example: r#"{"path": "wf.json"}"#,
        mutation: false,
    },
    VerbInfo {
        name: "workflow_run",
        summary: "Run a workflow over input files, directories or glob patterns",
        parameters: WORKFLOW_RUN_PARAMETERS,
        result: "workflow_run",
        example: r#"{"path": "wf.json", "inputs": ["structures/*.cif"], "jobs": 4}"#,
        mutation: false,
    },
];

/// The catalogue an assistant reads before its first call.
#[must_use]
pub fn catalogue() -> &'static [VerbInfo] {
    VERBS
}

/// Looks up a verb by name.
#[must_use]
pub fn lookup(name: &str) -> Option<&'static VerbInfo> {
    VERBS.iter().find(|verb| verb.name == name)
}

/// Every verb name, in catalogue order.
#[must_use]
pub fn names() -> Vec<&'static str> {
    VERBS.iter().map(|verb| verb.name).collect()
}
