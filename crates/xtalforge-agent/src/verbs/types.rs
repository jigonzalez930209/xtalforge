//! Verb and parameter types shared by the catalogue, the MCP server and the
//! bindings.

use serde::Serialize;

/// JSON type of one verb parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ParameterKind {
    /// A JSON string.
    String,
    /// A JSON integer.
    Integer,
    /// A JSON number.
    Number,
    /// A JSON boolean.
    Boolean,
    /// A JSON array.
    Array,
    /// A JSON object.
    Object,
}

/// One parameter of a verb.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct ParameterInfo {
    /// Argument name.
    pub name: &'static str,
    /// JSON type.
    pub kind: ParameterKind,
    /// Whether the caller must provide it.
    pub required: bool,
    /// One-line description.
    pub summary: &'static str,
    /// Default rendered as a JSON literal, when there is one.
    pub default: Option<&'static str>,
}

/// One verb an assistant can call.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct VerbInfo {
    /// Verb name as used by `run`.
    pub name: &'static str,
    /// One-line description.
    pub summary: &'static str,
    /// Parameters in call order.
    pub parameters: &'static [ParameterInfo],
    /// Result type returned by the verb.
    pub result: &'static str,
    /// A complete JSON argument example.
    pub example: &'static str,
    /// Whether the verb mutates the active document.
    pub mutation: bool,
}

pub(crate) const fn param(
    name: &'static str,
    kind: ParameterKind,
    required: bool,
    summary: &'static str,
    default: Option<&'static str>,
) -> ParameterInfo {
    ParameterInfo {
        name,
        kind,
        required,
        summary,
        default,
    }
}
