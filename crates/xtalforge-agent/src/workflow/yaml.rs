//! A line-oriented YAML reader for workflow documents (Task 17.2.1).
//!
//! JSON stays the canonical interchange format. This module is a deliberately
//! small, dependency-free YAML reader that turns the documented authoring
//! subset into a [`serde_json::Value`]; the normal serde deserialization and
//! the workflow validation rules then run unchanged, so a YAML document and
//! its JSON equivalent cannot validate differently.
//!
//! # Supported subset
//!
//! - Block mappings (`key: value`) and block sequences (`- item`), including
//!   compact `- name: value` mapping items and `- - item` sequences.
//! - Flow sequences and mappings on a single line (`[a, b]`, `{a: 1, b: [2]}`)
//!   with simple scalar keys.
//! - Plain, single-quoted and double-quoted scalars. Double quotes support
//!   `\"`, `\\`, `\/`, `\n`, `\t`, `\r`, `\b`, `\f` and `\uXXXX`; a doubled
//!   single quote (`''`) is an escaped quote.
//! - `#` comments (at the start of a line or after whitespace) and blank lines.
//! - `null`/`~`, `true`/`false` (and their capitalized forms), integers and
//!   finite decimal or exponent floats.
//! - One optional leading `---` document marker.
//!
//! Template strings such as `{{input}}` are kept as plain scalars even though
//! they start with a brace.
//!
//! # Rejected syntax
//!
//! Anchors (`&name`), aliases (`*name`), tags (`!tag`), directives (`%YAML`),
//! merge keys (`<<`), block scalars (`|` and `>`), multi-line plain scalars,
//! flow collections split over several lines and multiple documents are
//! rejected with a [`WorkflowError::Parse`] that names the line and column.
//! The reader does not try to be a general YAML implementation.

use serde_json::{Map, Value};

use super::WorkflowError;

mod flow;

use self::flow::{parse_inline, parse_key};

/// Parses the supported YAML subset into a JSON value.
///
/// # Errors
///
/// Returns [`WorkflowError::Parse`] with line and column context for syntax
/// outside the supported subset.
pub(crate) fn parse(text: &str) -> Result<Value, WorkflowError> {
    let lines = preprocess(text)?;
    if lines.is_empty() {
        return Ok(Value::Null);
    }
    let root_indent = lines.first().map_or(0, |line| line.indent);
    let mut parser = Parser { lines, position: 0 };
    let value = parser.parse_node(root_indent)?;
    if let Some(line) = parser.lines.get(parser.position) {
        return Err(syntax(
            line.number,
            line.indent + 1,
            format!("unexpected content `{}`", line.content),
        ));
    }
    Ok(value)
}

/// One significant source line with its indentation stripped.
#[derive(Debug, Clone)]
struct Line {
    /// One-based source line.
    number: usize,
    /// Width of the leading spaces.
    indent: usize,
    /// Content without indentation, comments or trailing whitespace.
    content: String,
}

fn syntax(line: usize, column: usize, message: impl Into<String>) -> WorkflowError {
    WorkflowError::Parse {
        reason: format!("YAML line {line}, column {column}: {}", message.into()),
    }
}

fn preprocess(text: &str) -> Result<Vec<Line>, WorkflowError> {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let mut lines = Vec::new();
    let mut has_content = false;
    let mut has_marker = false;
    for (index, raw) in text.split('\n').enumerate() {
        let number = index + 1;
        let raw = raw.strip_suffix('\r').unwrap_or(raw);
        let mut indent = 0;
        for character in raw.chars() {
            match character {
                ' ' => indent += 1,
                '\t' => {
                    return Err(syntax(
                        number,
                        indent + 1,
                        "tabs are not allowed in indentation",
                    ));
                }
                _ => break,
            }
        }
        let rest = raw.get(indent..).unwrap_or("");
        let content = strip_comment(rest).trim_end();
        if content.is_empty() {
            continue;
        }
        let trimmed = content.trim_start();
        if trimmed == "---" && !has_content && !has_marker {
            has_marker = true;
            continue;
        }
        if trimmed == "---"
            || trimmed == "..."
            || trimmed.starts_with("--- ")
            || trimmed.starts_with("... ")
        {
            return Err(syntax(
                number,
                indent + 1,
                "multi-document YAML is not supported",
            ));
        }
        if trimmed.starts_with('%') {
            return Err(syntax(
                number,
                indent + 1,
                "YAML directives are not supported",
            ));
        }
        lines.push(Line {
            number,
            indent,
            content: trimmed.to_owned(),
        });
        has_content = true;
    }
    Ok(lines)
}

fn strip_comment(text: &str) -> &str {
    let mut quote: Option<char> = None;
    let mut escaped = false;
    for (index, character) in text.char_indices() {
        if let Some(active) = quote {
            if active == '"' {
                if escaped {
                    escaped = false;
                } else if character == '\\' {
                    escaped = true;
                } else if character == '"' {
                    quote = None;
                }
            } else if character == '\'' {
                quote = None;
            }
            continue;
        }
        match character {
            '"' | '\'' => quote = Some(character),
            '#' if index == 0 || comment_boundary(text, index) => {
                return text.get(..index).unwrap_or(text);
            }
            _ => {}
        }
    }
    text
}

fn comment_boundary(text: &str, index: usize) -> bool {
    text.get(..index)
        .and_then(|prefix| prefix.chars().next_back())
        .is_some_and(|previous| previous == ' ' || previous == '\t')
}

struct Parser {
    lines: Vec<Line>,
    position: usize,
}

impl Parser {
    fn peek(&self) -> Option<&Line> {
        self.lines.get(self.position)
    }

    fn take(&mut self) -> Option<Line> {
        let line = self.lines.get(self.position)?.clone();
        self.position += 1;
        Some(line)
    }

    fn parse_node(&mut self, indent: usize) -> Result<Value, WorkflowError> {
        let Some(line) = self.peek() else {
            return Ok(Value::Null);
        };
        if line.indent < indent {
            return Ok(Value::Null);
        }
        if line.indent > indent {
            return Err(syntax(
                line.number,
                line.indent + 1,
                "unexpected indentation",
            ));
        }
        if sequence_rest(&line.content).is_some() {
            return self.parse_sequence(indent, None);
        }
        if split_key_value(&line.content).is_some() {
            return self.parse_mapping(indent, None);
        }
        let Some(line) = self.take() else {
            return Ok(Value::Null);
        };
        let value = parse_inline(&line.content, line.number, line.indent + 1)?;
        if let Some(next) = self.peek().filter(|next| next.indent > indent) {
            return Err(syntax(
                next.number,
                next.indent + 1,
                "multi-line plain scalars are not supported",
            ));
        }
        Ok(value)
    }

    fn parse_mapping(
        &mut self,
        indent: usize,
        first: Option<Line>,
    ) -> Result<Value, WorkflowError> {
        let mut map = Map::new();
        let mut first = first;
        loop {
            let line = if let Some(line) = first.take() {
                line
            } else {
                let Some(line) = self.peek() else {
                    break;
                };
                if line.indent < indent {
                    break;
                }
                if line.indent > indent {
                    return Err(syntax(
                        line.number,
                        line.indent + 1,
                        "unexpected indentation in a mapping",
                    ));
                }
                if sequence_rest(&line.content).is_some() {
                    break;
                }
                let Some(line) = self.take() else {
                    break;
                };
                line
            };
            let Some((key_text, value_text, value_offset)) = split_key_value(&line.content) else {
                return Err(syntax(line.number, indent + 1, "expected `key: value`"));
            };
            let key = parse_key(key_text, line.number, indent + 1)?;
            if map.contains_key(&key) {
                return Err(syntax(
                    line.number,
                    indent + 1,
                    format!("duplicate mapping key {key:?}"),
                ));
            }
            let value = if value_text.is_empty() {
                self.parse_child_or_null(indent)?
            } else {
                let column = indent
                    + 1
                    + line
                        .content
                        .get(..value_offset)
                        .map_or(0, |prefix| prefix.chars().count());
                let value = parse_inline(value_text, line.number, column)?;
                if let Some(next) = self.peek().filter(|next| next.indent > indent) {
                    return Err(syntax(
                        next.number,
                        next.indent + 1,
                        "unexpected indentation after a value",
                    ));
                }
                value
            };
            map.insert(key, value);
        }
        Ok(Value::Object(map))
    }

    fn parse_child_or_null(&mut self, indent: usize) -> Result<Value, WorkflowError> {
        let child = match self.peek() {
            Some(line) if line.indent > indent => Some(line.indent),
            Some(line) if line.indent == indent && sequence_rest(&line.content).is_some() => {
                Some(indent)
            }
            _ => None,
        };
        match child {
            Some(child) if child > indent => self.parse_node(child),
            Some(_) => self.parse_sequence(indent, None),
            None => Ok(Value::Null),
        }
    }

    fn parse_sequence(
        &mut self,
        indent: usize,
        first: Option<Line>,
    ) -> Result<Value, WorkflowError> {
        let mut items = Vec::new();
        let mut first = first;
        loop {
            let line = if let Some(line) = first.take() {
                line
            } else {
                let Some(line) = self.peek() else {
                    break;
                };
                if line.indent < indent {
                    break;
                }
                if line.indent > indent {
                    return Err(syntax(
                        line.number,
                        line.indent + 1,
                        "unexpected indentation in a sequence",
                    ));
                }
                if sequence_rest(&line.content).is_none() {
                    break;
                }
                let Some(line) = self.take() else {
                    break;
                };
                line
            };
            let Some((rest, spaces)) = sequence_rest(&line.content) else {
                return Err(syntax(line.number, indent + 1, "expected a sequence item"));
            };
            if rest.is_empty() {
                let child = self
                    .peek()
                    .and_then(|next| (next.indent > indent).then_some(next.indent));
                let value = match child {
                    Some(child) => self.parse_node(child)?,
                    None => Value::Null,
                };
                items.push(value);
                continue;
            }
            let item_indent = indent + 1 + spaces;
            if sequence_rest(rest).is_some() {
                let virtual_line = Line {
                    number: line.number,
                    indent: item_indent,
                    content: rest.to_owned(),
                };
                items.push(self.parse_sequence(item_indent, Some(virtual_line))?);
            } else if split_key_value(rest).is_some() {
                let virtual_line = Line {
                    number: line.number,
                    indent: item_indent,
                    content: rest.to_owned(),
                };
                items.push(self.parse_mapping(item_indent, Some(virtual_line))?);
            } else {
                let value = parse_inline(rest, line.number, item_indent + 1)?;
                if let Some(next) = self.peek().filter(|next| next.indent > indent) {
                    return Err(syntax(
                        next.number,
                        next.indent + 1,
                        "multi-line plain scalars are not supported",
                    ));
                }
                items.push(value);
            }
        }
        Ok(Value::Array(items))
    }
}

fn sequence_rest(content: &str) -> Option<(&str, usize)> {
    let rest = content.strip_prefix('-')?;
    if rest.is_empty() {
        return Some(("", 0));
    }
    if !rest.starts_with(' ') {
        return None;
    }
    let trimmed = rest.trim_start_matches(' ');
    let spaces = rest.len() - trimmed.len();
    Some((trimmed, spaces))
}

fn split_key_value(content: &str) -> Option<(&str, &str, usize)> {
    let mut quote: Option<char> = None;
    let mut escaped = false;
    let mut depth = 0usize;
    for (index, character) in content.char_indices() {
        if let Some(active) = quote {
            if active == '"' {
                if escaped {
                    escaped = false;
                } else if character == '\\' {
                    escaped = true;
                } else if character == '"' {
                    quote = None;
                }
            } else if character == '\'' {
                quote = None;
            }
            continue;
        }
        match character {
            '"' | '\'' => quote = Some(character),
            '[' | '{' => depth += 1,
            ']' | '}' => depth = depth.saturating_sub(1),
            ':' if depth == 0 => {
                let next = content
                    .get(index + 1..)
                    .and_then(|rest| rest.chars().next());
                let delimiter =
                    next.is_none() || next.is_some_and(|next| next == ' ' || next == '\t');
                if delimiter {
                    let key = content.get(..index)?;
                    let rest = content.get(index + 1..).unwrap_or("");
                    let value = rest.trim_start_matches(|character: char| {
                        character == ' ' || character == '\t'
                    });
                    let offset = index + 1 + (rest.len() - value.len());
                    return Some((key, value, offset));
                }
            }
            _ => {}
        }
    }
    None
}
