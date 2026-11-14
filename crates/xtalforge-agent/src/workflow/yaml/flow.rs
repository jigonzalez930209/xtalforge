//! Scalars, quoted strings, flow collections and inline values.
//!
//! This is the value-level half of the workflow YAML reader in the parent
//! module. It has no knowledge of indentation; the block parser hands it the
//! text of one value and receives a [`serde_json::Value`].

use serde_json::{Map, Number, Value};

use super::super::WorkflowError;
use super::syntax;

pub(super) fn parse_key(text: &str, line: usize, column: usize) -> Result<String, WorkflowError> {
    let text = text.trim();
    if text.is_empty() {
        return Err(syntax(line, column, "mapping key must not be empty"));
    }
    if matches!(text.chars().next(), Some('&' | '*' | '!')) {
        return Err(syntax(
            line,
            column,
            "YAML anchors, aliases and tags are not supported",
        ));
    }
    if text == "<<" {
        return Err(syntax(line, column, "YAML merge keys are not supported"));
    }
    if text.starts_with('"') || text.starts_with('\'') {
        let mut parser = FlowParser::new(text, line, column);
        let value = parser.parse_value()?;
        parser.skip_spaces();
        if parser.peek().is_some() {
            return Err(parser.error("unexpected text after a quoted key"));
        }
        match value {
            Value::String(key) if !key.is_empty() => Ok(key),
            Value::String(_) => Err(syntax(line, column, "mapping key must not be empty")),
            _ => Err(syntax(line, column, "mapping keys must be strings")),
        }
    } else {
        Ok(text.to_owned())
    }
}

pub(super) fn parse_inline(text: &str, line: usize, column: usize) -> Result<Value, WorkflowError> {
    let text = text.trim();
    if text.is_empty() {
        return Ok(Value::Null);
    }
    match text.chars().next() {
        Some('&' | '*' | '!') => {
            return Err(syntax(
                line,
                column,
                "YAML anchors, aliases and tags are not supported",
            ));
        }
        Some('|' | '>') => {
            return Err(syntax(
                line,
                column,
                "block scalars (`|` and `>`) are not supported",
            ));
        }
        _ => {}
    }
    if text.starts_with("{{") {
        return plain_value(text, line, column);
    }
    let flow = text.starts_with('[')
        || text.starts_with('{')
        || text.starts_with('"')
        || text.starts_with('\'');
    if flow {
        let mut parser = FlowParser::new(text, line, column);
        let value = parser.parse_value()?;
        parser.skip_spaces();
        if parser.peek().is_some() {
            return Err(parser.error("unexpected text after a value"));
        }
        Ok(value)
    } else {
        plain_value(text, line, column)
    }
}

fn plain_value(text: &str, line: usize, column: usize) -> Result<Value, WorkflowError> {
    let text = text.trim();
    match text {
        "null" | "Null" | "NULL" | "~" => return Ok(Value::Null),
        "true" | "True" | "TRUE" => return Ok(Value::Bool(true)),
        "false" | "False" | "FALSE" => return Ok(Value::Bool(false)),
        _ => {}
    }
    if let Ok(value) = text.parse::<i64>() {
        return Ok(Value::Number(Number::from(value)));
    }
    if let Ok(value) = text.parse::<u64>() {
        return Ok(Value::Number(Number::from(value)));
    }
    let numeric = text
        .chars()
        .next()
        .is_some_and(|first| first.is_ascii_digit() || matches!(first, '+' | '-' | '.'));
    if numeric && let Ok(value) = text.parse::<f64>() {
        let Some(number) = Number::from_f64(value) else {
            return Err(syntax(
                line,
                column,
                "non-finite floats cannot be represented in JSON",
            ));
        };
        return Ok(Value::Number(number));
    }
    Ok(Value::String(text.to_owned()))
}

struct FlowParser {
    characters: Vec<char>,
    position: usize,
    line: usize,
    column: usize,
}

impl FlowParser {
    fn new(text: &str, line: usize, column: usize) -> Self {
        Self {
            characters: text.chars().collect(),
            position: 0,
            line,
            column,
        }
    }

    fn peek(&self) -> Option<char> {
        self.characters.get(self.position).copied()
    }

    fn bump(&mut self) -> Option<char> {
        let character = self.peek()?;
        self.position += 1;
        Some(character)
    }

    fn skip_spaces(&mut self) {
        while matches!(self.peek(), Some(' ' | '\t')) {
            self.position += 1;
        }
    }

    fn error(&self, message: impl Into<String>) -> WorkflowError {
        self.error_at(self.position, message)
    }

    fn error_at(&self, position: usize, message: impl Into<String>) -> WorkflowError {
        syntax(self.line, self.column + position, message)
    }

    fn parse_value(&mut self) -> Result<Value, WorkflowError> {
        self.skip_spaces();
        match self.peek() {
            Some('[') => self.parse_sequence(),
            Some('{') => self.parse_mapping(),
            Some('"') => self.parse_double_quoted().map(Value::String),
            Some('\'') => self.parse_single_quoted().map(Value::String),
            Some('&' | '*' | '!') => {
                Err(self.error("YAML anchors, aliases and tags are not supported"))
            }
            Some('|' | '>') => Err(self.error("block scalars are not supported in flow values")),
            Some(',' | ']' | '}') => Err(self.error("expected a value")),
            Some(_) => self.parse_plain(),
            None => Err(self.error("expected a value")),
        }
    }

    fn parse_sequence(&mut self) -> Result<Value, WorkflowError> {
        self.bump();
        let mut items = Vec::new();
        loop {
            self.skip_spaces();
            match self.peek() {
                Some(']') => {
                    self.bump();
                    return Ok(Value::Array(items));
                }
                None => {
                    return Err(self.error(
                        "unterminated flow sequence; put the whole collection on one line",
                    ));
                }
                _ => {}
            }
            items.push(self.parse_value()?);
            self.skip_spaces();
            match self.peek() {
                Some(',') => {
                    self.bump();
                }
                Some(']') => {
                    self.bump();
                    return Ok(Value::Array(items));
                }
                Some(_) => {
                    return Err(self.error("expected `,` or `]` after a flow sequence item"));
                }
                None => {
                    return Err(self.error(
                        "unterminated flow sequence; put the whole collection on one line",
                    ));
                }
            }
        }
    }

    fn parse_mapping(&mut self) -> Result<Value, WorkflowError> {
        self.bump();
        let mut map = Map::new();
        loop {
            self.skip_spaces();
            match self.peek() {
                Some('}') => {
                    self.bump();
                    return Ok(Value::Object(map));
                }
                None => {
                    return Err(self
                        .error("unterminated flow mapping; put the whole collection on one line"));
                }
                _ => {}
            }
            let key_position = self.position;
            let key = self.parse_key()?;
            if map.contains_key(&key) {
                return Err(self.error_at(key_position, format!("duplicate mapping key {key:?}")));
            }
            self.skip_spaces();
            let value = match self.peek() {
                Some(',' | '}') => Value::Null,
                _ => self.parse_value()?,
            };
            map.insert(key, value);
            self.skip_spaces();
            match self.peek() {
                Some(',') => {
                    self.bump();
                }
                Some('}') => {
                    self.bump();
                    return Ok(Value::Object(map));
                }
                Some(_) => {
                    return Err(self.error("expected `,` or `}` after a flow mapping entry"));
                }
                None => {
                    return Err(self
                        .error("unterminated flow mapping; put the whole collection on one line"));
                }
            }
        }
    }

    fn parse_key(&mut self) -> Result<String, WorkflowError> {
        self.skip_spaces();
        match self.peek() {
            Some('"') => {
                let key = self.parse_double_quoted()?;
                self.skip_spaces();
                if self.bump() != Some(':') {
                    return Err(self.error("expected `:` after a mapping key"));
                }
                Ok(key)
            }
            Some('\'') => {
                let key = self.parse_single_quoted()?;
                self.skip_spaces();
                if self.bump() != Some(':') {
                    return Err(self.error("expected `:` after a mapping key"));
                }
                Ok(key)
            }
            Some('&' | '*' | '!') => {
                Err(self.error("YAML anchors, aliases and tags are not supported"))
            }
            Some(',' | '}') | None => Err(self.error("expected a mapping key")),
            Some(_) => {
                let start = self.position;
                while let Some(character) = self.peek() {
                    if character == ':' {
                        break;
                    }
                    if character == ',' || character == '}' {
                        return Err(self.error("expected `:` after a mapping key"));
                    }
                    self.position += 1;
                }
                let key: String = self
                    .characters
                    .iter()
                    .skip(start)
                    .take(self.position - start)
                    .collect();
                let key = key.trim().to_owned();
                if key.is_empty() {
                    return Err(self.error_at(start, "mapping key must not be empty"));
                }
                if key == "<<" {
                    return Err(self.error_at(start, "YAML merge keys are not supported"));
                }
                if self.bump() != Some(':') {
                    return Err(self.error_at(start, "expected `:` after a mapping key"));
                }
                Ok(key)
            }
        }
    }

    fn parse_plain(&mut self) -> Result<Value, WorkflowError> {
        let start = self.position;
        while let Some(character) = self.peek() {
            if matches!(character, ',' | ']' | '}') {
                break;
            }
            self.position += 1;
        }
        let raw: String = self
            .characters
            .iter()
            .skip(start)
            .take(self.position - start)
            .collect();
        let raw = raw.trim_end();
        if raw.is_empty() {
            return Err(self.error_at(start, "expected a value"));
        }
        plain_value(raw, self.line, self.column + start)
    }

    fn parse_single_quoted(&mut self) -> Result<String, WorkflowError> {
        let start = self.position;
        self.bump();
        let mut value = String::new();
        loop {
            match self.bump() {
                Some('\'') => {
                    if self.peek() == Some('\'') {
                        self.position += 1;
                        value.push('\'');
                    } else {
                        return Ok(value);
                    }
                }
                Some(character) => value.push(character),
                None => return Err(self.error_at(start, "unterminated single-quoted scalar")),
            }
        }
    }

    fn parse_double_quoted(&mut self) -> Result<String, WorkflowError> {
        let start = self.position;
        self.bump();
        let mut value = String::new();
        loop {
            match self.peek() {
                Some('"') => {
                    self.bump();
                    return Ok(value);
                }
                Some('\\') => {
                    self.bump();
                    let Some(escape) = self.bump() else {
                        return Err(self.error_at(start, "unterminated escape sequence"));
                    };
                    match escape {
                        '"' => value.push('"'),
                        '\\' => value.push('\\'),
                        '/' => value.push('/'),
                        'n' => value.push('\n'),
                        't' => value.push('\t'),
                        'r' => value.push('\r'),
                        'b' => value.push('\u{0008}'),
                        'f' => value.push('\u{000c}'),
                        'u' => {
                            let code = self.take_hex(4, start)?;
                            if (0xD800..=0xDFFF).contains(&code) {
                                return Err(self.error_at(
                                    start,
                                    "surrogate unicode escapes are not supported",
                                ));
                            }
                            let Some(character) = char::from_u32(code) else {
                                return Err(self.error_at(start, "invalid unicode escape"));
                            };
                            value.push(character);
                        }
                        other => {
                            return Err(self
                                .error_at(start, format!("unknown escape sequence `\\{other}`")));
                        }
                    }
                }
                Some(character) => {
                    self.bump();
                    value.push(character);
                }
                None => return Err(self.error_at(start, "unterminated double-quoted scalar")),
            }
        }
    }

    fn take_hex(&mut self, count: usize, start: usize) -> Result<u32, WorkflowError> {
        let mut code = 0u32;
        for _ in 0..count {
            let Some(character) = self.bump() else {
                return Err(self.error_at(start, "incomplete unicode escape"));
            };
            let Some(digit) = character.to_digit(16) else {
                return Err(self.error_at(start, "invalid unicode escape"));
            };
            code = code * 16 + digit;
        }
        Ok(code)
    }
}
