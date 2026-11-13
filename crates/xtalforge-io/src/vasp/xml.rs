//! Minimal XML tree for the vasprun reader (Task 13.2.4).
//!
//! Only the subset found in `vasprun.xml` is modelled: elements, attributes,
//! text, comments, processing instructions and the five predefined entities.
//! It is deliberately small and allocation-light; a malformed document returns
//! an `invalid_input` error instead of panicking.

#![allow(clippy::indexing_slicing)] // byte scanner, boundaries checked

use std::collections::BTreeMap;

use xtalforge_core::Error;

/// One XML element.
#[derive(Debug, Clone)]
pub(crate) struct Node {
    /// Element name.
    pub name: String,
    /// Attributes.
    pub attrs: BTreeMap<String, String>,
    /// Child elements.
    pub children: Vec<Node>,
    /// Concatenated text directly inside this element.
    pub text: String,
}

impl Node {
    fn new(name: &str) -> Self {
        Self {
            name: name.to_owned(),
            attrs: BTreeMap::new(),
            children: Vec::new(),
            text: String::new(),
        }
    }

    /// Attribute value.
    pub(crate) fn attr(&self, name: &str) -> Option<&str> {
        self.attrs.get(name).map(String::as_str)
    }

    /// Trimmed text content.
    pub(crate) fn trimmed(&self) -> &str {
        self.text.trim()
    }

    /// First direct child with a name.
    pub(crate) fn child(&self, name: &str) -> Option<&Node> {
        self.children.iter().find(|child| child.name == name)
    }

    /// Every direct child with a name.
    pub(crate) fn children_named<'a>(
        &'a self,
        name: &'a str,
    ) -> impl Iterator<Item = &'a Node> + 'a {
        self.children.iter().filter(move |child| child.name == name)
    }

    /// First descendant (depth first, document order).
    pub(crate) fn find(&self, name: &str) -> Option<&Node> {
        if self.name == name {
            return Some(self);
        }
        self.children.iter().find_map(|child| child.find(name))
    }

    /// Every descendant with a name, in document order.
    pub(crate) fn find_all<'a>(&'a self, name: &str, out: &mut Vec<&'a Node>) {
        if self.name == name {
            out.push(self);
        }
        for child in &self.children {
            child.find_all(name, out);
        }
    }

    /// Numeric text content.
    pub(crate) fn number(&self) -> Option<f64> {
        self.trimmed().parse().ok()
    }

    /// A `<varray name="...">` with one `<v>` row per line.
    pub(crate) fn varray(&self, name: &str) -> Option<Vec<Vec<f64>>> {
        let array = self
            .children
            .iter()
            .find(|child| child.name == "varray" && child.attr("name") == Some(name))?;
        let rows = array
            .children_named("v")
            .map(|row| {
                row.trimmed()
                    .split_whitespace()
                    .filter_map(|value| value.parse().ok())
                    .collect::<Vec<f64>>()
            })
            .collect();
        Some(rows)
    }
}

/// Parses a document into a root node with the name `#document`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unclosed tag, a mismatched closing
/// tag or a document with no element.
pub(crate) fn parse(text: &str) -> Result<Node, Error> {
    let bytes = text.as_bytes();
    let mut index = 0usize;
    let mut stack: Vec<Node> = vec![Node::new("#document")];
    while index < bytes.len() {
        if bytes[index] == b'<' {
            let rest = &text[index..];
            if rest.starts_with("<!--") {
                index = skip_until(text, index + 4, "-->")?;
                continue;
            }
            if rest.starts_with("<?") {
                index = skip_until(text, index + 2, "?>")?;
                continue;
            }
            if rest.starts_with("<!") {
                index = skip_until(text, index + 2, ">")?;
                continue;
            }
            let end = find_tag_end(text, index + 1)?;
            let inner = &text[index + 1..end];
            if let Some(close) = inner.strip_prefix('/') {
                let name = close.trim();
                let Some(open) = stack.pop() else {
                    return Err(Error::InvalidInput(
                        "XML closing tag without open tag".to_owned(),
                    ));
                };
                if open.name != name {
                    return Err(Error::InvalidInput(format!(
                        "XML closes <{}> with </{name}>",
                        open.name
                    )));
                }
                attach(&mut stack, open)?;
                index = end + 1;
                continue;
            }
            let self_closing = inner.trim_end().ends_with('/');
            let inner = inner.trim_end().trim_end_matches('/');
            let mut parts = inner.splitn(2, char::is_whitespace);
            let name = parts.next().unwrap_or("").trim();
            if name.is_empty() {
                return Err(Error::InvalidInput("XML tag has no name".to_owned()));
            }
            let mut node = Node::new(name);
            if let Some(rest) = parts.next() {
                parse_attributes(rest, &mut node)?;
            }
            if self_closing {
                attach(&mut stack, node)?;
            } else {
                stack.push(node);
            }
            index = end + 1;
            continue;
        }
        let start = index;
        while index < bytes.len() && bytes[index] != b'<' {
            index += 1;
        }
        let text = decode(&text[start..index]);
        if !text.trim().is_empty()
            && let Some(parent) = stack.last_mut()
        {
            if !parent.text.is_empty() {
                parent.text.push(' ');
            }
            parent.text.push_str(text.trim());
        }
    }
    if stack.len() != 1 {
        let name = stack.last().map_or("#document", |node| node.name.as_str());
        return Err(Error::InvalidInput(format!(
            "XML element <{name}> is not closed"
        )));
    }
    let document = stack
        .pop()
        .ok_or_else(|| Error::InvalidInput("XML has no root".to_owned()))?;
    if document.children.is_empty() {
        return Err(Error::InvalidInput("XML has no root element".to_owned()));
    }
    Ok(document)
}

fn attach(stack: &mut [Node], node: Node) -> Result<(), Error> {
    let Some(parent) = stack.last_mut() else {
        return Err(Error::InvalidInput(
            "XML element outside the document".to_owned(),
        ));
    };
    parent.children.push(node);
    Ok(())
}

fn find_tag_end(text: &str, from: usize) -> Result<usize, Error> {
    let bytes = text.as_bytes();
    let mut quote: Option<u8> = None;
    let mut index = from;
    while index < bytes.len() {
        let byte = bytes[index];
        if let Some(active) = quote {
            if byte == active {
                quote = None;
            }
        } else if byte == b'"' || byte == b'\'' {
            quote = Some(byte);
        } else if byte == b'>' {
            return Ok(index);
        }
        index += 1;
    }
    Err(Error::InvalidInput("XML tag is not closed".to_owned()))
}

fn skip_until(text: &str, from: usize, marker: &str) -> Result<usize, Error> {
    text.get(from..)
        .and_then(|rest| rest.find(marker))
        .map(|found| from + found + marker.len())
        .ok_or_else(|| Error::InvalidInput(format!("XML is missing {marker}")))
}

fn parse_attributes(text: &str, node: &mut Node) -> Result<(), Error> {
    let bytes = text.as_bytes();
    let mut index = 0usize;
    while index < bytes.len() {
        while index < bytes.len() && bytes[index].is_ascii_whitespace() {
            index += 1;
        }
        if index >= bytes.len() {
            break;
        }
        let name_start = index;
        while index < bytes.len() && !bytes[index].is_ascii_whitespace() && bytes[index] != b'=' {
            index += 1;
        }
        let name = &text[name_start..index];
        while index < bytes.len() && bytes[index].is_ascii_whitespace() {
            index += 1;
        }
        if index >= bytes.len() || bytes[index] != b'=' {
            return Err(Error::InvalidInput(format!(
                "XML attribute {name:?} has no value"
            )));
        }
        index += 1;
        while index < bytes.len() && bytes[index].is_ascii_whitespace() {
            index += 1;
        }
        let Some(&quote) = bytes.get(index) else {
            return Err(Error::InvalidInput(
                "XML attribute value is missing".to_owned(),
            ));
        };
        if quote != b'"' && quote != b'\'' {
            return Err(Error::InvalidInput(format!(
                "XML attribute {name:?} is not quoted"
            )));
        }
        index += 1;
        let value_start = index;
        while index < bytes.len() && bytes[index] != quote {
            index += 1;
        }
        if index >= bytes.len() {
            return Err(Error::InvalidInput(format!(
                "XML attribute {name:?} is not closed"
            )));
        }
        node.attrs
            .insert(name.to_owned(), decode(&text[value_start..index]));
        index += 1;
    }
    Ok(())
}

fn decode(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(start) = rest.find('&') {
        out.push_str(&rest[..start]);
        let tail = &rest[start..];
        let Some(end) = tail.find(';') else {
            out.push_str(tail);
            return out;
        };
        let entity = &tail[1..end];
        let replacement = match entity {
            "amp" => Some('&'),
            "lt" => Some('<'),
            "gt" => Some('>'),
            "quot" => Some('"'),
            "apos" => Some('\''),
            _ => None,
        };
        match replacement {
            Some(character) => out.push(character),
            None => out.push_str(&tail[..=end]),
        }
        rest = &tail[end + 1..];
    }
    out.push_str(rest);
    out
}
