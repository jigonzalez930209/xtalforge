//! Minimal streaming XML reader and tree for the engine outputs (Task 18.2).
//!
//! The reader walks one tag at a time, so a large `vasprun.xml` is never fully
//! materialised: the caller captures only the sections it needs with
//! [`read_start`]. The optional [`Node`] tree is built for small documents
//! such as the Quantum ESPRESSO schema output. Malformed markup returns a
//! typed [`Error::Parse`] instead of panicking.

use std::borrow::Cow;
use std::collections::BTreeMap;

use crate::Error;

use super::common::number;

/// One XML event.
pub(crate) enum Event<'a> {
    /// Opening tag. `empty` marks `<tag/>`.
    Start {
        /// Element name.
        name: &'a str,
        /// Decoded attributes in document order.
        attrs: Vec<(&'a str, String)>,
        /// Whether the tag closes itself.
        empty: bool,
    },
    /// Closing tag.
    End {
        /// Element name.
        name: &'a str,
    },
    /// Text between tags, with the predefined entities decoded.
    Text(Cow<'a, str>),
    /// End of input.
    Eof,
}

/// A streaming XML tag reader.
pub(crate) struct Reader<'a> {
    text: &'a str,
    index: usize,
}

impl<'a> Reader<'a> {
    /// Creates a reader over a document.
    pub(crate) fn new(text: &'a str) -> Self {
        Self { text, index: 0 }
    }

    /// Returns the next event.
    pub(crate) fn next(&mut self) -> Result<Event<'a>, Error> {
        loop {
            let Some(rest) = self.text.get(self.index..) else {
                return Ok(Event::Eof);
            };
            if rest.is_empty() {
                return Ok(Event::Eof);
            }
            let Some(offset) = rest.find('<') else {
                self.index = self.text.len();
                if rest.trim().is_empty() {
                    return Ok(Event::Eof);
                }
                return Ok(Event::Text(decode(rest)));
            };
            if offset > 0 {
                let Some(text) = rest.get(..offset) else {
                    return Err(malformed("XML text is not on a character boundary"));
                };
                self.index += offset;
                if text.trim().is_empty() {
                    continue;
                }
                return Ok(Event::Text(decode(text)));
            }
            if rest.starts_with("<!--") {
                self.skip("-->")?;
                continue;
            }
            if rest.starts_with("<![CDATA[") {
                let Some(body) = rest.get(9..) else {
                    return Err(malformed("CDATA section is not closed"));
                };
                let Some(end) = body.find("]]>") else {
                    return Err(malformed("CDATA section is not closed"));
                };
                let Some(text) = body.get(..end) else {
                    return Err(malformed("CDATA section is not on a character boundary"));
                };
                self.index += 9 + end + 3;
                return Ok(Event::Text(Cow::Borrowed(text)));
            }
            if rest.starts_with("<?") {
                self.skip("?>")?;
                continue;
            }
            if rest.starts_with("<!") {
                self.skip(">")?;
                continue;
            }
            let Some(end) = find_tag_end(rest) else {
                return Err(malformed("XML tag is not closed"));
            };
            let Some(inner) = rest.get(1..end) else {
                return Err(malformed("XML tag is not on a character boundary"));
            };
            self.index += end + 1;
            if let Some(close) = inner.strip_prefix('/') {
                let name = close.trim();
                if name.is_empty() {
                    return Err(malformed("XML closing tag has no name"));
                }
                return Ok(Event::End { name });
            }
            let empty = inner.trim_end().ends_with('/');
            let inner = inner.trim_end().trim_end_matches('/');
            let mut parts = inner.splitn(2, char::is_whitespace);
            let name = parts.next().unwrap_or("").trim();
            if name.is_empty() {
                return Err(malformed("XML tag has no name"));
            }
            let mut attrs = Vec::new();
            if let Some(rest) = parts.next() {
                parse_attributes(rest, &mut attrs)?;
            }
            return Ok(Event::Start { name, attrs, empty });
        }
    }

    fn skip(&mut self, marker: &str) -> Result<(), Error> {
        let Some(rest) = self.text.get(self.index..) else {
            return Err(malformed("XML comment is not closed"));
        };
        let Some(found) = rest.find(marker) else {
            return Err(Error::Parse(format!("XML is missing {marker}")));
        };
        self.index += found + marker.len();
        Ok(())
    }
}

/// A captured XML element.
#[derive(Debug, Clone)]
pub(crate) struct Node {
    /// Element name.
    pub(crate) name: String,
    /// Decoded attributes.
    pub(crate) attrs: BTreeMap<String, String>,
    /// Child elements.
    pub(crate) children: Vec<Node>,
    /// Concatenated direct text.
    pub(crate) text: String,
}

impl Node {
    fn new(name: &str, attrs: Vec<(&str, String)>) -> Self {
        Self {
            name: name.to_owned(),
            attrs: attrs
                .into_iter()
                .map(|(key, value)| (key.to_owned(), value))
                .collect(),
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

    /// Numeric text content.
    pub(crate) fn value(&self) -> Option<f64> {
        number(self.trimmed())
    }

    /// Numbers in the text content, in order.
    pub(crate) fn numbers(&self) -> Vec<f64> {
        self.text.split_whitespace().filter_map(number).collect()
    }

    /// First direct child with a name.
    pub(crate) fn child(&self, name: &str) -> Option<&Node> {
        self.children.iter().find(|child| child.name == name)
    }

    /// Every direct child with a name.
    pub(crate) fn children_named<'b>(
        &'b self,
        name: &'b str,
    ) -> impl Iterator<Item = &'b Node> + 'b {
        self.children.iter().filter(move |child| child.name == name)
    }

    /// First descendant with a name, depth first in document order.
    pub(crate) fn find(&self, name: &str) -> Option<&Node> {
        if self.name == name {
            return Some(self);
        }
        self.children.iter().find_map(|child| child.find(name))
    }

    /// Every descendant with a name, in document order.
    pub(crate) fn find_all<'b>(&'b self, name: &str, out: &mut Vec<&'b Node>) {
        if self.name == name {
            out.push(self);
        }
        for child in &self.children {
            child.find_all(name, out);
        }
    }

    /// A `<varray name="...">` with one `<v>` row per line.
    pub(crate) fn varray(&self, name: &str) -> Option<Vec<Vec<f64>>> {
        let array = self
            .children
            .iter()
            .find(|child| child.name == "varray" && child.attr("name") == Some(name))?;
        Some(
            array
                .children_named("v")
                .map(|row| {
                    row.trimmed()
                        .split_whitespace()
                        .filter_map(number)
                        .collect()
                })
                .collect(),
        )
    }
}

/// Reads the rest of the current element (whose start tag was just returned)
/// into a [`Node`], consuming its closing tag.
pub(crate) fn read_start(
    reader: &mut Reader<'_>,
    name: &str,
    attrs: Vec<(&str, String)>,
    empty: bool,
) -> Result<Node, Error> {
    let node = Node::new(name, attrs);
    if empty {
        return Ok(node);
    }
    let open = node.name.clone();
    let mut stack = vec![node];
    loop {
        if stack.len() > 512 {
            return Err(malformed("XML nesting is deeper than 512 elements"));
        }
        match reader.next()? {
            Event::Start { name, attrs, empty } => {
                let child = Node::new(name, attrs);
                if empty {
                    if let Some(parent) = stack.last_mut() {
                        parent.children.push(child);
                    }
                } else {
                    stack.push(child);
                }
            }
            Event::End { name } => {
                let Some(finished) = stack.pop() else {
                    return Err(malformed("XML closing tag without open tag"));
                };
                if finished.name != name {
                    return Err(Error::Parse(format!(
                        "XML closes <{}> with </{name}>",
                        finished.name
                    )));
                }
                if stack.is_empty() {
                    return Ok(finished);
                }
                if let Some(parent) = stack.last_mut() {
                    parent.children.push(finished);
                }
            }
            Event::Text(text) => {
                if let Some(current) = stack.last_mut() {
                    let trimmed = text.trim();
                    if !trimmed.is_empty() {
                        if !current.text.is_empty() {
                            current.text.push(' ');
                        }
                        current.text.push_str(trimmed);
                    }
                }
            }
            Event::Eof => {
                return Err(Error::Parse(format!("XML element <{open}> is not closed")));
            }
        }
    }
}

/// Reads the first element of a document into a [`Node`].
pub(crate) fn read_root(reader: &mut Reader<'_>) -> Result<Node, Error> {
    loop {
        match reader.next()? {
            Event::Start { name, attrs, empty } => return read_start(reader, name, attrs, empty),
            Event::Text(_) => continue,
            Event::End { name } => {
                return Err(Error::Parse(format!(
                    "XML closes </{name}> before the first element"
                )));
            }
            Event::Eof => return Err(malformed("XML document is empty")),
        }
    }
}

fn find_tag_end(text: &str) -> Option<usize> {
    let mut quote: Option<char> = None;
    for (index, character) in text.char_indices() {
        match quote {
            Some(active) if character == active => quote = None,
            Some(_) => {}
            None if character == '"' || character == '\'' => quote = Some(character),
            None if character == '>' => return Some(index),
            None => {}
        }
    }
    None
}

fn parse_attributes<'a>(text: &'a str, out: &mut Vec<(&'a str, String)>) -> Result<(), Error> {
    let mut rest = text;
    loop {
        rest = rest.trim_start();
        if rest.is_empty() {
            return Ok(());
        }
        let Some(equal) = rest.find('=') else {
            return Err(Error::Parse(format!(
                "XML attribute {:?} has no value",
                rest.trim()
            )));
        };
        let Some(name) = rest.get(..equal) else {
            return Err(malformed("XML attribute is not on a character boundary"));
        };
        let name = name.trim();
        if name.is_empty() || name.contains(char::is_whitespace) {
            return Err(Error::Parse(format!("XML attribute {name:?} has no name")));
        }
        let Some(after) = rest.get(equal + 1..) else {
            return Err(malformed("XML attribute value is missing"));
        };
        let after = after.trim_start();
        let Some(quote) = after.chars().next() else {
            return Err(malformed("XML attribute value is missing"));
        };
        if quote != '"' && quote != '\'' {
            return Err(Error::Parse(format!(
                "XML attribute {name:?} is not quoted"
            )));
        }
        let Some(body) = after.get(1..) else {
            return Err(malformed("XML attribute value is missing"));
        };
        let Some(end) = body.find(quote) else {
            return Err(Error::Parse(format!(
                "XML attribute {name:?} is not closed"
            )));
        };
        let Some(value) = body.get(..end) else {
            return Err(malformed("XML attribute is not on a character boundary"));
        };
        out.push((name, decode(value).into_owned()));
        rest = body.get(end + 1..).unwrap_or("");
    }
}

fn decode(text: &str) -> Cow<'_, str> {
    if !text.contains('&') {
        return Cow::Borrowed(text);
    }
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(start) = rest.find('&') {
        out.push_str(rest.get(..start).unwrap_or(""));
        let Some(tail) = rest.get(start..) else {
            break;
        };
        let Some(end) = tail.find(';') else {
            out.push_str(tail);
            return Cow::Owned(out);
        };
        let Some(entity) = tail.get(1..end) else {
            out.push_str(tail);
            return Cow::Owned(out);
        };
        match entity {
            "amp" => out.push('&'),
            "lt" => out.push('<'),
            "gt" => out.push('>'),
            "quot" => out.push('"'),
            "apos" => out.push('\''),
            _ => out.push_str(tail.get(..=end).unwrap_or(tail)),
        }
        rest = tail.get(end + 1..).unwrap_or("");
    }
    out.push_str(rest);
    Cow::Owned(out)
}

fn malformed(message: &str) -> Error {
    Error::Parse(format!("malformed XML: {message}"))
}
