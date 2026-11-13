//! The query half of the RCSR net search (Task 14.4.2).

#![allow(clippy::indexing_slicing)] // bounded local parse arrays

use std::collections::BTreeSet;

use xtalforge_core::Error;

use crate::netsearch::{NetFacts, TRANSITIVITY_FIELDS};

/// A search over the RCSR facts. An empty query matches every net.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct NetQuery {
    /// Case-insensitive substring of the name; empty matches everything.
    pub name: String,
    /// Coordination numbers that must all be present (or exactly present when
    /// [`NetQuery::exclusive`]).
    pub coordinations: BTreeSet<u8>,
    /// Match the coordination set exactly instead of as a subset.
    pub exclusive: bool,
    /// Inclusive space-group or plane-group number ranges.
    pub numbers: Vec<(u16, u16)>,
    /// Group symbols to match, normalized.
    pub symbols: Vec<String>,
    /// Transitivity `p q r s`; `None` matches anything.
    pub transitivity: [Option<u8>; TRANSITIVITY_FIELDS],
    /// Periodicity dimension, when restricted.
    pub periodicity: Option<u8>,
}

impl NetQuery {
    /// An empty query.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Whether the query matches every net.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.name.trim().is_empty()
            && self.coordinations.is_empty()
            && self.numbers.is_empty()
            && self.symbols.is_empty()
            && self.transitivity.iter().all(Option::is_none)
            && self.periodicity.is_none()
    }

    /// Restricts the name substring.
    #[must_use]
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = name.into();
        self
    }

    /// Requires every listed coordination number.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a zero coordination number.
    pub fn with_coordination(
        mut self,
        numbers: impl IntoIterator<Item = u8>,
    ) -> Result<Self, Error> {
        for number in numbers {
            if number == 0 {
                return Err(Error::InvalidInput(
                    "a coordination number must be positive".to_owned(),
                ));
            }
            self.coordinations.insert(number);
        }
        Ok(self)
    }

    /// Matches the coordination set exactly.
    #[must_use]
    pub fn exclusive(mut self) -> Self {
        self.exclusive = true;
        self
    }

    /// Requires one space-group or plane-group number.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a number outside 1 to 230.
    pub fn with_space_group(self, number: u16) -> Result<Self, Error> {
        self.with_space_group_range(number, number)
    }

    /// Requires a number in an inclusive range.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a range outside 1 to 230 or with
    /// `low > high`.
    pub fn with_space_group_range(mut self, low: u16, high: u16) -> Result<Self, Error> {
        if !(1..=230).contains(&low) || !(1..=230).contains(&high) || low > high {
            return Err(Error::InvalidInput(
                "space-group numbers are 1 to 230 and a range runs low to high".to_owned(),
            ));
        }
        self.numbers.push((low, high));
        Ok(self)
    }

    /// Requires a group symbol (normalized comparison).
    #[must_use]
    pub fn with_group_symbol(mut self, symbol: impl Into<String>) -> Self {
        self.symbols.push(normalize(&symbol.into()));
        self
    }

    /// Requires exact transitivity fields; `None` is a wildcard.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a zero count.
    pub fn with_transitivity(
        mut self,
        values: [Option<u8>; TRANSITIVITY_FIELDS],
    ) -> Result<Self, Error> {
        if values.contains(&Some(0)) {
            return Err(Error::InvalidInput(
                "a transitivity count must be positive".to_owned(),
            ));
        }
        self.transitivity = values;
        Ok(self)
    }

    /// Requires one periodicity dimension.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] unless the dimension is 2 or 3.
    pub fn with_periodicity(mut self, periodicity: u8) -> Result<Self, Error> {
        if !(2..=3).contains(&periodicity) {
            return Err(Error::InvalidInput(
                "a net is two- or three-periodic".to_owned(),
            ));
        }
        self.periodicity = Some(periodicity);
        Ok(self)
    }

    /// Reads a query from text, the way a search box would.
    ///
    /// `coordination` accepts `3,6` or `3,6-c`; `number` accepts `225`,
    /// `191,194` or `221-230`; `transitivity` accepts `1 1 * *`, `1,2` or the
    /// packed `11**`. Fields left off the end of the transitivity are
    /// wildcards.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] naming the field that cannot be read.
    pub fn parse(
        name: &str,
        coordination: &str,
        exclusive: bool,
        number: &str,
        transitivity: &str,
    ) -> Result<Self, Error> {
        let mut query = Self::new().with_name(name.trim());
        for token in tokens(coordination.replace("-c", " ").as_str()) {
            let value = token.parse::<u8>().map_err(|_| {
                Error::InvalidInput(format!(
                    "coordination {token:?} is not a number; write them like 3,6"
                ))
            })?;
            query = query.with_coordination([value])?;
        }
        if exclusive {
            query = query.exclusive();
        }
        for token in tokens(number) {
            let (low, high) = match token.split_once('-') {
                Some((low, high)) => (
                    low.parse::<u16>().map_err(|_| bad_number(&token))?,
                    high.parse::<u16>().map_err(|_| bad_number(&token))?,
                ),
                None => {
                    let value = token.parse::<u16>().map_err(|_| bad_number(&token))?;
                    (value, value)
                }
            };
            query = query.with_space_group_range(low, high)?;
        }
        query.with_transitivity(parse_transitivity(transitivity)?)
    }

    /// Whether one row satisfies the query.
    #[must_use]
    pub fn matches(&self, facts: &NetFacts) -> bool {
        let name = self.name.trim().to_lowercase();
        if !name.is_empty() && !facts.name.to_lowercase().contains(&name) {
            return false;
        }
        if let Some(periodicity) = self.periodicity
            && facts.periodicity != periodicity
        {
            return false;
        }
        if !self.coordinations.is_empty() {
            let present: BTreeSet<u8> = facts.coordination.iter().copied().collect();
            if self.exclusive {
                if present != self.coordinations {
                    return false;
                }
            } else if !self.coordinations.is_subset(&present) {
                return false;
            }
        }
        if !self.numbers.is_empty() {
            let Some(number) = facts.number else {
                return false;
            };
            if !self
                .numbers
                .iter()
                .any(|(low, high)| *low <= number && number <= *high)
            {
                return false;
            }
        }
        if !self.symbols.is_empty() {
            let group = normalize(&facts.group);
            if !self.symbols.contains(&group) {
                return false;
            }
        }
        self.transitivity
            .iter()
            .zip(facts.transitivity.iter())
            .all(|(wanted, known)| wanted.is_none() || wanted == known)
    }
}

pub(crate) fn normalize(text: &str) -> String {
    text.chars()
        .filter(|character| !character.is_whitespace() && *character != '_')
        .flat_map(char::to_lowercase)
        .collect()
}

fn tokens(text: &str) -> Vec<String> {
    text.split(|character: char| character.is_whitespace() || character == ',' || character == ';')
        .filter(|token| !token.is_empty())
        .map(str::to_owned)
        .collect()
}

fn bad_number(token: &str) -> Error {
    Error::InvalidInput(format!(
        "space group {token:?} is not a number or a range like 221-230"
    ))
}

fn parse_transitivity(text: &str) -> Result<[Option<u8>; TRANSITIVITY_FIELDS], Error> {
    let trimmed = text
        .trim()
        .trim_matches(|character| "[]()".contains(character))
        .trim();
    let mut fields = Vec::new();
    if trimmed.is_empty() {
        return Ok([None; TRANSITIVITY_FIELDS]);
    }
    let packed = tokens(trimmed);
    let parts: Vec<String> = if packed.len() == 1 {
        packed[0]
            .chars()
            .map(|character| character.to_string())
            .collect()
    } else {
        packed
    };
    if parts.len() > TRANSITIVITY_FIELDS {
        return Err(Error::InvalidInput(
            "transitivity has four fields at most (p q r s)".to_owned(),
        ));
    }
    for part in parts {
        if part == "*" {
            fields.push(None);
        } else if let Ok(value) = part.parse::<u8>() {
            if value == 0 {
                return Err(Error::InvalidInput(
                    "a transitivity count must be positive".to_owned(),
                ));
            }
            fields.push(Some(value));
        } else {
            return Err(Error::InvalidInput(format!(
                "transitivity {part:?} is not a count or *; write it like 1 2 * *"
            )));
        }
    }
    while fields.len() < TRANSITIVITY_FIELDS {
        fields.push(None);
    }
    let mut values = [None; TRANSITIVITY_FIELDS];
    for (slot, value) in values.iter_mut().zip(fields) {
        *slot = value;
    }
    Ok(values)
}
