//! Bond-perception rules, including pair-specific ranges (Task 14.3.2).

use std::collections::BTreeMap;

use serde_json::{Map, Value};

use crate::elements::{Element, element_by_symbol};
use crate::{Error, Structure};

/// Perception rules.
///
/// The global criterion is `minimum_distance <= d <= (r_i + r_j) * scale +
/// delta`, with `r` the Cordero covalent radii. A `pair_ranges` entry replaces
/// that criterion for one element pair, which is how a metal-oxygen contact is
/// allowed at a length no single scale can express. Metal-metal pairs are not
/// bonded unless [`allow_metal_metal`](Self::allow_metal_metal) is set, and the
/// dummy species `X` never bonds.
#[derive(Debug, Clone, PartialEq)]
pub struct BondRules {
    /// Multiplier on the sum of covalent radii.
    pub scale: f64,
    /// Constant added to the scaled radius sum, in Angstrom.
    pub delta: f64,
    /// Bonds shorter than this are ignored.
    pub minimum_distance: f64,
    /// Pair-specific `(minimum, maximum)` distance in Angstrom.
    ///
    /// The key is the two element symbols; either order names the same range.
    pub pair_ranges: BTreeMap<(String, String), (f64, f64)>,
    /// Whether two metal atoms may bond.
    pub allow_metal_metal: bool,
}

impl Default for BondRules {
    fn default() -> Self {
        Self {
            scale: 1.15,
            delta: 0.0,
            minimum_distance: 0.4,
            pair_ranges: BTreeMap::new(),
            allow_metal_metal: false,
        }
    }
}

impl BondRules {
    /// Creates the global rules with no pair-specific range.
    #[must_use]
    pub fn new(scale: f64, minimum_distance: f64) -> Self {
        Self {
            scale,
            minimum_distance,
            ..Self::default()
        }
    }

    /// Adds a pair-specific range and returns the rules.
    ///
    /// The pair is stored with its symbols sorted, so `("O", "Fe")` and
    /// `("Fe", "O")` name the same range.
    #[must_use]
    pub fn with_pair_range(
        mut self,
        first: impl Into<String>,
        second: impl Into<String>,
        minimum: f64,
        maximum: f64,
    ) -> Self {
        let first = first.into();
        let second = second.into();
        let key = sorted_pair(&first, &second);
        self.pair_ranges.insert(key, (minimum, maximum));
        self
    }

    /// Sets the metal-metal switch and returns the rules.
    #[must_use]
    pub fn with_metal_metal(mut self, allow: bool) -> Self {
        self.allow_metal_metal = allow;
        self
    }

    /// Checks every field and every pair-specific range.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite or negative value, a
    /// range whose maximum is below its minimum, or an unknown element symbol
    /// in a pair.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.scale.is_finite() || self.scale <= 0.0 {
            return Err(Error::InvalidInput(
                "bond scale must be positive and finite".to_owned(),
            ));
        }
        if !self.delta.is_finite() {
            return Err(Error::InvalidInput("bond delta must be finite".to_owned()));
        }
        if !self.minimum_distance.is_finite() || self.minimum_distance < 0.0 {
            return Err(Error::InvalidInput(
                "bond minimum distance must be non-negative and finite".to_owned(),
            ));
        }
        for ((first, second), (minimum, maximum)) in &self.pair_ranges {
            element_by_symbol(first)?;
            element_by_symbol(second)?;
            if !minimum.is_finite() || !maximum.is_finite() || *minimum < 0.0 {
                return Err(Error::InvalidInput(format!(
                    "bond range {first}-{second} must be finite and non-negative"
                )));
            }
            if maximum < minimum {
                return Err(Error::InvalidInput(format!(
                    "bond range {first}-{second} has its maximum below its minimum"
                )));
            }
        }
        Ok(())
    }

    /// Distance window `(minimum, maximum)` for one element pair in Angstrom.
    ///
    /// A stored pair range wins; otherwise the window is
    /// `(minimum_distance, (r_i + r_j) * scale + delta)`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an unknown element symbol.
    pub fn range(&self, first: &str, second: &str) -> Result<(f64, f64), Error> {
        Ok(self.range_of(element_by_symbol(first)?, element_by_symbol(second)?))
    }

    /// Whether a pair of elements may bond at all.
    ///
    /// Dummy species never bond, and two metals do not unless
    /// [`allow_metal_metal`](Self::allow_metal_metal) is set.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an unknown element symbol.
    pub fn allows(&self, first: &str, second: &str) -> Result<bool, Error> {
        Ok(self.allows_elements(element_by_symbol(first)?, element_by_symbol(second)?))
    }

    /// Search radius that cannot miss a bond among these symbols.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an unknown element symbol.
    pub fn max_cutoff(&self, symbols: &[String]) -> Result<f64, Error> {
        let mut widest: f64 = self
            .pair_ranges
            .values()
            .map(|(_, maximum)| *maximum)
            .fold(0.0, f64::max);
        let mut unique: Vec<&str> = symbols.iter().map(String::as_str).collect();
        unique.sort_unstable();
        unique.dedup();
        for (index, first) in unique.iter().enumerate() {
            for second in unique.iter().skip(index) {
                let (_, maximum) = self.range(first, second)?;
                widest = widest.max(maximum);
            }
        }
        Ok(widest)
    }

    /// Distance window for two resolved elements.
    pub(crate) fn range_of(&self, first: &Element, second: &Element) -> (f64, f64) {
        if let Some(range) = self.pair_range(&first.symbol, &second.symbol) {
            return range;
        }
        (
            self.minimum_distance,
            (first.covalent_radius + second.covalent_radius).mul_add(self.scale, self.delta),
        )
    }

    /// Whether two resolved elements may bond.
    pub(crate) fn allows_elements(&self, first: &Element, second: &Element) -> bool {
        if first.is_dummy() || second.is_dummy() {
            return false;
        }
        self.allow_metal_metal || !(first.is_metal && second.is_metal)
    }

    fn pair_range(&self, first: &str, second: &str) -> Option<(f64, f64)> {
        self.pair_ranges
            .iter()
            .find(|((a, b), _)| (a == first && b == second) || (a == second && b == first))
            .map(|(_, range)| *range)
    }

    /// JSON value stored under the structure's `bond_rules` metadata.
    #[must_use]
    pub fn to_json(&self) -> Value {
        let mut pairs = Map::new();
        for ((first, second), (minimum, maximum)) in &self.pair_ranges {
            pairs.insert(
                format!("{first}-{second}"),
                Value::Array(vec![Value::from(*minimum), Value::from(*maximum)]),
            );
        }
        let mut object = Map::new();
        object.insert("scale".to_owned(), Value::from(self.scale));
        object.insert("delta".to_owned(), Value::from(self.delta));
        object.insert(
            "minimum_distance".to_owned(),
            Value::from(self.minimum_distance),
        );
        object.insert("pair_ranges".to_owned(), Value::Object(pairs));
        object.insert(
            "allow_metal_metal".to_owned(),
            Value::from(self.allow_metal_metal),
        );
        Value::Object(object)
    }

    /// Reads rules from the JSON value stored under `bond_rules`.
    ///
    /// Missing keys take their default value, so a structure that stored the
    /// older two-field form keeps its scale and minimum distance.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-object value, an ill-typed
    /// field, a malformed pair key or a range that fails [`Self::validate`].
    pub fn from_json(value: &Value) -> Result<Self, Error> {
        let object = value.as_object().ok_or_else(|| {
            Error::InvalidInput("bond_rules metadata must be a JSON object".to_owned())
        })?;
        let default = Self::default();
        let allow_metal_metal = match object.get("allow_metal_metal") {
            Some(value) => value.as_bool().ok_or_else(|| {
                Error::InvalidInput("bond_rules allow_metal_metal must be a boolean".to_owned())
            })?,
            None => default.allow_metal_metal,
        };
        let mut rules = Self {
            scale: number(object, "scale", default.scale)?,
            delta: number(object, "delta", default.delta)?,
            minimum_distance: number(object, "minimum_distance", default.minimum_distance)?,
            pair_ranges: BTreeMap::new(),
            allow_metal_metal,
        };
        if let Some(ranges) = object.get("pair_ranges") {
            let ranges = ranges.as_object().ok_or_else(|| {
                Error::InvalidInput("bond_rules pair_ranges must be an object".to_owned())
            })?;
            for (key, value) in ranges {
                let (first, second) = key.split_once('-').ok_or_else(|| {
                    Error::InvalidInput(format!(
                        "bond_rules pair key {key:?} must be ELEMENT-ELEMENT"
                    ))
                })?;
                let window = value.as_array().and_then(|values| {
                    let [minimum, maximum] = values.as_slice() else {
                        return None;
                    };
                    Some((minimum.as_f64()?, maximum.as_f64()?))
                });
                let Some((minimum, maximum)) = window else {
                    return Err(Error::InvalidInput(format!(
                        "bond_rules pair {key:?} must be [minimum, maximum]"
                    )));
                };
                rules
                    .pair_ranges
                    .insert(sorted_pair(first, second), (minimum, maximum));
            }
        }
        rules.validate()?;
        Ok(rules)
    }

    /// Rules stored with a structure, or the defaults when none are stored.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Self::from_json`].
    pub fn from_structure(structure: &Structure) -> Result<Self, Error> {
        match structure.metadata().get("bond_rules") {
            Some(value) => Self::from_json(value),
            None => Ok(Self::default()),
        }
    }
}

fn number(object: &Map<String, Value>, key: &str, default: f64) -> Result<f64, Error> {
    match object.get(key) {
        Some(value) => value
            .as_f64()
            .ok_or_else(|| Error::InvalidInput(format!("bond_rules {key} must be a number"))),
        None => Ok(default),
    }
}

fn sorted_pair(first: &str, second: &str) -> (String, String) {
    if first <= second {
        (first.to_owned(), second.to_owned())
    } else {
        (second.to_owned(), first.to_owned())
    }
}
