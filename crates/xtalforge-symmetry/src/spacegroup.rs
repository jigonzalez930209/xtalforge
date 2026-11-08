//! Space-group lookup and operations (Task 2.1.2 and 2.2).

use moyo::data::{HallNumber, Setting, hall_symbol_entry, iter_wyckoff_positions_from_hall_number};
use xtalforge_core::{Error, SymOp};

/// One Hall setting of a space group.
#[derive(Debug, Clone, PartialEq)]
pub struct SpaceGroup {
    /// ITA number (1 to 230).
    pub number: u32,
    /// Hall number (1 to 530).
    pub hall_number: u32,
    /// Hermann-Mauguin short symbol.
    pub hm_short: String,
    /// Hermann-Mauguin full symbol.
    pub hm_full: String,
    /// Hall symbol.
    pub hall_symbol: String,
    /// Setting description from the moyo database.
    pub setting: String,
    /// Centering letter (`P`, `A`, `B`, `C`, `I`, `R` or `F`).
    pub centering: char,
    /// Crystal system.
    pub crystal_system: &'static str,
    /// Lattice system.
    pub lattice_system: &'static str,
    /// Operations of the conventional cell.
    pub operations: Vec<SymOp>,
}

impl SpaceGroup {
    /// Builds the space group for a Hall number.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Symmetry`] for a Hall number outside 1 to 530 or a
    /// Hall symbol moyo cannot parse.
    pub fn from_hall_number(hall_number: u32) -> Result<Self, Error> {
        let entry = hall_entry(hall_number)?;
        let symbol = moyo::data::HallSymbol::from_hall_number(hall_number as HallNumber)
            .ok_or_else(|| Error::Symmetry(format!("unknown Hall symbol for {hall_number}")))?;
        let base = symbol.traverse();
        let mut expanded = base.clone();
        for operation in &base {
            for translation in &symbol.centering_translations {
                expanded.push(
                    operation.clone()
                        * moyo::base::Operation::new(nalgebra::Matrix3::identity(), *translation),
                );
            }
        }
        let operations = crate::convert::operations_from_moyo(&expanded)?;
        Ok(Self {
            number: entry.number as u32,
            hall_number,
            hm_short: entry.hm_short.to_owned(),
            hm_full: entry.hm_full.to_owned(),
            hall_symbol: entry.hall_symbol.to_owned(),
            setting: entry.setting.to_owned(),
            centering: format!("{:?}", entry.centering)
                .chars()
                .next()
                .unwrap_or('P'),
            crystal_system: crystal_system(entry.number as u32),
            lattice_system: lattice_system(entry.number as u32),
            operations,
        })
    }

    /// Builds the default setting of a space-group number.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Symmetry`] for a number outside 1 to 230.
    pub fn from_number(number: u32) -> Result<Self, Error> {
        if !(1..=230).contains(&number) {
            return Err(Error::Symmetry(format!(
                "space-group number must be in 1..=230, got {number}"
            )));
        }
        let hall = (1..=530u32)
            .find(|hall| hall_entry(*hall).is_ok_and(|entry| entry.number as u32 == number))
            .ok_or_else(|| Error::Symmetry(format!("no Hall setting for number {number}")))?;
        Self::from_hall_number(hall)
    }

    /// Parses a Hall number, H-M symbol, Hall symbol or number.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Symmetry`] when no setting matches.
    pub fn parse(text: &str) -> Result<Self, Error> {
        let text = text.trim();
        if let Ok(number) = text.parse::<u32>() {
            return Self::from_number(number);
        }
        let (base, setting_hint) = match text.split_once(':') {
            Some((base, hint)) => (base, Some(normalize(hint))),
            None => (text, None),
        };
        let wanted = normalize(base);
        for hall in 1..=530u32 {
            let Ok(entry) = hall_entry(hall) else {
                continue;
            };
            let setting_text = normalize(entry.setting);
            let setting_matches = setting_hint.as_ref().is_none_or(|hint| {
                setting_text.starts_with(hint.as_str())
                    || (hint == "1" && setting_text.contains("origin choice 1"))
                    || (hint == "2" && setting_text.contains("origin choice 2"))
            });
            if setting_matches
                && (normalize(entry.hm_short) == wanted
                    || normalize(entry.hm_full) == wanted
                    || normalize(entry.hall_symbol) == wanted)
            {
                return Self::from_hall_number(hall);
            }
        }
        Err(Error::Symmetry(format!(
            "unknown space-group symbol or Hall symbol: {text:?}"
        )))
    }

    /// Number of operations in the conventional cell.
    #[must_use]
    pub fn order(&self) -> usize {
        self.operations.len()
    }

    /// Whether the group contains the inversion operation.
    #[must_use]
    pub fn is_centrosymmetric(&self) -> bool {
        self.operations.iter().any(|operation| {
            operation.rotation() == [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
                && operation
                    .translation()
                    .iter()
                    .all(|value| value.abs() < 1e-9)
        })
    }

    /// Whether the group is chiral (no improper operation).
    #[must_use]
    pub fn is_chiral(&self) -> bool {
        self.operations.iter().all(|operation| {
            let [[a, b, c], [d, e, f], [g, h, i]] = operation.rotation();
            a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g) == 1
        })
    }
}

/// One Wyckoff position of a Hall setting.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WyckoffPosition {
    /// Multiplicity in the conventional cell.
    pub multiplicity: usize,
    /// Wyckoff letter.
    pub letter: char,
    /// Site-symmetry symbol.
    pub site_symmetry: String,
    /// Coordinate triplet template.
    pub coordinates: String,
}

impl WyckoffPosition {
    /// Lists the Wyckoff positions of a Hall setting.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Symmetry`] for a Hall number outside 1 to 530.
    pub fn list(hall_number: u32) -> Result<Vec<Self>, Error> {
        hall_entry(hall_number)?;
        Ok(
            iter_wyckoff_positions_from_hall_number(hall_number as HallNumber)
                .map(|position| Self {
                    multiplicity: position.multiplicity,
                    letter: position.letter,
                    site_symmetry: position.site_symmetry.to_owned(),
                    coordinates: position.coordinates.to_owned(),
                })
                .collect(),
        )
    }
}

/// All 530 Hall settings, ordered by Hall number.
///
/// # Errors
///
/// Returns the first lookup error, which cannot happen for numbers in range.
pub fn all_settings() -> Result<Vec<SpaceGroup>, Error> {
    (1..=530).map(SpaceGroup::from_hall_number).collect()
}

/// Crystal system from the ITA number.
#[must_use]
pub fn crystal_system(number: u32) -> &'static str {
    match number {
        1..=2 => "triclinic",
        3..=15 => "monoclinic",
        16..=74 => "orthorhombic",
        75..=142 => "tetragonal",
        143..=167 => "trigonal",
        168..=194 => "hexagonal",
        _ => "cubic",
    }
}

/// Lattice system from the ITA number.
#[must_use]
pub fn lattice_system(number: u32) -> &'static str {
    match number {
        1..=2 => "triclinic",
        3..=15 => "monoclinic",
        16..=74 => "orthorhombic",
        75..=142 => "tetragonal",
        143..=194 => "hexagonal",
        _ => "cubic",
    }
}

/// Maps a [`SettingChoice`] to moyo's setting.
pub(crate) fn setting(choice: crate::SettingChoice) -> Setting {
    match choice {
        crate::SettingChoice::Spglib => Setting::Spglib,
        crate::SettingChoice::Standard => Setting::Standard,
        crate::SettingChoice::HallNumber(number) => Setting::HallNumber(number as HallNumber),
    }
}

fn hall_entry(hall_number: u32) -> Result<moyo::data::HallSymbolEntry, Error> {
    if !(1..=530).contains(&hall_number) {
        return Err(Error::Symmetry(format!(
            "Hall number must be in 1..=530, got {hall_number}"
        )));
    }
    hall_symbol_entry(hall_number as HallNumber)
        .ok_or_else(|| Error::Symmetry(format!("unknown Hall number {hall_number}")))
}

fn normalize(text: &str) -> String {
    text.chars()
        .filter(|character| !character.is_whitespace() && *character != '_')
        .flat_map(char::to_lowercase)
        .collect()
}
