//! Size limits: arithmetic estimates made before anything is allocated.
//!
//! Two levels: a soft limit warns, a hard limit refuses. Profiles are
//! `Standard` (8 GB class), `Generous` (every limit doubled), `Huge`
//! (million-atom work, 64 GB class) and `WarnOnly` (never refuses). The
//! numbers are the Crystal-Builder table (Task 1.3) plus the `Huge` preset
//! (Task 15.2.5).

use std::fmt;

use crate::Error;

/// Undo history sites allowed in whole-structure steps.
pub const UNDO_ATOMS: u64 = 2_000_000;
/// Most recent undo steps that are never dropped.
pub const UNDO_KEEP: u64 = 5;
/// Analytic Hessian memory estimate per atom, in MiB.
pub const HESSIAN_MIB_PER_ATOM: f64 = 13.0;

/// Limit profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Profile {
    /// 8 GB-class machine.
    #[default]
    Standard,
    /// Every limit doubled.
    Generous,
    /// Warn only; never refuses.
    WarnOnly,
    /// Million-atom work: soft and hard limits raised to the scale of a
    /// 64 GB-class machine.
    Huge,
}

impl Profile {
    /// Parses a profile name, falling back to Standard for unknown values.
    #[must_use]
    pub fn parse(name: &str) -> Self {
        match name {
            "generous" => Self::Generous,
            "warn_only" | "warn-only" => Self::WarnOnly,
            "huge" => Self::Huge,
            _ => Self::Standard,
        }
    }

    /// Stable profile name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Standard => "standard",
            Self::Generous => "generous",
            Self::WarnOnly => "warn_only",
            Self::Huge => "huge",
        }
    }

    /// Wording for the allowance a refusal mentions.
    const fn allowance(self) -> &'static str {
        match self {
            Self::Huge => "the huge limits preset allows",
            Self::Standard | Self::Generous | Self::WarnOnly => {
                "the Large structures setting allows"
            }
        }
    }
}

/// A size budget.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Budget {
    /// Atoms drawn for a display range.
    Drawn,
    /// Atoms in a supercell.
    Supercell,
    /// Points in a porosity grid.
    Grid,
    /// Atoms in a carbon build.
    Carbon,
    /// Atoms in a nanotube build.
    Nanotube,
    /// Atoms in a graphene build.
    Graphene,
    /// Atoms in a phonon supercell.
    Phonons,
    /// Atoms in an analytic Hessian.
    Hessian,
}

impl Budget {
    /// Stable budget name, the inverse of [`Budget::parse`].
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Drawn => "drawn",
            Self::Supercell => "supercell",
            Self::Grid => "grid",
            Self::Carbon => "carbon",
            Self::Nanotube => "nanotube",
            Self::Graphene => "graphene",
            Self::Phonons => "phonons",
            Self::Hessian => "hessian",
        }
    }

    /// Parses a budget name; unknown names are an error.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an unknown budget.
    pub fn parse(name: &str) -> Result<Self, Error> {
        match name {
            "drawn" => Ok(Self::Drawn),
            "supercell" => Ok(Self::Supercell),
            "grid" => Ok(Self::Grid),
            "carbon" => Ok(Self::Carbon),
            "nanotube" => Ok(Self::Nanotube),
            "graphene" => Ok(Self::Graphene),
            "phonons" => Ok(Self::Phonons),
            "hessian" => Ok(Self::Hessian),
            _ => Err(Error::InvalidInput(format!("unknown limit budget: {name}"))),
        }
    }

    fn table(self) -> (u64, u64) {
        match self {
            Self::Drawn => (100_000, 400_000),
            Self::Supercell => (50_000, 200_000),
            Self::Grid => (1_500_000, 4_000_000),
            Self::Carbon => (10_000, 30_000),
            Self::Nanotube => (50_000, 200_000),
            Self::Graphene => (50_000, 200_000),
            Self::Phonons => (2_000, 8_000),
            Self::Hessian => (100, 300),
        }
    }

    /// Million-atom table used by [`Profile::Huge`].
    ///
    /// The geometry budgets cross the million-atom line. The phonon and
    /// Hessian budgets stay far smaller on purpose: a million-atom analytic
    /// Hessian needs tens of terabytes, so a profile that claims to allow it
    /// would be lying. Insertion-code builds stop at a few million atoms
    /// because the bond graph, not the raw atom count, is the cost.
    fn huge_table(self) -> (u64, u64) {
        match self {
            Self::Drawn => (4_000_000, 20_000_000),
            Self::Supercell => (2_000_000, 20_000_000),
            Self::Grid => (32_000_000, 256_000_000),
            Self::Carbon => (1_000_000, 5_000_000),
            Self::Nanotube => (1_000_000, 5_000_000),
            Self::Graphene => (1_000_000, 5_000_000),
            Self::Phonons => (20_000, 100_000),
            Self::Hessian => (1_000, 3_000),
        }
    }
}

/// Limit level of a verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LimitLevel {
    /// Inside the soft limit.
    Ok,
    /// Over the soft limit but under the hard limit.
    Warn,
    /// Over the hard limit.
    Refuse,
}

/// Result of a size check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Verdict {
    /// Limit level.
    pub level: LimitLevel,
    /// Estimated count.
    pub estimate: u64,
    /// Human-readable sentence; empty when the estimate fits.
    pub sentence: String,
}

impl LimitLevel {
    /// Stable level name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Ok => "ok",
            Self::Warn => "warn",
            Self::Refuse => "refuse",
        }
    }
}

impl Verdict {
    /// Whether the operation is refused.
    #[must_use]
    pub fn refused(&self) -> bool {
        self.level == LimitLevel::Refuse
    }

    /// Whether the operation is warned or refused.
    #[must_use]
    pub fn warned(&self) -> bool {
        self.level != LimitLevel::Ok
    }
}

impl fmt::Display for Verdict {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.sentence.is_empty() {
            write!(formatter, "{}", self.estimate)
        } else {
            formatter.write_str(&self.sentence)
        }
    }
}

/// Limits model with a profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Limits {
    profile: Profile,
}

impl Limits {
    /// Creates limits with a profile.
    #[must_use]
    pub const fn new(profile: Profile) -> Self {
        Self { profile }
    }

    /// Profile used by these limits.
    #[must_use]
    pub const fn profile(&self) -> Profile {
        self.profile
    }

    /// `(soft, hard)` for a budget; `hard` is `None` under `WarnOnly`.
    #[must_use]
    pub fn bounds(&self, budget: Budget) -> (u64, Option<u64>) {
        let (soft, hard) = match self.profile {
            Profile::Huge => budget.huge_table(),
            Profile::Standard | Profile::Generous | Profile::WarnOnly => budget.table(),
        };
        match self.profile {
            Profile::Standard => (soft, Some(hard)),
            Profile::Huge => (soft, Some(hard)),
            Profile::Generous => (2 * soft, Some(2 * hard)),
            Profile::WarnOnly => (soft, None),
        }
    }

    /// Judges an estimate against a budget.
    #[must_use]
    pub fn judge(
        &self,
        budget: Budget,
        estimate: u64,
        what: &str,
        fits: impl Fn(u64) -> String,
    ) -> Verdict {
        let (soft, hard) = self.bounds(budget);
        if estimate <= soft {
            return Verdict {
                level: LimitLevel::Ok,
                estimate,
                sentence: String::new(),
            };
        }
        if let Some(hard) = hard
            && estimate > hard
        {
            return Verdict {
                level: LimitLevel::Refuse,
                estimate,
                sentence: format!(
                    "{} over the {} {}; {}",
                    what.replace("{}", &count(estimate)),
                    count(hard),
                    self.profile.allowance(),
                    fits(hard)
                ),
            };
        }
        Verdict {
            level: LimitLevel::Warn,
            estimate,
            sentence: format!(
                "{} over the {} that is comfortable; {}",
                what.replace("{}", &count(estimate)),
                count(soft),
                fits(soft)
            ),
        }
    }
}

/// Formats a count as a sentence says it: `52,000` or `1.4 M`.
#[must_use]
pub fn count(value: u64) -> String {
    if value >= 1_000_000 {
        format!("{:.1} M", value as f64 / 1e6)
    } else {
        let text = value.to_string();
        let mut out = String::new();
        for (index, character) in text.chars().rev().enumerate() {
            if index > 0 && index % 3 == 0 {
                out.push(',');
            }
            out.push(character);
        }
        out.chars().rev().collect()
    }
}

/// Bounds for a profile and budget name.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown budget.
pub fn limits_bounds(profile: &str, budget: &str) -> Result<(u64, Option<u64>), Error> {
    let limits = Limits::new(Profile::parse(profile));
    Ok(limits.bounds(Budget::parse(budget)?))
}

/// Generic limit check by profile, budget and estimate.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown budget.
pub fn limits_check(profile: &str, budget: &str, estimate: u64) -> Result<Verdict, Error> {
    let limits = Limits::new(Profile::parse(profile));
    let budget = Budget::parse(budget)?;
    Ok(limits.judge(
        budget,
        estimate,
        "the operation would have {} items",
        |limit| format!("up to {} items fit", count(limit)),
    ))
}

mod checks;

pub use checks::{
    check_drawn, check_grid, check_hessian, check_phonons, check_supercell, drawn_atoms,
    grid_points, grid_shape, multiplier, supercell_atoms,
};
