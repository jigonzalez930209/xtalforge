//! Size limits: arithmetic estimates made before anything is allocated.
//!
//! Two levels: a soft limit warns, a hard limit refuses. Profiles are
//! `Standard` (8 GB class), `Generous` (every limit doubled) and `WarnOnly`
//! (never refuses). The numbers are the Crystal-Builder table (Task 1.3).

use std::fmt;

use crate::Error;
use crate::lattice::Lattice;

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
}

impl Profile {
    /// Parses a profile name, falling back to Standard for unknown values.
    #[must_use]
    pub fn parse(name: &str) -> Self {
        match name {
            "generous" => Self::Generous,
            "warn_only" | "warn-only" => Self::WarnOnly,
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
        let (soft, hard) = budget.table();
        match self.profile {
            Profile::Standard => (soft, Some(hard)),
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
                    "{} over the {} the Large structures setting allows; {}",
                    what.replace("{}", &count(estimate)),
                    count(hard),
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

/// Atoms drawn for a display range: cell atoms times the range volume.
#[must_use]
pub fn drawn_atoms(n_atoms: u64, ranges: [(f64, f64); 3]) -> u64 {
    let volume = ranges
        .iter()
        .map(|(low, high)| (high - low).max(0.0))
        .product::<f64>();
    (n_atoms as f64 * volume).round() as u64
}

/// `|det P|` for a 3x3 transformation matrix.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a singular matrix or non-finite values.
pub fn multiplier(p: [[f64; 3]; 3]) -> Result<u64, Error> {
    if !p.iter().flatten().all(|value| value.is_finite()) {
        return Err(Error::InvalidInput(
            "transformation matrix contains non-finite values".to_owned(),
        ));
    }
    let [[a, b, c], [d, e, f], [g, h, i]] = p;
    let determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    if determinant.abs() < 1e-9 {
        return Err(Error::InvalidInput(
            "transformation matrix is singular".to_owned(),
        ));
    }
    Ok(determinant.abs().round() as u64)
}

/// Atoms in the cell a transformation builds from `n_atoms`.
///
/// # Errors
///
/// Returns the errors of [`multiplier`].
pub fn supercell_atoms(n_atoms: u64, p: [[f64; 3]; 3]) -> Result<u64, Error> {
    Ok(n_atoms.saturating_mul(multiplier(p)?))
}

/// Grid shape from cell lengths and a target spacing.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive or non-finite spacing.
pub fn grid_shape(lattice: &Lattice, spacing: f64) -> Result<[u32; 3], Error> {
    if !spacing.is_finite() || spacing <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "grid spacing must be positive and finite, got {spacing}"
        )));
    }
    let spacing = spacing.max(0.05);
    Ok(lattice
        .lengths()
        .map(|length| (length / spacing).round().max(8.0) as u32))
}

/// Points in a porosity grid.
///
/// # Errors
///
/// Returns the errors of [`grid_shape`].
pub fn grid_points(lattice: &Lattice, spacing: f64) -> Result<u64, Error> {
    let [a, b, c] = grid_shape(lattice, spacing)?;
    Ok(u64::from(a) * u64::from(b) * u64::from(c))
}

fn cube_that_fits(per_cell: f64, limit: u64, cells: &str) -> String {
    let k = ((limit as f64 / per_cell.max(1.0)).powf(1.0 / 3.0) + 1e-9).floor() as u64;
    if k < 1 {
        format!("one cell alone is {} atoms", count(per_cell.round() as u64))
    } else {
        format!("{k} x {k} x {k}{cells} is the largest that fits")
    }
}

/// Checks the atoms a display range would draw.
#[must_use]
pub fn check_drawn(limits: &Limits, n_atoms: u64, ranges: [(f64, f64); 3]) -> Verdict {
    let estimate = drawn_atoms(n_atoms, ranges);
    limits.judge(
        Budget::Drawn,
        estimate,
        "the picture would draw {} atoms",
        |limit| cube_that_fits(n_atoms as f64, limit, " cells"),
    )
}

/// Checks a supercell transformation.
///
/// # Errors
///
/// Returns the errors of [`supercell_atoms`].
pub fn check_supercell(limits: &Limits, n_atoms: u64, p: [[f64; 3]; 3]) -> Result<Verdict, Error> {
    let estimate = supercell_atoms(n_atoms, p)?;
    Ok(limits.judge(
        Budget::Supercell,
        estimate,
        "the new cell would have {} atoms",
        |limit| cube_that_fits(n_atoms as f64, limit, ""),
    ))
}

/// Checks a porosity grid.
///
/// # Errors
///
/// Returns the errors of [`grid_points`].
pub fn check_grid(limits: &Limits, lattice: &Lattice, spacing: f64) -> Result<Verdict, Error> {
    let estimate = grid_points(lattice, spacing)?;
    Ok(limits.judge(
        Budget::Grid,
        estimate,
        "the grid would have {} points",
        |limit| {
            let mut step = (spacing / 0.05).ceil() * 0.05;
            step = step.max(0.05);
            while grid_points(lattice, step).unwrap_or(u64::MAX) > limit && step < 5.0 {
                step += 0.05;
            }
            format!(
                "at {step:.2} A it would have {}",
                grid_points(lattice, step)
                    .map(count)
                    .unwrap_or_else(|_| "an invalid count".to_owned())
            )
        },
    ))
}

/// Checks a phonon supercell.
///
/// # Errors
///
/// Returns the errors of [`supercell_atoms`].
pub fn check_phonons(limits: &Limits, n_atoms: u64, p: [[f64; 3]; 3]) -> Result<Verdict, Error> {
    let estimate = supercell_atoms(n_atoms, p)?;
    Ok(limits.judge(
        Budget::Phonons,
        estimate,
        "the phonon supercell would have {} atoms",
        |limit| {
            format!(
                "{}, or Gamma only needs none",
                cube_that_fits(n_atoms as f64, limit, "")
            )
        },
    ))
}

/// Checks an analytic Hessian over `n_atoms` free atoms.
#[must_use]
pub fn check_hessian(limits: &Limits, n_atoms: u64) -> Verdict {
    limits.judge(
        Budget::Hessian,
        n_atoms,
        &format!(
            "an analytic Hessian over {{}} atoms would hold about {:.1} GB",
            n_atoms as f64 * HESSIAN_MIB_PER_ATOM / 1024.0
        ),
        |limit| {
            format!(
                "up to {} atoms fit, and finite differences need no more than one evaluation does",
                count(limit)
            )
        },
    )
}
