//! The `spacegroup` command family.

use clap::Subcommand;

use xtalforge::symmetry::{SpaceGroup, SpaceGroupFilter, catalogue};

use crate::{Cli, CliError};

/// Prints one space group and its operations.
///
/// # Errors
///
/// Returns [`CliError`] when the symbol or number is unknown.
pub(crate) fn show(cli: &Cli, text: &str) -> Result<(), CliError> {
    let group = SpaceGroup::parse(text).map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::json!({
                "number": group.number,
                "hall_number": group.hall_number,
                "hm_short": group.hm_short,
                "hm_full": group.hm_full,
                "hall_symbol": group.hall_symbol,
                "setting": group.setting,
                "centering": group.centering.to_string(),
                "crystal_system": group.crystal_system,
                "lattice_system": group.lattice_system,
                "point_group": group.point_group,
                "schoenflies": group.schoenflies,
                "bravais": group.bravais,
                "laue_class": group.laue_class,
                "centrosymmetric": group.is_centrosymmetric(),
                "chiral": group.is_chiral(),
                "operations": group.operations.iter().map(|operation| operation.triplet()).collect::<Vec<_>>(),
            })
        );
    } else {
        println!(
            "{} (No. {}, Hall {}, setting {})",
            group.hm_short, group.number, group.hall_number, group.setting
        );
        println!(
            "point group {} / {} / {} ({})",
            group.point_group, group.schoenflies, group.bravais, group.laue_class
        );
        println!(
            "system {} / {} centering {} operations {}",
            group.crystal_system,
            group.lattice_system,
            group.centering,
            group.operations.len()
        );
        for operation in &group.operations {
            println!("  {}", operation.triplet());
        }
    }
    Ok(())
}

/// Lists catalogue entries matching the filters.
///
/// # Errors
///
/// Returns [`CliError`] when the catalogue cannot be built.
pub(crate) fn list(
    cli: &Cli,
    crystal_system: Option<&str>,
    centering: Option<char>,
    point_group: Option<&str>,
) -> Result<(), CliError> {
    let groups = catalogue(SpaceGroupFilter {
        crystal_system,
        centering,
        point_group,
    })
    .map_err(CliError::core)?;
    if cli.json {
        let items: Vec<_> = groups
            .iter()
            .map(|group| {
                serde_json::json!({
                    "number": group.number,
                    "hall_number": group.hall_number,
                    "hm_short": group.hm_short,
                    "setting": group.setting,
                    "centering": group.centering.to_string(),
                    "point_group": group.point_group,
                })
            })
            .collect();
        println!("{}", serde_json::json!(items));
    } else {
        for group in groups {
            println!(
                "{:>3} {:>3} {:<14} {:<10} {:>1} {}",
                group.number,
                group.hall_number,
                group.hm_short,
                group.point_group,
                group.centering,
                group.setting
            );
        }
    }
    Ok(())
}

/// Searches the catalogue by number, H-M symbol or Hall symbol substring.
///
/// # Errors
///
/// Returns [`CliError`] when the catalogue cannot be built.
pub(crate) fn search(cli: &Cli, text: &str) -> Result<(), CliError> {
    let wanted = normalize(text);
    let groups: Vec<SpaceGroup> = catalogue(SpaceGroupFilter::default())
        .map_err(CliError::core)?
        .into_iter()
        .filter(|group| {
            normalize(&group.hm_short).contains(&wanted)
                || normalize(&group.hall_symbol).contains(&wanted)
                || group.number.to_string() == wanted
        })
        .collect();
    if cli.json {
        println!(
            "{}",
            serde_json::json!(
                groups
                    .iter()
                    .map(|group| serde_json::json!({
                        "number": group.number,
                        "hall_number": group.hall_number,
                        "hm_short": group.hm_short,
                        "hall_symbol": group.hall_symbol,
                    }))
                    .collect::<Vec<_>>()
            )
        );
    } else {
        for group in groups {
            println!(
                "{:>3} {:>3} {:<14} {}",
                group.number, group.hall_number, group.hm_short, group.hall_symbol
            );
        }
    }
    Ok(())
}

fn normalize(text: &str) -> String {
    text.chars()
        .filter(|character| !character.is_whitespace() && *character != '_')
        .flat_map(char::to_lowercase)
        .collect()
}

#[derive(Debug, Subcommand)]
pub(crate) enum SpacegroupAction {
    /// Show one space group by symbol, Hall symbol or number.
    Show {
        /// H-M symbol, Hall symbol or number.
        symbol: String,
    },
    /// List catalogue entries.
    List {
        /// Filter by crystal system.
        #[arg(long)]
        system: Option<String>,
        /// Filter by centering letter.
        #[arg(long)]
        centering: Option<char>,
        /// Filter by point-group symbol.
        #[arg(long)]
        point_group: Option<String>,
    },
    /// Search by symbol or number substring.
    Search {
        /// Text to search for.
        text: String,
    },
}
