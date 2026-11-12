//! MOF structure checker (Task 10.1.1).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure, element_info};

use crate::platform::util::{bond_adjacency, components, p1_and_bonds};

/// Severity of a MOF check finding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MofSeverity {
    /// Useful context.
    Info,
    /// A quality problem that does not invalidate the framework.
    Warning,
    /// The structure cannot be used without a fix.
    Error,
}

impl MofSeverity {
    /// Lower-case label used by the surfaces.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Info => "info",
            Self::Warning => "warning",
            Self::Error => "error",
        }
    }
}

/// One problem found in a MOF.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MofFinding {
    /// Stable code, for example `overlapping_atoms`.
    pub code: String,
    /// Severity.
    pub severity: MofSeverity,
    /// Human-readable message.
    pub message: String,
    /// P1 atom indices involved.
    pub atoms: Vec<usize>,
    /// Suggested agent verb, when one can fix it.
    pub fix: Option<String>,
}

/// The complete MOF check.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MofCheckReport {
    /// Findings sorted by severity and code.
    pub findings: Vec<MofFinding>,
    /// Number of findings per severity.
    pub counts: BTreeMap<String, usize>,
    /// Whether no warning or error was found.
    pub clean: bool,
}

/// Runs the MOF quality checks on a structure.
///
/// The checks are geometric and do not require the reference tables of
/// MOFChecker. They cover overlaps, unbound atoms and fragments, impossible
/// valences, disorder, metal coordination and charge balance.
///
/// # Errors
///
/// Returns the errors of the P1 expansion, bond perception, element lookup or
/// oxidation-state estimation.
pub fn check_mof(structure: &Structure) -> Result<MofCheckReport, Error> {
    let (p1, graph) = p1_and_bonds(structure)?;
    let adjacency = bond_adjacency(&graph, p1.sites().len());
    let mut findings = Vec::new();
    let mut add = |code: &str,
                   severity: MofSeverity,
                   message: String,
                   atoms: Vec<usize>,
                   fix: Option<&str>| {
        findings.push(MofFinding {
            code: code.to_owned(),
            severity,
            message,
            atoms,
            fix: fix.map(str::to_owned),
        });
    };
    // Overlapping atoms.
    if let Ok(pairs) = xtalforge_core::neighbors(&p1, 0.4) {
        for pair in pairs {
            if pair.i <= pair.j {
                add(
                    "overlapping_atoms",
                    MofSeverity::Error,
                    format!(
                        "atoms {0} and {1} are {2:.3} A apart",
                        pair.i, pair.j, pair.distance
                    ),
                    vec![pair.i, pair.j],
                    Some("features edit --remove"),
                );
            }
        }
    }
    // Unbound atoms and fragments.
    let groups = components(&adjacency);
    let main = groups
        .iter()
        .enumerate()
        .max_by_key(|(_, group)| group.len())
        .map(|(index, _)| index)
        .unwrap_or(0);
    for (index, group) in groups.iter().enumerate() {
        if group.len() == 1 {
            if let Some(&atom) = group.first() {
                add(
                    "unbound_atom",
                    MofSeverity::Warning,
                    format!("atom {atom} has no bond"),
                    vec![atom],
                    Some("features edit --remove"),
                );
            }
        } else if index != main {
            let counts = crate::platform::util::element_counts(&p1, group)?;
            add(
                "unbound_fragment",
                MofSeverity::Warning,
                format!(
                    "fragment {} with {} atoms is not bonded to the framework",
                    crate::platform::util::formula(&counts),
                    group.len()
                ),
                group.clone(),
                Some("mof solvent-remove"),
            );
        }
    }
    // Impossible valences.
    let mut order_sums = vec![0.0_f64; p1.sites().len()];
    for bond in &graph.bonds {
        if let Some(value) = order_sums.get_mut(bond.i) {
            *value += bond.order;
        }
        if bond.j != bond.i
            && let Some(value) = order_sums.get_mut(bond.j)
        {
            *value += bond.order;
        }
    }
    for (atom, &total) in order_sums.iter().enumerate() {
        let Some(site) = p1.sites().get(atom) else {
            continue;
        };
        let limit = match site.element() {
            "C" => 4.2,
            "H" | "F" | "Cl" | "Br" | "I" => 1.2,
            "O" => 2.2,
            "N" => 4.2,
            "B" => 4.2,
            "P" => 6.2,
            "S" => 6.2,
            _ => continue,
        };
        if total > limit {
            add(
                if site.element() == "C" && total > 4.2 {
                    "hypervalent_carbon"
                } else {
                    "impossible_valence"
                },
                MofSeverity::Error,
                format!(
                    "{} atom {atom} has total bond order {total:.2} (limit {limit:.1})",
                    site.element()
                ),
                vec![atom],
                Some("features edit --replace"),
            );
        }
        if site.element() == "C"
            && total < 3.2
            && adjacency.get(atom).is_some_and(|list| list.len() < 3)
        {
            add(
                "missing_hydrogen",
                MofSeverity::Warning,
                format!("carbon {atom} has only {total:.2} bond order; a hydrogen may be missing"),
                vec![atom],
                Some("hydrogens"),
            );
        }
    }
    // Disorder.
    for (index, site) in structure.sites().iter().enumerate() {
        if site.is_partial() || site.disorder_group().is_some() {
            add(
                "disordered_site",
                MofSeverity::Warning,
                format!(
                    "site {index} ({}) has occupancy {:.3}{}",
                    site.element(),
                    site.occupancy(),
                    site.disorder_group()
                        .map_or(String::new(), |group| format!(" in group {group}"))
                ),
                vec![index],
                Some("mof disorder-resolve"),
            );
        }
    }
    // Metal coordination.
    for (index, site) in p1.sites().iter().enumerate() {
        let info = match element_info(site.element()) {
            Ok(info) => info,
            Err(_) => continue,
        };
        if !info.is_metal {
            continue;
        }
        let degree = adjacency.get(index).map(Vec::len).unwrap_or(0);
        if degree < 4 {
            add(
                "missing_metal_coordination",
                MofSeverity::Warning,
                format!(
                    "{} atom {index} has only {degree} bonds; coordination is likely incomplete",
                    site.element()
                ),
                vec![index],
                Some("mof open-metal-sites"),
            );
        }
    }
    // Charge balance only makes sense for a metal-organic framework. A
    // covalent network without metals is reported as not assessed.
    let has_metal = p1
        .sites()
        .iter()
        .any(|site| element_info(site.element()).is_ok_and(|info| info.is_metal));
    if !has_metal {
        add(
            "charge_not_assessed",
            MofSeverity::Info,
            "no metal atoms; framework charge was not assessed".to_owned(),
            Vec::new(),
            None,
        );
    } else {
        match xtalforge_analysis::oxidation_states(&p1, &xtalforge_analysis::BvsOptions::default())
        {
            Ok(report) if !report.neutral => add(
                "charge_imbalance",
                MofSeverity::Error,
                format!(
                    "the framework carries a net charge of {:.3}",
                    report.total_charge
                ),
                Vec::new(),
                Some("mof charge-balance"),
            ),
            Ok(_) => {}
            Err(error) => add(
                "charge_unassigned",
                MofSeverity::Info,
                format!("oxidation states could not be assigned: {error}"),
                Vec::new(),
                None,
            ),
        }
    }
    findings.sort_by(|first, second| {
        second
            .severity
            .cmp(&first.severity)
            .then_with(|| first.code.cmp(&second.code))
            .then_with(|| first.atoms.cmp(&second.atoms))
    });
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for finding in &findings {
        *counts
            .entry(finding.severity.as_str().to_owned())
            .or_insert(0) += 1;
    }
    for label in ["info", "warning", "error"] {
        counts.entry(label.to_owned()).or_insert(0);
    }
    let clean = findings
        .iter()
        .all(|finding| finding.severity == MofSeverity::Info);
    Ok(MofCheckReport {
        findings,
        counts,
        clean,
    })
}
