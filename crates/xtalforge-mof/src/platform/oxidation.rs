//! Metal oxidation states and framework charge balance (Task 10.1.4).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure, element_info};

use crate::platform::util::{bond_adjacency, p1_and_bonds};

/// One metal oxidation-state assignment.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MetalState {
    /// Element symbol.
    pub element: String,
    /// Estimated oxidation state.
    pub state: i32,
    /// P1 atom indices with this element and state.
    pub atoms: Vec<usize>,
    /// Mean bond-valence sum of those atoms.
    pub bvs: f64,
}

/// Framework charge-balance report.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ChargeBalanceReport {
    /// Metal states, sorted by element.
    pub metals: Vec<MetalState>,
    /// Charge contributed by deprotonated carboxylate and similar linkers.
    pub linker_charge: f64,
    /// Charge of free counter-ions.
    pub counter_ion_charge: f64,
    /// Sum of the framework charges.
    pub framework_charge: f64,
    /// Total charge of the cell.
    pub total_charge: f64,
    /// Whether the cell is neutral within `1e-3`.
    pub neutral: bool,
}

/// Estimates metal oxidation states and the framework charge.
///
/// The metal states come from the bond-valence MAP search of Task 9.3.3. The
/// linker charge counts deprotonated carboxylate groups: a carbon bonded to
/// two oxygens, at least one of which is bonded to a metal and not to
/// hydrogen, contributes `-1`.
///
/// # Errors
///
/// Returns the errors of the P1 expansion, bond perception or oxidation-state
/// estimation.
pub fn charge_balance(structure: &Structure) -> Result<ChargeBalanceReport, Error> {
    let (p1, graph) = p1_and_bonds(structure)?;
    let adjacency = bond_adjacency(&graph, p1.sites().len());
    let oxidation =
        xtalforge_analysis::oxidation_states(&p1, &xtalforge_analysis::BvsOptions::default())?;
    let mut metals: BTreeMap<String, MetalState> = BTreeMap::new();
    for site in &oxidation.sites {
        let info = element_info(&site.element)?;
        if !info.is_metal {
            continue;
        }
        let entry = metals.entry(site.element.clone()).or_insert(MetalState {
            element: site.element.clone(),
            state: site.state,
            atoms: Vec::new(),
            bvs: 0.0,
        });
        entry.atoms.push(site.index);
        entry.bvs += site.bvs;
    }
    for entry in metals.values_mut() {
        entry.bvs /= entry.atoms.len().max(1) as f64;
    }
    let metal_charge: f64 = metals
        .values()
        .map(|entry| f64::from(entry.state) * entry.atoms.len() as f64)
        .sum();
    let mut linker_charge = 0.0;
    for (carbon, site) in p1.sites().iter().enumerate() {
        if site.element() != "C" {
            continue;
        }
        let neighbours = adjacency.get(carbon).map(Vec::as_slice).unwrap_or(&[]);
        let oxygens: Vec<usize> = neighbours
            .iter()
            .copied()
            .filter(|index| {
                p1.sites()
                    .get(*index)
                    .is_some_and(|site| site.element() == "O")
            })
            .collect();
        if oxygens.len() < 2 {
            continue;
        }
        let deprotonated = oxygens.iter().any(|oxygen| {
            let bonds = adjacency.get(*oxygen).map(Vec::as_slice).unwrap_or(&[]);
            let has_hydrogen = bonds.iter().any(|index| {
                p1.sites()
                    .get(*index)
                    .is_some_and(|site| site.element() == "H")
            });
            let has_metal = bonds.iter().any(|index| {
                p1.sites().get(*index).is_some_and(|site| {
                    element_info(site.element()).is_ok_and(|info| info.is_metal)
                })
            });
            !has_hydrogen && has_metal
        });
        if deprotonated {
            linker_charge -= 1.0;
        }
    }
    let counter_ion_charge = 0.0;
    let framework_charge = metal_charge + linker_charge;
    let total_charge = framework_charge + counter_ion_charge;
    Ok(ChargeBalanceReport {
        metals: metals.into_values().collect(),
        linker_charge,
        counter_ion_charge,
        framework_charge,
        total_charge,
        neutral: total_charge.abs() < 1e-3,
    })
}
