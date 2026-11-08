//! Perception helpers for the UFF typer: aromatic rings, terminal-atom
//! refinement and UFF bond orders.
//!
//! Ported from Crystal-Builder `xtal/ff/uff/typer.py`. Aromaticity is a
//! planar five- or six-membered ring of C/N/O/S; terminal atoms are refined
//! against the bond length their type implies, because a structure without
//! hydrogens has no angles at a terminal carbon. Bond orders are a documented
//! subset of the reference core: aromatic bonds are 1.5, otherwise the
//! minimum pi capacity of the two types decides, and an amide C-N bond is
//! UFF's 1.41.

use crate::Error;
use crate::uff::geometry::{AtomView, Geometry, LocalBond};
use crate::uff::terms::natural_bond_length;
use crate::uff::typer::{AtomType, Confidence};

const REFINABLE: &[(&str, &[&str])] = &[
    ("C", &["C_3", "C_2", "C_1"]),
    ("N", &["N_3", "N_2", "N_1"]),
    ("O", &["O_3", "O_2"]),
    ("S", &["S_3+2", "S_2"]),
];

fn pi_capacity(name: &str) -> u8 {
    match name {
        "C_2" | "N_2" | "O_2" | "S_2" | "B_2" => 1,
        "C_1" | "N_1" => 2,
        _ => 0,
    }
}

fn implied_length(candidate: &str, partner: &str) -> Result<f64, Error> {
    let order = 1.0 + f64::from(pi_capacity(candidate).min(pi_capacity(partner)));
    natural_bond_length(candidate, partner, order)
}

pub(crate) fn refine_terminal(
    types: &mut [AtomType],
    view: &AtomView,
    geometry: &Geometry,
) -> Result<(), Error> {
    let seeds: Vec<String> = types.iter().map(|atom| atom.name.clone()).collect();
    for index in 0..types.len() {
        let Some(current) = types.get(index) else {
            continue;
        };
        if current.overridden || geometry.coordination(index) != 1 {
            continue;
        }
        let Some(element) = view.elements.get(index) else {
            continue;
        };
        let Some((_, candidates)) = REFINABLE
            .iter()
            .find(|(symbol, _)| *symbol == element.as_str())
        else {
            continue;
        };
        let Some(neighbor) = geometry.partners(index).first() else {
            continue;
        };
        let Some(partner) = seeds.get(neighbor.index) else {
            continue;
        };
        let Some(distance) = geometry.distances(view, index).first().copied() else {
            continue;
        };
        let mut scored: Vec<(f64, &str)> = candidates
            .iter()
            .map(|candidate| {
                implied_length(candidate, partner)
                    .map(|length| ((distance - length).abs(), *candidate))
            })
            .collect::<Result<_, _>>()?;
        scored.sort_by(|first, second| first.0.total_cmp(&second.0));
        let Some((_, best)) = scored.first().copied() else {
            continue;
        };
        let Some((_, runner_up)) = scored.get(1).copied() else {
            continue;
        };
        let best_length = implied_length(best, partner)?;
        let runner_length = implied_length(runner_up, partner)?;
        if let Some(slot) = types.get_mut(index) {
            *slot = AtomType::new(
                best,
                Confidence::Likely,
                format!(
                    "one bond to {partner} at {distance:.2} A, which suits {best} ({best_length:.2} A) better than {runner_up} ({runner_length:.2} A)"
                ),
            );
        }
    }
    Ok(())
}

pub(crate) fn aromatic_atoms(rings: &[Vec<usize>], view: &AtomView) -> Vec<bool> {
    let mut aromatic = vec![false; view.n_atoms()];
    for ring in rings {
        if !matches!(ring.len(), 5 | 6) {
            continue;
        }
        if !ring.iter().all(|index| {
            view.elements
                .get(*index)
                .is_some_and(|element| matches!(element.as_str(), "C" | "N" | "O" | "S"))
        }) {
            continue;
        }
        if !ring_is_planar(ring, view) {
            continue;
        }
        for index in ring {
            if let Some(slot) = aromatic.get_mut(*index) {
                *slot = true;
            }
        }
    }
    aromatic
}

fn ring_is_planar(ring: &[usize], view: &AtomView) -> bool {
    let Some(first) = ring
        .first()
        .and_then(|index| view.positions.get(*index))
        .copied()
    else {
        return false;
    };
    let Some(second) = ring
        .get(1)
        .and_then(|index| view.positions.get(*index))
        .copied()
    else {
        return false;
    };
    let first_vector = crate::uff::math::sub(second, first);
    let mut normal = None;
    for index in ring.iter().skip(2) {
        let Some(point) = view.positions.get(*index) else {
            continue;
        };
        let candidate = crate::uff::math::cross(first_vector, crate::uff::math::sub(*point, first));
        if crate::uff::math::norm(candidate) > 0.1 {
            normal = Some(candidate);
            break;
        }
    }
    let Some(normal) = normal else {
        return false;
    };
    let length = crate::uff::math::norm(normal);
    ring.iter().all(|index| {
        view.positions.get(*index).is_some_and(|point| {
            crate::uff::math::dot(normal, crate::uff::math::sub(*point, first)).abs() / length
                <= 0.1
        })
    })
}

pub(crate) fn assign_bond_orders(
    bonds: &[LocalBond],
    types: &[AtomType],
    aromatic: &[bool],
    aromatic_rings: &[Vec<usize>],
) -> Vec<f64> {
    let mut orders: Vec<f64> = bonds
        .iter()
        .map(|bond| {
            if bond.order != 1.0 {
                return bond.order;
            }
            let ring_bond = aromatic.get(bond.i).copied().unwrap_or(false)
                && aromatic.get(bond.j).copied().unwrap_or(false)
                && aromatic_rings
                    .iter()
                    .any(|ring| ring.contains(&bond.i) && ring.contains(&bond.j));
            if ring_bond {
                return 1.5;
            }
            let capacity = types
                .get(bond.i)
                .map(|atom| pi_capacity(&atom.name))
                .unwrap_or(0)
                .min(
                    types
                        .get(bond.j)
                        .map(|atom| pi_capacity(&atom.name))
                        .unwrap_or(0),
                );
            1.0 + f64::from(capacity)
        })
        .collect();
    let mut carbonyl_carbons = Vec::new();
    for (index, bond) in bonds.iter().enumerate() {
        if orders.get(index).copied() != Some(2.0) {
            continue;
        }
        let names = types
            .get(bond.i)
            .zip(types.get(bond.j))
            .map(|(first, second)| (first.name.as_str(), second.name.as_str()));
        if let Some((first, second)) = names
            && ((first == "C_2" && second == "O_2") || (first == "O_2" && second == "C_2"))
        {
            carbonyl_carbons.push(if first == "C_2" { bond.i } else { bond.j });
        }
    }
    for (index, bond) in bonds.iter().enumerate() {
        if bond.order != 1.0 || orders.get(index).copied() == Some(1.5) {
            continue;
        }
        let Some(first) = types.get(bond.i) else {
            continue;
        };
        let Some(second) = types.get(bond.j) else {
            continue;
        };
        let pair_has_carbon = first.name == "C_2" || second.name == "C_2";
        let pair_has_nitrogen = matches!(first.name.as_str(), "N_2" | "N_3" | "N_R")
            || matches!(second.name.as_str(), "N_2" | "N_3" | "N_R");
        if pair_has_carbon
            && pair_has_nitrogen
            && (carbonyl_carbons.contains(&bond.i) || carbonyl_carbons.contains(&bond.j))
            && let Some(slot) = orders.get_mut(index)
        {
            *slot = 1.41;
        }
    }
    orders
}
