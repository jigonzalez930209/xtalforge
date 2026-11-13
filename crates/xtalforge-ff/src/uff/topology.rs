//! Topology enumeration for UFF: bonds, angles, torsions, inversions and the
//! periodic van der Waals pair list (Task 4.2.4).
//!
//! Ported from Crystal-Builder `xtal/ff/uff/calculator.py`. Everything here
//! is built once for a fixed bond graph; positions move under it. The pair
//! list is built once at the requested cutoff without a neighbour-list skin
//! because the [`crate::Calculator`] interface evaluates whole structures
//! rather than a sequence of optimisation steps; [`crate::UffModel`] keeps
//! the list for callers that do iterate.

use std::collections::{BTreeMap, HashMap};

use xtalforge_core::{Structure, neighbors};

use crate::Error;
use crate::uff::geometry::{AtomView, Geometry, LocalBond, PairKey, exclusions, pair_key};
use crate::uff::inversion::{InversionTerm, inversion_parameters};
use crate::uff::math::{add, merge_virial};
use crate::uff::params::get;
use crate::uff::terms::{
    AngleForm, AngleTerm, BondTerm, angle_force_constant, bond_force_constant,
    general_coefficients, natural_bond_length, vdw_pair,
};
use crate::uff::torsion::{TorsionTerm, torsion_parameters};
use crate::uff::typer::Typing;
use crate::uff::vdw::VanDerWaalsTerm;

#[derive(Debug, Clone)]
pub(crate) struct Topology {
    pub bonds: BondTerm,
    pub angles: AngleTerm,
    pub torsions: TorsionTerm,
    pub inversions: InversionTerm,
    pub vdw: VanDerWaalsTerm,
}

impl Topology {
    pub(crate) fn build(
        structure: &Structure,
        view: &AtomView,
        bonds: &[LocalBond],
        geometry: &Geometry,
        typing: &Typing,
        vdw_enabled: bool,
        vdw_cutoff: f64,
    ) -> Result<Self, Error> {
        Ok(Self {
            bonds: build_bonds(bonds, typing)?,
            angles: build_angles(geometry, typing)?,
            torsions: build_torsions(bonds, geometry, typing)?,
            inversions: build_inversions(geometry, typing)?,
            vdw: if vdw_enabled {
                build_vdw(structure, view, bonds, geometry, typing, vdw_cutoff)?
            } else {
                VanDerWaalsTerm::default()
            },
        })
    }

    pub(crate) fn energy_and_gradient(
        &self,
        positions: &[[f64; 3]],
        matrix: [[f64; 3]; 3],
        virial: Option<&mut [[f64; 3]; 3]>,
        parallel: bool,
    ) -> Result<(f64, Vec<[f64; 3]>), Error> {
        let mut gradient = vec![[0.0; 3]; positions.len()];
        let mut virial_local = [[0.0; 3]; 3];
        let mut energy = 0.0;
        let (part, part_gradient) = self.bonds.energy_and_gradient_with(
            positions,
            matrix,
            Some(&mut virial_local),
            parallel,
        )?;
        energy += part;
        add_gradient(&mut gradient, &part_gradient)?;
        let (part, part_gradient) = self.angles.energy_and_gradient_with(
            positions,
            matrix,
            Some(&mut virial_local),
            parallel,
        )?;
        energy += part;
        add_gradient(&mut gradient, &part_gradient)?;
        let (part, part_gradient) = self.torsions.energy_and_gradient_with(
            positions,
            matrix,
            Some(&mut virial_local),
            parallel,
        )?;
        energy += part;
        add_gradient(&mut gradient, &part_gradient)?;
        let (part, part_gradient) = self.inversions.energy_and_gradient_with(
            positions,
            matrix,
            Some(&mut virial_local),
            parallel,
        )?;
        energy += part;
        add_gradient(&mut gradient, &part_gradient)?;
        let (part, part_gradient) = self.vdw.energy_and_gradient_with(
            positions,
            matrix,
            Some(&mut virial_local),
            parallel,
        )?;
        energy += part;
        add_gradient(&mut gradient, &part_gradient)?;
        if let Some(target) = virial {
            merge_virial(target, virial_local);
        }
        Ok((energy, gradient))
    }
}

fn add_gradient(target: &mut [[f64; 3]], source: &[[f64; 3]]) -> Result<(), Error> {
    if target.len() != source.len() {
        return Err(Error::Topology(
            "term gradient length does not match the atom count".to_owned(),
        ));
    }
    for (slot, value) in target.iter_mut().zip(source.iter()) {
        *slot = add(*slot, *value);
    }
    Ok(())
}

fn type_name(typing: &Typing, index: usize) -> Result<&str, Error> {
    typing
        .types
        .get(index)
        .map(|atom| atom.name.as_str())
        .ok_or_else(|| Error::Topology(format!("no UFF type for atom {index}")))
}

fn order_of(typing: &Typing, bond: usize, fallback: f64) -> f64 {
    typing.bond_orders.get(bond).copied().unwrap_or(fallback)
}

fn build_bonds(bonds: &[LocalBond], typing: &Typing) -> Result<BondTerm, Error> {
    let mut term = BondTerm::default();
    for (index, bond) in bonds.iter().enumerate() {
        let first = type_name(typing, bond.i)?;
        let second = type_name(typing, bond.j)?;
        let order = order_of(typing, index, bond.order);
        let r0 = natural_bond_length(first, second, order)?;
        term.i.push(bond.i);
        term.j.push(bond.j);
        term.shift.push(bond.image);
        term.r0.push(r0);
        term.k.push(bond_force_constant(first, second, r0)?);
    }
    Ok(term)
}

type AngleKey = (String, String, String, u64, u64, u64);

fn build_angles(geometry: &Geometry, typing: &Typing) -> Result<AngleTerm, Error> {
    let mut term = AngleTerm::default();
    let mut cache: HashMap<AngleKey, f64> = HashMap::new();
    for centre in 0..geometry.neighbors.len() {
        let partners = geometry.partners(centre);
        if partners.len() < 2 {
            continue;
        }
        let centre_type = type_name(typing, centre)?;
        let theta0 = get(centre_type)?.theta0;
        let form = AngleForm::for_theta0(theta0);
        let (c0, c1, c2) = general_coefficients(theta0);
        for (first, partner_a) in partners.iter().enumerate() {
            for partner_b in partners.iter().skip(first + 1) {
                let type_a = type_name(typing, partner_a.index)?;
                let type_b = type_name(typing, partner_b.index)?;
                let order_a = order_of(typing, partner_a.bond, 1.0);
                let order_b = order_of(typing, partner_b.bond, 1.0);
                let key = (
                    type_a.to_owned(),
                    centre_type.to_owned(),
                    type_b.to_owned(),
                    order_a.to_bits(),
                    order_b.to_bits(),
                    theta0.to_bits(),
                );
                let force = if let Some(value) = cache.get(&key) {
                    *value
                } else {
                    let value = angle_force_constant(
                        type_a,
                        centre_type,
                        type_b,
                        theta0,
                        order_a,
                        order_b,
                    )?;
                    cache.insert(key, value);
                    value
                };
                term.i.push(partner_a.index);
                term.j.push(centre);
                term.k.push(partner_b.index);
                term.shift_i.push(partner_a.shift);
                term.shift_k.push(partner_b.shift);
                term.force.push(force);
                term.c0.push(c0);
                term.c1.push(c1);
                term.c2.push(c2);
                term.form.push(form);
            }
        }
    }
    Ok(term)
}

fn build_torsions(
    bonds: &[LocalBond],
    geometry: &Geometry,
    typing: &Typing,
) -> Result<TorsionTerm, Error> {
    let mut term = TorsionTerm::default();
    for (index, bond) in bonds.iter().enumerate() {
        let near = geometry.partners(bond.i);
        let far = geometry.partners(bond.j);
        let parameters = torsion_parameters(
            type_name(typing, bond.i)?,
            type_name(typing, bond.j)?,
            order_of(typing, index, bond.order),
            near.len(),
            far.len(),
        )?;
        if parameters.barrier == 0.0 {
            continue;
        }
        let scale =
            1.0 / ((near.len().saturating_sub(1) * far.len().saturating_sub(1)).max(1) as f64);
        for neighbor_i in near {
            if neighbor_i.index == bond.j && neighbor_i.shift == bond.image {
                continue;
            }
            for neighbor_l in far {
                let shift_l = [
                    bond.image[0] + neighbor_l.shift[0],
                    bond.image[1] + neighbor_l.shift[1],
                    bond.image[2] + neighbor_l.shift[2],
                ];
                if neighbor_l.index == bond.i && shift_l == [0, 0, 0] {
                    continue;
                }
                if neighbor_i.index == neighbor_l.index && neighbor_i.shift == shift_l {
                    continue;
                }
                term.i.push(neighbor_i.index);
                term.j.push(bond.i);
                term.k.push(bond.j);
                term.l.push(neighbor_l.index);
                term.shift_i.push(neighbor_i.shift);
                term.shift_k.push(bond.image);
                term.shift_l.push(shift_l);
                term.barrier.push(parameters.barrier * scale);
                term.periodicity.push(parameters.periodicity);
                term.cos_n_phi0.push(parameters.cos_n_phi0);
            }
        }
    }
    Ok(term)
}

fn build_inversions(geometry: &Geometry, typing: &Typing) -> Result<InversionTerm, Error> {
    let mut term = InversionTerm::default();
    for centre in 0..geometry.neighbors.len() {
        let partners = geometry.partners(centre);
        if partners.len() != 3 {
            continue;
        }
        let mut neighbour_types = Vec::with_capacity(3);
        for partner in partners {
            neighbour_types.push(type_name(typing, partner.index)?);
        }
        let parameters = inversion_parameters(type_name(typing, centre)?, &neighbour_types);
        if parameters.force == 0.0 {
            continue;
        }
        for out in 0..3 {
            let plane: Vec<_> = partners
                .iter()
                .enumerate()
                .filter(|(position, _)| *position != out)
                .map(|(_, partner)| *partner)
                .collect();
            let (Some(partner_j), Some(partner_k)) =
                (plane.first().copied(), plane.get(1).copied())
            else {
                continue;
            };
            let Some(partner_l) = partners.get(out) else {
                continue;
            };
            term.centre.push(centre);
            term.j.push(partner_j.index);
            term.k.push(partner_k.index);
            term.l.push(partner_l.index);
            term.shift_j.push(partner_j.shift);
            term.shift_k.push(partner_k.shift);
            term.shift_l.push(partner_l.shift);
            term.force.push(parameters.force);
            term.c0.push(parameters.c0);
            term.c1.push(parameters.c1);
            term.c2.push(parameters.c2);
        }
    }
    Ok(term)
}

fn build_vdw(
    structure: &Structure,
    view: &AtomView,
    bonds: &[LocalBond],
    geometry: &Geometry,
    typing: &Typing,
    cutoff: f64,
) -> Result<VanDerWaalsTerm, Error> {
    let excluded: BTreeMap<PairKey, ()> = exclusions(bonds, geometry);
    let mut term = VanDerWaalsTerm {
        cutoff,
        ..VanDerWaalsTerm::default()
    };
    for neighbor in neighbors(structure, cutoff)? {
        let (Some(atom_i), Some(atom_j)) = (view.local(neighbor.i), view.local(neighbor.j)) else {
            continue;
        };
        let key = pair_key(atom_i, atom_j, neighbor.image);
        if excluded.contains_key(&key) {
            continue;
        }
        let (x, depth) = vdw_pair(type_name(typing, atom_i)?, type_name(typing, atom_j)?)?;
        let ratio6 = (x / cutoff).powi(6);
        term.i.push(atom_i);
        term.j.push(atom_j);
        term.shift.push(neighbor.image);
        term.x.push(x);
        term.d.push(depth);
        term.shift_energy
            .push(depth * (ratio6 * ratio6 - 2.0 * ratio6));
    }
    Ok(term)
}
