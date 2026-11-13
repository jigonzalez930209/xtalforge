//! Merging coincident sites (Task 14.3.3 support).

use crate::spatial::{FractionalGrid, perpendicular_widths};
use crate::{Error, ExplicitBond, Site, Structure, StructureBuilder, TopologyBond};

/// Merges same-element P1 atoms that sit on the same position.
///
/// Two atoms are duplicates when their periodic distance is at most
/// `tolerance` Angstrom. The lowest-numbered atom of each group is kept with
/// its site properties (occupancy, displacement parameters, charge, disorder
/// group and props); the explicit and topology bonds of the removed atoms are
/// moved onto the keeper and their lattice images are recomputed from the
/// keeper's position. A different element at the same position is disorder and
/// is left alone.
///
/// The result is P1: a duplicate is judged on the whole cell, and the operation
/// is offered next to `FindSymmetry`, which re-imposes a group afterwards. A
/// structure with no duplicates is returned unchanged.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive or non-finite tolerance,
/// or the errors of the P1 expansion and the rebuild.
pub fn merge_duplicates(structure: &Structure, tolerance: f64) -> Result<Structure, Error> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "merge tolerance must be positive and finite, got {tolerance}"
        )));
    }
    let p1 = structure.p1()?;
    if p1.atoms.len() < 2 {
        return Ok(structure.clone());
    }
    let widths = perpendicular_widths(&p1.lattice)?;
    let mut grid = FractionalGrid::new(widths, tolerance);
    let mut keeper: Vec<Option<usize>> = vec![None; p1.atoms.len()];
    let mut merged = false;
    for (index, atom) in p1.atoms.iter().enumerate() {
        let mut best: Option<usize> = None;
        grid.for_each_candidate(&atom.element, atom.frac, |candidate| {
            let Some(other) = p1.atoms.get(candidate) else {
                return true;
            };
            if p1.lattice.distance(other.frac, atom.frac) <= tolerance
                && best.is_none_or(|value| candidate < value)
            {
                best = Some(candidate);
            }
            true
        });
        if let Some(candidate) = best {
            keeper[index] = Some(candidate);
            merged = true;
        } else {
            grid.insert(&atom.element, atom.frac, index);
        }
    }
    if !merged {
        return Ok(structure.clone());
    }
    let mut mapping: Vec<Option<usize>> = vec![None; p1.atoms.len()];
    let mut sites: Vec<Site> = Vec::new();
    let mut kept_frac: Vec<[f64; 3]> = Vec::new();
    for (index, atom) in p1.atoms.iter().enumerate() {
        match keeper[index] {
            Some(candidate) => mapping[index] = mapping.get(candidate).copied().flatten(),
            None => {
                let target = sites.len();
                mapping[index] = Some(target);
                sites.push(copied_site(structure, atom, target)?);
                kept_frac.push(atom.frac);
            }
        }
    }
    let explicit: Vec<ExplicitBond> = structure
        .explicit_bonds()
        .iter()
        .filter_map(|bond| remap_explicit(bond, &p1.atoms, &mapping, &kept_frac))
        .collect();
    let topology: Vec<TopologyBond> = structure
        .topology_bonds()
        .iter()
        .filter_map(|bond| remap_topology(bond, &p1.atoms, &mapping, &kept_frac))
        .collect();
    let mut builder = StructureBuilder::new(structure.name(), p1.lattice);
    for site in sites {
        builder = builder.site(site);
    }
    for bond in explicit {
        builder = builder.explicit_bond(bond);
    }
    for bond in topology {
        builder = builder.topology_bond(bond);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

fn copied_site(structure: &Structure, atom: &crate::P1Atom, index: usize) -> Result<Site, Error> {
    let mut site = Site::new(
        format!("{}{}", atom.element, index + 1),
        &atom.element,
        atom.frac,
    )?
    .with_occupancy(atom.occupancy)?;
    if let Some(parent) = structure.sites().get(atom.site_idx) {
        if let Some(u_iso) = parent.u_iso() {
            site = site.with_u_iso(u_iso)?;
        }
        if let Some(u_aniso) = parent.u_aniso() {
            site = site.with_u_aniso(u_aniso)?;
        }
        if let Some(charge) = parent.charge() {
            site = site.with_charge(charge)?;
        }
        if let Some(group) = parent.disorder_group() {
            site = site.with_disorder_group(group);
        }
        for (key, value) in parent.props() {
            site = site.with_prop(key.clone(), value.clone());
        }
    }
    Ok(site)
}

fn remap_explicit(
    bond: &ExplicitBond,
    atoms: &[crate::P1Atom],
    mapping: &[Option<usize>],
    kept_frac: &[[f64; 3]],
) -> Option<ExplicitBond> {
    let (i, j, image) = remapped(bond.i, bond.j, bond.image, atoms, mapping, kept_frac)?;
    if i == j && image == [0, 0, 0] {
        return None;
    }
    Some(ExplicitBond {
        i,
        j,
        image,
        ..*bond
    })
}

fn remap_topology(
    bond: &TopologyBond,
    atoms: &[crate::P1Atom],
    mapping: &[Option<usize>],
    kept_frac: &[[f64; 3]],
) -> Option<TopologyBond> {
    let (i, j, image) = remapped(bond.i, bond.j, bond.image, atoms, mapping, kept_frac)?;
    if i == j && image == [0, 0, 0] {
        return None;
    }
    Some(TopologyBond { i, j, image })
}

fn remapped(
    i: usize,
    j: usize,
    image: [i32; 3],
    atoms: &[crate::P1Atom],
    mapping: &[Option<usize>],
    kept_frac: &[[f64; 3]],
) -> Option<(usize, usize, [i32; 3])> {
    let new_i = mapping.get(i).copied().flatten()?;
    let new_j = mapping.get(j).copied().flatten()?;
    let old_i = atoms.get(i)?;
    let old_j = atoms.get(j)?;
    let new_i_frac = kept_frac.get(new_i)?;
    let new_j_frac = kept_frac.get(new_j)?;
    let image = [
        round(
            image[0],
            old_j.frac[0] - old_i.frac[0],
            new_j_frac[0] - new_i_frac[0],
        ),
        round(
            image[1],
            old_j.frac[1] - old_i.frac[1],
            new_j_frac[1] - new_i_frac[1],
        ),
        round(
            image[2],
            old_j.frac[2] - old_i.frac[2],
            new_j_frac[2] - new_i_frac[2],
        ),
    ];
    Some((new_i, new_j, image))
}

fn round(image: i32, old_delta: f64, new_delta: f64) -> i32 {
    (f64::from(image) + old_delta - new_delta).round() as i32
}
