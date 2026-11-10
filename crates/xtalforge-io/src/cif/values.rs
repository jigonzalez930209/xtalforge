//! Helpers of the CIF reader: numbers, symmetry codes and explicit-bond mapping.

use std::collections::{BTreeMap, BTreeSet};

use xtalforge_core::{Error, ExplicitBond, P1Cell, Site, SymOp};

use super::CifBlock;

pub(super) fn number(text: Option<&str>) -> Option<f64> {
    let text = text?.trim();
    if text.is_empty() || text == "?" || text == "." {
        return None;
    }
    let text = text.split('(').next().unwrap_or(text);
    text.parse().ok()
}

pub(super) fn missing(tag: &str) -> Error {
    Error::InvalidInput(format!("CIF is missing {tag}"))
}

pub(super) fn oxidation_numbers(block: &CifBlock) -> BTreeMap<String, f64> {
    let mut out = BTreeMap::new();
    if let Some(loop_) = block.loop_with("_atom_type_oxidation_number") {
        for row in 0..loop_.rows.len() {
            let (Some(symbol), Some(value)) = (
                loop_.value(row, "_atom_type_symbol"),
                number(loop_.value(row, "_atom_type_oxidation_number")),
            ) else {
                continue;
            };
            out.insert(symbol.to_owned(), value);
        }
    }
    out
}

pub(super) fn parse_image(text: &str) -> Option<[i32; 3]> {
    let values = text
        .split(',')
        .map(|value| value.trim().parse::<i32>().ok())
        .collect::<Option<Vec<_>>>()?;
    let [x, y, z] = <[i32; 3]>::try_from(values).ok()?;
    Some([x, y, z])
}

pub(super) fn parse_symmetry_code(text: &str) -> Option<(usize, [i32; 3])> {
    let (operation, code) = text.trim().split_once('_')?;
    let operation = operation.parse::<usize>().ok()?.checked_sub(1)?;
    let code = code.parse::<u32>().ok()?;
    let digits = [code / 100, (code / 10) % 10, code % 10];
    Some((operation, digits.map(|digit| digit as i32 - 5)))
}

type RawBond = (usize, usize, usize, [i32; 3], f64, bool);

pub(super) fn map_explicit(
    raw: &[RawBond],
    sites: &[Site],
    operations: &[SymOp],
    p1: &P1Cell,
) -> Vec<ExplicitBond> {
    let identity = [SymOp::identity()];
    let operations = if operations.is_empty() {
        identity.as_slice()
    } else {
        operations
    };
    let mut found = BTreeMap::new();
    for &(site_i, site_j, operation, image, order, remove) in raw {
        let (Some(near), Some(far), Some(partner)) = (
            sites.get(site_i),
            sites.get(site_j),
            operations.get(operation),
        ) else {
            continue;
        };
        let start = near.frac();
        let far = add(partner.apply(far.frac()), image.map(f64::from));
        for symmetry in operations {
            let here_raw = symmetry.apply(start);
            let (here, tau) = wrap_with_tau(here_raw);
            let there_raw = add(symmetry.apply(far), tau);
            let (there, _) = wrap_with_tau(there_raw);
            let (Some(a), Some(b)) = (find_atom(p1, here), find_atom(p1, there)) else {
                continue;
            };
            let (Some(atom_a), Some(atom_b)) = (p1.atoms.get(a), p1.atoms.get(b)) else {
                continue;
            };
            let bond_image = [
                round(there_raw[0] - atom_b.frac[0] - (here_raw[0] - atom_a.frac[0])),
                round(there_raw[1] - atom_b.frac[1] - (here_raw[1] - atom_a.frac[1])),
                round(there_raw[2] - atom_b.frac[2] - (here_raw[2] - atom_a.frac[2])),
            ];
            let (i, j, image) = canonical(a, b, bond_image);
            found.insert(
                (i, j, image, remove),
                ExplicitBond {
                    i,
                    j,
                    image,
                    order,
                    remove,
                },
            );
        }
    }
    let mut unique = BTreeMap::new();
    for bond in found.into_values() {
        if let Some(key) = contact_key(p1, &bond) {
            unique.entry(key).or_insert(bond);
        }
    }
    unique.into_values().collect()
}

pub(super) fn contact_key(p1: &P1Cell, bond: &ExplicitBond) -> Option<(usize, usize, [i64; 3])> {
    let first = p1.atoms.get(bond.i)?.frac;
    let second = p1.atoms.get(bond.j)?.frac;
    let mut delta = [
        second[0] + f64::from(bond.image[0]) - first[0],
        second[1] + f64::from(bond.image[1]) - first[1],
        second[2] + f64::from(bond.image[2]) - first[2],
    ];
    for value in &mut delta {
        *value -= value.round();
    }
    Some((
        bond.i,
        bond.j,
        delta.map(|value| (value * 1e6).round() as i64),
    ))
}

pub(super) fn wrap_with_tau(raw: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let floor = raw.map(f64::floor);
    (
        [raw[0] - floor[0], raw[1] - floor[1], raw[2] - floor[2]],
        [-floor[0], -floor[1], -floor[2]],
    )
}

pub(super) fn find_atom(p1: &P1Cell, frac: [f64; 3]) -> Option<usize> {
    p1.atoms
        .iter()
        .position(|atom| p1.lattice.distance(atom.frac, frac) < 1e-3)
}

pub(super) fn add(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] + second[0],
        first[1] + second[1],
        first[2] + second[2],
    ]
}

pub(super) fn round(value: f64) -> i32 {
    value.round() as i32
}

pub(super) fn canonical(i: usize, j: usize, image: [i32; 3]) -> (usize, usize, [i32; 3]) {
    if i > j
        || (i == j
            && (image[0] < 0
                || (image[0] == 0 && image[1] < 0)
                || (image[0] == 0 && image[1] == 0 && image[2] < 0)))
    {
        (j, i, [-image[0], -image[1], -image[2]])
    } else {
        (i, j, image)
    }
}

pub(super) fn same_operations(first: &[SymOp], second: &[SymOp]) -> bool {
    fn keys(operations: &[SymOp]) -> BTreeSet<([[i32; 3]; 3], [i64; 3])> {
        operations
            .iter()
            .map(|operation| {
                (
                    operation.rotation(),
                    operation
                        .translation()
                        .map(|value| (value.rem_euclid(1.0) * 1e8).round() as i64),
                )
            })
            .collect()
    }
    keys(first) == keys(second)
}
