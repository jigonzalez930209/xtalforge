//! Bond-angle targets for the reduced distance-geometry embedder.
//!
//! Each pair of bonds sharing an atom gets a 1-3 distance target from the law
//! of cosines with the ideal angle of the shared atom: 180 degrees for a
//! triple bond or two cumulated double bonds, 120 degrees for an aromatic or
//! double-bonded atom, and 109.47 degrees otherwise. Rings of three and four
//! atoms keep only their bond terms, since their angles cannot reach the ideal.

use std::collections::BTreeSet;

use crate::smiles::SmilesMolecule;

/// Tetrahedral angle, degrees.
const TETRAHEDRAL: f64 = 109.471_220_634_490_7;

/// Ideal angle at one atom, degrees.
fn ideal_angle(graph: &SmilesMolecule, atom: usize) -> f64 {
    let mut doubles = 0;
    let mut planar = graph.atoms.get(atom).is_some_and(|value| value.aromatic);
    for bond in &graph.bonds {
        if bond.a != atom && bond.b != atom {
            continue;
        }
        if bond.order >= 2.5 {
            return 180.0;
        }
        if (bond.order - 2.0).abs() < 1e-6 {
            doubles += 1;
            planar = true;
        }
        if (bond.order - 1.5).abs() < 1e-6 {
            planar = true;
        }
    }
    if doubles >= 2 {
        180.0
    } else if planar {
        120.0
    } else {
        TETRAHEDRAL
    }
}

/// 1-3 distance targets `(first, second, distance)` for every bond angle.
pub(super) fn angle_targets(
    graph: &SmilesMolecule,
    bonds: &[(usize, usize, f64)],
) -> Vec<(usize, usize, f64)> {
    let bonded: BTreeSet<(usize, usize)> = bonds
        .iter()
        .map(|&(a, b, _)| (a.min(b), a.max(b)))
        .collect();
    let mut targets = Vec::new();
    let mut seen = BTreeSet::new();
    for center in 0..graph.atoms.len() {
        let arms: Vec<(usize, f64)> = bonds
            .iter()
            .filter_map(|&(a, b, length)| {
                if a == center {
                    Some((b, length))
                } else if b == center {
                    Some((a, length))
                } else {
                    None
                }
            })
            .collect();
        let cosine = ideal_angle(graph, center).to_radians().cos();
        for (index, &(first, first_length)) in arms.iter().enumerate() {
            for &(second, second_length) in arms.iter().skip(index + 1) {
                let key = (first.min(second), first.max(second));
                if first == second || bonded.contains(&key) || in_four_ring(&bonded, center, key) {
                    continue;
                }
                if !seen.insert(key) {
                    continue;
                }
                let squared = first_length * first_length + second_length * second_length
                    - 2.0 * first_length * second_length * cosine;
                targets.push((key.0, key.1, squared.max(0.0).sqrt()));
            }
        }
    }
    targets
}

/// Whether `first`-`center`-`second` closes a four-membered ring.
fn in_four_ring(
    bonded: &BTreeSet<(usize, usize)>,
    center: usize,
    (first, second): (usize, usize),
) -> bool {
    let linked = |a: usize, b: usize| bonded.contains(&(a.min(b), a.max(b)));
    bonded.iter().any(|&(a, b)| {
        [a, b].into_iter().any(|fourth| {
            fourth != center
                && fourth != first
                && fourth != second
                && linked(fourth, first)
                && linked(fourth, second)
        })
    })
}
