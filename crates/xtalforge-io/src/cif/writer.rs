//! CIF writer (Task 1.15.4): sites, displacement parameters, disorder,
//! charges and the program's explicit-bond loop.

use serde_json::Value;
use xtalforge_core::{Error, Provenance, Structure};

use crate::provenance::CIF_PROVENANCE_TAG;

/// Writes a structure as a CIF data block.
///
/// A structure carrying a [`Provenance`] record in its metadata also writes
/// the `_xtalforge_provenance` tag as a semicolon text block.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] if a site has no usable label or the
/// attached provenance record is invalid.
pub fn write_cif(structure: &Structure) -> Result<String, Error> {
    let mut out = String::new();
    out.push_str(&format!("data_{}\n", sanitize_name(structure.name())));
    if let Some(record) = Provenance::from_structure(structure)? {
        out.push_str(CIF_PROVENANCE_TAG);
        out.push('\n');
        out.push(';');
        out.push_str(&record.to_json()?);
        out.push_str("\n;\n");
    }
    let [a, b, c, alpha, beta, gamma] = structure.lattice().parameters();
    for (tag, value) in [
        ("_cell_length_a", a),
        ("_cell_length_b", b),
        ("_cell_length_c", c),
        ("_cell_angle_alpha", alpha),
        ("_cell_angle_beta", beta),
        ("_cell_angle_gamma", gamma),
    ] {
        out.push_str(&format!("{tag} {value:.6}\n"));
    }
    out.push_str(&format!("_space_group_name_Hall '{}'\n", structure.hall()));
    if !structure.operations().is_empty() {
        out.push_str("loop_\n_symmetry_equiv_pos_as_xyz\n");
        for operation in structure.operations() {
            out.push_str(&format!("'{}'\n", operation.triplet()));
        }
    }

    let has_uiso = structure.sites().iter().any(|site| site.u_iso().is_some());
    let has_disorder = structure
        .sites()
        .iter()
        .any(|site| site.disorder_group().is_some());
    let has_assembly = structure
        .sites()
        .iter()
        .any(|site| site.props().contains_key("disorder_assembly"));
    let has_charge = structure.sites().iter().any(|site| site.charge().is_some());
    let mut tags = vec![
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "_atom_site_occupancy",
    ];
    if has_uiso {
        tags.push("_atom_site_U_iso_or_equiv");
    }
    if has_disorder {
        tags.push("_atom_site_disorder_group");
    }
    if has_assembly {
        tags.push("_atom_site_disorder_assembly");
    }
    if has_charge {
        tags.push("_atom_site_charge");
    }
    out.push_str("loop_\n");
    for tag in &tags {
        out.push_str(tag);
        out.push('\n');
    }
    for site in structure.sites() {
        let [x, y, z] = site.frac();
        let mut row = vec![
            site.label().to_owned(),
            site.element().to_owned(),
            format_number(x),
            format_number(y),
            format_number(z),
            format_number(site.occupancy()),
        ];
        if has_uiso {
            row.push(site.u_iso().map_or("?".to_owned(), format_number));
        }
        if has_disorder {
            row.push(
                site.disorder_group()
                    .map_or(".".to_owned(), |group| group.to_string()),
            );
        }
        if has_assembly {
            row.push(
                site.props()
                    .get("disorder_assembly")
                    .and_then(Value::as_str)
                    .unwrap_or(".")
                    .to_owned(),
            );
        }
        if has_charge {
            row.push(site.charge().map_or("?".to_owned(), format_number));
        }
        out.push_str(&row.join(" "));
        out.push('\n');
    }

    let aniso: Vec<_> = structure
        .sites()
        .iter()
        .filter(|site| site.u_aniso().is_some())
        .collect();
    if !aniso.is_empty() {
        out.push_str("loop_\n_atom_site_aniso_label\n_atom_site_aniso_U_11\n");
        out.push_str("_atom_site_aniso_U_22\n_atom_site_aniso_U_33\n_atom_site_aniso_U_12\n");
        out.push_str("_atom_site_aniso_U_13\n_atom_site_aniso_U_23\n");
        for site in aniso {
            let Some([u11, u22, u33, u12, u13, u23]) = site.u_aniso() else {
                continue;
            };
            out.push_str(&format!(
                "{} {} {} {} {} {} {}\n",
                site.label(),
                format_number(u11),
                format_number(u22),
                format_number(u33),
                format_number(u12),
                format_number(u13),
                format_number(u23)
            ));
        }
    }

    if !structure.explicit_bonds().is_empty() {
        if structure.is_p1() {
            out.push_str("loop_\n_xtal_bond_kind\n_geom_bond_atom_site_label_1\n");
            out.push_str("_geom_bond_atom_site_label_2\n_xtal_bond_image\n");
            out.push_str("_xtal_bond_order\n");
            for bond in structure.explicit_bonds() {
                let first = structure
                    .sites()
                    .get(bond.i)
                    .map_or("?", |site| site.label());
                let second = structure
                    .sites()
                    .get(bond.j)
                    .map_or("?", |site| site.label());
                let kind = if bond.remove {
                    "suppressed"
                } else {
                    "explicit"
                };
                out.push_str(&format!(
                    "{kind} {first} {second} {},{},{} {}\n",
                    bond.image[0],
                    bond.image[1],
                    bond.image[2],
                    format_number(bond.order)
                ));
            }
        } else {
            out.push_str(&format!(
                "_xtal_bond_p1_count {}\n",
                structure.p1()?.n_atoms()
            ));
            out.push_str("loop_\n_xtal_bond_kind\n_xtal_bond_atom_1\n");
            out.push_str("_xtal_bond_atom_2\n_xtal_bond_image\n_xtal_bond_order\n");
            for bond in structure.explicit_bonds() {
                let kind = if bond.remove {
                    "suppressed"
                } else {
                    "explicit"
                };
                out.push_str(&format!(
                    "{kind} {} {} {},{},{} {}\n",
                    bond.i,
                    bond.j,
                    bond.image[0],
                    bond.image[1],
                    bond.image[2],
                    format_number(bond.order)
                ));
            }
        }
    }
    Ok(out)
}

fn format_number(value: f64) -> String {
    let text = format!("{value:.6}");
    text.trim_end_matches('0').trim_end_matches('.').to_owned()
}

fn sanitize_name(name: &str) -> String {
    name.chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || character == '_' {
                character
            } else {
                '_'
            }
        })
        .collect()
}
