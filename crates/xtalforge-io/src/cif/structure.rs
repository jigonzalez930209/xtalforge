//! CIF structure extraction (Task 1.15.4).

use std::collections::BTreeMap;

use xtalforge_core::{Error, ExplicitBond, Lattice, Site, Structure, StructureBuilder, SymOp};

use super::CifBlock;

pub(super) fn structure_from_block(block: &CifBlock) -> Result<Structure, Error> {
    let a = number(block.item("_cell_length_a")).ok_or_else(|| missing("_cell_length_a"))?;
    let b = number(block.item("_cell_length_b")).ok_or_else(|| missing("_cell_length_b"))?;
    let c = number(block.item("_cell_length_c")).ok_or_else(|| missing("_cell_length_c"))?;
    let alpha =
        number(block.item("_cell_angle_alpha")).ok_or_else(|| missing("_cell_angle_alpha"))?;
    let beta = number(block.item("_cell_angle_beta")).ok_or_else(|| missing("_cell_angle_beta"))?;
    let gamma =
        number(block.item("_cell_angle_gamma")).ok_or_else(|| missing("_cell_angle_gamma"))?;
    let lattice = Lattice::from_parameters(a, b, c, alpha, beta, gamma)?;
    let hall = block
        .item("_space_group_name_Hall")
        .or_else(|| block.item("_symmetry_space_group_name_Hall"))
        .or_else(|| block.item("_symmetry_space_group_name_H-M"))
        .or_else(|| block.item("_space_group_name_H-M_alt"))
        .unwrap_or("P 1")
        .trim_matches('\'')
        .to_owned();
    let mut operations = Vec::new();
    for loop_ in &block.loops {
        let column = loop_.tags.iter().position(|tag| {
            tag == "_symmetry_equiv_pos_as_xyz" || tag == "_space_group_symop_operation_xyz"
        });
        if let Some(column) = column {
            for row in &loop_.rows {
                if let Some(text) = row.get(column) {
                    operations.push(SymOp::parse(text)?);
                }
            }
        }
    }
    let atoms = block
        .loop_with("_atom_site_fract_x")
        .ok_or_else(|| Error::InvalidInput("CIF has no _atom_site loop".to_owned()))?;
    let oxidation = oxidation_numbers(block);
    let mut sites = Vec::new();
    let mut occupancy_clamped = 0u32;
    for row in 0..atoms.rows.len() {
        let label = atoms
            .value(row, "_atom_site_label")
            .or_else(|| atoms.value(row, "_atom_site_type_symbol"))
            .unwrap_or("X");
        let element = atoms.value(row, "_atom_site_type_symbol").unwrap_or(label);
        let x = number(atoms.value(row, "_atom_site_fract_x"))
            .ok_or_else(|| missing("_atom_site_fract_x"))?;
        let y = number(atoms.value(row, "_atom_site_fract_y"))
            .ok_or_else(|| missing("_atom_site_fract_y"))?;
        let z = number(atoms.value(row, "_atom_site_fract_z"))
            .ok_or_else(|| missing("_atom_site_fract_z"))?;
        let raw_occupancy = number(atoms.value(row, "_atom_site_occupancy")).unwrap_or(1.0);
        let occupancy = raw_occupancy.clamp(0.0, 1.0);
        if (occupancy - raw_occupancy).abs() > 1e-9 {
            occupancy_clamped += 1;
        }
        let mut site = Site::new(label, element, [x, y, z])?.with_occupancy(occupancy)?;
        if let Some(u_iso) = number(atoms.value(row, "_atom_site_U_iso_or_equiv")) {
            site = site.with_u_iso(u_iso)?;
        }
        if let Some(charge) = number(atoms.value(row, "_atom_site_charge")) {
            site = site.with_charge(charge)?;
        } else if let Some(charge) = oxidation.get(element) {
            site = site.with_charge(*charge)?;
        }
        if let Some(group) = atoms.value(row, "_atom_site_disorder_group") {
            let group = group.trim();
            if group != "." && group != "?" {
                match group.parse::<u32>() {
                    Ok(number) => site = site.with_disorder_group(number),
                    Err(_) => {
                        site = site.with_prop(
                            "disorder_group",
                            serde_json::Value::String(group.to_owned()),
                        );
                    }
                }
            }
        }
        if let Some(assembly) = atoms.value(row, "_atom_site_disorder_assembly") {
            let assembly = assembly.trim();
            if assembly != "." && assembly != "?" {
                site = site.with_prop(
                    "disorder_assembly",
                    serde_json::Value::String(assembly.to_owned()),
                );
            }
        }
        sites.push(site);
    }
    if let Some(aniso) = block.loop_with("_atom_site_aniso_U_11") {
        let has_label = aniso.column("_atom_site_aniso_label").is_some();
        for row in 0..aniso.rows.len() {
            let target = if has_label {
                aniso
                    .value(row, "_atom_site_aniso_label")
                    .and_then(|label| sites.iter().position(|site| site.label() == label))
            } else if row < sites.len() {
                Some(row)
            } else {
                None
            };
            let Some(target) = target else {
                continue;
            };
            let values = [
                number(aniso.value(row, "_atom_site_aniso_U_11")),
                number(aniso.value(row, "_atom_site_aniso_U_22")),
                number(aniso.value(row, "_atom_site_aniso_U_33")),
                number(aniso.value(row, "_atom_site_aniso_U_12")),
                number(aniso.value(row, "_atom_site_aniso_U_13")),
                number(aniso.value(row, "_atom_site_aniso_U_23")),
            ];
            if let [
                Some(u11),
                Some(u22),
                Some(u33),
                Some(u12),
                Some(u13),
                Some(u23),
            ] = values
                && let Some(site) = sites.get_mut(target)
            {
                *site = site.clone().with_u_aniso([u11, u22, u33, u12, u13, u23])?;
            }
        }
    }
    let mut raw_bonds = Vec::new();
    if let Some(bonds) = block.loop_with("_xtal_bond_kind") {
        for row in 0..bonds.rows.len() {
            let (Some(first), Some(second)) = (
                bonds.value(row, "_geom_bond_atom_site_label_1"),
                bonds.value(row, "_geom_bond_atom_site_label_2"),
            ) else {
                continue;
            };
            let (Some(i), Some(j)) = (
                sites.iter().position(|site| site.label() == first),
                sites.iter().position(|site| site.label() == second),
            ) else {
                continue;
            };
            let (operation, coded_image) = bonds
                .value(row, "_geom_bond_site_symmetry_2")
                .and_then(parse_symmetry_code)
                .unwrap_or((0, [0, 0, 0]));
            let image = bonds
                .value(row, "_xtal_bond_image")
                .and_then(parse_image)
                .unwrap_or(coded_image);
            let order = number(bonds.value(row, "_xtal_bond_order")).unwrap_or(1.0);
            let remove = bonds.value(row, "_xtal_bond_kind") == Some("suppressed");
            raw_bonds.push((i, j, operation, image, order, remove));
        }
    }
    let mut source = "operation loop";
    if operations.is_empty() && hall != "P 1" && hall != "P1" {
        match xtalforge_symmetry::SpaceGroup::parse(&hall) {
            Ok(group) => {
                operations = group.operations;
                source = "Hall symbol";
            }
            Err(_) => {
                source = "identity";
            }
        }
    }
    let has_operations = !operations.is_empty();
    // A P1 cell names its atoms one-to-one with the site loop. In a higher
    // symmetry the same pair needs the operation and tau bookkeeping of
    // `map_explicit_bond`; until that lands, keep the bonds in memory only
    // and do not invent a wrong periodic image.
    let explicit = if operations.is_empty() {
        raw_bonds
            .iter()
            .map(|&(i, j, _, image, order, remove)| ExplicitBond {
                i,
                j,
                image,
                order,
                remove,
            })
            .collect()
    } else {
        Vec::new()
    };
    let mut builder = StructureBuilder::new(block.name.clone(), lattice).hall(hall.clone());
    for operation in operations {
        builder = builder.operation(operation);
    }
    if !has_operations && hall != "P 1" && hall != "P1" {
        builder = builder.metadata(
            "symmetry_warning",
            serde_json::Value::String(
                "no symmetry operation loop and the Hall symbol was not recognised; P1 expansion uses the identity".to_owned(),
            ),
        );
    } else if has_operations {
        builder = builder.metadata(
            "symmetry_source",
            serde_json::Value::String(source.to_owned()),
        );
    }
    for bond in explicit {
        builder = builder.explicit_bond(bond);
    }
    for site in sites {
        builder = builder.site(site);
    }
    builder = builder.metadata("cif_block", serde_json::Value::String(block.name.clone()));
    if occupancy_clamped > 0 {
        builder = builder.metadata(
            "occupancy_clamped",
            serde_json::Value::from(occupancy_clamped),
        );
    }
    builder.build()
}

fn number(text: Option<&str>) -> Option<f64> {
    let text = text?.trim();
    if text.is_empty() || text == "?" || text == "." {
        return None;
    }
    let text = text.split('(').next().unwrap_or(text);
    text.parse().ok()
}

fn missing(tag: &str) -> Error {
    Error::InvalidInput(format!("CIF is missing {tag}"))
}

fn oxidation_numbers(block: &CifBlock) -> BTreeMap<String, f64> {
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

fn parse_image(text: &str) -> Option<[i32; 3]> {
    let values = text
        .split(',')
        .map(|value| value.trim().parse::<i32>().ok())
        .collect::<Option<Vec<_>>>()?;
    let [x, y, z] = <[i32; 3]>::try_from(values).ok()?;
    Some([x, y, z])
}

fn parse_symmetry_code(text: &str) -> Option<(usize, [i32; 3])> {
    let (operation, code) = text.trim().split_once('_')?;
    let operation = operation.parse::<usize>().ok()?.checked_sub(1)?;
    let code = code.parse::<u32>().ok()?;
    let digits = [code / 100, (code / 10) % 10, code % 10];
    Some((operation, digits.map(|digit| digit as i32 - 5)))
}
