//! CIF structure extraction (Task 1.15.4).

use xtalforge_core::{Error, ExplicitBond, Lattice, Site, Structure, StructureBuilder, SymOp};

use super::CifBlock;
use super::coincident;
use super::values::{
    map_explicit, missing, number, oxidation_numbers, parse_image, parse_symmetry_code,
    same_operations,
};

pub(crate) fn structure_from_block(block: &CifBlock) -> Result<Structure, Error> {
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
        .or_else(|| block.loop_with("_atom_site_Cartn_x"))
        .ok_or_else(|| Error::InvalidInput("CIF has no _atom_site loop".to_owned()))?;
    let cartesian = atoms.column("_atom_site_fract_x").is_none();
    let first_model = atoms
        .value(0, "_atom_site_pdbx_PDB_model_num")
        .filter(|model| *model != "." && *model != "?")
        .map(str::to_owned);
    let oxidation = oxidation_numbers(block);
    let mut sites = Vec::new();
    let mut occupancy_clamped = 0u32;
    for row in 0..atoms.rows.len() {
        if let (Some(first), Some(model)) = (
            first_model.as_deref(),
            atoms
                .value(row, "_atom_site_pdbx_PDB_model_num")
                .filter(|model| *model != "." && *model != "?"),
        ) && model != first
        {
            continue;
        }
        let label = atoms
            .value(row, "_atom_site_label")
            .or_else(|| atoms.value(row, "_atom_site_auth_atom_id"))
            .or_else(|| atoms.value(row, "_atom_site_label_atom_id"))
            .or_else(|| atoms.value(row, "_atom_site_type_symbol"))
            .unwrap_or("X");
        let element = atoms.value(row, "_atom_site_type_symbol").unwrap_or(label);
        let (x, y, z, tags) = if cartesian {
            (
                number(atoms.value(row, "_atom_site_Cartn_x")),
                number(atoms.value(row, "_atom_site_Cartn_y")),
                number(atoms.value(row, "_atom_site_Cartn_z")),
                (
                    "_atom_site_Cartn_x",
                    "_atom_site_Cartn_y",
                    "_atom_site_Cartn_z",
                ),
            )
        } else {
            (
                number(atoms.value(row, "_atom_site_fract_x")),
                number(atoms.value(row, "_atom_site_fract_y")),
                number(atoms.value(row, "_atom_site_fract_z")),
                (
                    "_atom_site_fract_x",
                    "_atom_site_fract_y",
                    "_atom_site_fract_z",
                ),
            )
        };
        let x = x.ok_or_else(|| missing(tags.0))?;
        let y = y.ok_or_else(|| missing(tags.1))?;
        let z = z.ok_or_else(|| missing(tags.2))?;
        let frac = if cartesian {
            lattice.to_frac([x, y, z])
        } else {
            [x, y, z]
        };
        let raw_occupancy = number(atoms.value(row, "_atom_site_occupancy")).unwrap_or(1.0);
        let occupancy = raw_occupancy.clamp(0.0, 1.0);
        if (occupancy - raw_occupancy).abs() > 1e-9 {
            occupancy_clamped += 1;
        }
        let mut site = Site::new(label, element, frac)?.with_occupancy(occupancy)?;
        if let Some(u_iso) = number(atoms.value(row, "_atom_site_U_iso_or_equiv")) {
            site = site.with_u_iso(u_iso)?;
        } else if let Some(b_iso) = number(atoms.value(row, "_atom_site_B_iso_or_equiv")) {
            let u_iso = b_iso / (8.0 * std::f64::consts::PI * std::f64::consts::PI);
            site = site.with_u_iso(u_iso.max(0.0))?;
        }
        if let Some(charge) = number(atoms.value(row, "_atom_site_charge")) {
            site = site.with_charge(charge)?;
        } else if let Some(charge) = oxidation.get(element) {
            site = site.with_charge(*charge)?;
        }
        for (prop, tags) in [
            (
                "chain_id",
                ["_atom_site_auth_asym_id", "_atom_site_label_asym_id"],
            ),
            (
                "residue_name",
                ["_atom_site_auth_comp_id", "_atom_site_label_comp_id"],
            ),
            (
                "residue_number",
                ["_atom_site_auth_seq_id", "_atom_site_label_seq_id"],
            ),
        ] {
            if let Some(value) = tags
                .iter()
                .find_map(|tag| atoms.value(row, tag))
                .filter(|value| *value != "." && *value != "?")
            {
                site = site.with_prop(prop, serde_json::Value::String(value.to_owned()));
            }
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
    let indexed_count = block
        .item("_xtal_bond_p1_count")
        .and_then(|value| value.trim().parse::<usize>().ok());
    let mut raw_indices = Vec::new();
    if let Some(bonds) = block.loop_with("_xtal_bond_atom_1") {
        for row in 0..bonds.rows.len() {
            let (Some(i), Some(j)) = (
                bonds
                    .value(row, "_xtal_bond_atom_1")
                    .and_then(|value| value.parse::<usize>().ok()),
                bonds
                    .value(row, "_xtal_bond_atom_2")
                    .and_then(|value| value.parse::<usize>().ok()),
            ) else {
                continue;
            };
            let image = bonds
                .value(row, "_xtal_bond_image")
                .and_then(parse_image)
                .unwrap_or([0, 0, 0]);
            let order = number(bonds.value(row, "_xtal_bond_order")).unwrap_or(1.0);
            let remove = bonds.value(row, "_xtal_bond_kind") == Some("suppressed");
            raw_indices.push((i, j, image, order, remove));
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
    let operations_count = operations.len();
    let operation_warning = if has_operations && hall != "P 1" && hall != "P1" {
        match xtalforge_symmetry::SpaceGroup::parse(&hall) {
            Ok(group) if group.order() != operations_count => Some(format!(
                "operation loop has {operations_count} operations but {hall:?} (Hall {}) has {}",
                group.hall_symbol,
                group.order()
            )),
            Ok(group) if !same_operations(&group.operations, &operations) => Some(format!(
                "operation loop does not match {hall:?} (Hall {})",
                group.hall_symbol
            )),
            Ok(_) => None,
            Err(_) => Some(format!(
                "operation loop is present but the space-group symbol {hall:?} was not recognised"
            )),
        }
    } else {
        None
    };
    let needs_p1 = indexed_count.is_some() || !raw_bonds.is_empty();
    let base_p1 = if needs_p1 {
        let mut base = StructureBuilder::new(block.name.clone(), lattice).hall(hall.clone());
        for operation in &operations {
            base = base.operation(*operation);
        }
        for site in &sites {
            base = base.site(site.clone());
        }
        Some(base.build()?.p1()?)
    } else {
        None
    };
    let explicit = if let Some(count) = indexed_count {
        let p1 = base_p1.as_ref().ok_or_else(|| {
            Error::InvalidInput("CIF has no P1 cell for explicit bonds".to_owned())
        })?;
        if p1.n_atoms() != count {
            return Err(Error::InvalidInput(format!(
                "CIF explicit bonds declare {count} P1 atoms but the expansion has {}",
                p1.n_atoms()
            )));
        }
        raw_indices
            .into_iter()
            .map(|(i, j, image, order, remove)| ExplicitBond {
                i,
                j,
                image,
                order,
                remove,
            })
            .collect()
    } else if raw_bonds.is_empty() {
        Vec::new()
    } else if operations.is_empty() {
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
        let p1 = base_p1.as_ref().ok_or_else(|| {
            Error::InvalidInput("CIF has no P1 cell for explicit bonds".to_owned())
        })?;
        map_explicit(&raw_bonds, &sites, &operations, p1)
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
        if let Some(warning) = operation_warning {
            builder = builder.metadata("symmetry_warning", serde_json::Value::String(warning));
        }
    }
    for bond in explicit {
        builder = builder.explicit_bond(bond);
    }
    for site in sites {
        builder = builder.site(site);
    }
    builder = builder.metadata("cif_block", serde_json::Value::String(block.name.clone()));
    if let Some(formula) = block
        .item("_chemical_formula_sum")
        .or_else(|| block.item("_chemical_formula_structural"))
    {
        builder = builder.metadata(
            "chemical_formula",
            serde_json::Value::String(formula.trim_matches('\'').to_owned()),
        );
    }
    if occupancy_clamped > 0 {
        builder = builder.metadata(
            "occupancy_clamped",
            serde_json::Value::from(occupancy_clamped),
        );
    }
    let structure = builder.clone().build()?;
    match coincident::warning(&structure) {
        Some(message) => builder
            .metadata(
                "coincident_atoms_warning",
                serde_json::Value::String(message),
            )
            .build(),
        None => Ok(structure),
    }
}
