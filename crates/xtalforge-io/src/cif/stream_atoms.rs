//! Atom-site row parsing for the streaming CIF reader (Task 15.2.2).

use serde_json::Value;
use xtalforge_core::{Error, Lattice, Site};

use super::stream::loop_value;
use super::tag_matches;
use super::values::{missing, number};

pub(super) struct AtomColumns {
    pub(super) cartesian: bool,
    label: [Option<usize>; 4],
    element: Option<usize>,
    x: Option<usize>,
    y: Option<usize>,
    z: Option<usize>,
    occupancy: Option<usize>,
    u_iso: Option<usize>,
    b_iso: Option<usize>,
    charge: Option<usize>,
    chain: [Option<usize>; 2],
    residue_name: [Option<usize>; 2],
    residue_number: [Option<usize>; 2],
    disorder_group: Option<usize>,
    disorder_assembly: Option<usize>,
    model: Option<usize>,
}

impl AtomColumns {
    pub(super) fn new(tags: &[String]) -> Self {
        let find = |name: &str| {
            tags.iter()
                .position(|candidate| tag_matches(candidate, name))
        };
        let cartesian = find("_atom_site_fract_x").is_none();
        let (x, y, z) = if cartesian {
            (
                find("_atom_site_Cartn_x"),
                find("_atom_site_Cartn_y"),
                find("_atom_site_Cartn_z"),
            )
        } else {
            (
                find("_atom_site_fract_x"),
                find("_atom_site_fract_y"),
                find("_atom_site_fract_z"),
            )
        };
        Self {
            cartesian,
            label: [
                find("_atom_site_label"),
                find("_atom_site_auth_atom_id"),
                find("_atom_site_label_atom_id"),
                find("_atom_site_type_symbol"),
            ],
            element: find("_atom_site_type_symbol"),
            x,
            y,
            z,
            occupancy: find("_atom_site_occupancy"),
            u_iso: find("_atom_site_U_iso_or_equiv"),
            b_iso: find("_atom_site_B_iso_or_equiv"),
            charge: find("_atom_site_charge"),
            chain: [
                find("_atom_site_auth_asym_id"),
                find("_atom_site_label_asym_id"),
            ],
            residue_name: [
                find("_atom_site_auth_comp_id"),
                find("_atom_site_label_comp_id"),
            ],
            residue_number: [
                find("_atom_site_auth_seq_id"),
                find("_atom_site_label_seq_id"),
            ],
            disorder_group: find("_atom_site_disorder_group"),
            disorder_assembly: find("_atom_site_disorder_assembly"),
            model: find("_atom_site_pdbx_PDB_model_num"),
        }
    }
}

pub(super) fn parse_atom(
    row: &[String],
    row_index: usize,
    columns: &AtomColumns,
    lattice: Option<&Lattice>,
    first_model: &mut Option<String>,
    occupancy_clamped: &mut u32,
) -> Result<Option<Site>, Error> {
    let value = |column: Option<usize>| column.and_then(|index| row.get(index)).map(String::as_str);
    let model = value(columns.model).filter(|model| *model != "." && *model != "?");
    if row_index == 0 {
        *first_model = model.map(str::to_owned);
    } else if let (Some(first), Some(model)) = (first_model.as_deref(), model)
        && first != model
    {
        return Ok(None);
    }
    let label = columns
        .label
        .iter()
        .find_map(|column| value(*column))
        .unwrap_or("X");
    let element = value(columns.element).unwrap_or(label);
    let (x_tag, y_tag, z_tag) = if columns.cartesian {
        (
            "_atom_site_Cartn_x",
            "_atom_site_Cartn_y",
            "_atom_site_Cartn_z",
        )
    } else {
        (
            "_atom_site_fract_x",
            "_atom_site_fract_y",
            "_atom_site_fract_z",
        )
    };
    let x = number(value(columns.x)).ok_or_else(|| missing(x_tag))?;
    let y = number(value(columns.y)).ok_or_else(|| missing(y_tag))?;
    let z = number(value(columns.z)).ok_or_else(|| missing(z_tag))?;
    let frac = if columns.cartesian {
        let lattice = lattice.ok_or_else(|| {
            Error::InvalidInput(
                "the streaming CIF reader needs the cell items before a Cartesian atom loop; use read_cif"
                    .to_owned(),
            )
        })?;
        lattice.to_frac([x, y, z])
    } else {
        [x, y, z]
    };
    let raw_occupancy = number(value(columns.occupancy)).unwrap_or(1.0);
    let occupancy = raw_occupancy.clamp(0.0, 1.0);
    if (occupancy - raw_occupancy).abs() > 1e-9 {
        *occupancy_clamped += 1;
    }
    let mut site = Site::new(label, element, frac)?.with_occupancy(occupancy)?;
    if let Some(u_iso) = number(value(columns.u_iso)) {
        site = site.with_u_iso(u_iso)?;
    } else if let Some(b_iso) = number(value(columns.b_iso)) {
        let u_iso = b_iso / (8.0 * std::f64::consts::PI * std::f64::consts::PI);
        site = site.with_u_iso(u_iso.max(0.0))?;
    }
    if let Some(charge) = number(value(columns.charge)) {
        site = site.with_charge(charge)?;
    }
    for (prop, keys) in [
        ("chain_id", columns.chain),
        ("residue_name", columns.residue_name),
        ("residue_number", columns.residue_number),
    ] {
        if let Some(value) = keys
            .iter()
            .find_map(|column| value(*column))
            .filter(|value| *value != "." && *value != "?")
        {
            site = site.with_prop(prop, Value::String(value.to_owned()));
        }
    }
    if let Some(group) = value(columns.disorder_group) {
        let group = group.trim();
        if group != "." && group != "?" {
            match group.parse::<u32>() {
                Ok(number) => site = site.with_disorder_group(number),
                Err(_) => {
                    site = site.with_prop("disorder_group", Value::String(group.to_owned()));
                }
            }
        }
    }
    if let Some(assembly) = value(columns.disorder_assembly) {
        let assembly = assembly.trim();
        if assembly != "." && assembly != "?" {
            site = site.with_prop("disorder_assembly", Value::String(assembly.to_owned()));
        }
    }
    Ok(Some(site))
}

pub(super) fn apply_aniso(
    sites: &mut [Site],
    tags: &[String],
    rows: &[Vec<String>],
) -> Result<(), Error> {
    if tags.is_empty() || rows.is_empty() {
        return Ok(());
    }
    let has_label = tags
        .iter()
        .any(|tag| tag_matches(tag, "_atom_site_aniso_label"));
    for (row_index, row) in rows.iter().enumerate() {
        let target = if has_label {
            loop_value(tags, row, "_atom_site_aniso_label")
                .and_then(|label| sites.iter().position(|site| site.label() == label))
        } else if row_index < sites.len() {
            Some(row_index)
        } else {
            None
        };
        let Some(target) = target else {
            continue;
        };
        let values = [
            number(loop_value(tags, row, "_atom_site_aniso_U_11")),
            number(loop_value(tags, row, "_atom_site_aniso_U_22")),
            number(loop_value(tags, row, "_atom_site_aniso_U_33")),
            number(loop_value(tags, row, "_atom_site_aniso_U_12")),
            number(loop_value(tags, row, "_atom_site_aniso_U_13")),
            number(loop_value(tags, row, "_atom_site_aniso_U_23")),
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
    Ok(())
}
