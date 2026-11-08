//! Element-specific UFF typing rules.
//!
//! Ported from Crystal-Builder `xtal/ff/uff/typer.py`. These are the
//! hand-written rules for the common elements; everything else goes through
//! `crate::uff::coordination::by_coordination`, which reads the shape a UFF
//! type name declares. Every rule receives the periodic bond geometry and
//! whether the atom belongs to a planar aromatic ring.

use crate::Error;
use crate::uff::geometry::{AtomView, Geometry};
use crate::uff::typer::{AtomType, Confidence};

const SIMPLE_TYPES: &[(&str, &str)] = &[
    ("Rb", "Rb"),
    ("Cs", "Cs"),
    ("Fr", "Fr"),
    ("F", "F_"),
    ("Cl", "Cl"),
    ("Br", "Br"),
    ("I", "I_"),
    ("At", "At"),
    ("Si", "Si3"),
    ("Ge", "Ge3"),
    ("Sn", "Sn3"),
    ("He", "He4+4"),
    ("Ne", "Ne4+4"),
    ("Ar", "Ar4+4"),
    ("Kr", "Kr4+4"),
    ("Xe", "Xe4+4"),
    ("Rn", "Rn4+4"),
];

const ZEOLITE_FORMERS: &[&str] = &["Si", "Al", "P", "Ge", "B"];
const SQUARE_PLANAR_ANGLE: f64 = 145.0;

pub(crate) fn simple_type(element: &str) -> Option<&'static str> {
    SIMPLE_TYPES
        .iter()
        .find(|(symbol, _)| *symbol == element)
        .map(|(_, name)| *name)
}

fn is_metal(element: &str) -> bool {
    xtalforge_core::elements::element_by_symbol(element).is_ok_and(|value| value.is_metal)
}

pub(crate) fn framework_node(view: &AtomView, geometry: &Geometry, index: usize) -> bool {
    let Some(element) = view.elements.get(index) else {
        return false;
    };
    if !is_metal(element) {
        return false;
    }
    let partners = geometry.partners(index);
    if partners.len() < 2 {
        return false;
    }
    if partners.iter().any(|partner| {
        view.elements
            .get(partner.index)
            .is_some_and(|value| is_metal(value))
    }) {
        return false;
    }
    for partner in partners {
        let Some(symbol) = view.elements.get(partner.index) else {
            continue;
        };
        if !matches!(symbol.as_str(), "O" | "N") {
            continue;
        }
        if geometry.partners(partner.index).iter().any(|next| {
            view.elements
                .get(next.index)
                .is_some_and(|value| value == "C")
        }) {
            return true;
        }
    }
    false
}

fn bridges_framework_metals(view: &AtomView, geometry: &Geometry, index: usize) -> bool {
    let partners = geometry.partners(index);
    if partners.is_empty() {
        return false;
    }
    if !partners.iter().all(|partner| {
        view.elements
            .get(partner.index)
            .is_some_and(|value| is_metal(value))
    }) {
        return false;
    }
    partners
        .iter()
        .any(|partner| framework_node(view, geometry, partner.index))
}

pub(crate) fn shape_character(geometry: &Geometry, view: &AtomView, index: usize) -> &'static str {
    if geometry.coordination(index) != 4 {
        return "";
    }
    if geometry.max_angle(view, index) > SQUARE_PLANAR_ANGLE {
        "4"
    } else {
        "3"
    }
}

type Rule = fn(usize, &AtomView, &Geometry, bool) -> Result<AtomType, Error>;

fn hydrogen(
    index: usize,
    view: &AtomView,
    geometry: &Geometry,
    _: bool,
) -> Result<AtomType, Error> {
    let partners = geometry.partners(index);
    if partners.len() == 2
        && partners.iter().all(|partner| {
            view.elements
                .get(partner.index)
                .is_some_and(|value| value == "B")
        })
    {
        return Ok(AtomType::new(
            "H_b",
            Confidence::Certain,
            "bridging two borons",
        ));
    }
    if partners.len() > 1 {
        return Ok(AtomType::new(
            "H_",
            Confidence::Uncertain,
            format!(
                "{} neighbours; hydrogen should have one, so the bond criteria may be too generous",
                partners.len()
            ),
        ));
    }
    Ok(AtomType::new("H_", Confidence::Certain, ""))
}

fn boron(index: usize, _: &AtomView, geometry: &Geometry, _: bool) -> Result<AtomType, Error> {
    let count = geometry.coordination(index);
    if count >= 4 {
        Ok(AtomType::new("B_3", Confidence::Certain, "four neighbours"))
    } else {
        Ok(AtomType::new(
            "B_2",
            Confidence::Likely,
            format!("{count} neighbours, taken as trigonal"),
        ))
    }
}

fn carbon(
    index: usize,
    view: &AtomView,
    geometry: &Geometry,
    aromatic: bool,
) -> Result<AtomType, Error> {
    let count = geometry.coordination(index);
    if count >= 4 {
        return Ok(AtomType::new("C_3", Confidence::Certain, "four neighbours"));
    }
    if aromatic {
        return Ok(AtomType::new(
            "C_R",
            Confidence::Certain,
            "in a flat aromatic ring",
        ));
    }
    if count == 3 {
        let total = geometry.angle_sum(view, index);
        if total >= 345.0 {
            return Ok(AtomType::new(
                "C_2",
                Confidence::Certain,
                format!("planar, angles sum to {total:.0}"),
            ));
        }
        return Ok(AtomType::new(
            "C_3",
            Confidence::Likely,
            format!("pyramidal, angles sum to {total:.0}; a hydrogen is probably missing"),
        ));
    }
    if count == 2 {
        let angle = geometry.max_angle(view, index);
        if angle > 155.0 {
            return Ok(AtomType::new(
                "C_1",
                Confidence::Certain,
                format!("linear, {angle:.0}"),
            ));
        }
        if angle >= 114.0 {
            return Ok(AtomType::new(
                "C_2",
                Confidence::Likely,
                format!("{angle:.0} between its two neighbours"),
            ));
        }
        return Ok(AtomType::new(
            "C_3",
            Confidence::Likely,
            format!("{angle:.0} between its two neighbours; hydrogens are probably missing"),
        ));
    }
    Ok(AtomType::new(
        "C_3",
        Confidence::Uncertain,
        format!("only {count} neighbour(s); nothing in the geometry says what this carbon is"),
    ))
}

fn nitrogen(
    index: usize,
    view: &AtomView,
    geometry: &Geometry,
    aromatic: bool,
) -> Result<AtomType, Error> {
    let count = geometry.coordination(index);
    if count >= 4 {
        return Ok(AtomType::new("N_3", Confidence::Certain, "four neighbours"));
    }
    if aromatic {
        return Ok(AtomType::new(
            "N_R",
            Confidence::Certain,
            "in a flat aromatic ring",
        ));
    }
    if count == 3 {
        let total = geometry.angle_sum(view, index);
        if total >= 345.0 {
            return Ok(AtomType::new(
                "N_2",
                Confidence::Certain,
                format!("planar, angles sum to {total:.0}"),
            ));
        }
        return Ok(AtomType::new(
            "N_3",
            Confidence::Certain,
            format!("pyramidal, angles sum to {total:.0}"),
        ));
    }
    if count == 2 {
        let angle = geometry.max_angle(view, index);
        if angle > 155.0 {
            return Ok(AtomType::new(
                "N_1",
                Confidence::Certain,
                format!("linear, {angle:.0}"),
            ));
        }
        return Ok(AtomType::new(
            "N_2",
            Confidence::Likely,
            format!("bent at {angle:.0}, taken as sp2"),
        ));
    }
    if count == 1 {
        let distance = geometry
            .distances(view, index)
            .first()
            .copied()
            .unwrap_or(0.0);
        if distance < 1.20 {
            return Ok(AtomType::new(
                "N_1",
                Confidence::Likely,
                format!("one neighbour at {distance:.2} A -- a triple bond"),
            ));
        }
        return Ok(AtomType::new(
            "N_2",
            Confidence::Likely,
            format!("one neighbour at {distance:.2} A"),
        ));
    }
    Ok(AtomType::new("N_3", Confidence::Uncertain, "no neighbours"))
}

fn oxygen(
    index: usize,
    view: &AtomView,
    geometry: &Geometry,
    aromatic: bool,
) -> Result<AtomType, Error> {
    let count = geometry.coordination(index);
    if aromatic {
        return Ok(AtomType::new(
            "O_R",
            Confidence::Certain,
            "in a flat aromatic ring",
        ));
    }
    if count >= 3 {
        if bridges_framework_metals(view, geometry, index) {
            return Ok(AtomType::new(
                "O_3_f",
                Confidence::Certain,
                format!("bridges {count} framework metals -- the oxide at the centre of a node"),
            ));
        }
        return Ok(AtomType::new(
            "O_3",
            Confidence::Likely,
            format!("{count} neighbours; bridging or over-bonded"),
        ));
    }
    if count == 2 {
        let elements: Vec<&str> = geometry
            .partners(index)
            .iter()
            .filter_map(|partner| view.elements.get(partner.index).map(String::as_str))
            .collect();
        let angle = geometry.max_angle(view, index);
        if elements.len() == 2
            && elements
                .iter()
                .all(|element| ZEOLITE_FORMERS.contains(element))
            && angle >= 130.0
        {
            return Ok(AtomType::new(
                "O_3_z",
                Confidence::Certain,
                format!(
                    "bridges {} and {} at {angle:.0} -- a framework oxygen",
                    elements.first().copied().unwrap_or("?"),
                    elements.get(1).copied().unwrap_or("?")
                ),
            ));
        }
        // UFF4MOF's O_2_z is reachable by override only: the reference
        // measured that applying it to carboxylate oxygens relaxes MOF-5
        // further from experiment, so the code does not claim that
        // environment on its own.
        return Ok(AtomType::new(
            "O_3",
            Confidence::Certain,
            format!("two neighbours at {angle:.0}"),
        ));
    }
    if count == 1 {
        let distance = geometry
            .distances(view, index)
            .first()
            .copied()
            .unwrap_or(0.0);
        let partner = geometry
            .partners(index)
            .first()
            .and_then(|neighbor| view.elements.get(neighbor.index))
            .map_or("?", String::as_str);
        if distance < 1.30 {
            return Ok(AtomType::new(
                "O_2",
                Confidence::Certain,
                format!("one short bond to {partner} at {distance:.2} A -- a double bond"),
            ));
        }
        return Ok(AtomType::new(
            "O_3",
            Confidence::Likely,
            format!("one bond to {partner} at {distance:.2} A, too long to be a double bond"),
        ));
    }
    Ok(AtomType::new("O_3", Confidence::Uncertain, "no neighbours"))
}

fn phosphorus(index: usize, _: &AtomView, geometry: &Geometry, _: bool) -> Result<AtomType, Error> {
    let count = geometry.coordination(index);
    if count >= 4 {
        Ok(AtomType::new(
            "P_3+5",
            Confidence::Certain,
            format!("{count} neighbours"),
        ))
    } else {
        Ok(AtomType::new(
            "P_3+3",
            Confidence::Likely,
            format!("{count} neighbours"),
        ))
    }
}

fn sulfur(
    index: usize,
    view: &AtomView,
    geometry: &Geometry,
    aromatic: bool,
) -> Result<AtomType, Error> {
    let count = geometry.coordination(index);
    if aromatic {
        return Ok(AtomType::new(
            "S_R",
            Confidence::Certain,
            "in a flat aromatic ring",
        ));
    }
    if count >= 4 {
        return Ok(AtomType::new(
            "S_3+6",
            Confidence::Certain,
            format!("{count} neighbours"),
        ));
    }
    if count == 3 {
        return Ok(AtomType::new(
            "S_3+4",
            Confidence::Certain,
            "three neighbours",
        ));
    }
    if count == 2 {
        return Ok(AtomType::new(
            "S_3+2",
            Confidence::Certain,
            "two neighbours",
        ));
    }
    if count == 1 {
        let distance = geometry
            .distances(view, index)
            .first()
            .copied()
            .unwrap_or(0.0);
        if distance < 1.70 {
            return Ok(AtomType::new(
                "S_2",
                Confidence::Certain,
                format!("one short bond at {distance:.2} A"),
            ));
        }
        return Ok(AtomType::new(
            "S_3+2",
            Confidence::Likely,
            format!("one bond at {distance:.2} A"),
        ));
    }
    Ok(AtomType::new(
        "S_3+2",
        Confidence::Uncertain,
        "no neighbours",
    ))
}

pub(crate) fn rule_for(element: &str) -> Option<Rule> {
    match element {
        "H" => Some(hydrogen),
        "B" => Some(boron),
        "C" => Some(carbon),
        "N" => Some(nitrogen),
        "O" => Some(oxygen),
        "P" => Some(phosphorus),
        "S" => Some(sulfur),
        _ => None,
    }
}
