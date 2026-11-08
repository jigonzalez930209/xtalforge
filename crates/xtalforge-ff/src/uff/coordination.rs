//! The general UFF typing rule: let the type table decide from coordination.
//!
//! Ported from Crystal-Builder `xtal/ff/uff/typer.py`. A UFF type name states
//! the coordination it was fitted for, so choosing the row whose geometry
//! matches what the atom actually has needs no per-element code. UFF4MOF's
//! fitted rows are only reachable for a framework node.

use crate::Error;
use crate::uff::geometry::{AtomView, Geometry};
use crate::uff::params::{ParameterSet, UffParams, types_for_set};
use crate::uff::rules::{framework_node, shape_character};
use crate::uff::typer::{AtomType, Confidence};

pub(crate) fn by_coordination(
    index: usize,
    view: &AtomView,
    geometry: &Geometry,
    parameter_set: ParameterSet,
) -> Result<AtomType, Error> {
    let element = view
        .elements
        .get(index)
        .ok_or_else(|| Error::Typing(format!("atom index {index} is out of range")))?;
    let count = geometry.coordination(index);
    let mut candidates = types_for_set(element, parameter_set)?;
    if candidates.is_empty() {
        return Err(Error::UnknownElement(element.clone()));
    }
    if parameter_set == ParameterSet::Uff4Mof && !framework_node(view, geometry, index) {
        let plain: Vec<&UffParams> = candidates
            .iter()
            .copied()
            .filter(|params| !params.fitted)
            .collect();
        if !plain.is_empty() {
            candidates = plain;
        }
    }
    if candidates.len() == 1 {
        let only = candidates
            .first()
            .copied()
            .ok_or_else(|| Error::Typing("empty candidate list".to_owned()))?;
        if only.coordination().is_none() || only.coordination() == Some(count) || count == 0 {
            return Ok(AtomType::new(
                only.name.as_str(),
                Confidence::Certain,
                "the only type UFF gives this element",
            ));
        }
        return Ok(AtomType::new(
            only.name.as_str(),
            Confidence::Uncertain,
            format!(
                "the only type UFF gives this element, fitted for {}-coordinate and this one has {count}",
                only.coordination().unwrap_or(0)
            ),
        ));
    }
    let exact: Vec<&UffParams> = candidates
        .iter()
        .copied()
        .filter(|params| params.coordination() == Some(count))
        .collect();
    if !exact.is_empty() {
        let shape = shape_character(geometry, view, index);
        let mut chosen = exact;
        if !shape.is_empty() {
            let shaped: Vec<&UffParams> = chosen
                .iter()
                .copied()
                .filter(|params| params.geometry == shape)
                .collect();
            if !shaped.is_empty() {
                chosen = shaped;
            }
        }
        if chosen.len() > 1 {
            let fitted: Vec<&UffParams> = chosen
                .iter()
                .copied()
                .filter(|params| params.fitted)
                .collect();
            if !fitted.is_empty() {
                chosen = fitted;
            }
        }
        let best = chosen
            .first()
            .copied()
            .ok_or_else(|| Error::Typing("empty exact candidate list".to_owned()))?;
        let reason = format!(
            "{count} neighbours{}{}",
            if shape.is_empty() {
                String::new()
            } else {
                format!(", {shape}")
            },
            if best.fitted {
                ", in a framework node"
            } else {
                ""
            }
        );
        if chosen.len() == 1 {
            return Ok(AtomType::new(
                best.name.as_str(),
                Confidence::Likely,
                reason,
            ));
        }
        let rest = chosen
            .iter()
            .skip(1)
            .map(|params| params.name.as_str())
            .collect::<Vec<_>>()
            .join(", ");
        return Ok(AtomType::new(
            best.name.as_str(),
            Confidence::Uncertain,
            format!(
                "{count} neighbours, and {element} has more than one type fitted for that: {} was taken over {rest}",
                best.name
            ),
        ));
    }
    let best = candidates
        .iter()
        .copied()
        .min_by_key(|params| (params.coordination().unwrap_or(0) as isize - count as isize).abs())
        .ok_or_else(|| Error::UnknownElement(element.clone()))?;
    Ok(AtomType::new(
        best.name.as_str(),
        Confidence::Uncertain,
        format!(
            "{count} neighbours, which no {element} type in UFF was fitted for; {} ({}-coordinate) is the closest",
            best.name,
            best.coordination().unwrap_or(0)
        ),
    ))
}
