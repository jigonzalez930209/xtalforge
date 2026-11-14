//! VASP INCAR, KPOINTS, POSCAR and POTCAR-spec writer (Task 18.1.2).
//!
//! The POSCAR is the P1 expansion. MAGMOM and LDAU are read from the site
//! properties `magmom`, `ldau_u`, `ldau_orbital` and `ldau_j`, falling back to
//! the same keys in the structure metadata and finally to the preset settings
//! (`magmom`, `ldau_u`, `ldau_l`, `ldau_j`). The fallback keeps the input
//! valid when no magnetic or Hubbard data exists; it never invents licensed
//! POTCAR content. Elements missing from the POTCAR label table fall back to
//! their symbol, which the caller can override with `potcar.labels`.

use serde_json::Value;
use xtalforge_core::{P1Cell, Structure};

use crate::Error;
use crate::input_set::settings::Settings;
use crate::input_set::util;

mod incar;

use incar::incar_text;

const INCAR: &str = "INCAR";
const KPOINTS: &str = "KPOINTS";
const POSCAR: &str = "POSCAR";
const POTCAR: &str = "POTCAR.spec";

pub(crate) fn build(
    structure: &Structure,
    settings: &Settings,
) -> Result<Vec<super::InputFile>, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::Job(
            "VASP needs a structure with at least one atom".to_owned(),
        ));
    }
    let symbols = util::species(&p1.atoms);
    let files = vec![
        super::InputFile {
            name: INCAR.to_owned(),
            text: incar_text(structure, &p1, settings)?,
        },
        super::InputFile {
            name: KPOINTS.to_owned(),
            text: kpoints_text(structure, settings)?,
        },
        super::InputFile {
            name: POSCAR.to_owned(),
            text: poscar_text(structure, &p1, &symbols, settings)?,
        },
        super::InputFile {
            name: POTCAR.to_owned(),
            text: potcar_text(&symbols, settings)?,
        },
    ];
    Ok(files)
}

fn kpoints_text(structure: &Structure, settings: &Settings) -> Result<String, Error> {
    let mode = settings.text("kpoints.mode")?;
    match mode {
        "density" => {
            let spacing = settings.float("kpoints.spacing")?;
            if !spacing.is_finite() || spacing <= 0.0 {
                return Err(Error::Job(
                    "kpoints.spacing must be positive and finite".to_owned(),
                ));
            }
            let reciprocal = structure.lattice().reciprocal(true)?;
            let [row_a, row_b, row_c] = reciprocal.matrix();
            let grid = [
                subdivisions(norm(row_a), spacing),
                subdivisions(norm(row_b), spacing),
                subdivisions(norm(row_c), spacing),
            ];
            let gamma = settings.flag("kpoints.gamma")?;
            let shift = if gamma { "0 0 0" } else { "0.5 0.5 0.5" };
            let scheme = if gamma { "Gamma" } else { "Monkhorst-Pack" };
            Ok(format!(
                "Automatic mesh from a target spacing of {} 1/Angstrom\n0\n{scheme}\n{} {} {}\n{shift}\n",
                vasp_float(spacing),
                grid[0],
                grid[1],
                grid[2]
            ))
        }
        "mesh" => {
            let mesh: [i64; 3] = settings
                .integer_array("kpoints.mesh")?
                .try_into()
                .map_err(|_| Error::Job("kpoints.mesh must have three entries".to_owned()))?;
            if mesh.iter().any(|value| *value < 1) {
                return Err(Error::Job(
                    "kpoints.mesh entries must be at least one".to_owned(),
                ));
            }
            let gamma = settings.flag("kpoints.gamma")?;
            let shift = if gamma { "0 0 0" } else { "0.5 0.5 0.5" };
            let scheme = if gamma { "Gamma" } else { "Monkhorst-Pack" };
            Ok(format!(
                "Explicit k-mesh\n0\n{scheme}\n{} {} {}\n{shift}\n",
                mesh[0], mesh[1], mesh[2]
            ))
        }
        "path" => band_kpoints(structure, settings),
        other => Err(Error::Job(format!(
            "unknown kpoints.mode {other:?}; use density, mesh or path"
        ))),
    }
}

fn band_kpoints(structure: &Structure, settings: &Settings) -> Result<String, Error> {
    let path = match xtalforge_symmetry::band_path_for_structure(
        structure,
        xtalforge_symmetry::SymmetryOptions::default(),
    ) {
        Ok(path) => path,
        Err(_) => xtalforge_symmetry::band_path(structure.lattice())?,
    };
    if path.path.is_empty() {
        return Err(Error::Job("the band path has no segments".to_owned()));
    }
    let points = settings.usize("kpoints.points")?;
    if points == 0 {
        return Err(Error::Job("kpoints.points must be at least one".to_owned()));
    }
    let mut out = format!(
        "Band path (Bravais {})\n{points}\nLine-mode\nReciprocal\n",
        path.bravais
    );
    for (start, end) in &path.path {
        let first = path.special_points.get(start).ok_or_else(|| {
            Error::Job(format!("band path point {start} is missing from the table"))
        })?;
        let second = path.special_points.get(end).ok_or_else(|| {
            Error::Job(format!("band path point {end} is missing from the table"))
        })?;
        out.push_str(&format!(
            "{} {} {} 1 ! {start}\n",
            point_coordinate(*first, 0),
            point_coordinate(*first, 1),
            point_coordinate(*first, 2)
        ));
        out.push_str(&format!(
            "{} {} {} 1 ! {end}\n",
            point_coordinate(*second, 0),
            point_coordinate(*second, 1),
            point_coordinate(*second, 2)
        ));
    }
    Ok(out)
}

fn poscar_text(
    structure: &Structure,
    p1: &P1Cell,
    symbols: &[String],
    settings: &Settings,
) -> Result<String, Error> {
    let mut counts: Vec<usize> = symbols.iter().map(|_| 0).collect();
    for atom in &p1.atoms {
        let index = symbols
            .iter()
            .position(|symbol| symbol == &atom.element)
            .ok_or_else(|| Error::Job(format!("unknown species {}", atom.element)))?;
        if let Some(count) = counts.get_mut(index) {
            *count += 1;
        }
    }
    let mut out = String::new();
    out.push_str(&clean_line(structure.name(), "xtalforge"));
    out.push('\n');
    out.push_str("1.0\n");
    for row in structure.lattice().matrix() {
        for value in row {
            out.push_str(&format!("{value:.10} "));
        }
        out.pop();
        out.push('\n');
    }
    out.push_str(&symbols.join(" "));
    out.push('\n');
    out.push_str(
        &counts
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(" "),
    );
    out.push('\n');
    let selective = settings.flag("selective")?;
    if selective {
        out.push_str("Selective dynamics\n");
    }
    out.push_str("Direct\n");
    for atom in &p1.atoms {
        out.push_str(&format!(
            "{:.10} {:.10} {:.10}",
            atom.frac[0], atom.frac[1], atom.frac[2]
        ));
        if selective {
            out.push_str(" T T T");
        }
        out.push('\n');
    }
    Ok(out)
}

fn potcar_text(symbols: &[String], settings: &Settings) -> Result<String, Error> {
    let functional = util::token(settings.text("potcar.functional")?, "POTCAR functional")?;
    let labels = settings.table("potcar.labels")?;
    let mut out = String::from("# POTCAR names only; provide the licensed PAW files separately.\n");
    for symbol in symbols {
        let label = labels
            .get(symbol)
            .and_then(Value::as_str)
            .unwrap_or(symbol.as_str());
        let label = util::token(label, "POTCAR label")?;
        out.push_str(&format!("PAW_{functional} {label}\n"));
    }
    Ok(out)
}

pub(super) fn vasp_float(value: f64) -> String {
    if value == 0.0 {
        return "0".to_owned();
    }
    if value.abs() < 1e-8 || value.abs() >= 1e8 {
        return format!("{value:.6e}");
    }
    let text = format!("{value:.10}");
    let trimmed = text.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() {
        "0".to_owned()
    } else {
        trimmed.to_owned()
    }
}

fn clean_line(value: &str, fallback: &str) -> String {
    let cleaned: String = value
        .chars()
        .map(|character| {
            if character.is_control() {
                ' '
            } else if matches!(character, '/' | '\\') {
                '_'
            } else {
                character
            }
        })
        .collect();
    let cleaned = cleaned.trim();
    if cleaned.is_empty() {
        fallback.to_owned()
    } else {
        cleaned.to_owned()
    }
}

fn norm(vector: [f64; 3]) -> f64 {
    let [x, y, z] = vector;
    (x * x + y * y + z * z).sqrt()
}

fn subdivisions(length: f64, spacing: f64) -> i64 {
    let ratio = length / spacing;
    if !ratio.is_finite() || ratio <= 0.0 {
        return 1;
    }
    let rounded = (ratio - 1e-12).ceil();
    if rounded < 1.0 { 1 } else { rounded as i64 }
}

fn point_coordinate(point: [f64; 3], index: usize) -> f64 {
    point.get(index).copied().unwrap_or(0.0)
}
