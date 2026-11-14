//! Quantum ESPRESSO `pw.x` input writer (Task 18.1.3).
//!
//! Pseudopotential names are SSSP names only. `bands` and `nscf` use the
//! symmetry band path (Task 2.7) through `K_POINTS crystal_b`; every other
//! preset uses an automatic mesh from a target spacing. Elements missing from
//! the `pseudopotentials` map fall back to `<symbol>.UPF`; override that map
//! to use another SSSP name.

use xtalforge_core::{P1Cell, Structure};
use xtalforge_symmetry::SymmetryOptions;

use crate::Error;
use crate::input_set::settings::Settings;
use crate::input_set::util;

const INPUT: &str = "pw.in";

pub(crate) fn build(
    structure: &Structure,
    settings: &Settings,
) -> Result<Vec<super::InputFile>, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::Job(
            "Quantum ESPRESSO needs a structure with at least one atom".to_owned(),
        ));
    }
    Ok(vec![super::InputFile {
        name: INPUT.to_owned(),
        text: input_text(structure, &p1, settings)?,
    }])
}

fn input_text(structure: &Structure, p1: &P1Cell, settings: &Settings) -> Result<String, Error> {
    let calculation = settings.text("calculation")?;
    let valid = matches!(calculation, "scf" | "nscf" | "bands" | "relax" | "vc-relax");
    if !valid {
        return Err(Error::Job(format!(
            "unknown pw.x calculation {calculation:?}"
        )));
    }
    let prefix = util::identifier(settings.text("prefix")?, "xtalforge");
    let pseudo_dir = util::relative_path(settings.text("pseudo_dir")?, "pseudo_dir")?;
    let outdir = util::relative_path(settings.text("outdir")?, "outdir")?;
    let symbols = util::species(&p1.atoms);
    let pseudos = settings.table("pseudopotentials")?;
    let mut out = String::from("&CONTROL\n");
    out.push_str(&format!("  calculation = '{calculation}'\n"));
    out.push_str(&format!("  prefix = '{prefix}'\n"));
    out.push_str(&format!("  pseudo_dir = '{pseudo_dir}'\n"));
    out.push_str(&format!("  outdir = '{outdir}'\n"));
    out.push_str(&format!(
        "  tprnfor = {}\n  tstress = {}\n",
        fortran(settings.flag("tprnfor")?),
        fortran(settings.flag("tstress")?)
    ));
    out.push_str("/\n&SYSTEM\n  ibrav = 0\n");
    out.push_str(&format!("  nat = {}\n", p1.atoms.len()));
    out.push_str(&format!("  ntyp = {}\n", symbols.len()));
    out.push_str(&format!(
        "  ecutwfc = {}\n",
        float(settings.float("ecutwfc")?)
    ));
    out.push_str(&format!(
        "  ecutrho = {}\n",
        float(settings.float("ecutrho")?)
    ));
    let occupations = settings.text("occupations")?;
    out.push_str(&format!(
        "  occupations = '{}'\n",
        util::token(occupations, "occupations")?
    ));
    if occupations == "smearing" {
        out.push_str(&format!(
            "  smearing = '{}'\n",
            util::token(settings.text("smearing")?, "smearing")?
        ));
        out.push_str(&format!(
            "  degauss = {}\n",
            float(settings.float("degauss")?)
        ));
    }
    if let Some(nbnd) = settings.raw("nbnd").and_then(|value| value.as_i64())
        && nbnd > 0
    {
        out.push_str(&format!("  nbnd = {nbnd}\n"));
    }
    out.push_str("/\n&ELECTRONS\n");
    out.push_str(&format!(
        "  conv_thr = {}\n  mixing_beta = {}\n/\n",
        float(settings.float("conv_thr")?),
        float(settings.float("mixing_beta")?)
    ));
    if matches!(calculation, "relax" | "vc-relax") {
        out.push_str("&IONS\n");
        out.push_str(&format!(
            "  ion_dynamics = '{}'\n",
            util::token(settings.text("ion_dynamics")?, "ion_dynamics")?
        ));
        out.push_str("/\n");
    }
    if calculation == "vc-relax" {
        out.push_str("&CELL\n");
        out.push_str(&format!(
            "  cell_dynamics = '{}'\n",
            util::token(settings.text("cell_dynamics")?, "cell_dynamics")?
        ));
        out.push_str("/\n");
    }
    out.push_str("ATOMIC_SPECIES\n");
    for symbol in &symbols {
        let mass = xtalforge_core::elements::element_by_symbol(symbol)?.mass;
        let pseudo = pseudos
            .get(symbol)
            .and_then(|value| value.as_str())
            .map_or_else(|| format!("{symbol}.UPF"), str::to_owned);
        let pseudo = util::token(&pseudo, "pseudopotential")?;
        out.push_str(&format!("  {symbol} {mass:.4} {pseudo}\n"));
    }
    out.push_str("CELL_PARAMETERS angstrom\n");
    for row in structure.lattice().matrix() {
        for value in row {
            out.push_str(&format!(" {value:.10}"));
        }
        out.push('\n');
    }
    out.push_str("ATOMIC_POSITIONS crystal\n");
    for atom in &p1.atoms {
        out.push_str(&format!(
            "  {} {:.10} {:.10} {:.10}\n",
            atom.element, atom.frac[0], atom.frac[1], atom.frac[2]
        ));
    }
    out.push_str(&kpoints_text(structure, settings)?);
    Ok(out)
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
            Ok(format!(
                "K_POINTS automatic\n  {} {} {} 0 0 0\n",
                grid[0], grid[1], grid[2]
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
            Ok(format!(
                "K_POINTS automatic\n  {} {} {} 0 0 0\n",
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
    let path =
        match xtalforge_symmetry::band_path_for_structure(structure, SymmetryOptions::default()) {
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
    let mut vertices: Vec<(&str, [f64; 3])> = Vec::new();
    for (start, end) in &path.path {
        let first = path.special_points.get(start).ok_or_else(|| {
            Error::Job(format!("band path point {start} is missing from the table"))
        })?;
        let second = path.special_points.get(end).ok_or_else(|| {
            Error::Job(format!("band path point {end} is missing from the table"))
        })?;
        if vertices.last().map(|(_, point)| point) != Some(first) {
            vertices.push((start.as_str(), *first));
        }
        vertices.push((end.as_str(), *second));
    }
    let mut out = format!("K_POINTS crystal_b\n  {}\n", vertices.len());
    for (label, point) in vertices {
        out.push_str(&format!(
            "  {} {} {} {} ! {label}\n",
            point[0], point[1], point[2], points
        ));
    }
    Ok(out)
}

fn fortran(value: bool) -> &'static str {
    if value { ".true." } else { ".false." }
}

fn float(value: f64) -> String {
    if value == 0.0 {
        return "0.0".to_owned();
    }
    if value.abs() < 1e-8 || value.abs() >= 1e8 {
        return format!("{value:.6e}");
    }
    let text = format!("{value:.10}");
    let trimmed = text.trim_end_matches('0');
    let trimmed = trimmed.strip_suffix('.').unwrap_or(trimmed);
    trimmed.to_owned()
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
