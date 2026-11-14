//! CP2K GPW input writer (Task 18.1.4).
//!
//! Basis sets and pseudopotentials use the MOLOPT/GTH names stored in the
//! preset data. The cell is written as `ABC` plus `ALPHA_BETA_GAMMA`, the
//! coordinates as fractional `SCALED` positions, and DFT-D3 is enabled by the
//! presets that request it. An element missing from `potentials` falls back to
//! `GTH-PBE-q<valence>` from the element table.

use serde_json::Value;
use xtalforge_core::{P1Cell, Structure};

use crate::Error;
use crate::input_set::settings::Settings;
use crate::input_set::util;

const INPUT: &str = "cp2k.inp";

pub(crate) fn build(
    structure: &Structure,
    settings: &Settings,
) -> Result<Vec<super::InputFile>, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::Job(
            "CP2K needs a structure with at least one atom".to_owned(),
        ));
    }
    Ok(vec![super::InputFile {
        name: INPUT.to_owned(),
        text: input_text(structure, &p1, settings)?,
    }])
}

fn input_text(structure: &Structure, p1: &P1Cell, settings: &Settings) -> Result<String, Error> {
    let run_type = settings.text("run_type")?;
    if !matches!(run_type, "ENERGY" | "GEO_OPT" | "CELL_OPT") {
        return Err(Error::Job(format!(
            "unknown CP2K run_type {run_type:?}; use ENERGY, GEO_OPT or CELL_OPT"
        )));
    }
    let project = util::identifier(settings.text("project")?, "xtalforge");
    let [a, b, c, alpha, beta, gamma] = structure.lattice().parameters();
    let mut out = String::from("&GLOBAL\n");
    out.push_str(&format!("  PROJECT {project}\n"));
    out.push_str(&format!("  RUN_TYPE {run_type}\n"));
    out.push_str(&format!(
        "  PRINT_LEVEL {}\n&END GLOBAL\n",
        util::token(settings.text("print_level")?, "print level")?
    ));
    out.push_str("&FORCE_EVAL\n  METHOD Quickstep\n  &DFT\n");
    out.push_str(&format!(
        "    BASIS_SET_FILE_NAME {}\n",
        util::token(settings.text("basis_set_file")?, "basis set file")?
    ));
    out.push_str(&format!(
        "    POTENTIAL_FILE_NAME {}\n",
        util::token(settings.text("potential_file")?, "potential file")?
    ));
    out.push_str("    &MGRID\n");
    out.push_str(&format!(
        "      CUTOFF {}\n      NGRIDS {}\n      REL_CUTOFF {}\n    &END MGRID\n",
        float(settings.float("cutoff")?),
        settings.integer("ngrids")?,
        float(settings.float("rel_cutoff")?)
    ));
    out.push_str("    &QS\n");
    out.push_str(&format!(
        "      EPS_DEFAULT {}\n    &END QS\n",
        scientific(settings.float("eps_default")?)
    ));
    out.push_str("    &SCF\n      SCF_GUESS ATOMIC\n");
    out.push_str(&format!(
        "      EPS_SCF {}\n      MAX_SCF {}\n",
        scientific(settings.float("eps_scf")?),
        settings.integer("max_scf")?
    ));
    out.push_str("      &OT\n        MINIMIZER DIIS\n");
    out.push_str("        PRECONDITIONER FULL_SINGLE_INVERSE\n      &END OT\n");
    out.push_str("      &OUTER_SCF\n");
    out.push_str(&format!(
        "        EPS_SCF {}\n        MAX_SCF {}\n      &END OUTER_SCF\n",
        scientific(settings.float("eps_scf")?),
        settings.integer("outer_max_scf")?
    ));
    out.push_str("    &END SCF\n    &XC\n");
    out.push_str(&format!(
        "      &XC_FUNCTIONAL {}\n      &END XC_FUNCTIONAL\n",
        util::token(
            settings.text("reference_functional")?,
            "exchange functional"
        )?
    ));
    if settings.flag("dft_d3")? {
        out.push_str("      &VDW_POTENTIAL\n        POTENTIAL_TYPE PAIR_POTENTIAL\n");
        out.push_str("        &PAIR_POTENTIAL\n          TYPE DFTD3\n");
        out.push_str("          PARAMETER_FILE_NAME dftd3.dat\n");
        out.push_str(&format!(
            "          REFERENCE_FUNCTIONAL {}\n        &END PAIR_POTENTIAL\n",
            util::token(
                settings.text("reference_functional")?,
                "reference functional"
            )?
        ));
        out.push_str("      &END VDW_POTENTIAL\n");
    }
    out.push_str("    &END XC\n  &END DFT\n  &SUBSYS\n    &CELL\n");
    out.push_str(&format!(
        "      ABC {} {} {}\n      ALPHA_BETA_GAMMA {} {} {}\n",
        float(a),
        float(b),
        float(c),
        float(alpha),
        float(beta),
        float(gamma)
    ));
    out.push_str("      PERIODIC XYZ\n    &END CELL\n    &COORD\n      SCALED\n");
    for atom in &p1.atoms {
        out.push_str(&format!(
            "      {} {:.10} {:.10} {:.10}\n",
            atom.element, atom.frac[0], atom.frac[1], atom.frac[2]
        ));
    }
    out.push_str("    &END COORD\n");
    let symbols = util::species(&p1.atoms);
    let basis = util::token(settings.text("basis_set")?, "basis set")?;
    let potentials = settings.table("potentials")?;
    for symbol in &symbols {
        let potential = potentials
            .get(symbol)
            .and_then(Value::as_str)
            .map_or_else(|| fallback_potential(symbol), str::to_owned);
        let potential = util::token(&potential, "GTH potential")?;
        out.push_str(&format!(
            "    &KIND {symbol}\n      BASIS_SET {basis}\n      POTENTIAL {potential}\n    &END KIND\n"
        ));
    }
    out.push_str("  &END SUBSYS\n&END FORCE_EVAL\n");
    if matches!(run_type, "GEO_OPT" | "CELL_OPT") {
        let block = if run_type == "CELL_OPT" {
            "CELL_OPT"
        } else {
            "GEO_OPT"
        };
        out.push_str(&format!("&MOTION\n  &{block}\n"));
        out.push_str(&format!(
            "    OPTIMIZER {}\n    MAX_ITER {}\n  &END {block}\n&END MOTION\n",
            util::token(settings.text("optimizer")?, "optimizer")?,
            settings.integer("max_iter")?
        ));
    }
    Ok(out)
}

fn fallback_potential(symbol: &str) -> String {
    match xtalforge_core::element_info(symbol)
        .ok()
        .and_then(|info| info.valence)
    {
        Some(valence) => format!("GTH-PBE-q{valence}"),
        None => "GTH-PBE".to_owned(),
    }
}

fn float(value: f64) -> String {
    if value == 0.0 {
        return "0.0".to_owned();
    }
    let text = format!("{value:.8}");
    let trimmed = text.trim_end_matches('0');
    let trimmed = trimmed.strip_suffix('.').unwrap_or(trimmed);
    trimmed.to_owned()
}

fn scientific(value: f64) -> String {
    format!("{value:.1E}")
}
