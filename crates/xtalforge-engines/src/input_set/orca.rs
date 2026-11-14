//! ORCA simple-input adapter (Task 18.1.6).
//!
//! The text is produced by the existing `xtalforge_ff` ORCA writer, so the
//! preset registry exposes that export without changing a single byte.

use xtalforge_core::Structure;
use xtalforge_ff::OrcaOptions;

use crate::Error;
use crate::input_set::settings::Settings;

const INPUT: &str = "orca.inp";

pub(crate) fn build(
    structure: &Structure,
    settings: &Settings,
) -> Result<Vec<super::InputFile>, Error> {
    let charge = settings.integer("charge")?;
    let charge = i32::try_from(charge)
        .map_err(|_| Error::Job("ORCA charge does not fit a 32-bit integer".to_owned()))?;
    let multiplicity = settings.usize("multiplicity")?;
    let multiplicity = u32::try_from(multiplicity)
        .map_err(|_| Error::Job("ORCA multiplicity does not fit a 32-bit integer".to_owned()))?;
    let options = OrcaOptions {
        method: settings.text("method")?.to_owned(),
        basis: settings.text("basis")?.to_owned(),
        charge,
        multiplicity,
        solvent: settings.optional_text("solvent")?.map(str::to_owned),
    };
    let text = xtalforge_ff::external::write_orca_input(structure, &options)
        .map_err(|error| Error::Job(format!("ORCA input: {error}")))?;
    Ok(vec![super::InputFile {
        name: INPUT.to_owned(),
        text,
    }])
}
