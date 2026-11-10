//! Picture verb.

use serde_json::Value;

use crate::answers::Answer;
use crate::render::{RenderOptions, render_png, render_svg};
use crate::run::{document_argument, integer_or, optional_string};
use crate::{Error, Session};

pub(crate) fn render(session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let name = session.resolve(document_argument(arguments)?)?;
    let format = optional_string(arguments, "format")?
        .unwrap_or("svg")
        .trim()
        .to_ascii_lowercase();
    let width = integer_or(arguments, "width", 800)?;
    let height = integer_or(arguments, "height", 600)?;
    let width = u32::try_from(width).map_err(|_| Error::Argument {
        name: "width".to_owned(),
        reason: "does not fit in 32 bits".to_owned(),
    })?;
    let height = u32::try_from(height).map_err(|_| Error::Argument {
        name: "height".to_owned(),
        reason: "does not fit in 32 bits".to_owned(),
    })?;
    let path = optional_string(arguments, "path")?.map(str::to_owned);
    let options = RenderOptions {
        width,
        height,
        ..Default::default()
    };
    let document = session.document(Some(&name))?;
    {
        let n_atoms = document.structure().p1()?.atoms.len() as u64;
        let verdict = xtalforge_core::limits::check_drawn(
            &session.options().limits,
            n_atoms,
            [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)],
        );
        crate::ops::enforce(verdict)?;
    }
    let bytes = match format.as_str() {
        "svg" => render_svg(document.structure(), &options)?.into_bytes(),
        "png" => render_png(document.structure(), &options)?,
        other => {
            return Err(Error::Argument {
                name: "format".to_owned(),
                reason: format!("unknown picture format {other:?}; use svg or png"),
            });
        }
    };
    if let Some(path) = &path {
        std::fs::write(path, &bytes).map_err(|error| Error::Io(format!("{path}: {error}")))?;
    }
    let written = path.is_some();
    let encoded = bytes.len();
    Ok(Answer::Picture {
        document: name,
        format,
        width,
        height,
        path,
        bytes: if written { None } else { Some(encoded) },
    })
}
