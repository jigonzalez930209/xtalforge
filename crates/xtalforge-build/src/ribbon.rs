//! Periodic honeycomb nanoribbon (Task 6.8.3).
//!
//! The ribbon is a graphene strip with vacuum added along the second lattice
//! vector. It is named by its repeat counts rather than by an armchair or
//! zigzag label: the edge character follows from how the strip is cut, and the
//! builder only promises a honeycomb strip periodic along `a`.

use xtalforge_core::{Error, Lattice, Structure, StructureBuilder};

use crate::carbon::graphene;

/// Builds a periodic honeycomb strip with vacuum on both sides.
///
/// `width` and `length` are graphene repeats along `b` and `a`; the strip has
/// `2 * width * length` carbon atoms. The `b` vector is lengthened by
/// `vacuum` Angstrom while every fractional coordinate is kept, so the strip
/// is unchanged and the empty space around it grows.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for zero width or length, a negative or
/// non-finite vacuum, or a strip whose cell cannot be rebuilt.
pub fn nanoribbon(width: u32, length: u32, vacuum: f64) -> Result<Structure, Error> {
    if width == 0 || length == 0 {
        return Err(Error::InvalidInput(
            "nanoribbon width and length must be positive".to_owned(),
        ));
    }
    if !vacuum.is_finite() || vacuum < 0.0 {
        return Err(Error::InvalidInput(
            "nanoribbon vacuum must be non-negative and finite".to_owned(),
        ));
    }
    let sheet = graphene(length, width, 1, 3.35, 0.0)?;
    let matrix = sheet.lattice().matrix();
    let b = matrix[1];
    let b_length = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
    if !b_length.is_finite() || b_length <= 0.0 {
        return Err(Error::InvalidInput(
            "the nanoribbon second cell vector is degenerate".to_owned(),
        ));
    }
    let scale = (b_length + vacuum) / b_length;
    let enlarged = [b[0] * scale, b[1] * scale, b[2] * scale];
    let lattice = Lattice::new([matrix[0], enlarged, matrix[2]])?;
    let mut builder = StructureBuilder::new("nanoribbon", lattice);
    for site in sheet.sites() {
        builder = builder.site(site.clone());
    }
    builder.build()
}
