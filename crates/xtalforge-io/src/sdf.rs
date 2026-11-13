//! SDF and MDL MOL reader and writer (Task 13.2.5).
//!
//! V2000 is read and written. V3000 is read as a documented subset: the
//! `COUNTS`, `ATOM` and `BOND` blocks, atom-map numbers and `CHG=` attributes.
//! Files without a cell become a molecule in a padded P1 box. Data fields
//! after `M  END` are kept in the `sdf_data` metadata object and written back.

mod read;
mod write;

pub use read::{read_sdf, read_sdf_all};
pub use write::write_sdf;
