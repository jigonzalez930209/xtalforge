//! The global limit profile flag.

use clap::ValueEnum;

#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum Limits {
    Standard,
    Generous,
    WarnOnly,
}

impl Limits {
    pub(crate) fn profile(self) -> xtalforge::core::Profile {
        match self {
            Self::Standard => xtalforge::core::Profile::Standard,
            Self::Generous => xtalforge::core::Profile::Generous,
            Self::WarnOnly => xtalforge::core::Profile::WarnOnly,
        }
    }
}
