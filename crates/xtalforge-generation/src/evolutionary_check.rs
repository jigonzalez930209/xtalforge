//! Option validation for the evolutionary search.

use xtalforge_core::Error;

use crate::evolutionary::EvolutionaryOptions;

pub(crate) fn validate(options: &EvolutionaryOptions) -> Result<(), Error> {
    if options.population < 2 {
        return Err(Error::InvalidInput(
            "the population must be at least two".to_owned(),
        ));
    }
    if options.composition.is_empty() {
        return Err(Error::InvalidInput(
            "the composition needs at least one element".to_owned(),
        ));
    }
    if options.space_groups.is_empty() {
        return Err(Error::InvalidInput(
            "space_groups must list at least one group".to_owned(),
        ));
    }
    if !(0.0..1.0).contains(&options.survival) || options.survival == 0.0 {
        return Err(Error::InvalidInput("survival must be in (0, 1)".to_owned()));
    }
    if !options.mutation_rate.is_finite() || !(0.0..=1.0).contains(&options.mutation_rate) {
        return Err(Error::InvalidInput(
            "mutation_rate must be in [0, 1]".to_owned(),
        ));
    }
    if options.mutation_sigma < 0.0 || !options.mutation_sigma.is_finite() {
        return Err(Error::InvalidInput(
            "mutation_sigma must be finite and non-negative".to_owned(),
        ));
    }
    if options.strain_sigma < 0.0 || !options.strain_sigma.is_finite() {
        return Err(Error::InvalidInput(
            "strain_sigma must be finite and non-negative".to_owned(),
        ));
    }
    if options.max_steps == 0 {
        return Err(Error::InvalidInput(
            "max_steps must be at least one".to_owned(),
        ));
    }
    Ok(())
}
