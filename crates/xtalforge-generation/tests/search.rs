//! Random and evolutionary search tests (Tasks 20.2.1 and 20.2.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_generation::{
    EvolutionaryOptions, RandomSearchOptions, evolutionary_search, random_search,
};

fn composition() -> Vec<(String, u32)> {
    vec![("Na".to_owned(), 4), ("Cl".to_owned(), 4)]
}

#[test]
fn random_search_ranks_by_energy_per_atom() {
    let options = RandomSearchOptions {
        composition: composition(),
        space_groups: vec![225],
        n_structures: 3,
        seed: 21,
        volume_per_atom: 22.0,
        min_distance: 1.2,
        max_attempts: 64,
        relax: false,
        relax_cell: false,
        max_steps: 10,
        deduplicate: false,
        matcher: Default::default(),
        minimizer: Default::default(),
    };
    let report = random_search(&options).expect("search");
    assert_eq!(report.generated, 3);
    assert_eq!(report.relaxed, 0);
    assert_eq!(report.hits.len(), 3);
    for pair in report.hits.windows(2) {
        assert!(pair[0].energy_per_atom <= pair[1].energy_per_atom);
    }
    assert!(report.best().is_some());
    assert_eq!(report.calculator, "uff");
}

#[test]
fn random_search_removes_duplicates() {
    let options = RandomSearchOptions {
        composition: composition(),
        space_groups: vec![225],
        n_structures: 3,
        seed: 3,
        volume_per_atom: 22.0,
        min_distance: 1.2,
        max_attempts: 64,
        relax: false,
        relax_cell: false,
        max_steps: 10,
        deduplicate: true,
        matcher: Default::default(),
        minimizer: Default::default(),
    };
    let report = random_search(&options).expect("search");
    assert!(report.duplicates_removed >= 2);
    assert_eq!(report.hits.len(), 1);
}

#[test]
fn random_search_validates_options() {
    let mut options = RandomSearchOptions::new(composition(), vec![225], 1, 1);
    options.n_structures = 0;
    assert!(random_search(&options).is_err());
    options = RandomSearchOptions::new(composition(), Vec::new(), 1, 1);
    assert!(random_search(&options).is_err());
    options = RandomSearchOptions::new(composition(), vec![0], 1, 1);
    assert!(random_search(&options).is_err());
    options = RandomSearchOptions::new(Vec::new(), vec![225], 1, 1);
    assert!(random_search(&options).is_err());
}

#[test]
fn evolutionary_search_reports_generations() {
    let options = EvolutionaryOptions {
        population: 2,
        generations: 1,
        max_steps: 3,
        ..EvolutionaryOptions::new(composition(), vec![225], 5)
    };
    let report = evolutionary_search(&options).expect("search");
    assert_eq!(report.generations.len(), 2);
    assert_eq!(report.generations[0].generation, 0);
    assert_eq!(report.generations[1].generation, 1);
    assert!(report.relaxations >= 2);
    assert!(!report.hits.is_empty());
    for pair in report.hits.windows(2) {
        assert!(pair[0].energy_per_atom <= pair[1].energy_per_atom);
    }
}

#[test]
fn evolutionary_search_validates_options() {
    let mut options = EvolutionaryOptions::new(composition(), vec![225], 1);
    options.population = 1;
    assert!(evolutionary_search(&options).is_err());
    options = EvolutionaryOptions::new(composition(), Vec::new(), 1);
    assert!(evolutionary_search(&options).is_err());
    options = EvolutionaryOptions::new(composition(), vec![225], 1);
    options.survival = 0.0;
    assert!(evolutionary_search(&options).is_err());
    options = EvolutionaryOptions::new(composition(), vec![225], 1);
    options.mutation_rate = 2.0;
    assert!(evolutionary_search(&options).is_err());
}
