//! Integration tests for RCSR net search by space group and transitivity
//! (Task 14.4.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{NetQuery, RcsrFactsIndex, plane_group_number, space_group_number};

#[test]
fn the_facts_index_holds_the_catalogue() {
    let index = RcsrFactsIndex::from_embedded().unwrap();
    assert_eq!(index.len(), 2931);
    let pcu = index.facts("pcu").unwrap();
    assert_eq!(pcu.periodicity, 3);
    assert_eq!(pcu.coordination, vec![6]);
    assert_eq!(pcu.number, Some(221));
    assert_eq!(pcu.transitivity, [Some(1), Some(1), Some(1), Some(1)]);
    assert!(pcu.summary().contains("Pm-3m (221)"));
    assert!(pcu.summary().contains("[1 1 1 1]"));
}

#[test]
fn plane_groups_answer_with_their_own_number() {
    let index = RcsrFactsIndex::from_embedded().unwrap();
    let hcb = index.facts("hcb").unwrap();
    assert_eq!(hcb.number, Some(17));
    assert_eq!(space_group_number("p6mm"), Some(183));
    assert_eq!(plane_group_number("p6mm"), Some(17));
    assert_eq!(plane_group_number("p6m"), Some(17));
    assert_eq!(plane_group_number("nonsense"), None);
    assert_eq!(space_group_number("Pm-3m"), Some(221));
    assert_eq!(space_group_number("nonsense"), None);
}

#[test]
fn search_by_space_group_and_transitivity() {
    let index = RcsrFactsIndex::from_embedded().unwrap();
    let by_number = NetQuery::new().with_space_group(221).unwrap();
    let names = index.search_names(&by_number);
    assert!(names.contains(&"pcu"));
    assert!(
        names
            .iter()
            .all(|name| index.facts(name).unwrap().number == Some(221))
    );

    let by_symbol = NetQuery::new().with_group_symbol("Fd-3m:2");
    assert!(index.search_names(&by_symbol).contains(&"dia"));

    let by_transitivity = NetQuery::new()
        .with_transitivity([Some(1), Some(1), Some(1), Some(1)])
        .unwrap();
    let names = index.search_names(&by_transitivity);
    assert!(names.contains(&"pcu"));
    assert!(names.contains(&"dia"));
    assert!(names.iter().all(|name| {
        index.facts(name).unwrap().transitivity == [Some(1), Some(1), Some(1), Some(1)]
    }));

    let unknown_face = NetQuery::new()
        .with_transitivity([None, None, Some(250), None])
        .unwrap();
    assert!(index.search(&unknown_face).is_empty());
}

#[test]
fn query_text_parses_ranges_lists_and_packed_transitivity() {
    let index = RcsrFactsIndex::from_embedded().unwrap();
    let query = NetQuery::parse("", "3,4", false, "221-230", "2 1 * *").unwrap();
    let names = index.search_names(&query);
    assert!(names.contains(&"tbo"));

    let packed = NetQuery::parse("", "", false, "", "11**").unwrap();
    assert_eq!(packed.transitivity, [Some(1), Some(1), None, None]);
    let exclusive = NetQuery::parse("", "4-c", true, "", "").unwrap();
    assert!(exclusive.exclusive);
    assert!(index.search_names(&exclusive).contains(&"dia"));
}

#[test]
fn query_errors_name_the_field() {
    assert!(NetQuery::parse("", "three", false, "", "").is_err());
    assert!(NetQuery::parse("", "", false, "0-230", "").is_err());
    assert!(NetQuery::parse("", "", false, "10-5", "").is_err());
    assert!(NetQuery::parse("", "", false, "", "1 2 3 4 5").is_err());
    assert!(NetQuery::new().with_coordination([0]).is_err());
    assert!(NetQuery::new().with_space_group(0).is_err());
    assert!(
        NetQuery::new()
            .with_transitivity([Some(0), None, None, None])
            .is_err()
    );
    assert!(NetQuery::new().with_periodicity(1).is_err());
    assert!(NetQuery::new().with_periodicity(2).is_ok());
}
