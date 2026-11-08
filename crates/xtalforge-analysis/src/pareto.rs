//! Pareto front for combined objectives (Task 7.9.2).

/// Indices of the non-dominated points, minimising both coordinates.
///
/// A point is dominated when another point is no worse in both objectives and
/// strictly better in at least one. The result is sorted by the first
/// objective and then the second.
#[must_use]
pub fn pareto_front(points: &[(f64, f64)]) -> Vec<usize> {
    let mut front = Vec::new();
    for (index, (first, second)) in points.iter().enumerate() {
        let dominated = points.iter().enumerate().any(|(other, (x, y))| {
            other != index && x <= first && y <= second && (x < first || y < second)
        });
        if !dominated {
            front.push(index);
        }
    }
    front.sort_by(|a, b| {
        let first = points
            .get(*a)
            .copied()
            .unwrap_or((f64::INFINITY, f64::INFINITY));
        let second = points
            .get(*b)
            .copied()
            .unwrap_or((f64::INFINITY, f64::INFINITY));
        first
            .0
            .partial_cmp(&second.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                first
                    .1
                    .partial_cmp(&second.1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    front
}
