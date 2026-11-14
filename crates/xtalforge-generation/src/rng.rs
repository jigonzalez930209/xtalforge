//! Deterministic SplitMix64 generator and sampling helpers.
//!
//! Every generator in this crate takes an explicit `seed` (AGENTS.md rule 6)
//! and draws from this state, so the same seed and options produce the same
//! structures on every platform and thread count.

/// Small deterministic generator.
pub(crate) struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub(crate) fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform value in `[0, 1)`.
    pub(crate) fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }

    /// Uniform value in `[low, high)`.
    pub(crate) fn range(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * self.uniform()
    }

    /// Uniform integer in `[0, n)`; `n` must be positive.
    pub(crate) fn below(&mut self, n: usize) -> usize {
        if n <= 1 {
            return 0;
        }
        (self.next_u64() % n as u64) as usize
    }

    /// In-place Fisher-Yates shuffle.
    pub(crate) fn shuffle<T>(&mut self, values: &mut [T]) {
        for index in (1..values.len()).rev() {
            let other = self.below(index + 1);
            values.swap(index, other);
        }
    }

    /// A uniform unit quaternion `(w, x, y, z)`.
    ///
    /// A point uniformly distributed in the unit ball supplies the vector
    /// part and `sqrt(1 - |v|^2)` the scalar part; a random sign covers both
    /// hemispheres. This is uniform on the 4-sphere because the ball volume
    /// element projects onto the spherical surface element.
    pub(crate) fn quaternion(&mut self) -> [f64; 4] {
        loop {
            let x = self.range(-1.0, 1.0);
            let y = self.range(-1.0, 1.0);
            let z = self.range(-1.0, 1.0);
            let s = x * x + y * y + z * z;
            if s >= 1.0 {
                continue;
            }
            let scalar = (1.0 - s).sqrt();
            return if self.next_u64() & 1 == 0 {
                [scalar, x, y, z]
            } else {
                [-scalar, x, y, z]
            };
        }
    }
}
