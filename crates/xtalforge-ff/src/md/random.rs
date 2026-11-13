//! A small deterministic random-number generator for molecular dynamics.
//!
//! SplitMix64 supplies the uniform bits; a normal deviate is the Irwin-Hall
//! sum of twelve uniforms minus six.  That sum uses only arithmetic, so the
//! stream is identical on every platform and the far tails are bounded at
//! `+-6`, which is harmless for thermostat noise and barostat proposals and is
//! documented in the module help of [`crate::md`].

/// Seeded SplitMix64 generator.
#[derive(Debug, Clone)]
pub(crate) struct SeededRng {
    state: u64,
}

impl SeededRng {
    /// A generator seeded with `seed`.
    pub(crate) fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Next 64 random bits.
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut value = self.state;
        value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        value ^ (value >> 31)
    }

    /// Uniform deviate in `(0, 1)`.
    pub(crate) fn uniform(&mut self) -> f64 {
        let bits = self.next_u64() >> 11;
        (bits as f64 + 0.5) / 9_007_199_254_740_992.0
    }

    /// Approximate standard normal deviate in `[-6, 6]`.
    pub(crate) fn normal(&mut self) -> f64 {
        let mut total = 0.0;
        for _ in 0..12 {
            total += self.uniform();
        }
        total - 6.0
    }
}
