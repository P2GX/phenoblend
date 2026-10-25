//! Simulating the null distribution by weighted sampling.
//!
//! The exact test treats every urn term as equally likely to be among the terms an individual
//! covers. Here each term has a weight (for example its annotation frequency), and a
//! simulated individual covers `m` urn terms drawn without replacement, each draw made with
//! probability proportional to weight among the terms not yet drawn. The statistic of each
//! simulated individual is the number of drawn terms that are shared.
//!
//! With equal weights this is the hypergeometric distribution of
//! [`crate::hypergeometric`], estimated by sampling, which is how the simulation is tested.
//!
//! Everything here works on plain numbers: terms are positions `0..n` in a slice, and
//! whether a term is shared is a `bool`. The caller decides what the positions mean.
use rand::{Rng, RngExt};

/// Draws `m` distinct positions from `0..weights.len()`, one after another, each with
/// probability proportional to its weight among the positions not yet drawn.
///
/// Uses the Efraimidis–Spirakis keys: give position `i` the key `ln(u_i) / w_i` with `u_i`
/// uniform on `(0, 1)`, and keep the `m` largest keys. This gives exactly the sequential
/// sampling described above, in one pass.
///
/// All weights must be finite and positive (checked by [`simulate_hits`]), and `m` must not
/// exceed `weights.len()`. The positions are returned in no particular order.
/// 
/// 




pub fn weighted_sample(weights: &[f64], m: usize, rng: &mut impl Rng) -> Vec<usize> {
    debug_assert!(m <= weights.len());
    debug_assert!(weights.iter().all(|w| w.is_finite() && *w > 0.0));
    if m == weights.len() {
        return (0..m).collect();
    }
    if m == 0 {
        return Vec::new();
    }
    // ln(u) < 0, so a larger weight gives a key closer to zero, i.e. a larger key.
    let mut keyed: Vec<(f64, usize)> = weights
        .iter()
        .enumerate()
        .map(|(i, w)| (rng.random_range(f64::EPSILON..1.0).ln() / w, i))
        .collect();
    // Move the m largest keys to the front; their order does not matter.
    keyed.select_nth_unstable_by(m - 1, |a, b| b.0.total_cmp(&a.0));
    keyed[..m].iter().map(|&(_, i)| i).collect()
}

/// Simulates `n_sim` individuals and returns the histogram of their statistic:
/// `counts[x]` is the number of simulated individuals with exactly `x` shared terms among the
/// `m` they cover, for `x` in `0..=m`.
///
/// `weights[i]` is the weight of urn term `i` and `shared[i]` says whether it is a shared term.
///
/// # Errors
/// Returns an error if `weights` and `shared` differ in length, if `m` exceeds the number of
/// terms, if `n_sim` is zero, or if a weight is not finite and positive.
pub fn simulate_hits(
    weights: &[f64],
    shared: &[bool],
    m: usize,
    n_sim: usize,
    rng: &mut impl Rng,
) -> Result<Vec<usize>, String> {
    if weights.len() != shared.len() {
        return Err(format!(
            "{} weights but {} shared flags",
            weights.len(),
            shared.len()
        ));
    }
    if m > weights.len() {
        return Err(format!("cannot draw {m} terms from an urn of {}", weights.len()));
    }
    if n_sim == 0 {
        return Err("the number of simulations must be at least 1".to_string());
    }
    if let Some(bad) = weights.iter().find(|w| !(w.is_finite() && **w > 0.0)) {
        return Err(format!("weights must be finite and positive, found {bad}"));
    }

    let mut counts = vec![0usize; m + 1];
    for _ in 0..n_sim {
        let hits = weighted_sample(weights, m, rng)
            .into_iter()
            .filter(|&i| shared[i])
            .count();
        counts[hits] += 1;
    }
    Ok(counts)
}

/// The empirical one-sided p-value: the fraction of simulated individuals with at least
/// `observed` shared terms, counting the observed individual itself among the simulated ones.
///
/// That is `(tail + 1) / (n_sim + 1)`, where `tail` is the number of simulated individuals
/// with `observed` or more shared terms. It is never zero (the smallest value is
/// `1 / (n_sim + 1)`), and it is a valid p-value for a simulated null.
///
/// `observed` must be at most `counts.len() − 1`.
pub fn empirical_p_value(counts: &[usize], observed: usize) -> f64 {
    let n_sim: usize = counts.iter().sum();
    let tail: usize = counts[observed..].iter().sum();
    (tail as f64 + 1.0) / (n_sim as f64 + 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enrichment::hypergeometric::hypergeometric_pmf;
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use rstest::rstest;

    fn rng(seed: u64) -> StdRng {
        StdRng::seed_from_u64(seed)
    }

    /// Four standard errors of a proportion `p` estimated from `n` trials: a tolerance that a
    /// correct simulation with a fixed seed essentially never exceeds.
    fn tolerance(p: f64, n: usize) -> f64 {
        4.0 * (p * (1.0 - p) / n as f64).sqrt() + 1e-9
    }

    // ---------------------------------------------------------------------------------
    // weighted_sample
    // ---------------------------------------------------------------------------------

    #[rstest]
    #[case(10, 0)]
    #[case(10, 3)]
    #[case(10, 9)]
    #[case(10, 10)]
    #[case(1, 1)]
    fn a_sample_has_m_distinct_positions_in_range(#[case] n: usize, #[case] m: usize) {
        let weights: Vec<f64> = (1..=n).map(|i| i as f64).collect();
        let mut sample = weighted_sample(&weights, m, &mut rng(1));
        assert_eq!(sample.len(), m);
        sample.sort_unstable();
        sample.dedup();
        assert_eq!(sample.len(), m, "positions must be distinct");
        assert!(sample.iter().all(|&i| i < n));
    }

    #[test]
    fn the_same_seed_gives_the_same_sample() {
        let weights = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let a = weighted_sample(&weights, 3, &mut rng(42));
        let b = weighted_sample(&weights, 3, &mut rng(42));
        assert_eq!(a, b);
    }

    #[test]
    fn with_equal_weights_every_position_is_equally_likely() {
        let (n, m, trials) = (10, 3, 40_000);
        let weights = vec![1.0; n];
        let mut rng = rng(7);
        let mut seen = vec![0usize; n];
        for _ in 0..trials {
            for i in weighted_sample(&weights, m, &mut rng) {
                seen[i] += 1;
            }
        }
        let expected = m as f64 / n as f64; // each position is in 3 of the 10
        for (i, c) in seen.iter().enumerate() {
            let p = *c as f64 / trials as f64;
            assert!((p - expected).abs() < tolerance(expected, trials), "position {i}: {p}");
        }
    }

    #[test]
    fn a_heavy_term_is_drawn_in_proportion_to_its_weight() {
        // Weights 10, 1, 1, 1, 1; one draw: P(position 0) = 10/14.
        let weights = [10.0, 1.0, 1.0, 1.0, 1.0];
        let trials = 40_000;
        let mut rng = rng(11);
        let hits = (0..trials)
            .filter(|_| weighted_sample(&weights, 1, &mut rng) == [0])
            .count();
        let expected = 10.0 / 14.0;
        assert!((hits as f64 / trials as f64 - expected).abs() < tolerance(expected, trials));
    }

    #[test]
    fn two_draws_are_sequential_not_independent() {
        // Weights 10, 1, 1, 1, 1 and two draws. Position 0 is missed only if it is passed over
        // twice: 4/14 for the first draw, then 3/13 of the remaining weight for the second.
        // So P(position 0 is in the sample) = 1 − (4/14)(3/13) ≈ 0.934.
        let weights = [10.0, 1.0, 1.0, 1.0, 1.0];
        let trials = 40_000;
        let mut rng = rng(13);
        let hits = (0..trials)
            .filter(|_| weighted_sample(&weights, 2, &mut rng).contains(&0))
            .count();
        let expected = 1.0 - (4.0 / 14.0) * (3.0 / 13.0);
        assert!((hits as f64 / trials as f64 - expected).abs() < tolerance(expected, trials));
    }

    // ---------------------------------------------------------------------------------
    // simulate_hits
    // ---------------------------------------------------------------------------------

    /// The simulation with equal weights is the hypergeometric distribution, sampled.
    #[rstest]
    #[case(30, 8, 10)]
    #[case(20, 6, 8)] // the example on the documentation page
    #[case(15, 3, 12)]
    fn equal_weights_reproduce_the_hypergeometric_distribution(
        #[case] n: usize,
        #[case] k: usize,
        #[case] m: usize,
    ) {
        let n_sim = 100_000;
        let weights = vec![1.0; n];
        let shared: Vec<bool> = (0..n).map(|i| i < k).collect();
        let counts = simulate_hits(&weights, &shared, m, n_sim, &mut rng(2024)).unwrap();
        let exact = hypergeometric_pmf(n, k, m);
        for x in 0..=m {
            let estimate = counts[x] as f64 / n_sim as f64;
            assert!(
                (estimate - exact[x]).abs() < tolerance(exact[x], n_sim),
                "x={x}: simulated {estimate} but exact {}",
                exact[x]
            );
        }
    }

    #[test]
    fn heavier_shared_terms_give_more_shared_hits() {
        let (n, k, m, n_sim) = (30, 8, 10, 20_000);
        let shared: Vec<bool> = (0..n).map(|i| i < k).collect();
        let mean = |weights: &[f64]| {
            let counts = simulate_hits(weights, &shared, m, n_sim, &mut rng(5)).unwrap();
            counts.iter().enumerate().map(|(x, c)| x * c).sum::<usize>() as f64 / n_sim as f64
        };
        let equal = mean(&vec![1.0; n]);
        let heavy_shared: Vec<f64> = (0..n).map(|i| if i < k { 4.0 } else { 1.0 }).collect();
        let heavy_unshared: Vec<f64> = (0..n).map(|i| if i < k { 0.25 } else { 1.0 }).collect();
        assert!(mean(&heavy_shared) > equal + 0.5);
        assert!(mean(&heavy_unshared) < equal - 0.5);
    }

    #[test]
    fn the_same_seed_gives_the_same_histogram() {
        let weights = [1.0, 0.5, 0.25, 1.0, 0.75, 1.0];
        let shared = [true, true, false, false, true, false];
        let a = simulate_hits(&weights, &shared, 3, 1_000, &mut rng(99)).unwrap();
        let b = simulate_hits(&weights, &shared, 3, 1_000, &mut rng(99)).unwrap();
        assert_eq!(a, b);
        assert_eq!(a.iter().sum::<usize>(), 1_000);
    }

    #[test]
    fn hits_cannot_exceed_the_number_of_shared_terms() {
        // Only one shared term: no simulated individual can have two.
        let shared = [true, false, false, false, false, false];
        let counts = simulate_hits(&[1.0; 6], &shared, 4, 5_000, &mut rng(3)).unwrap();
        assert_eq!(counts.len(), 5);
        assert!(counts[2..].iter().all(|c| *c == 0));
    }

    #[rstest]
    #[case(&[1.0, 1.0], &[true], 1, 10, "shared flags")]
    #[case(&[1.0, 1.0], &[true, false], 3, 10, "cannot draw")]
    #[case(&[1.0, 1.0], &[true, false], 1, 0, "at least 1")]
    #[case(&[1.0, 0.0], &[true, false], 1, 10, "positive")]
    #[case(&[1.0, -2.0], &[true, false], 1, 10, "positive")]
    #[case(&[1.0, f64::NAN], &[true, false], 1, 10, "positive")]
    #[case(&[1.0, f64::INFINITY], &[true, false], 1, 10, "positive")]
    fn invalid_input_is_an_error(
        #[case] weights: &[f64],
        #[case] shared: &[bool],
        #[case] m: usize,
        #[case] n_sim: usize,
        #[case] message: &str,
    ) {
        let err = simulate_hits(weights, shared, m, n_sim, &mut rng(0)).unwrap_err();
        assert!(err.contains(message), "{err}");
    }

    // ---------------------------------------------------------------------------------
    // empirical_p_value
    // ---------------------------------------------------------------------------------

    #[rstest]
    #[case(0, 11.0 / 11.0)] // every simulated individual has 0 or more
    #[case(1, 6.0 / 11.0)] // 3 + 2 simulated individuals have 1 or more, plus the observed one
    #[case(2, 3.0 / 11.0)] // 2 simulated individuals have 2 or more
    fn empirical_p_value_counts_the_tail_plus_one(#[case] observed: usize, #[case] expected: f64) {
        // 10 simulated individuals: 5 with no shared terms, 3 with one, 2 with two.
        let counts = [5, 3, 2];
        assert!((empirical_p_value(&counts, observed) - expected).abs() < 1e-12);
    }

    #[test]
    fn the_empirical_p_value_is_never_zero() {
        // Nobody in the simulation reached the observed value.
        let counts = [10_000, 0, 0, 0];
        let p = empirical_p_value(&counts, 3);
        assert!((p - 1.0 / 10_001.0).abs() < 1e-15);
        assert!(p > 0.0);
    }
}
