//! The hypergeometric distribution.
//!
//! Draw `m` items without replacement from an urn of `n` items, `k` of which are "successes".
//! The number of successes drawn, `X`, has a hypergeometric distribution. Everything here is
//! pure arithmetic on small integers; it knows nothing about ontologies or genes.
//!
//! Example: an urn of 10 items with 4 successes, and 5 draws. There are C(10,5) = 252 equally
//! likely draws. Exactly 3 successes can occur in C(4,3)·C(6,2) = 60 of them, and 4 successes
//! in C(4,4)·C(6,1) = 6, so `P(X ≥ 3) = 66/252`.

/// The smallest and largest value `X` can take: `(max(0, m + k − n), min(m, k))`.
///
/// If the urn has fewer failures than draws (`m > n − k`), some draws must be successes, which
/// raises the lower end. If there are fewer successes than draws, the upper end is `k`.
pub fn support(n: usize, k: usize, m: usize) -> (usize, usize) {
    ((m + k).saturating_sub(n), m.min(k))
}

/// `P(X ≥ x)`, summed from the upper tail so that tiny values keep their precision.
///
/// Computing `1 − P(X < x)` instead would lose all precision for tails below about 1e-16 and
/// could even give a negative number.
///
/// `x` must be at most `pmf.len() − 1`.
pub fn upper_tail(pmf: &[f64], x: usize) -> f64 {
    pmf[x..].iter().sum::<f64>().clamp(0.0, 1.0)
}

/// `pmf[x] = P(X = x)` for `x` in `0..=m`, where `X` is the number of successes in `m` draws
/// without replacement from `n` items of which `k` are successes. Requires `k ≤ n` and
/// `m ≤ n`. Values outside [`support`] are exactly zero.
///
/// Computed from the ratio of neighbouring probabilities, starting at the mode with weight 1
/// and normalizing at the end. Every step multiplies by a ratio of small integers, so the
/// relative error stays tiny even far in the tails. (A table of log factorials accumulates an
/// error of the order of 1e-12 for a few hundred items.)
pub fn hypergeometric_pmf(n: usize, k: usize, m: usize) -> Vec<f64> {
    debug_assert!(k <= n && m <= n);
    let (low, high) = support(n, k, m);
    let mode = (((m + 1) * (k + 1)) / (n + 2)).clamp(low, high);

    // weights[x] is proportional to P(X = x)
    let mut weights = vec![0.0_f64; m + 1];
    weights[mode] = 1.0;
    // P(x+1) / P(x) = (k-x)(m-x) / ((x+1)(n-k-m+x+1))
    for x in mode..high {
        let up = ((k - x) * (m - x)) as f64;
        let down = ((x + 1) * (n + x + 1 - k - m)) as f64;
        weights[x + 1] = weights[x] * up / down;
    }
    // P(x-1) / P(x) = x(n-k-m+x) / ((k-x+1)(m-x+1))
    for x in (low + 1..=mode).rev() {
        let up = (x * (n + x - k - m)) as f64;
        let down = ((k - x + 1) * (m - x + 1)) as f64;
        weights[x - 1] = weights[x] * up / down;
    }

    let total: f64 = weights.iter().sum();
    weights.iter().map(|w| w / total).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[rstest]
    #[case(10, 4, 5, 0, 4)]
    #[case(10, 4, 8, 2, 4)] // 8 draws from 6 failures: at least 2 successes
    #[case(10, 4, 2, 0, 2)] // fewer draws than successes
    #[case(10, 0, 5, 0, 0)] // no successes at all
    #[case(10, 10, 5, 5, 5)] // everything is a success
    #[case(10, 4, 10, 4, 4)] // the whole urn is drawn
    fn support_of_the_distribution(
        #[case] n: usize,
        #[case] k: usize,
        #[case] m: usize,
        #[case] low: usize,
        #[case] high: usize,
    ) {
        assert_eq!(support(n, k, m), (low, high));
    }

    /// N=10 items, K=4 successes, m=5 draws: C(10,5) = 252 equally likely draws.
    #[test]
    fn matches_hand_computation() {
        let pmf = hypergeometric_pmf(10, 4, 5);
        let expected = [
            6.0 / 252.0,   // x=0: C(4,0)·C(6,5) = 1·6
            60.0 / 252.0,  // x=1: C(4,1)·C(6,4) = 4·15
            120.0 / 252.0, // x=2: C(4,2)·C(6,3) = 6·20
            60.0 / 252.0,  // x=3: C(4,3)·C(6,2) = 4·15
            6.0 / 252.0,   // x=4: C(4,4)·C(6,1) = 1·6
            0.0,           // x=5: more successes than exist
        ];
        for (x, e) in expected.iter().enumerate() {
            assert!(close(pmf[x], *e), "pmf[{x}] = {} but expected {e}", pmf[x]);
        }
        assert!(close(upper_tail(&pmf, 3), 66.0 / 252.0));
        assert!(close(upper_tail(&pmf, 4), 6.0 / 252.0));
        assert!(close(upper_tail(&pmf, 0), 1.0));
    }

    #[rstest]
    #[case(10, 4, 5)]
    #[case(10, 0, 5)]
    #[case(10, 10, 5)]
    #[case(10, 4, 10)]
    #[case(10, 4, 0)]
    #[case(1, 1, 1)]
    #[case(0, 0, 0)]
    #[case(500, 120, 80)]
    fn probabilities_sum_to_one(#[case] n: usize, #[case] k: usize, #[case] m: usize) {
        let total: f64 = hypergeometric_pmf(n, k, m).iter().sum();
        assert!(close(total, 1.0), "sum = {total}");
    }

    /// C(n, r) with exact integer arithmetic; fine for n up to about 100.
    fn binomial(n: u128, r: u128) -> u128 {
        (0..r).fold(1, |acc, i| acc * (n - i) / (i + 1))
    }

    /// An independent calculation: C(K,x)·C(N-K,m-x) / C(N,m), from exact integers.
    #[rstest]
    #[case(10, 4, 5)]
    #[case(40, 7, 12)]
    #[case(60, 30, 30)]
    #[case(60, 3, 45)]
    #[case(25, 25, 10)]
    #[case(25, 0, 10)]
    fn agrees_with_exact_integer_arithmetic(#[case] n: u128, #[case] k: u128, #[case] m: u128) {
        let pmf = hypergeometric_pmf(n as usize, k as usize, m as usize);
        for x in 0..=m {
            let exact = if x > k || m - x > n - k {
                0.0
            } else {
                (binomial(k, x) * binomial(n - k, m - x)) as f64 / binomial(n, m) as f64
            };
            assert!(
                (pmf[x as usize] - exact).abs() <= 1e-12 * exact.max(1e-300) + 1e-15,
                "x={x}: {} vs {exact}",
                pmf[x as usize]
            );
        }
    }

    #[test]
    fn tiny_tails_keep_their_precision() {
        // All 40 shared out of 400, and 40 drawn: P(X = 40) is astronomically small but must
        // be neither zero nor negative.
        let pmf = hypergeometric_pmf(400, 40, 40);
        let p = upper_tail(&pmf, 40);
        assert!(p > 0.0 && p < 1e-40, "{p}");
    }

    #[test]
    fn the_whole_tail_is_one() {
        let pmf = hypergeometric_pmf(30, 9, 12);
        assert!(close(upper_tail(&pmf, 0), 1.0));
    }
}
