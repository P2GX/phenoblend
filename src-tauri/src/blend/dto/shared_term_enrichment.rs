/// Result of the exact enrichment test for one set of genes: does the individual have more
/// terms in the overlap of the genes than expected by chance?
///
/// See `blend::shared_enrichment` for the model. In short: the *urn* is every term annotated
/// to at least one of the genes (`n_urn`), the *shared* terms are those all the genes
/// account for (`n_shared`), and the individual *covers* an urn term if they have it or a
/// more specific one. The individual covers `n_drawn` urn terms, `n_observed` of them
/// shared; under the null the covered terms are a random subset of the urn.

use serde::{Deserialize, Serialize};

use crate::blend::dto::GeneKey;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SharedTermEnrichment {
    /// The genes tested, in sorted order.
    pub genes: Vec<GeneKey>,
    /// `N`: terms annotated to at least one of the genes.
    pub n_urn: usize,
    /// `K`: terms shared by all of the genes.
    pub n_shared: usize,
    /// `m`: urn terms the individual covers (the number of draws).
    pub n_drawn: usize,
    /// `k`: shared terms the individual covers (the statistic).
    pub n_observed: usize,
    /// `pmf[x] = P(X = x)` under the null, for `x` in `0..=n_drawn`. Zero outside the support.
    pub pmf: Vec<f64>,
    /// One-sided p-value, `P(X ≥ n_observed)`.
    pub p_value: f64,
    /// Mid-p value, `P(X > n_observed) + ½·P(X = n_observed)`. Under the null its mean is ½,
    /// which makes it better behaved than `p_value` when comparing many discrete tests.
    pub p_mid: f64,
    /// The expected number of shared terms among the `n_drawn` covered ones under the null,
    /// `n_drawn · n_shared / n_urn`. Compare with `n_observed`.
    pub expected_shared: f64,
    /// True if the null has a single possible outcome (no shared terms, nothing covered, or
    /// the whole urn covered). The test carries no information then and `p_value` is 1.
    pub degenerate: bool,
     /// The frequency-weighted simulation for the same genes, if it was run.
    pub simulated: Option<SimulatedEnrichment>,
}


/// Serializes a `u64` as a decimal string.
///
/// JSON numbers are read as 64-bit floats by JavaScript, which keeps only 53 bits of integer
/// precision, so a 64-bit seed sent as a number would arrive in the frontend with its last
/// digits changed and could no longer be used to reproduce a result.
mod u64_as_string {
    use serde::{de::Error, Deserialize, Deserializer, Serializer};

    pub fn serialize<S: Serializer>(value: &u64, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&value.to_string())
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(deserializer: D) -> Result<u64, D::Error> {
        let text = String::deserialize(deserializer)?;
        text.parse().map_err(Error::custom)
    }
}

/// Result of the frequency-weighted simulation (see [`monte_carlo_enrichment`]).
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SimulatedEnrichment {
    /// Number of simulated individuals.
    pub n_sim: usize,
    /// Seed of the random number generator, so the result can be reproduced.
    #[serde(with = "u64_as_string")]
    pub seed: u64,
    /// `pmf[x]`: the fraction of simulated individuals with `x` shared terms among the
    /// `n_drawn` they cover, for `x` in `0..=n_drawn`.
    pub pmf: Vec<f64>,
    /// Empirical one-sided p-value, `(tail + 1) / (n_sim + 1)`: never zero, and at least
    /// `1 / (n_sim + 1)`.
    pub p_value: f64,
    /// Number of urn terms with a weight below 1, that is, with a reported frequency below
    /// 100%. If this is 0 the weights are all equal and the simulation estimates the exact test.
    pub n_downweighted_terms: usize,
    /// The expected number of shared terms among the covered ones under the weighted null
    /// (the mean of the simulated distribution). If it is far below the exact test's
    /// `expected_shared`, the weights make shared terms less likely to be covered.
    pub expected_shared: f64,
    /// Mean weight of the shared urn terms, or `None` if there are none.
    pub mean_weight_shared: Option<f64>,
    /// Mean weight of the urn terms that are not shared, or `None` if there are none.
    pub mean_weight_unshared: Option<f64>,
}