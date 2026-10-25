//! Enrichment of shared terms: does an individual have more terms in the overlap of the
//! diagnosed genes than expected by chance?
//!
//! The test is built on the same structures as the upset plot ([`UpsetPlotInput`]) and uses
//! the same definition of "shared" ([`shared_terms`]), so a result and the plot's column for
//! the same genes describe the same terms.
//!
//! # Model
//!
//! For a set of two or more genes:
//!
//! * the **urn** `U` is the set of terms annotated to at least one of the genes (explicit
//!   annotations only, no ancestors), `N = |U|`;
//! * the **shared** terms `S ⊆ U` are the terms shared by *all* the genes under the true-path
//!   rule, `K = |S|` (see [`shared_terms`]);
//! * the individual **covers** a term if the individual has it or a more specific term, that
//!   is, if the term equals or is an ancestor of an observed term;
//! * `m` is the number of urn terms the individual covers (the draw), and `k` the number of
//!   shared terms the individual covers (the statistic).
//!
//! Under the null, the `m` covered terms are a uniform random subset of the urn, so `k` follows
//! a hypergeometric distribution with parameters `N`, `K` and `m`. The alternative is that the
//! individual's terms fall into the shared part more often, as expected if the phenotype
//! reflects all the genes. The p-value is one-sided, `P(X ≥ k)`.
//!
//! Counting `k` and `m` over urn terms (not over the individual's own terms) keeps the
//! statistic and the null in the same space: `k ≤ min(m, K)` always holds, however many of the
//! individual's terms fall under a single shared annotation.
//!
//! # Simulation
//!
//! The exact test treats every urn term as equally likely to be covered. The simulated test
//! gives each urn term a weight, the mean of the annotation frequencies reported for it in the
//! diseases of the genes tested (1.0 when none is reported, by the convention of the HPO
//! project), and draws `m` urn terms without replacement with probability proportional to the
//! weight. The statistic of each simulated individual is the number of shared terms drawn, and
//! the p-value is the empirical one, `(tail + 1) / (n_sim + 1)`. With no frequency
//! information all weights are 1 and the simulation estimates the exact test.
use std::collections::BTreeSet;

use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::blend::dto::shared_term_enrichment::{SharedTermEnrichment, SimulatedEnrichment};
use crate::blend::dto::{
    GeneKey, GeneTerms, TermKey, UpsetPlotInput,
};
use crate::blend::upset_plot::{shared_terms, MAX_GENES_PER_COLUMN};
use crate::enrichment::hypergeometric::{hypergeometric_pmf, support, upper_tail};
use crate::enrichment::montecarlo::{empirical_p_value, simulate_hits};

pub(crate) const DEFAULT_N_SIM: usize = 10_000 as usize;

/// The mean of a distribution given as `pmf[x] = P(X = x)`.
fn expected_value(pmf: &[f64]) -> f64 {
    pmf.iter().enumerate().map(|(x, p)| x as f64 * p).sum()
}

/// The mean of the weights of the urn terms whose shared flag equals `shared`, if there are any.
fn mean_weight(weights: &[f64], flags: &[bool], shared: bool) -> Option<f64> {
    let selected: Vec<f64> =
        weights.iter().zip(flags).filter(|(_, f)| **f == shared).map(|(w, _)| *w).collect();
    if selected.is_empty() {
        None
    } else {
        Some(selected.iter().sum::<f64>() / selected.len() as f64)
    }
}

/// What both tests need: the genes, the urn, the shared terms and what the individual covers.
struct Setup<'a> {
    genes: Vec<&'a GeneKey>,
    gene_terms: Vec<&'a GeneTerms>,
    /// The urn in a fixed (sorted) order.
    urn: Vec<&'a TermKey>,
    /// `shared_flags[i]` says whether `urn[i]` is a shared term.
    shared_flags: Vec<bool>,
    /// `K`
    n_shared: usize,
    /// `m`
    n_drawn: usize,
    /// `k`
    n_observed: usize,
}

impl<'a> Setup<'a> {
    fn new(input: &'a UpsetPlotInput, genes: &[&str]) -> Result<Self, String> {
        // Sorted and distinct, so the order of the arguments never matters.
        let names: BTreeSet<&str> = genes.iter().copied().collect();
        if names.len() != genes.len() {
            return Err(format!("duplicate gene in {genes:?}"));
        }
        if names.len() < 2 {
            return Err(format!("need at least two genes, got {}", names.len()));
        }
        let mut selected: Vec<(&GeneKey, &GeneTerms)> = Vec::with_capacity(names.len());
        for name in &names {
            let entry = input
                .gene_annotations
                .get_key_value(*name)
                .ok_or_else(|| format!("gene {name} has no annotations in the input"))?;
            selected.push(entry);
        }
        let gene_terms: Vec<&GeneTerms> = selected.iter().map(|(_, gt)| *gt).collect();
        let indices: Vec<usize> = (0..gene_terms.len()).collect();

        let urn_set: BTreeSet<&TermKey> =
            gene_terms.iter().flat_map(|g| g.annotated.iter()).collect();
        let urn: Vec<&TermKey> = urn_set.into_iter().collect();
        let shared = shared_terms(&indices, &gene_terms);
        let shared_flags: Vec<bool> = urn.iter().map(|t| shared.contains(*t)).collect();
        debug_assert_eq!(shared.len(), shared_flags.iter().filter(|f| **f).count());

        // The individual covers an urn term if it equals or is an ancestor of an observed term,
        // that is, if it is in the union of the closures of the observed terms.
        let covered: BTreeSet<&TermKey> = input.proband_closures.values().flatten().collect();
        let n_drawn = urn.iter().filter(|t| covered.contains(*t)).count();
        let n_observed = urn
            .iter()
            .zip(&shared_flags)
            .filter(|(t, is_shared)| **is_shared && covered.contains(*t))
            .count();

        Ok(Self {
            genes: selected.iter().map(|(g, _)| *g).collect(),
            gene_terms,
            urn,
            n_shared: shared.len(),
            shared_flags,
            n_drawn,
            n_observed,
        })
    }
}

/// Runs the exact (hypergeometric) enrichment test for the named genes. The `simulated` field
/// of the result is `None`.
///
/// Two or more genes are required; three or four are tested for the terms they all share.
///
/// # Errors
/// Returns an error if fewer than two distinct genes are given or if a gene is not in
/// `input.gene_annotations`.
pub fn exact_enrichment(input: &UpsetPlotInput, genes: &[&str]) -> Result<SharedTermEnrichment, String> {
    let s = Setup::new(input, genes)?;
    let n_urn = s.urn.len();

    let pmf = hypergeometric_pmf(n_urn, s.n_shared, s.n_drawn);
    let expected_shared = expected_value(&pmf);
    let p_value = upper_tail(&pmf, s.n_observed);
    let p_mid = (p_value - 0.5 * pmf[s.n_observed]).clamp(0.0, 1.0);
    let (low, high) = support(n_urn, s.n_shared, s.n_drawn);

    Ok(SharedTermEnrichment {
        genes: s.genes.iter().map(|g| (*g).clone()).collect(),
        n_urn,
        n_shared: s.n_shared,
        n_drawn: s.n_drawn,
        n_observed: s.n_observed,
        pmf,
        p_value,
        p_mid,
        expected_shared,
        degenerate: low == high,
        simulated: None,
    })
}

/// The weight of each urn term: the mean of all the frequencies reported for it in the diseases
/// of the genes tested, or 1.0 if none is reported.
///
/// The mean is over diseases, pooled across the genes (not a mean of per-gene means), so a
/// term annotated in three diseases counts three times.
fn urn_weights(urn: &[&TermKey], gene_terms: &[&GeneTerms]) -> Vec<f64> {
    urn.iter()
        .map(|term| {
            let freqs: Vec<f64> = gene_terms
                .iter()
                .filter_map(|g| g.frequencies.get(*term))
                .flatten()
                .copied()
                .collect();
            if freqs.is_empty() {
                1.0
            } else {
                freqs.iter().sum::<f64>() / freqs.len() as f64
            }
        })
        .collect()
}

/// Runs the frequency-weighted simulation for the named genes: `n_sim` simulated individuals
/// each cover `m` urn terms drawn with probability proportional to the term's weight, and the
/// empirical p-value is the fraction of them with at least as many shared terms as observed
/// (see [`urn_weights`] for the weights). The result is reproducible for a given `seed`.
///
/// # Errors
/// Returns an error for an invalid gene list (as [`exact_enrichment`]), for `n_sim = 0`, or if
/// a weight is not finite and positive.
pub fn monte_carlo_enrichment(
    input: &UpsetPlotInput,
    genes: &[&str],
    n_sim: usize,
    seed: u64,
) -> Result<SimulatedEnrichment, String> {
    let s = Setup::new(input, genes)?;
    let weights = urn_weights(&s.urn, &s.gene_terms);
    let mut rng = StdRng::seed_from_u64(seed);
    let counts = simulate_hits(&weights, &s.shared_flags, s.n_drawn, n_sim, &mut rng)?;
    let pmf: Vec<f64> = counts.iter().map(|&c| c as f64 / n_sim as f64).collect();
    Ok(SimulatedEnrichment {
        n_sim,
        seed,
        expected_shared: expected_value(&pmf),
        mean_weight_shared: mean_weight(&weights, &s.shared_flags, true),
        mean_weight_unshared: mean_weight(&weights, &s.shared_flags, false),
        pmf,
        p_value: empirical_p_value(&counts, s.n_observed),
        n_downweighted_terms: weights.iter().filter(|w| **w < 1.0).count(),
    })
}

/// A seed for the simulation, derived from the individual and the genes tested, so the same
/// analysis always gives the same result, whatever the order of the genes.
///
/// Uses FNV-1a over the identifier and the sorted gene symbols. (The standard library's hashers
/// may change between Rust versions, which would change every result.)
pub fn seed_for(individual_id: &str, genes: &[&str]) -> u64 {
    const PRIME: u64 = 0x0000_0100_0000_01b3;
    let mut sorted = genes.to_vec();
    sorted.sort_unstable();

    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |bytes: &[u8]| {
        for byte in bytes {
            hash ^= u64::from(*byte);
            hash = hash.wrapping_mul(PRIME);
        }
        // a separator that cannot occur in a symbol, so ("AB", "C") differs from ("A", "BC")
        hash ^= 0xff;
        hash = hash.wrapping_mul(PRIME);
    };
    feed(individual_id.as_bytes());
    for gene in sorted {
        feed(gene.as_bytes());
    }
    hash
}

/// The exact test and the simulation for one set of genes.
fn enrich(
    input: &UpsetPlotInput,
    genes: &[&str],
    n_sim: usize,
    individual_id: &str,
) -> Result<SharedTermEnrichment, String> {
    let mut result = exact_enrichment(input, genes)?;
    let seed = seed_for(individual_id, genes);
    result.simulated = Some(monte_carlo_enrichment(input, genes, n_sim, seed)?);
    Ok(result)
}

/// Runs both tests for every pair of genes in the input, and, for three or four genes, also for
/// all of them together (the terms every gene shares).
///
/// The results are in a fixed order: the pairs in alphabetical order of their gene symbols,
/// then the test for all genes. With fewer than two genes there is nothing to test and the
/// list is empty. With more than four genes only the pairs are tested, as for the plot's
/// columns. Each simulation uses `n_sim` simulated individuals and a seed derived from
/// `individual_id` and the genes, so repeated runs give identical results.
///
/// # Errors
/// Propagates any error from [`exact_enrichment`] and [`monte_carlo_enrichment`].
pub fn all_enrichments(
    input: &UpsetPlotInput,
    n_sim: usize,
    individual_id: &str,
) -> Result<Vec<SharedTermEnrichment>, String> {
    let genes: Vec<&str> = input.gene_annotations.keys().map(String::as_str).collect();
    let mut results = Vec::new();
    for (i, a) in genes.iter().enumerate() {
        for b in &genes[i + 1..] {
            results.push(enrich(input, &[*a, *b], n_sim, individual_id)?);
        }
    }
    if (3..=MAX_GENES_PER_COLUMN).contains(&genes.len()) {
        results.push(enrich(input, &genes, n_sim, individual_id)?);
    }
    Ok(results)
}



#[cfg(test)]
mod tests {
    use super::*;
    use crate::blend::disease_gene_entity::GeneDiseaseEntity;
    use crate::blend::upset_plot::{plot_data, plot_input};
    use crate::tests::hpo_fixture::{
        entity, hpo, tids, FOCAL, SEIZURE, SHORT_STATURE, TONIC,
    };
    use ontolius::ontology::csr::FullCsrOntology;
    use rstest::rstest;
    use std::sync::Arc;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    fn input(
        hpo: &FullCsrOntology,
        entities: &[GeneDiseaseEntity],
        observed: &[&str],
    ) -> UpsetPlotInput {
        plot_input(entities, &tids(observed), hpo)
    }

    fn run(
        hpo: &FullCsrOntology,
        entities: &[GeneDiseaseEntity],
        observed: &[&str],
        genes: &[&str],
    ) -> SharedTermEnrichment {
        exact_enrichment(&input(hpo, entities, observed), genes).expect("valid input")
    }

    // A = {focal, short stature}, B = {seizure, short stature}.
    //   urn    = {SEIZURE, FOCAL, SHORT_STATURE}, N = 3
    //   shared = {SEIZURE, SHORT_STATURE}, K = 2 (SEIZURE stands for the pair FOCAL / SEIZURE)
    fn pair() -> [GeneDiseaseEntity; 2] {
        [entity("A", &[FOCAL, SHORT_STATURE]), entity("B", &[SEIZURE, SHORT_STATURE])]
    }

    // ---------------------------------------------------------------------------------
    // The counts come from the same definitions as the plot
    // ---------------------------------------------------------------------------------

    #[rstest]
    fn counts_and_probabilities_for_a_small_urn(hpo: Arc<FullCsrOntology>) {
        // The individual has a tonic seizure and short stature. They cover SEIZURE (an ancestor
        // of the tonic seizure) and SHORT_STATURE: both urn terms are shared.
        let r = run(&hpo, &pair(), &[TONIC, SHORT_STATURE], &["A", "B"]);
        assert_eq!((r.n_urn, r.n_shared, r.n_drawn, r.n_observed), (3, 2, 2, 2));
        // Hypergeometric(N=3, K=2, m=2): there are C(3,2) = 3 equally likely picks of two
        // urn terms. One of them contains both shared terms and two contain only one, so
        // P(X=2) = 1/3 and P(X=1) = 2/3.
        assert!(close(r.pmf[1], 2.0 / 3.0) && close(r.pmf[2], 1.0 / 3.0));
        assert!(close(r.p_value, 1.0 / 3.0));
        assert!(close(r.p_mid, 1.0 / 6.0));
        assert!(!r.degenerate);
    }

    #[rstest]
    fn the_statistic_agrees_with_the_upset_column(hpo: Arc<FullCsrOntology>) {
        let observed = [TONIC, SHORT_STATURE];
        let data = plot_data(&input(&hpo, &pair(), &observed)).expect("valid input");
        let column = data.columns.iter().find(|c| c.row_indices == [0, 1]).expect("A+B column");
        let r = run(&hpo, &pair(), &observed, &["A", "B"]);
        assert_eq!(r.n_shared, column.n_annotated_terms);
        assert_eq!(r.n_observed, column.n_overlapping_terms);
    }

    #[rstest]
    fn several_observed_terms_under_one_shared_term_count_once(hpo: Arc<FullCsrOntology>) {
        // Both genes have SEIZURE; the individual has two kinds of seizure. The urn term SEIZURE
        // is covered once, so k can never exceed K.
        let entities = [entity("A", &[SEIZURE, SHORT_STATURE]), entity("B", &[SEIZURE])];
        let r = run(&hpo, &entities, &[FOCAL, TONIC], &["A", "B"]);
        assert_eq!((r.n_urn, r.n_shared, r.n_drawn, r.n_observed), (2, 1, 1, 1));
        assert!(r.n_observed <= r.n_shared.min(r.n_drawn));
        assert!((0.0..=1.0).contains(&r.p_value));
    }

    // ---------------------------------------------------------------------------------
    // p-values by hand
    // ---------------------------------------------------------------------------------

    // A = {focal, tonic, short stature}, B = {focal, short stature}:
    //   urn = {FOCAL, TONIC, SHORT_STATURE}, N = 3;  shared = {FOCAL, SHORT_STATURE}, K = 2.
    // Only the individual changes between the cases, so the p-values can be compared.
    #[rstest]
    #[case(&[TONIC], 1, 0, 1.0)] // covers only the unshared tonic seizure
    #[case(&[FOCAL], 1, 1, 2.0 / 3.0)]
    #[case(&[FOCAL, SHORT_STATURE], 2, 2, 1.0 / 3.0)] // both shared terms: the smallest p-value
    #[case(&[FOCAL, TONIC, SHORT_STATURE], 3, 2, 1.0)] // covers the whole urn: no information
    fn p_value_for_a_fixed_urn(
        hpo: Arc<FullCsrOntology>,
        #[case] observed: &[&str],
        #[case] drawn: usize,
        #[case] hits: usize,
        #[case] p: f64,
    ) {
        let entities = [
            entity("A", &[FOCAL, TONIC, SHORT_STATURE]),
            entity("B", &[FOCAL, SHORT_STATURE]),
        ];
        let r = run(&hpo, &entities, observed, &["A", "B"]);
        assert_eq!((r.n_urn, r.n_shared), (3, 2));
        assert_eq!((r.n_drawn, r.n_observed), (drawn, hits));
        assert!(close(r.p_value, p), "p = {} but expected {p}", r.p_value);
    }

    // ---------------------------------------------------------------------------------
    // Degenerate tests carry no information
    // ---------------------------------------------------------------------------------

    #[rstest]
    #[case(&[("A", &[SHORT_STATURE][..]), ("B", &[FOCAL][..])], &[SHORT_STATURE, FOCAL])] // nothing shared
    #[case(&[("A", &[FOCAL, SHORT_STATURE][..]), ("B", &[FOCAL, SHORT_STATURE][..])], &[])] // nothing observed
    #[case(&[("A", &[FOCAL][..]), ("B", &[FOCAL][..])], &[FOCAL])] // whole urn covered
    fn degenerate_tests_report_p_value_one(
        hpo: Arc<FullCsrOntology>,
        #[case] genes: &[(&str, &[&str])],
        #[case] observed: &[&str],
    ) {
        let entities: Vec<_> = genes.iter().map(|(g, ts)| entity(g, ts)).collect();
        let r = run(&hpo, &entities, observed, &["A", "B"]);
        assert!(r.degenerate);
        assert!(close(r.p_value, 1.0));
    }

    // ---------------------------------------------------------------------------------
    // More than two genes, ordering, and bad input
    // ---------------------------------------------------------------------------------

    #[rstest]
    fn three_genes_are_tested_for_the_terms_they_all_share(hpo: Arc<FullCsrOntology>) {
        let entities = [
            entity("A", &[FOCAL, SHORT_STATURE]),
            entity("B", &[FOCAL, SHORT_STATURE]),
            entity("C", &[FOCAL, TONIC]),
        ];
        let r = run(&hpo, &entities, &[FOCAL, SHORT_STATURE], &["A", "B", "C"]);
        assert_eq!(r.genes, ["A", "B", "C"]);
        assert_eq!(r.n_urn, 3); // FOCAL, TONIC, SHORT_STATURE
        assert_eq!(r.n_shared, 1); // only FOCAL is shared by all three
        assert_eq!(r.n_drawn, 2); // the individual covers FOCAL and SHORT_STATURE
        assert_eq!(r.n_observed, 1);
    }

    #[rstest]
    fn the_order_of_the_genes_does_not_matter(hpo: Arc<FullCsrOntology>) {
        let a = run(&hpo, &pair(), &[TONIC, SHORT_STATURE], &["A", "B"]);
        let b = run(&hpo, &pair(), &[TONIC, SHORT_STATURE], &["B", "A"]);
        assert_eq!(a.genes, b.genes);
        assert!(close(a.p_value, b.p_value));
    }

    #[rstest]
    #[case(&["A"], "at least two")]
    #[case(&[], "at least two")]
    #[case(&["A", "A"], "duplicate")]
    #[case(&["A", "NOPE"], "NOPE")]
    fn invalid_gene_lists_are_errors(
        hpo: Arc<FullCsrOntology>,
        #[case] genes: &[&str],
        #[case] message: &str,
    ) {
        let err = exact_enrichment(&input(&hpo, &pair(), &[TONIC]), genes).unwrap_err();
        assert!(err.contains(message), "{err}");
    }

    #[rstest]
    fn mid_p_never_exceeds_p(hpo: Arc<FullCsrOntology>) {
        let r = run(&hpo, &pair(), &[TONIC, SHORT_STATURE], &["A", "B"]);
        assert!(r.p_mid <= r.p_value);
        assert!((0.0..=1.0).contains(&r.p_mid));
    }

    // ---------------------------------------------------------------------------------
    // Which gene sets are tested
    // ---------------------------------------------------------------------------------

    fn same_terms(genes: &[&str]) -> Vec<GeneDiseaseEntity> {
        genes.iter().map(|g| entity(g, &[FOCAL, SHORT_STATURE])).collect()
    }

    #[rstest]
    #[case(&["A"], 0)] // nothing to compare
    #[case(&["A", "B"], 1)]
    #[case(&["A", "B", "C"], 4)] // three pairs and the three genes together
    #[case(&["A", "B", "C", "D"], 7)] // six pairs and the four genes together
    #[case(&["A", "B", "C", "D", "E"], 10)] // ten pairs; five genes are not combined
    fn number_of_tests(hpo: Arc<FullCsrOntology>, #[case] genes: &[&str], #[case] expected: usize) {
        let data = input(&hpo, &same_terms(genes), &[FOCAL]);
        let results = all_enrichments(&data, 100, "id").unwrap();
        assert_eq!(results.len(), expected);
    }

    #[rstest]
    fn tests_come_in_a_fixed_order(hpo: Arc<FullCsrOntology>) {
        // The genes are given out of order; the results do not depend on that.
        let data = input(&hpo, &same_terms(&["C", "A", "B"]), &[FOCAL]);
        let results = all_enrichments(&data, 100, "id").unwrap();
        let sets: Vec<Vec<&str>> =
            results.iter().map(|r| r.genes.iter().map(String::as_str).collect()).collect();
        assert_eq!(sets, [vec!["A", "B"], vec!["A", "C"], vec!["B", "C"], vec!["A", "B", "C"]]);
    }

    // ---------------------------------------------------------------------------------
    // Weights
    // ---------------------------------------------------------------------------------

    fn gene_terms_with(freqs: &[(&str, &[f64])]) -> GeneTerms {
        GeneTerms {
            frequencies: freqs.iter().map(|(t, f)| (t.to_string(), f.to_vec())).collect(),
            ..Default::default()
        }
    }

    #[test]
    fn a_weight_is_the_mean_of_the_frequencies_across_the_diseases_of_all_genes() {
        let a = gene_terms_with(&[(FOCAL, &[0.4])]);
        let b = gene_terms_with(&[(FOCAL, &[0.8, 1.0])]);
        let (focal, tonic) = (FOCAL.to_string(), TONIC.to_string());
        let weights = urn_weights(&[&focal, &tonic], &[&a, &b]);
        assert!(close(weights[0], (0.4 + 0.8 + 1.0) / 3.0)); // pooled over three diseases
        assert!(close(weights[1], 1.0)); // no frequency reported: 100%
    }

    // ---------------------------------------------------------------------------------
    // The simulation
    // ---------------------------------------------------------------------------------

    #[rstest]
    fn the_exact_test_does_not_run_the_simulation(hpo: Arc<FullCsrOntology>) {
        let r = run(&hpo, &pair(), &[TONIC, SHORT_STATURE], &["A", "B"]);
        assert!(r.simulated.is_none());
    }

    #[rstest]
    fn without_frequencies_the_simulation_estimates_the_exact_test(hpo: Arc<FullCsrOntology>) {
        let data = input(&hpo, &pair(), &[TONIC, SHORT_STATURE]);
        let exact = exact_enrichment(&data, &["A", "B"]).unwrap();
        let sim = monte_carlo_enrichment(&data, &["A", "B"], 100_000, 1).unwrap();
        assert_eq!(sim.n_downweighted_terms, 0);
        assert_eq!(sim.pmf.len(), exact.pmf.len());
        for (x, e) in exact.pmf.iter().enumerate() {
            let tolerance = 4.0 * (e * (1.0 - e) / 100_000.0).sqrt() + 1e-9;
            assert!((sim.pmf[x] - e).abs() < tolerance, "x={x}: {} vs {e}", sim.pmf[x]);
        }
        assert!((sim.p_value - exact.p_value).abs() < 0.01, "{} vs {}", sim.p_value, exact.p_value);
    }

    #[rstest]
    fn a_rarely_reported_unshared_term_makes_shared_hits_expected(hpo: Arc<FullCsrOntology>) {
        // A = {focal, tonic, short stature}, B = {focal, short stature}. The urn is
        // {FOCAL, TONIC, SHORT_STATURE} and the shared terms are {FOCAL, SHORT_STATURE}. The
        // individual covers focal seizure and short stature: m = 2, k = 2, exact p = 1/3.
        let entities = [
            entity("A", &[FOCAL, TONIC, SHORT_STATURE]),
            entity("B", &[FOCAL, SHORT_STATURE]),
        ];
        let mut data = input(&hpo, &entities, &[FOCAL, SHORT_STATURE]);
        let equal = monte_carlo_enrichment(&data, &["A", "B"], 20_000, 3).unwrap();
        assert!((equal.p_value - 1.0 / 3.0).abs() < 0.02, "{}", equal.p_value);

        // If the unshared term (tonic seizure) is reported in only 1% of patients, a patient
        // is very unlikely to have it, so two shared terms are exactly what to expect.
        data.gene_annotations
            .get_mut("A")
            .unwrap()
            .frequencies
            .insert(TONIC.to_string(), vec![0.01]);
        let weighted = monte_carlo_enrichment(&data, &["A", "B"], 20_000, 3).unwrap();
        assert_eq!(weighted.n_downweighted_terms, 1);
        assert!(weighted.p_value > 0.9, "{}", weighted.p_value);
    }

    #[rstest]
    fn the_same_seed_gives_the_same_simulation_whatever_the_gene_order(hpo: Arc<FullCsrOntology>) {
        let data = input(&hpo, &pair(), &[TONIC, SHORT_STATURE]);
        let a = monte_carlo_enrichment(&data, &["A", "B"], 1_000, 9).unwrap();
        let b = monte_carlo_enrichment(&data, &["B", "A"], 1_000, 9).unwrap();
        assert_eq!(a.pmf, b.pmf);
        assert_eq!(a.p_value, b.p_value);
    }

    #[rstest]
    fn a_degenerate_test_has_a_simulated_p_value_of_one(hpo: Arc<FullCsrOntology>) {
        // Nothing observed, so the individual covers no urn terms.
        let entities = [entity("A", &[FOCAL, SHORT_STATURE]), entity("B", &[FOCAL, SHORT_STATURE])];
        let sim = monte_carlo_enrichment(&input(&hpo, &entities, &[]), &["A", "B"], 100, 1).unwrap();
        assert!(close(sim.p_value, 1.0));
    }

    #[rstest]
    fn zero_simulations_is_an_error(hpo: Arc<FullCsrOntology>) {
        let data = input(&hpo, &pair(), &[TONIC]);
        let err = monte_carlo_enrichment(&data, &["A", "B"], 0, 1).unwrap_err();
        assert!(err.contains("at least 1"), "{err}");
    }

    #[rstest]
    fn a_zero_frequency_in_the_input_is_an_error_not_a_silently_smaller_draw(
        hpo: Arc<FullCsrOntology>,
    ) {
        let mut data = input(&hpo, &pair(), &[TONIC]);
        data.gene_annotations
            .get_mut("A")
            .unwrap()
            .frequencies
            .insert(FOCAL.to_string(), vec![0.0]);
        let err = monte_carlo_enrichment(&data, &["A", "B"], 100, 1).unwrap_err();
        assert!(err.contains("positive"), "{err}");
    }

    // ---------------------------------------------------------------------------------
    // The seed
    // ---------------------------------------------------------------------------------

    #[test]
    fn the_seed_does_not_depend_on_the_order_of_the_genes() {
        assert_eq!(
            seed_for("patient-1", &["PTPN11", "NF1"]),
            seed_for("patient-1", &["NF1", "PTPN11"])
        );
    }

    #[test]
    fn the_seed_is_stable_across_versions() {
        // FNV-1a values computed independently of this code. If these change, every
        // simulated p-value would change with them.
        assert_eq!(seed_for("patient-1", &["NF1", "PTPN11"]), 1_470_397_748_979_433_974);
        assert_eq!(seed_for("patient-2", &["NF1", "PTPN11"]), 17_369_675_848_618_041_635);
    }

    #[test]
    fn different_genes_and_symbol_boundaries_give_different_seeds() {
        assert_ne!(seed_for("x", &["A", "B"]), seed_for("x", &["A", "C"]));
        assert_ne!(seed_for("x", &["AB", "C"]), seed_for("x", &["A", "BC"]));
    }

    #[rstest]
    fn every_result_carries_a_reproducible_simulation(hpo: Arc<FullCsrOntology>) {
        let entities: Vec<_> = ["A", "B", "C"].iter().map(|g| entity(g, &[FOCAL, SHORT_STATURE])).collect();
        let data = input(&hpo, &entities, &[FOCAL]);
        let first = all_enrichments(&data, 500, "patient-1").unwrap();
        let second = all_enrichments(&data, 500, "patient-1").unwrap();
        assert_eq!(first.len(), 4); // three pairs and the three genes together
        assert!(first.iter().all(|r| r.simulated.is_some()));
        let p_values = |results: &[SharedTermEnrichment]| {
            results.iter().map(|r| r.simulated.as_ref().unwrap().p_value).collect::<Vec<_>>()
        };
        assert_eq!(p_values(&first), p_values(&second));
        // the seed is the one derived from the individual and the genes
        assert_eq!(first[0].simulated.as_ref().unwrap().seed, seed_for("patient-1", &["A", "B"]));
    }

    // ---------------------------------------------------------------------------------
    // Diagnostics: why a weighted result differs from the exact one
    // ---------------------------------------------------------------------------------

    #[rstest]
    fn the_expected_number_of_shared_terms_is_draws_times_shared_over_urn(hpo: Arc<FullCsrOntology>) {
        // N = 3, K = 2, m = 2: the expectation is 2 * 2 / 3
        let r = run(&hpo, &pair(), &[TONIC, SHORT_STATURE], &["A", "B"]);
        assert!(close(r.expected_shared, 4.0 / 3.0), "{}", r.expected_shared);
    }

    #[rstest]
    fn without_frequencies_the_weights_are_all_one(hpo: Arc<FullCsrOntology>) {
        let data = input(&hpo, &pair(), &[TONIC, SHORT_STATURE]);
        let sim = monte_carlo_enrichment(&data, &["A", "B"], 20_000, 1).unwrap();
        assert_eq!(sim.mean_weight_shared, Some(1.0));
        assert_eq!(sim.mean_weight_unshared, Some(1.0));
        let exact = exact_enrichment(&data, &["A", "B"]).unwrap();
        assert!((sim.expected_shared - exact.expected_shared).abs() < 0.02);
    }

    #[rstest]
    fn the_diagnostics_show_why_a_weighted_result_differs(hpo: Arc<FullCsrOntology>) {
        // Urn {FOCAL, TONIC, SHORT_STATURE}; shared {FOCAL, SHORT_STATURE}; FOCAL is reported in
        // only 20% of patients.
        let entities = [
            entity("A", &[FOCAL, TONIC, SHORT_STATURE]),
            entity("B", &[FOCAL, SHORT_STATURE]),
        ];
        let mut data = input(&hpo, &entities, &[FOCAL, SHORT_STATURE]);
        data.gene_annotations.get_mut("A").unwrap().frequencies.insert(FOCAL.to_string(), vec![0.2]);
        let exact = exact_enrichment(&data, &["A", "B"]).unwrap();
        let sim = monte_carlo_enrichment(&data, &["A", "B"], 20_000, 4).unwrap();
        // the shared terms weigh 0.2 (FOCAL) and 1.0 (SHORT_STATURE); the other term weighs 1.0
        assert!(close(sim.mean_weight_shared.unwrap(), 0.6));
        assert_eq!(sim.mean_weight_unshared, Some(1.0));
        assert!(sim.expected_shared < exact.expected_shared);
    }

    #[rstest]
    fn a_group_with_no_terms_has_no_mean_weight(hpo: Arc<FullCsrOntology>) {
        // Every urn term is shared, so there are no unshared terms to average.
        let entities = [entity("A", &[FOCAL, SHORT_STATURE]), entity("B", &[FOCAL, SHORT_STATURE])];
        let data = input(&hpo, &entities, &[FOCAL]);
        let sim = monte_carlo_enrichment(&data, &["A", "B"], 100, 1).unwrap();
        assert_eq!(sim.mean_weight_unshared, None);
        assert!(sim.mean_weight_shared.is_some());
    }
}
