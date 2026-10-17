//! Data behind the spread plot: for each organ system (top-level HPO category), the share of
//! HPO terms that fall under it, for the proband ("Ppkt") and for every combination of the
//! genotypes' disease models (A, B, A+B, ...).
//!
//! Besides the numbers the graphic reads, the payload carries the terms behind every number
//! (`ppktTerms`, `geneTerms`) and each series' denominator and uncategorized terms
//! (`seriesSummaries`), so that any value can be checked: `value * nTerms = number of terms`.

use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use ontolius::ontology::csr::FullCsrOntology;
use ontolius::ontology::{OntologyTerms, TaxonomyQuery, TaxonomyWalk};
use ontolius::term::MinimalTerm;
use ontolius::TermId;

use crate::blend::disease_gene_entity::GeneDiseaseEntity;
use crate::blend::dto::{
    assign_to_categories, CategoryAssignment, SeriesSummary, SeriesTerms, SpreadPlotCategory,
    SpreadPlotPayload,
};
use crate::blend::upset_plot::{simple_term, UpsetData};
use crate::model::proband::Proband;
use crate::model::simple_term::SimpleOntologyTerm;

/// The plot has one column for every non-empty subset of the genotypes (2^n - 1).
const MAX_GENOTYPES: usize = 10;
/// A category whose value is below this in every column is left out of the plot.
const EPSILON: f64 = 1e-8;

/// Fetches the top-level HPO organ-system categories (direct children of HP:0000118),
/// sorted by term ID so that the result does not depend on the ontology's iteration order.
fn top_level_categories(hpo: &FullCsrOntology) -> Result<Vec<TermId>, String> {
    let phenotype_root: TermId = "HP:0000118".parse().map_err(|_| "Could not map root term".to_string())?;
    let mut categories: Vec<TermId> = hpo.iter_child_ids(&phenotype_root).cloned().collect();
    categories.sort_by_key(|t| t.to_string());
    Ok(categories)
}

/// Computes every non-empty subset (of every size) of the given genotype symbols.
fn genotype_power_set(genotypes: &[String]) -> Vec<Vec<String>> {
    let mut power_set = Vec::new();
    for size in 1..=genotypes.len() {
        for combo in genotypes.iter().cloned().combinations(size) {
            power_set.push(combo);
        }
    }
    power_set
}

/// How the union of a genotype combination's explicit disease annotations falls into the
/// categories.
fn combo_assignment(
    combo: &[String],
    entity_lookup: &HashMap<String, &GeneDiseaseEntity>,
    categories: &[TermId],
    is_ancestor_of: impl Fn(&TermId, &TermId) -> bool,
) -> CategoryAssignment {
    let union: Vec<TermId> = combo
        .iter()
        .filter_map(|g| entity_lookup.get(g))
        .flat_map(|e| e.disease_hpo_ids.iter().cloned())
        .collect();
    assign_to_categories(&union, categories, is_ancestor_of)
}

/// The category assignment of every combination in the power set, keyed by the
/// combination's canonical "+"-joined label.
fn build_combo_column_data(
    genes_power_set: &[Vec<String>],
    entity_lookup: &HashMap<String, &GeneDiseaseEntity>,
    categories: &[TermId],
    is_ancestor_of: impl Fn(&TermId, &TermId) -> bool,
) -> HashMap<String, CategoryAssignment> {
    genes_power_set
        .iter()
        .map(|combo| {
            let assignment = combo_assignment(combo, entity_lookup, categories, &is_ancestor_of);
            (combo.join("+"), assignment)
        })
        .collect()
}

/// Scores each individual gene by the dot product of its category fractions with the
/// patient's category fractions, then returns the gene symbols ranked descending by that
/// score (ties keep the input order). A gene without any annotations scores 0.
fn rank_genes_by_similarity(
    genotypes: &[String],
    columns: &HashMap<String, CategoryAssignment>,
    ppkt_fractions: &HashMap<TermId, f64>,
    categories: &[TermId],
) -> Vec<String> {
    let mut scored: Vec<(String, f64)> = genotypes
        .iter()
        .map(|gene| {
            let score = match columns.get(gene) {
                Some(column) if column.n_terms > 0 => {
                    let fractions = column.fractions();
                    categories
                        .iter()
                        .map(|c| {
                            fractions.get(c).copied().unwrap_or(0.0)
                                * ppkt_fractions.get(c).copied().unwrap_or(0.0)
                        })
                        .sum::<f64>()
                }
                _ => 0.0,
            };
            (gene.clone(), score)
        })
        .collect();
    scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    scored.into_iter().map(|(g, _)| g).collect()
}

/// Orders every power-set combination label according to the single-gene
/// similarity ranking, with each combo's internal gene order kept canonical
/// (matching original input index positions).
fn order_combination_labels(
    ranked_genes: &[String],
    gene_pos_map: &HashMap<String, usize>,
) -> Vec<String> {
    let mut ordered_labels = Vec::new();

    for size in 1..=ranked_genes.len() {
        for mut combo in ranked_genes.iter().cloned().combinations(size) {
            combo.sort_by_key(|g| gene_pos_map.get(g).cloned().unwrap_or(0));
            ordered_labels.push(combo.join("+"));
        }
    }

    ordered_labels
}

/// Sorts categories descending by their patient-profile ("Ppkt") fraction; ties are broken
/// by term ID so the order is the same on every run.
fn sort_categories_by_ppkt_value(
    categories: &[TermId],
    ppkt_fractions: &HashMap<TermId, f64>,
) -> Vec<TermId> {
    let mut sorted = categories.to_vec();
    sorted.sort_by(|a, b| {
        let a_val = ppkt_fractions.get(a).copied().unwrap_or(0.0);
        let b_val = ppkt_fractions.get(b).copied().unwrap_or(0.0);
        b_val
            .partial_cmp(&a_val)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.to_string().cmp(&b.to_string()))
    });
    sorted
}

/// Builds the final per-category rows, in the given sorted category order: the values the
/// graphic reads, plus the terms behind them. A category that is empty everywhere is skipped.
fn build_final_categories(
    sorted_categories: &[TermId],
    ordered_labels: &[String],
    ppkt: &CategoryAssignment,
    columns: &HashMap<String, CategoryAssignment>,
    name_of: impl Fn(&TermId) -> String,
    to_term: impl Fn(&TermId) -> SimpleOntologyTerm,
) -> Vec<SpreadPlotCategory> {
    let ppkt_fractions = ppkt.fractions();
    let column_fractions: HashMap<String, HashMap<TermId, f64>> = columns
        .iter()
        .map(|(series, a)| (series.clone(), a.fractions()))
        .collect();

    let mut final_categories = Vec::with_capacity(sorted_categories.len());
    for cat in sorted_categories {
        let ppkt_value = ppkt_fractions.get(cat).copied().unwrap_or(0.0);

        // None (JSON null) for a combination whose genes have no annotations at all
        let gene_values: Vec<Option<f64>> = ordered_labels
            .iter()
            .map(|series| {
                let column = columns.get(series)?;
                if column.n_terms == 0 {
                    return None;
                }
                column_fractions.get(series.as_str())?.get(cat).copied()
            })
            .collect();

        let max_gene_value = gene_values.iter().flatten().copied().fold(0.0, f64::max);
        if ppkt_value < EPSILON && max_gene_value < EPSILON {
            continue;
        }

        let terms_of = |a: &CategoryAssignment| -> Vec<SimpleOntologyTerm> {
            a.by_category
                .get(cat)
                .map(|ts| ts.iter().map(&to_term).collect())
                .unwrap_or_default()
        };
        final_categories.push(SpreadPlotCategory {
            id: cat.to_string(),
            name: name_of(cat),
            alias: None,
            ppkt_value,
            gene_values,
            ppkt_terms: terms_of(ppkt),
            gene_terms: ordered_labels
                .iter()
                .map(|series| SeriesTerms {
                    series: series.clone(),
                    terms: columns.get(series).map(&terms_of).unwrap_or_default(),
                })
                .collect(),
        });
    }
    final_categories
}

/// Builds the spread plot payload. The ontology enters only through the three closures, so
/// the logic can be tested without one: `is_ancestor_of(category, term)`, `name_of(category)`
/// and `to_term(term)` (the displayable term).
///
/// # Errors
/// Fails if there are no observed terms or no genotypes, if there are more than
/// [`MAX_GENOTYPES`], if a genotype is listed twice or has no gene-disease entity, or if the
/// result is internally inconsistent (see [`SpreadPlotPayload::validate`]).
pub fn spread_payload(
    observed: &[TermId],
    genotypes: &[String],
    entities: &[GeneDiseaseEntity],
    categories: &[TermId],
    is_ancestor_of: impl Fn(&TermId, &TermId) -> bool,
    name_of: impl Fn(&TermId) -> String,
    to_term: impl Fn(&TermId) -> SimpleOntologyTerm,
) -> Result<SpreadPlotPayload, String> {
    if observed.is_empty() || genotypes.is_empty() {
        return Err("Phenotypes and genotypes must not be empty.".to_string());
    }
    if genotypes.len() > MAX_GENOTYPES {
        return Err(format!(
            "At most {MAX_GENOTYPES} genotypes are supported ({} given).",
            genotypes.len()
        ));
    }
    let mut seen = HashSet::new();
    if let Some(g) = genotypes.iter().find(|g| !seen.insert(g.as_str())) {
        return Err(format!("Genotype {g} is listed more than once."));
    }
    let entity_lookup: HashMap<String, &GeneDiseaseEntity> =
        entities.iter().map(|e| (e.gene_symbol.clone(), e)).collect();
    if let Some(g) = genotypes.iter().find(|g| !entity_lookup.contains_key(*g)) {
        return Err(format!("No gene-disease entity for genotype {g}."));
    }

    let genes_power_set = genotype_power_set(genotypes);
    let gene_pos_map: HashMap<String, usize> = genotypes
        .iter()
        .enumerate()
        .map(|(i, g)| (g.clone(), i))
        .collect();

    let ppkt = assign_to_categories(observed, categories, &is_ancestor_of);
    let ppkt_fractions = ppkt.fractions();
    let columns = build_combo_column_data(&genes_power_set, &entity_lookup, categories, &is_ancestor_of);

    let ranked_genes = rank_genes_by_similarity(genotypes, &columns, &ppkt_fractions, categories);
    let ordered_labels = order_combination_labels(&ranked_genes, &gene_pos_map);
    if let Some(label) = ordered_labels.iter().find(|l| !columns.contains_key(*l)) {
        return Err(format!("Internal error: no data for series {label}."));
    }
    let sorted_categories = sort_categories_by_ppkt_value(categories, &ppkt_fractions);
    let final_categories =
        build_final_categories(&sorted_categories, &ordered_labels, &ppkt, &columns, &name_of, &to_term);

    let series_summaries: Vec<SeriesSummary> = ordered_labels
        .iter()
        .map(|series| {
            let column = &columns[series];
            SeriesSummary {
                series: series.clone(),
                n_terms: column.n_terms,
                uncategorized: column.uncategorized.iter().map(&to_term).collect(),
            }
        })
        .collect();

    let mut series_labels = vec!["Ppkt".to_string()];
    series_labels.extend(ordered_labels);

    let payload = SpreadPlotPayload {
        series_labels,
        categories: final_categories,
        ppkt_uncategorized: ppkt.uncategorized.iter().map(&to_term).collect(),
        n_ppkt_terms: ppkt.n_terms,
        series_summaries,
    };
    payload.validate()?;
    Ok(payload)
}

/// Builds the spread payload from the observed terms with the real ontology.
pub fn build_spread_payload(
    observed: &[TermId],
    genotypes: &[String],
    entities: &[GeneDiseaseEntity],
    hpo: &FullCsrOntology,
) -> Result<SpreadPlotPayload, String> {
    let categories = top_level_categories(hpo)?;
    spread_payload(
        observed,
        genotypes,
        entities,
        &categories,
        |c, t| hpo.is_ancestor_of(c, t),
        |t| hpo.term_by_id(t).map(|x| x.name().to_string()).unwrap_or_else(|| t.to_string()),
        |t| simple_term(hpo, t),
    )
}

/// Same, for callers that hold a `Proband`.
pub fn get_spread_plot_payload(
    proband: &Proband,
    genotypes: &[String],
    gd_entities: &[GeneDiseaseEntity],
    hpo: &FullCsrOntology,
) -> Result<SpreadPlotPayload, String> {
    build_spread_payload(&proband.observed_hpos, genotypes, gd_entities, hpo)
}

impl SpreadPlotPayload {
    /// Checks that the parallel parts line up and that every value matches the terms listed
    /// for it: `value * denominator = number of terms`.
    pub fn validate(&self) -> Result<(), String> {
        let series: Vec<&str> = self.series_labels.iter().skip(1).map(String::as_str).collect();

        if !self.series_summaries.is_empty() {
            let names: Vec<&str> = self.series_summaries.iter().map(|s| s.series.as_str()).collect();
            if names != series {
                return Err("series_summaries do not match series_labels".to_string());
            }
        }
        for c in &self.categories {
            if c.gene_values.len() != series.len() {
                return Err(format!("{}: {} values for {} series", c.id, c.gene_values.len(), series.len()));
            }
            if self.n_ppkt_terms > 0 {
                let expected = c.ppkt_value * self.n_ppkt_terms as f64;
                if (expected - c.ppkt_terms.len() as f64).abs() > 1e-9 {
                    return Err(format!(
                        "{}: ppktValue {} of {} terms but {} terms listed",
                        c.id, c.ppkt_value, self.n_ppkt_terms, c.ppkt_terms.len()
                    ));
                }
            }
            if c.gene_terms.is_empty() {
                continue;
            }
            let names: Vec<&str> = c.gene_terms.iter().map(|s| s.series.as_str()).collect();
            if names != series || self.series_summaries.len() != series.len() {
                return Err(format!("{}: gene_terms do not match series_labels", c.id));
            }
            for ((listed, value), summary) in c.gene_terms.iter().zip(&c.gene_values).zip(&self.series_summaries) {
                match value {
                    Some(v) => {
                        if (v * summary.n_terms as f64 - listed.terms.len() as f64).abs() > 1e-9 {
                            return Err(format!(
                                "{} / {}: value {v} of {} terms but {} terms listed",
                                c.id, listed.series, summary.n_terms, listed.terms.len()
                            ));
                        }
                    }
                    None => {
                        if summary.n_terms != 0 || !listed.terms.is_empty() {
                            return Err(format!("{} / {}: no value but terms present", c.id, listed.series));
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

/// Cross-checks the spread plot against the UpSet data and the gene list both were built from:
/// the same genes, the same number of observed terms, and every combination drawn in the UpSet
/// plot is also a column of the spread plot (same "+"-joined label, so pass the genotypes
/// sorted alphabetically, as the UpSet plot sorts them).
pub fn check_plots_agree(
    gene_symbols: &[String],
    upset: &UpsetData,
    spread: &SpreadPlotPayload,
) -> Result<(), String> {
    let mut expected: Vec<&str> = gene_symbols.iter().map(String::as_str).collect();
    expected.sort();
    expected.dedup();
    let mut upset_genes: Vec<&str> = upset.genes.iter().map(|g| g.gene_symbol.as_str()).collect();
    upset_genes.sort();
    if upset_genes != expected {
        return Err(format!("UpSet genes {upset_genes:?} differ from {expected:?}"));
    }
    if upset.n_proband_terms != spread.n_ppkt_terms {
        return Err(format!(
            "UpSet counts {} observed terms but the spread plot {}",
            upset.n_proband_terms, spread.n_ppkt_terms
        ));
    }
    let series: HashSet<&str> = spread.series_labels.iter().skip(1).map(String::as_str).collect();
    for i in &upset.intersections {
        let label = i.gene_symbol_list.join("+");
        if !series.contains(label.as_str()) {
            return Err(format!("UpSet combination {label} is missing from the spread plot"));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::{fixture, rstest};
    use std::collections::HashMap;

    // ---------------------------------------------------------------------
    // genotype_power_set
    // ---------------------------------------------------------------------

    #[fixture]
    fn three_genes() -> Vec<String> {
        vec!["GENE_A".to_string(), "GENE_B".to_string(), "GENE_C".to_string()]
    }

    #[rstest]
    fn power_set_of_empty_input_is_empty() {
        let result = genotype_power_set(&[]);
        assert!(result.is_empty());
    }

    #[rstest]
    fn power_set_of_single_gene_yields_one_combo() {
        let genes = vec!["GENE_A".to_string()];
        let result = genotype_power_set(&genes);
        assert_eq!(result, vec![vec!["GENE_A".to_string()]]);
    }

    #[rstest]
    fn power_set_has_correct_total_count(three_genes: Vec<String>) {
        // For n=3: C(3,1) + C(3,2) + C(3,3) = 3 + 3 + 1 = 7
        let result = genotype_power_set(&three_genes);
        assert_eq!(result.len(), 7);
    }

    #[rstest]
    fn power_set_contains_full_combination(three_genes: Vec<String>) {
        let result = genotype_power_set(&three_genes);
        assert!(result.contains(&three_genes));
    }

    #[rstest]
    fn power_set_contains_each_single_gene_alone(three_genes: Vec<String>) {
        let result = genotype_power_set(&three_genes);
        for gene in &three_genes {
            assert!(result.contains(&vec![gene.clone()]));
        }
    }

    #[rstest]
    fn power_set_contains_no_duplicate_combinations(three_genes: Vec<String>) {
        let result = genotype_power_set(&three_genes);
        let unique: std::collections::HashSet<_> = result.iter().cloned().collect();
        assert_eq!(unique.len(), result.len());
    }

    // ---------------------------------------------------------------------
    // order_combination_labels
    // ---------------------------------------------------------------------

    #[fixture]
    fn gene_pos_map() -> HashMap<String, usize> {
        // Original input order: GENE_A=0, GENE_B=1, GENE_C=2
        [
            ("GENE_A".to_string(), 0),
            ("GENE_B".to_string(), 1),
            ("GENE_C".to_string(), 2),
        ]
        .into_iter()
        .collect()
    }

    #[rstest]
    fn single_gene_ranking_yields_single_labels(gene_pos_map: HashMap<String, usize>) {
        let ranked = vec!["GENE_A".to_string()];
        let result = order_combination_labels(&ranked, &gene_pos_map);
        assert_eq!(result, vec!["GENE_A".to_string()]);
    }

    #[rstest]
    fn combo_labels_use_canonical_original_order_not_rank_order(
        gene_pos_map: HashMap<String, usize>,
    ) {
        // Ranked order puts GENE_C first (e.g. it scored highest), but the
        // combo label should still reflect original input order (A, B, C),
        // not the ranking order.
        let ranked = vec!["GENE_C".to_string(), "GENE_A".to_string(), "GENE_B".to_string()];
        let result = order_combination_labels(&ranked, &gene_pos_map);

        // The full 3-gene combination must be labeled in original order.
        assert!(result.contains(&"GENE_A+GENE_B+GENE_C".to_string()));
        // It must NOT appear in ranked order.
        assert!(!result.contains(&"GENE_C+GENE_A+GENE_B".to_string()));
    }

    #[rstest]
    fn total_label_count_matches_power_set_size(gene_pos_map: HashMap<String, usize>) {
        let ranked = vec!["GENE_A".to_string(), "GENE_B".to_string(), "GENE_C".to_string()];
        let result = order_combination_labels(&ranked, &gene_pos_map);
        // Same combinatorial count as genotype_power_set: 7 for n=3
        assert_eq!(result.len(), 7);
    }

    #[rstest]
    fn unknown_gene_defaults_to_position_zero(gene_pos_map: HashMap<String, usize>) {
        // A gene missing from gene_pos_map should sort as if at position 0,
        // per the `.unwrap_or(0)` fallback in the implementation.
        let ranked = vec!["GENE_B".to_string(), "UNKNOWN_GENE".to_string()];
        let result = order_combination_labels(&ranked, &gene_pos_map);

        // UNKNOWN_GENE (falls back to pos 0) should sort before GENE_B (pos 1)
        assert!(result.contains(&"UNKNOWN_GENE+GENE_B".to_string()));
    }
}

#[cfg(test)]
mod payload_tests {
    use super::*;
    use crate::blend::upset_plot::upset_data;

    fn tid(s: &str) -> TermId { s.parse().unwrap() }

    const EYE: &str = "HP:0000478";
    const NERV: &str = "HP:0000707";
    const EAR: &str = "HP:0000598";

    // ancestors among the categories: Zonular -> eye; Seizure -> nervous
    fn is_anc(c: &TermId, t: &TermId) -> bool {
        let under: &[&str] = match t.to_string().as_str() {
            "HP:0010920" => &[EYE],
            "HP:0001250" => &[NERV],
            _ => &[],
        };
        under.contains(&c.to_string().as_str())
    }
    fn name_of(t: &TermId) -> String { format!("name of {t}") }
    fn to_term(t: &TermId) -> SimpleOntologyTerm {
        SimpleOntologyTerm::from_term_id(t, format!("label of {t}"))
    }

    fn entity(symbol: &str, ids: &[&str]) -> GeneDiseaseEntity {
        GeneDiseaseEntity {
            gene_symbol: symbol.to_string(),
            disease_hpo_ids: ids.iter().map(|s| tid(s)).collect(),
            ncbi_gene_id: "0".to_string(),
            disease_list: vec![],
        }
    }

    fn payload_with(
        observed: &[&str],
        genotypes: &[&str],
        ents: &[GeneDiseaseEntity],
        categories: &[&str],
    ) -> Result<SpreadPlotPayload, String> {
        let obs: Vec<TermId> = observed.iter().map(|s| tid(s)).collect();
        let gts: Vec<String> = genotypes.iter().map(|s| s.to_string()).collect();
        let cats: Vec<TermId> = categories.iter().map(|s| tid(s)).collect();
        spread_payload(&obs, &gts, ents, &cats, is_anc, name_of, to_term)
    }
    fn payload(observed: &[&str], genotypes: &[&str], ents: &[GeneDiseaseEntity]) -> Result<SpreadPlotPayload, String> {
        payload_with(observed, genotypes, ents, &[EYE, NERV])
    }

    fn ids(v: &[SimpleOntologyTerm]) -> Vec<String> { v.iter().map(|t| t.term_id_as_string()).collect() }
    fn assert_values(actual: &[Option<f64>], expected: &[Option<f64>]) {
        assert_eq!(actual.len(), expected.len(), "{actual:?} vs {expected:?}");
        for (a, e) in actual.iter().zip(expected) {
            match (a, e) {
                (Some(a), Some(e)) => assert!((a - e).abs() < 1e-12, "{actual:?} vs {expected:?}"),
                (None, None) => {}
                _ => panic!("{actual:?} vs {expected:?}"),
            }
        }
    }
    fn two_genes() -> [GeneDiseaseEntity; 2] {
        [entity("A", &["HP:0010920", "HP:0001250"]), entity("B", &["HP:0001250", "HP:0009999"])]
    }

    #[test]
    fn values_terms_and_denominators_belong_together() {
        let p = payload(&["HP:0010920", "HP:0001250"], &["A", "B"], &two_genes()).unwrap();

        assert_eq!(p.series_labels, ["Ppkt", "A", "B", "A+B"]);
        assert_eq!(p.n_ppkt_terms, 2);
        assert!(p.ppkt_uncategorized.is_empty());
        let order: Vec<&str> = p.categories.iter().map(|c| c.id.as_str()).collect();
        assert_eq!(order, [EYE, NERV]); // tie on the Ppkt value, broken by term ID

        let eye = &p.categories[0];
        assert_eq!(eye.ppkt_value, 0.5);
        assert_eq!(ids(&eye.ppkt_terms), ["HP:0010920"]);
        assert_values(&eye.gene_values, &[Some(0.5), Some(0.0), Some(1.0 / 3.0)]);
        let per_series: Vec<(&str, Vec<String>)> =
            eye.gene_terms.iter().map(|s| (s.series.as_str(), ids(&s.terms))).collect();
        assert_eq!(
            per_series,
            [("A", vec!["HP:0010920".to_string()]), ("B", vec![]), ("A+B", vec!["HP:0010920".to_string()])]
        );

        let nerv = &p.categories[1];
        assert_values(&nerv.gene_values, &[Some(0.5), Some(0.5), Some(1.0 / 3.0)]);

        // denominators, and the term that falls under no category
        let summary: Vec<(&str, usize, Vec<String>)> = p
            .series_summaries
            .iter()
            .map(|s| (s.series.as_str(), s.n_terms, ids(&s.uncategorized)))
            .collect();
        assert_eq!(
            summary,
            [
                ("A", 2, vec![]),
                ("B", 2, vec!["HP:0009999".to_string()]),
                ("A+B", 3, vec!["HP:0009999".to_string()]),
            ]
        );
    }

    #[test]
    fn series_follow_the_similarity_ranking_but_labels_stay_canonical() {
        let ents = [entity("A", &["HP:0010920", "HP:0001250"]), entity("C", &["HP:0001250"])];
        let p = payload(&["HP:0001250"], &["A", "C"], &ents).unwrap();
        // C is nervous-system only, so it resembles the patient more than A does
        assert_eq!(p.series_labels, ["Ppkt", "C", "A", "A+C"]);
        let order: Vec<&str> = p.categories.iter().map(|c| c.id.as_str()).collect();
        assert_eq!(order, [NERV, EYE]); // sorted by the Ppkt value
    }

    #[test]
    fn inputs_that_cannot_be_plotted_are_rejected() {
        let ents = two_genes();
        let ok = ["HP:0010920"];
        assert!(payload(&[], &["A"], &ents).is_err());
        assert!(payload(&ok, &[], &ents).is_err());
        let err = payload(&ok, &["A", "NOPE"], &ents).unwrap_err();
        assert!(err.contains("NOPE"), "{err}");
        let err = payload(&ok, &["A", "A"], &ents).unwrap_err();
        assert!(err.contains("more than once"), "{err}");
        let many: Vec<String> = (0..=MAX_GENOTYPES).map(|i| format!("G{i}")).collect();
        let many_ents: Vec<GeneDiseaseEntity> = many.iter().map(|g| entity(g, &["HP:0001250"])).collect();
        let many_refs: Vec<&str> = many.iter().map(String::as_str).collect();
        assert!(payload(&ok, &many_refs, &many_ents).is_err());
    }

    #[test]
    fn a_gene_without_annotations_gives_null_values_and_the_payload_still_round_trips() {
        let ents = [entity("D", &[])];
        let p = payload(&["HP:0001250"], &["D"], &ents).unwrap();
        let nerv = &p.categories[0];
        assert_eq!(nerv.gene_values, [None]);
        assert_eq!(p.series_summaries[0].n_terms, 0);

        let v = serde_json::to_value(&p).unwrap();
        assert!(v["categories"][0]["geneValues"][0].is_null()); // what the graphic received before, too
        let back: SpreadPlotPayload = serde_json::from_value(v).unwrap();
        back.validate().unwrap();
    }

    #[test]
    fn nan_could_not_be_read_back_which_is_why_values_are_optional() {
        let v = serde_json::to_value(f64::NAN).unwrap();
        assert!(v.is_null());
        assert!(serde_json::from_value::<f64>(v).is_err());
    }

    #[test]
    fn categories_that_are_empty_everywhere_are_left_out() {
        let p = payload_with(&["HP:0010920", "HP:0001250"], &["A", "B"], &two_genes(), &[EYE, NERV, EAR]).unwrap();
        let order: Vec<&str> = p.categories.iter().map(|c| c.id.as_str()).collect();
        assert_eq!(order, [EYE, NERV]);
    }

    #[test]
    fn output_does_not_depend_on_entity_order() {
        let [a, b] = two_genes();
        let fwd = payload(&["HP:0010920", "HP:0001250"], &["A", "B"], &[a.clone(), b.clone()]).unwrap();
        let rev = payload(&["HP:0010920", "HP:0001250"], &["A", "B"], &[b, a]).unwrap();
        assert_eq!(serde_json::to_value(&fwd).unwrap(), serde_json::to_value(&rev).unwrap());
    }

    #[test]
    fn validate_catches_values_that_do_not_match_their_terms() {
        let good = payload(&["HP:0010920", "HP:0001250"], &["A", "B"], &two_genes()).unwrap();
        good.validate().unwrap();

        let mut wrong_value = good.clone();
        wrong_value.categories[0].gene_values[0] = Some(0.99);
        assert!(wrong_value.validate().is_err());

        let mut missing_term = good.clone();
        missing_term.categories[0].gene_terms[0].terms.clear();
        assert!(missing_term.validate().is_err());

        let mut wrong_denominator = good.clone();
        wrong_denominator.n_ppkt_terms = 5;
        assert!(wrong_denominator.validate().is_err());

        let mut misaligned = good.clone();
        misaligned.categories[0].gene_terms.swap(0, 1);
        assert!(misaligned.validate().is_err());
    }

    #[test]
    fn plots_agree_or_the_difference_is_named() {
        let ents = two_genes();
        let obs: Vec<TermId> = ["HP:0010920", "HP:0001250"].iter().map(|s| tid(s)).collect();
        let upset = upset_data(&obs, &ents, |t| vec![t.clone()], to_term);
        let genes = vec!["A".to_string(), "B".to_string()];

        let spread = payload(&["HP:0010920", "HP:0001250"], &["A", "B"], &ents).unwrap();
        check_plots_agree(&genes, &upset, &spread).unwrap();

        // genotypes passed in a different order give "B+A", which the UpSet plot calls "A+B"
        let unsorted = payload(&["HP:0010920", "HP:0001250"], &["B", "A"], &ents).unwrap();
        let err = check_plots_agree(&genes, &upset, &unsorted).unwrap_err();
        assert!(err.contains("A+B"), "{err}");

        // a different number of observed terms
        let one: Vec<TermId> = vec![tid("HP:0010920")];
        let upset_one = upset_data(&one, &ents, |t| vec![t.clone()], to_term);
        assert!(check_plots_agree(&genes, &upset_one, &spread).is_err());

        // a different gene list
        assert!(check_plots_agree(&["A".to_string()], &upset, &spread).is_err());
    }
}
