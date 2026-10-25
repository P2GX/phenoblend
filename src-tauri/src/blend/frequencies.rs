//! Annotation frequencies of a gene's disease terms, for the frequency-weighted test.
//!
//! HPOA gives, for some annotations, the proportion of individuals with the disease who have
//! the feature. By the convention of the HPO project, an annotation without an explicit
//! frequency is taken to mean 100%, so every annotated term has a frequency in `(0, 1]`.
use std::collections::{HashMap, HashSet};

use ontolius::TermId;

use crate::hpoa::disease_model::{GeneDiseaseAssociation, SimpleDiseaseModel};



/// The frequency of `term` in `model`: the reported value, or 1.0 if none is reported.
///
/// A reported value that is not a number or not positive cannot be a frequency of an annotated
/// term (a 0% frequency means the term is *excluded*, which is stored separately), so it is
/// treated like a missing one. Values above 1 are capped at 1.
pub fn frequency_in_model(model: &SimpleDiseaseModel, term: &TermId) -> f64 {
    model
        .term_frequencies
        .get(term)
        .copied()
        .filter(|f| f.is_finite() && *f > 0.0)
        .map_or(1.0, |f| f.min(1.0))
}

/// For every term annotated to at least one of the gene's diseases, the frequency in each
/// disease that annotates it (one value per disease, in no particular order).
///
/// `gdas` holds all associations of one gene. A disease that appears more than once (for
/// example under two association types) counts once, and associations without a disease
/// model are skipped.
pub fn collect_frequencies(gdas: &[GeneDiseaseAssociation]) -> HashMap<TermId, Vec<f64>> {
    let mut frequencies: HashMap<TermId, Vec<f64>> = HashMap::new();
    let mut seen: HashSet<&TermId> = HashSet::new();
    for gda in gdas {
        let Some(model) = gda.disease_model.as_ref() else {
            continue;
        };
        if !seen.insert(&gda.disease_id) {
            continue; // this disease was already counted
        }
        for term in &model.observed_hpo_ids {
            frequencies
                .entry(term.clone())
                .or_default()
                .push(frequency_in_model(model, term));
        }
    }
    frequencies
}

#[cfg(test)]
mod tests {
    use crate::hpoa::disease_model::GeneDiseaseAssociation;

use super::*;
    use rstest::rstest;

    fn t(s: &str) -> TermId {
        s.parse().unwrap()
    }

    /// A disease annotating the given terms; `freqs` are the explicit frequencies.
    fn gda(disease: &str, terms: &[&str], freqs: &[(&str, f64)]) -> GeneDiseaseAssociation {
        let tid = "OMIM:text".parse().unwrap();
        GeneDiseaseAssociation {
            gene_symbol: "G".to_string(),
            disease_id: t(disease),
            disease_model: Some(SimpleDiseaseModel {
                observed_hpo_ids: terms.iter().map(|x| t(x)).collect(),
                term_frequencies: freqs.iter().map(|(x, f)| (t(x), *f)).collect(),
                omim_disease_id: tid,
                omim_disease_name: "text disease label".to_string(),
                excluded_hpo_ids: HashSet::default(),
            }),
            ncbi_gene_id: String::default(),
            association_type: String::default(),
            source: String::default(),
        }
    }

    #[test]
    fn explicit_frequency_is_recorded() {
        let f = collect_frequencies(&[gda("OMIM:1", &["HP:1"], &[("HP:1", 0.4)])]);
        assert_eq!(f[&t("HP:1")], [0.4]);
    }

    #[test]
    fn a_missing_frequency_means_one_hundred_percent() {
        let f = collect_frequencies(&[gda("OMIM:1", &["HP:1", "HP:2"], &[("HP:1", 0.4)])]);
        assert_eq!(f[&t("HP:2")], [1.0]);
    }

    #[test]
    fn each_disease_contributes_one_value_per_term() {
        let f = collect_frequencies(&[
            gda("OMIM:1", &["HP:1"], &[("HP:1", 0.4)]),
            gda("OMIM:2", &["HP:1", "HP:2"], &[]),
        ]);
        let mut hp1 = f[&t("HP:1")].clone();
        hp1.sort_by(f64::total_cmp);
        assert_eq!(hp1, [0.4, 1.0]);
        assert_eq!(f[&t("HP:2")], [1.0]); // only the second disease annotates it
    }

    #[test]
    fn a_disease_listed_twice_counts_once() {
        let f = collect_frequencies(&[
            gda("OMIM:1", &["HP:1"], &[("HP:1", 0.4)]),
            gda("OMIM:1", &["HP:1"], &[("HP:1", 0.4)]),
        ]);
        assert_eq!(f[&t("HP:1")], [0.4]);
    }

    #[test]
    fn associations_without_a_disease_model_are_skipped() {
        let no_model = GeneDiseaseAssociation {
            gene_symbol: "G".to_string(),
            disease_id: t("OMIM:3"),
            disease_model: None,
            ncbi_gene_id: String::default(),
            association_type: String::default(),
            source: String::default(),
        };
        let f = collect_frequencies(&[no_model, gda("OMIM:1", &["HP:1"], &[])]);
        assert_eq!(f.len(), 1);
    }

    #[test]
    fn terms_that_are_only_in_the_frequency_map_are_not_annotations() {
        // HP:9 has a frequency entry but is not among the observed (annotated) terms.
        let f = collect_frequencies(&[gda("OMIM:1", &["HP:1"], &[("HP:9", 0.2)])]);
        assert!(!f.contains_key(&t("HP:9")));
    }

    #[rstest]
    #[case(0.25, 0.25)]
    #[case(1.0, 1.0)]
    #[case(1.5, 1.0)] // capped
    #[case(0.0, 1.0)] // a 0% frequency cannot belong to an annotated term
    #[case(-0.1, 1.0)]
    #[case(f64::NAN, 1.0)]
    #[case(f64::INFINITY, 1.0)]
    fn unusable_frequencies_are_treated_as_missing(#[case] reported: f64, #[case] expected: f64) {
        let g = gda("OMIM:1", &["HP:1"], &[("HP:1", reported)]);
        let f = collect_frequencies(&[g]);
        assert_eq!(f[&t("HP:1")], [expected]);
    }

    #[test]
    fn no_associations_gives_an_empty_map() {
        assert!(collect_frequencies(&[]).is_empty());
    }
}
