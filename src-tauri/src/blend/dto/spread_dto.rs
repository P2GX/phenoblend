use std::collections::HashMap;

use serde::{Deserialize,Serialize};
use ontolius::TermId;
use crate::model::simple_term::SimpleOntologyTerm;


#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SpreadPlotCategory {
    pub id: String,
    pub name: String,
    /// Optional field mapping to `alias?: string`
    pub alias: Option<String>,
    /// Patient/Phenopacket fraction (`ppktValue`)
    pub ppkt_value: f64,
    /// Array of values for each gene combo/series (`geneValues`)
    pub gene_values: Vec<f64>,
    /// The proband's observed terms under this category (the term itself
    /// or a descendant), sorted by term ID, intended to show details of spread matches in Notebook
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub ppkt_terms: Vec<SimpleOntologyTerm>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SpreadPlotPayload {
    /// e.g., ["Ppkt", "MBNL1", "DMPK", "MBNL1+DMPK"]
    pub series_labels: Vec<String>,
    pub categories: Vec<SpreadPlotCategory>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub ppkt_uncategorized: Vec<SimpleOntologyTerm>,
    /// Total count of distinct observed terms (the denominator of `ppktValue`).
    #[serde(default)]
    pub n_ppkt_terms: usize,
}



#[derive(Debug, Clone, Default)]
pub struct PpktCategoryAssignment {
    pub by_category: HashMap<TermId, Vec<TermId>>,
    pub uncategorized: Vec<TermId>,
    pub n_terms: usize,
}
 
pub fn assign_to_categories(
    observed: &[TermId],
    categories: &[TermId],
    is_ancestor_of: impl Fn(&TermId, &TermId) -> bool,
) -> PpktCategoryAssignment {
    let mut terms: Vec<TermId> = observed.to_vec();
    terms.sort_by_key(|t| t.to_string());
    terms.dedup();
 
    let mut by_category: HashMap<TermId, Vec<TermId>> =
        categories.iter().map(|c| (c.clone(), Vec::new())).collect();
    let mut uncategorized = Vec::new();
    for t in &terms {
        let mut placed = false;
        for c in categories {
            if c == t || is_ancestor_of(c, t) {
                by_category.entry(c.clone()).or_default().push(t.clone());
                placed = true;
            }
        }
        if !placed {
            uncategorized.push(t.clone());
        }
    }
    PpktCategoryAssignment { by_category, uncategorized, n_terms: terms.len() }
}
 
impl PpktCategoryAssignment {
    pub fn fractions(&self) -> HashMap<TermId, f64> {
        let total = self.n_terms as f64;
        self.by_category
            .iter()
            .map(|(c, ts)| {
                let f = if self.n_terms == 0 { 0.0 } else { ts.len() as f64 / total };
                (c.clone(), f)
            })
            .collect()
    }
}



#[cfg(test)]
mod tests {
    use ontolius::TermId;

use super::*;
 
    fn tid(s: &str) -> TermId { s.parse().unwrap() }
    fn ids(v: &[TermId]) -> Vec<String> { v.iter().map(|t| t.to_string()).collect() }
 
    const EYE: &str = "HP:0000478";
    const NERV: &str = "HP:0000707";
 
    // ancestors among the categories: Zonular -> eye; Seizure -> nervous; HP:0005555 -> both
    fn is_anc(c: &TermId, t: &TermId) -> bool {
        let under: &[&str] = match t.to_string().as_str() {
            "HP:0010920" => &[EYE],
            "HP:0001250" => &[NERV],
            "HP:0005555" => &[EYE, NERV],
            _ => &[],
        };
        under.contains(&c.to_string().as_str())
    }
 
    fn assign(observed: &[&str]) -> PpktCategoryAssignment {
        let obs: Vec<TermId> = observed.iter().map(|s| tid(s)).collect();
        assign_to_categories(&obs, &[tid(EYE), tid(NERV)], is_anc)
    }
 
    #[test]
    fn terms_are_listed_per_category_and_uncategorized_ones_are_kept() {
        let a = assign(&["HP:0010920", "HP:0001250", "HP:0009999", "HP:0005555"]);
        assert_eq!(ids(&a.by_category[&tid(EYE)]), ["HP:0005555", "HP:0010920"]);
        assert_eq!(ids(&a.by_category[&tid(NERV)]), ["HP:0001250", "HP:0005555"]);
        assert_eq!(ids(&a.uncategorized), ["HP:0009999"]);
        assert_eq!(a.n_terms, 4);
        let f = a.fractions();
        assert_eq!(f[&tid(EYE)], 0.5);
        assert_eq!(f[&tid(NERV)], 0.5);
    }
 
    #[test]
    fn a_category_term_counts_under_itself() {
        let a = assign(&[EYE]);
        assert_eq!(ids(&a.by_category[&tid(EYE)]), [EYE]);
        assert!(a.uncategorized.is_empty());
    }
 
    #[test]
    fn duplicates_are_counted_once() {
        let a = assign(&["HP:0001250", "HP:0001250"]);
        assert_eq!(a.n_terms, 1);
        assert_eq!(a.fractions()[&tid(NERV)], 1.0);
    }
 
    #[test]
    fn empty_input_gives_zero_fractions_not_nan() {
        let a = assign(&[]);
        assert_eq!(a.n_terms, 0);
        assert!(a.fractions().values().all(|f| *f == 0.0));
        assert_eq!(a.by_category.len(), 2); // every category is present, with no terms
    }
 
    #[test]
    fn every_term_is_either_categorized_or_listed_as_uncategorized() {
        let observed = ["HP:0010920", "HP:0001250", "HP:0009999", "HP:0005555"];
        let a = assign(&observed);
        let mut seen: Vec<String> = a.by_category.values().flatten().map(|t| t.to_string()).collect();
        seen.extend(ids(&a.uncategorized));
        seen.sort();
        seen.dedup();
        assert_eq!(seen.len(), observed.len());
    }
 
    #[test]
    fn new_fields_are_additive_in_json() {
        let c = SpreadPlotCategory {
            id: "x".into(), name: "X".into(), alias: None, ppkt_value: 0.5,
            gene_values: vec![0.1],
            ppkt_terms: vec![SimpleOntologyTerm::new("HP:0010920", "Zonular cataract").unwrap()],
        };
        let v = serde_json::to_value(&c).unwrap();
        assert_eq!(v["ppktValue"], 0.5);                       // existing field unchanged
        assert_eq!(v["ppktTerms"][0]["termId"], "HP:0010920"); // new detail
        // JSON written before the new fields existed still loads
        let old = r#"{"id":"x","name":"X","alias":null,"ppktValue":0.5,"geneValues":[0.1]}"#;
        let back: SpreadPlotCategory = serde_json::from_str(old).unwrap();
        assert!(back.ppkt_terms.is_empty());
    }
}