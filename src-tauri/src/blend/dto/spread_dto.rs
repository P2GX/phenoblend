use std::collections::HashMap;

use ontolius::TermId;
use serde::{Deserialize, Serialize};

use crate::model::simple_term::SimpleOntologyTerm;

/// The terms behind one column (gene combination) of one category.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SeriesTerms {
    /// Series label, e.g. "NF1+PTPN11"; the same string as in `series_labels`.
    pub series: String,
    /// The combination's HPO terms (union of its genes' explicit annotations) that fall under
    /// this category, sorted by term ID. The category's `geneValues` entry for this series is
    /// their count divided by the series' `nTerms` (see `SeriesSummary`).
    pub terms: Vec<SimpleOntologyTerm>,
}

/// Per-series facts needed to read the values: the denominator and what fell outside.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SeriesSummary {
    pub series: String,
    /// Distinct HPO terms in the combination's union: the denominator of its values.
    pub n_terms: usize,
    /// Terms of the union that fall under none of the categories.
    pub uncategorized: Vec<SimpleOntologyTerm>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SpreadPlotCategory {
    pub id: String,
    pub name: String,
    /// Optional field mapping to `alias?: string`
    pub alias: Option<String>,
    /// Patient/Phenopacket fraction (`ppktValue`)
    pub ppkt_value: f64,
    /// One value per gene combination/series (`geneValues`), in `series_labels[1..]` order.
    /// `None` (JSON `null`) for a combination whose genes have no annotations at all.
    /// An `Option` (not NaN) so that saved analyses can be loaded again: `null` cannot be
    /// read back into an `f64`.
    pub gene_values: Vec<Option<f64>>,
    /// Debugging detail: the proband's observed terms under this category (the term itself
    /// or a descendant), sorted by term ID.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub ppkt_terms: Vec<SimpleOntologyTerm>,
    /// Debugging detail: for every series, the terms behind its `geneValues` entry.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub gene_terms: Vec<SeriesTerms>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SpreadPlotPayload {
    /// e.g., ["Ppkt", "MBNL1", "DMPK", "MBNL1+DMPK"]
    pub series_labels: Vec<String>,
    pub categories: Vec<SpreadPlotCategory>,
    /// Debugging detail: observed terms that fall under none of the categories.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub ppkt_uncategorized: Vec<SimpleOntologyTerm>,
    /// Debugging detail: number of distinct observed terms (the denominator of `ppktValue`).
    #[serde(default)]
    pub n_ppkt_terms: usize,
    /// Debugging detail: denominator and uncategorized terms for each gene-combination series.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub series_summaries: Vec<SeriesSummary>,
}

/// How a set of HPO terms falls into the top-level categories (organ systems). Used for the
/// proband's observed terms and for each gene combination's union of annotations.
#[derive(Debug, Clone, Default)]
pub struct CategoryAssignment {
    /// Category -> the terms under it (the term itself or a descendant), sorted by term ID.
    /// Every category is present, possibly with an empty list.
    pub by_category: HashMap<TermId, Vec<TermId>>,
    /// Terms that fall under none of the categories, sorted by term ID.
    pub uncategorized: Vec<TermId>,
    /// Number of distinct terms (the denominator of the fractions).
    pub n_terms: usize,
}

/// Splits `terms` (de-duplicated) among `categories`. A term under several categories is
/// listed under each of them, so the fractions can sum to more than 1.
pub fn assign_to_categories(
    terms: &[TermId],
    categories: &[TermId],
    is_ancestor_of: impl Fn(&TermId, &TermId) -> bool,
) -> CategoryAssignment {
    let mut terms: Vec<TermId> = terms.to_vec();
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
    CategoryAssignment { by_category, uncategorized, n_terms: terms.len() }
}

impl CategoryAssignment {
    /// Fraction of the terms under each category; all 0.0 when there are no terms.
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