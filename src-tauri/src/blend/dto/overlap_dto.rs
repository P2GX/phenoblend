//! Overlap
//! Test the overlap of HPO terms of a phenopacket with the distribution of HPO terms in the disease model.
//! 
//! 
use std::{collections::{BTreeMap, HashMap}, sync::Arc};

use serde::{Deserialize,Serialize};
use ontolius::{TermId, ontology::{OntologyTerms, csr::FullCsrOntology}, term::MinimalTerm};
use crate::{blend::dto::constants::ROOT_PHENOTYPE_ID, model::simple_term::SimpleOntologyTerm};


/// How a gene-disease entity relates to one proband term.
#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "camelCase")]
pub enum MatchType {
    /// The disease is annotated to exactly the proband term.
    Exact,
    /// The disease is annotated to one or more descendants of the proband term
    /// (true path rule, e.g. proband "Cataract", disease "Nuclear cataract").
    Descendant,
    /// No perfect match; `score` is a partial-similarity score. 
    Partial,
    /// Only the root term is shared. Score=0.0 means match at root (no match!)
    NoMatch,
}

impl MatchType {
    /// True if the disease annotation accounts for the proband term
    /// (`Exact` or `Descendant`). `Partial` and `NoMatch` do not count.
    pub fn is_match(self) -> bool {
        matches!(self, MatchType::Exact | MatchType::Descendant)
    }

     pub fn is_partial_match(self) -> bool {
        matches!(self, MatchType::Partial)
    }
}


/// The outcome of matching one proband term against one gene-disease entity.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct GeneMatch {
    pub match_type: MatchType,
    /// 1.0 for `Exact` and `Descendant`; in `[0.0, 1.0]` for `Partial`.
    pub score: f64,
    /// The single term responsible for the match: the identical disease term for `Exact`;
    /// the closest matching descendant for `Descendant` (ties broken by term ID);
    /// the best shared term for `Partial`; `HP:0000118` for `None`.
    pub matched_term: SimpleOntologyTerm,
}


impl GeneMatch {
    pub fn exact(term: SimpleOntologyTerm) -> Self {
        Self { match_type: MatchType::Exact, score: 1.0, matched_term: term }
    }

    pub fn descendant(term: SimpleOntologyTerm) -> Self {
        Self { match_type: MatchType::Descendant, score: 1.0, matched_term: term }
    }

    /// Builds a partial match, or a `None` match when the shared term is the root.
    ///
    /// # Errors
    /// Returns an error if `score` is not finite or lies outside `[0, 1]`, or if the
    /// term is the root but the score is not 0 (a contradiction).
    pub fn partial(score: f64, term: SimpleOntologyTerm) -> Result<Self, String> {
        if !score.is_finite() || !(0.0..=1.0).contains(&score) {
            return Err(format!("partial match score out of range: {score}"));
        }
        if term.term_id_as_string() == ROOT_PHENOTYPE_ID {
            if score != 0.0 {
                return Err(format!("root term {ROOT_PHENOTYPE_ID} with non-zero score {score}"));
            }
            return Ok(Self { match_type: MatchType::NoMatch, score: 0.0, matched_term: term });
        }
        Ok(Self { match_type: MatchType::Partial, score, matched_term: term })
    }

    pub fn is_match(&self) -> bool {
        self.match_type.is_match()
    }

    pub fn is_partial_match(&self) -> bool {
        self.match_type.is_partial_match()
    }
}


/// All match results for one observed proband term, across all gene entities.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct TermMatches {
    pub term: SimpleOntologyTerm,
    /// Gene symbol -> match
    pub matches: BTreeMap<String, GeneMatch>,
}

impl TermMatches {
    pub fn new(term: SimpleOntologyTerm) -> Self {
        Self { term, matches: BTreeMap::new() }
    }

    pub fn insert(&mut self, gene_symbol: &str, m: GeneMatch) {
        self.matches.insert(gene_symbol.to_string(), m);
    }

    /// The wire format your overlap graphic already reads (`hpoId`, `hpoName`, `scores`),
    /// so the frontend contract doesn't change.
    pub fn to_overlap_item(&self) -> OverlapPlotItem {
        OverlapPlotItem {
            hpo_id: self.term.term_id_as_string(),
            hpo_name: self.term.term_label.clone(),
            scores: self.matches.iter().map(|(g, m)| (g.clone(), m.score)).collect(),
        }
    }
}

/// Represents the scores for one HPO term in the presence matrix
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct OverlapPlotItem {
    pub hpo_id: String,
    pub hpo_name: String,
    /// A map matching gene symbols to their calculated matrix match scores [0.0 - 1.0]
    pub scores: HashMap<String, f64>,
}

impl OverlapPlotItem {
    pub fn new(hpo_id: &TermId, hpo: Arc<FullCsrOntology>) -> Self {
        let term = hpo.term_by_id(hpo_id);
        let label = match term {
            Some(t) => t.name().to_string(),
            None => "n/a".to_string() // should never happen
        };
        println!("OverlapPlotItem::new {}-{}", label, hpo_id);
        Self {
            hpo_id: hpo_id.to_string(),
            hpo_name: label,
            scores: HashMap::default()
        }
    }

    pub fn add_perfect_match(&mut self, entity: &str) {
        self.scores.insert(entity.to_string(), 1.0);
    }

    pub fn add_match(&mut self, entity: &str, score: f64) {
        self.scores.insert(entity.to_string(), score);
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct OverlapPlotPayload {
    /// Ordered array of gene symbols or diseases following column sorting block rules
    pub entities: Vec<String>,
    /// Ordered rows following phenotype grouping tiers
    pub columns: Vec<OverlapPlotItem>,
}


impl OverlapPlotPayload {
    pub fn n_columns(&self) -> usize {
        self.columns.len()
    }

   
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct OverlapData {
     /// Ordered array of gene symbols or diseases following column sorting block rules
    pub entities: Vec<String>,
    /// Ordered rows following phenotype grouping tiers
    pub rows: Vec<TermMatches>,
     /// `n_matched_terms[i]` is the number of rows (proband terms) that have a
    /// match for `entities[i]`. Same length and order as `entities`.
    pub n_matched_terms: Vec<usize>,
     /// `n_matched_terms[i]` is the number of rows (proband terms) that have a
    /// partial match for `entities[i]`.
    pub n_partially_matched_terms: Vec<usize>
}




impl OverlapData {
    pub fn to_payload(&self) -> OverlapPlotPayload {
        let overlap_list: Vec<OverlapPlotItem> = self.rows.iter()
            .map(|entry| entry.to_overlap_item())
            .collect();
        OverlapPlotPayload {
            entities: self.entities.clone(),
            columns: overlap_list
        }
    }

     /// Derives the overlap graphic's payload from the stored match records.
    ///
    /// `entities` is the sorted, de-duplicated list of gene symbols; `columns` has one
    /// entry per record, in the same order.
    pub fn from_term_matches<I>(gene_symbols: I, term_matches: &[TermMatches]) -> Self
    where
        I: IntoIterator<Item = String>,
    {
        let mut entities: Vec<String> = gene_symbols.into_iter().collect();
        entities.sort();
        entities.dedup();
        let mut n_matched_terms = vec![0usize; entities.len()];
        let mut n_partially_matched_terms = vec![0usize; entities.len()];
        for tm in term_matches {
            for (count, gene) in n_matched_terms.iter_mut().zip(&entities) {
                if tm.matches.get(gene).is_some_and(GeneMatch::is_match) {
                    *count += 1;
                } 
            }
            for (count, gene) in n_partially_matched_terms.iter_mut().zip(&entities) {
                if tm.matches.get(gene).is_some_and(GeneMatch::is_partial_match) {
                    *count += 1;
                }
            }
        }
        Self {
            entities,
            rows: term_matches.to_vec(),
            n_matched_terms,
            n_partially_matched_terms
        }
    }
}
