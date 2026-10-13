use std::{collections::HashMap, sync::{Arc, LazyLock}};
use ontolius::{TermId, ontology::{OntologyTerms, csr::FullCsrOntology}, term::MinimalTerm};
use serde::{Serialize, Deserialize};

use std::collections::BTreeMap;

use crate::{blend::dto::MatchType::NoMatch, model::simple_term::SimpleOntologyTerm};


/// Root of the HPO phenotype branch; sharing only this term means no meaningful match.
pub const ROOT_PHENOTYPE_ID: &str = "HP:0000118";
pub static ROOT_PHENOTYPE_TERM_ID: std::sync::LazyLock<TermId> = std::sync::LazyLock::new(|| {
    let tid: TermId  = ROOT_PHENOTYPE_ID.parse().unwrap();
    tid
});
pub static SIMPLE_ROOT_TERM: LazyLock<SimpleOntologyTerm> = LazyLock::new(|| {
    SimpleOntologyTerm::new(ROOT_PHENOTYPE_ID, "Phenotypic abnormality").unwrap()
});


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

    pub fn no_match() -> Self {
        Self {
            match_type: NoMatch,
            score: 0.0,
            matched_term: SIMPLE_ROOT_TERM.clone(),
        }
    }

    pub fn matched_term(&self) -> SimpleOntologyTerm {
        self.matched_term.clone()
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

    pub fn proband_hpo_term(&self) -> SimpleOntologyTerm {
        self.term.clone()
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
        Self {
            entities,
            columns: term_matches.iter().map(TermMatches::to_overlap_item).collect(),
        }
    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetPlotPayload {
    /// The unique gene symbols being evaluated (e.g., ["MBNL1", "DMPK"])
    pub genes: Vec<String>,
    
    /// List of combinations, where each combination is a subset of gene symbols
    pub combinations: Vec<Vec<String>>,
    
    /// Total annotated HPO term intersection sizes per combination (Top-Right background bars)
    pub combination_annotated: Vec<u32>,
    
    /// Subsets of those terms overlapping with patient phenotypes (Top-Right foreground bars)
    pub combination_observed: Vec<u32>,
    
    /// Total annotated HPO term sizes per individual gene (Bottom-Left background bars)
    pub gene_annotated: Vec<u32>,
    
    /// Total terms per individual gene overlapping with patient phenotypes (Bottom-Left foreground bars)
    pub gene_observed: Vec<u32>,
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
    /// Array of values for each gene combo/series (`geneValues`)
    pub gene_values: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SpreadPlotPayload {
    /// e.g., ["Ppkt", "MBNL1", "DMPK", "MBNL1+DMPK"]
    pub series_labels: Vec<String>,
    pub categories: Vec<SpreadPlotCategory>,
}




#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct CohortAnalysis {
    /// Source of truth: how each observed proband term matches each gene entity.
    pub term_matches: Vec<TermMatches>,
    /// All annotated terms per gene entity (needed later for the UpSet "annotated" bars).
    pub gene_annotations: BTreeMap<String, Vec<SimpleOntologyTerm>>,
    pub overlap: OverlapPlotPayload,
    // pub upset: UpsetPlotPayload,   // step 2
    // pub spread: SpreadPlotPayload, // step 2
}