
use ontolius::TermId;
use serde::{Serialize, Deserialize};


use crate::{ model::simple_term::SimpleOntologyTerm, util::serde_util::serialize_term_id_vec};







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


/* 
#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetGene {
    /// Gene symbol of gene associated with a disease
    gene_symbol: String,
    /// HPO terms observed in proband and also associated with gene (explicit only)
    #[serde(serialize_with = "serialize_term_id_vec")]
    proband_overlapped_tids: Vec<TermId>,
    /// HPO Terms annotation to the gene but not reported in the proband (explicit only)
    #[serde(serialize_with = "serialize_term_id_vec")]
    remaining_annotated_tids: Vec<TermId>
}*/

/// One **horizontal** bar: all explicit annotations of one disease/gene entity, including
/// terms it shares with the others.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetGene {
    pub gene_symbol: String,
    #[serde(flatten)]
    pub terms: TermSplit,
}

/*
#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetIntersection {
    /// Gene symbol or symbols for intersection (on vertical of upset plot)
    gene_symbol_list: Vec<String>,
     /// HPO terms observed in proband and also associated with gene (explicit only)
    #[serde(serialize_with = "serialize_term_id_vec")]
    proband_overlapped_tids: Vec<TermId>,
    /// HPO Terms annotation to the gene but not reported in the proband (explicit only)
    #[serde(serialize_with = "serialize_term_id_vec")]
    remaining_annotated_tids: Vec<TermId>
} */


/// `gene_symbol_list` and to no other entity. Each annotation appears in at most one
/// intersection.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetIntersection {
    pub gene_symbol_list: Vec<String>,
    #[serde(flatten)]
    pub terms: TermSplit,
}

/// A disease annotation that the proband has, with the proband term that accounts for it.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct MatchedTerm {
    /// The disease's own (explicit) annotation, as shown in the bar.
    pub term: SimpleOntologyTerm,
    /// The proband term that accounts for it: the same term if the proband has it, otherwise
    /// a descendant of it (the smallest term ID if there are several). Observed terms are
    /// non-redundant, so an exact match is never listed together with descendants.
    pub matched_by: SimpleOntologyTerm,
}

/// A group of disease annotations split by whether the proband has them.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct TermSplit {
    /// Green: annotations the proband has (exactly, or through a descendant).
    pub proband_overlapped: Vec<MatchedTerm>,
    /// Purple: annotations the proband does not have.
    pub remaining_annotated: Vec<SimpleOntologyTerm>,
}

impl TermSplit {
    /// Height of the green part.
    pub fn n_observed(&self) -> u32 {
        self.proband_overlapped.len() as u32
    }

    /// Height of the whole bar (purple + green).
    pub fn n_annotated(&self) -> u32 {
        (self.proband_overlapped.len() + self.remaining_annotated.len()) as u32
    }
}



/// One **vertical** bar: the annotations that belong to exactly the entities in


/// Everything behind one UpSet plot.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetData {
    /// Horizontal bars, sorted by gene symbol.
    pub genes: Vec<UpsetGene>,
    /// Vertical bars, sorted by number of entities and then alphabetically.
    pub intersections: Vec<UpsetIntersection>,
    /// Proband terms that no disease annotation accounts for (no disease term equals the
    /// term or is one of its ancestors). They are not drawn in the plot.
    pub proband_only: Vec<SimpleOntologyTerm>,
    /// Number of distinct observed terms of the proband.
    pub n_proband_terms: usize,
}