
use ontolius::TermId;
use serde::{Serialize, Deserialize};


use crate::{ util::serde_util::serialize_term_id_vec};







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
}

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
}



