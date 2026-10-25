use std::collections::BTreeMap;
use serde::{Deserialize, Serialize};
use crate::{blend::dto::{UpsetPlotData, overlap_dto::{OverlapData, OverlapPlotPayload, TermMatches}, shared_term_enrichment::SharedTermEnrichment, spread_dto::SpreadPlotPayload}, model::simple_term::SimpleOntologyTerm};



#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct PhenoBlendAnalysis {
    /// Source of truth: how each observed proband term matches each gene entity.
    pub term_matches: Vec<TermMatches>,
    /// All annotated terms per gene entity (needed later for the UpSet "annotated" bars).
    pub gene_annotations: BTreeMap<String, Vec<SimpleOntologyTerm>>,
    pub overlap: OverlapData,
    pub upset: UpsetPlotData,   
    pub spread: SpreadPlotPayload, 
    /// Exact enrichment of shared terms: every pair of genes, plus all the genes together for
    /// three or four. See [`SharedTermEnrichment`].
    pub enrichment: Vec<SharedTermEnrichment>,
}


impl PhenoBlendAnalysis {
    pub fn new(
        term_matches: Vec<TermMatches>,
        annotations: BTreeMap<String, Vec<SimpleOntologyTerm>>,
        overlap: OverlapData,
        upset: UpsetPlotData,
        spread: SpreadPlotPayload,
        enrichment: Vec<SharedTermEnrichment>
    ) -> Self {
            Self { 
                term_matches, 
                gene_annotations: annotations, 
                overlap, 
                upset,
                spread,
                enrichment
            }
    } 

    pub fn upset_payload(&self) -> Result<UpsetPlotData, String> {
    Ok(self.upset.clone())
}

    pub fn overlap_payload(&self) -> Result<OverlapPlotPayload, String> {
        Ok(self.overlap.to_payload())
    }

}

