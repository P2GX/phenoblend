use std::collections::BTreeMap;
use serde::{Deserialize, Serialize};
use crate::{blend::{dto::{overlap_dto::{OverlapData, OverlapPlotPayload, TermMatches}, spread_dto::SpreadPlotPayload, upset_dto::UpsetPlotPayload}, upset_plot::UpsetData}, model::simple_term::SimpleOntologyTerm};



#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct PhenoBlendAnalysis {
    /// Source of truth: how each observed proband term matches each gene entity.
    pub term_matches: Vec<TermMatches>,
    /// All annotated terms per gene entity (needed later for the UpSet "annotated" bars).
    pub gene_annotations: BTreeMap<String, Vec<SimpleOntologyTerm>>,
    pub overlap: OverlapData,
    pub upset: UpsetData,   
    pub spread: SpreadPlotPayload, 
}


impl PhenoBlendAnalysis {
    pub fn new(
        term_matches: Vec<TermMatches>,
        annotations: BTreeMap<String, Vec<SimpleOntologyTerm>>,
        overlap: OverlapData,
        upset: UpsetData,
        spread: SpreadPlotPayload) -> Self {
            Self { 
                term_matches, 
                gene_annotations: annotations, 
                overlap, 
                upset,
                spread
            }
    } 

    pub fn upset_payload(&self) -> Result<UpsetPlotPayload, String> {
        self.upset.build_upset_payload()
    }

    pub fn overlap_payload(&self) -> Result<OverlapPlotPayload, String> {
        Ok(self.overlap.to_payload())
    }

}

