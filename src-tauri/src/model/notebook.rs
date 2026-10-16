//! Notebook
//! Keep track of analysis run
//! Intended to be displayed as a summary for the user in the GUI


use uuid::Uuid;

use crate::{blend::{disease_gene_entity::GeneDiseaseEntity, dto::{OverlapData, OverlapPlotPayload, PhenoBlendAnalysis, SpreadPlotPayload, UpsetPlotPayload}, upset_plot::UpsetData}, model::simple_term::SimpleOntologyTerm};
use std::time::{SystemTime, UNIX_EPOCH};
use serde::Serialize;
use tauri::{AppHandle, Emitter};


pub const NOTEBOOK_EVENT: &str = "notebook-updated";


#[derive(Serialize, Clone, Debug)]
#[serde(rename_all = "camelCase")]
pub struct Notebook {
    pub id: String,
    pub timestamp: u64,
    pub observed_hpo: Vec<SimpleOntologyTerm>,
    pub data_input_history: String,
    pub gene_disease_list: Vec<GeneDiseaseEntity>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub analysis: Option<PhenoBlendAnalysis>
}


pub fn emit_notebook(app: &AppHandle, notebook: &Notebook) -> Result<(), String> {
    app.emit(NOTEBOOK_EVENT, notebook).map_err(|e| e.to_string())
}

impl Notebook {
    fn new(observed_hpo: Vec<SimpleOntologyTerm>, data_input_history: String) -> Self {
        let id = Uuid::new_v4().to_string();
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("Time went backwards")
            .as_secs();
        Self {
            id,
            timestamp,
            observed_hpo,
            data_input_history,
            gene_disease_list: Vec::default(),
            analysis: None
        }
    }

    pub fn from_ppkt(observed_hpo: Vec<SimpleOntologyTerm>, fname: impl Into<String>) -> Self {
        let data_input_history = format!("HPO data from phenopacket '{}'.", fname.into());
        Notebook::new(observed_hpo, data_input_history)
    }

    pub fn from_ner(observed_hpo: Vec<SimpleOntologyTerm>) -> Self {
        let data_input_history = "HPO data entered manually.".to_string();
        Notebook::new(observed_hpo, data_input_history)
    }

    pub fn upset_payload(&self) -> Result<UpsetPlotPayload, String> {
           let analysis = self.analysis.as_ref().ok_or("Analysis not initialized")?;
           analysis.upset_payload()
    }

     pub fn upset_data(&self) -> Result<UpsetData, String> {
           let analysis = self.analysis.as_ref().ok_or("Analysis not initialized")?;
           Ok(analysis.upset.clone())
    }

    pub fn overlap_data(&self) -> Result<OverlapData, String> {
        let analysis = self.analysis.as_ref().ok_or("Analysis not initialized")?;
        Ok(analysis.overlap.clone())
    }

    pub fn overlap_payload(&self) -> Result<OverlapPlotPayload, String> {
        let analysis = self.analysis.as_ref().ok_or("Analysis not initialized")?;
        let overlap_data = &analysis.overlap;
        Ok(overlap_data.to_payload())
    }

    pub fn spread_payload(&self) -> Result<SpreadPlotPayload, String> {
           let analysis = self.analysis.as_ref().ok_or("Analysis not initialized")?;
           Ok(analysis.spread.clone())
    }

 
}