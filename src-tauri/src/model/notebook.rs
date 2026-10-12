//! Notebook
//! Keep track of analysis run
//! Intended to be displayed as a summary for the user in the GUI


use uuid::Uuid;

use crate::{blend::disease_gene_entity::GeneDiseaseEntity, model::simple_term::SimpleOntologyTerm};
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
}


pub fn emit_notebook(app: &AppHandle, notebook: &Notebook) -> Result<(), String> {
    eprintln!("emit notebook no={:?}", notebook);
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
            gene_disease_list: Vec::default()
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

 
}