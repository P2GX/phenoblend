use std::collections::{HashMap, HashSet};
use serde::{Deserialize, Serialize};
use ontolius::{Identified, TermId};
use crate::util::serde_util::*;

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SimpleDiseaseModel {
    #[serde(
        serialize_with = "serialize_term_id", 
        deserialize_with = "parse_term_id"
    )]
    pub omim_disease_id: TermId,
    pub omim_disease_name: String,
     #[serde(
        serialize_with = "serialize_term_id_set",
        deserialize_with = "deserialize_term_id_set")]
    pub observed_hpo_ids: HashSet<TermId>,
    #[serde(
        serialize_with = "serialize_term_id_freq_map",
        deserialize_with = "deserialize_term_id_freq_map"
    )]
    pub term_frequencies: HashMap<TermId, f64>,
     #[serde(
        serialize_with = "serialize_term_id_set",
        deserialize_with = "deserialize_term_id_set")]
    pub excluded_hpo_ids: HashSet<TermId>
}


impl SimpleDiseaseModel {
    pub fn new(omim_disease_id: TermId,
                omim_disease_name: String,
                observed_hpo_ids: HashSet<TermId>,
                term_frequencies: HashMap<TermId, f64>,
                excluded_hpo_ids: HashSet<TermId>)
    -> Self {
        Self {
            omim_disease_id,
            omim_disease_name,
            observed_hpo_ids,
            term_frequencies,
            excluded_hpo_ids
        }
    }

    /// Fall back in case we cannot find the disease model, which should never happen
    pub fn fallback_from_id(omim_disease_id: TermId) -> Self {
        let msg = format!("Could not find disease model for '{}'.", omim_disease_id.to_string());
        Self { 
            omim_disease_id, 
            omim_disease_name: msg, 
            observed_hpo_ids: HashSet::default(), 
            term_frequencies: HashMap::default(),
            excluded_hpo_ids: HashSet::default() 
        }
    }

    /// Terms whose frequency meets or exceeds `min_freq`.
    /// This is the function your widget will call, live, per slider move.
    pub fn terms_at_or_above(&self, min_freq: f64) -> Vec<&TermId> {
        self.term_frequencies
            .iter()
            .filter(|(_, &freq)| freq >= min_freq)
            .map(|(term, _)| term)
            .collect()
    }

    pub fn omim_id(&self) -> String {
        self.omim_disease_id.identifier().to_string()
    }
}


#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct GeneDiseaseAssociation {
    pub ncbi_gene_id: String,
    pub gene_symbol: String,
    pub association_type: String,
    #[serde(deserialize_with = "parse_term_id", serialize_with="serialize_term_id")]
    pub disease_id: TermId, 
    pub disease_model: Option<SimpleDiseaseModel>,
    pub source: String,
}



