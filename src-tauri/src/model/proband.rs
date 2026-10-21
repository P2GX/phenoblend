use std::collections::BTreeMap;

use ontolius::TermId;

use serde::{Serialize};
use crate::blend::dto::Severity;
use crate::model::simple_term::SimpleOntologyTerm;
use crate::util::serde_util::serialize_term_id_vec;

#[derive(Serialize, Debug, Clone)]
pub struct Proband {
    pub id: String,
    pub gene_list: Vec<String>,
    pub disease_list: Vec<SimpleOntologyTerm>,
    #[serde(serialize_with = "serialize_term_id_vec")]
    pub observed_hpos: Vec<TermId>,
    /// HPO term id (string) -> severity recorded for that feature in the phenopacket. Generally only a small minority of terms have severity
    #[serde(default)]
    pub severities: BTreeMap<String, Severity>,
}




impl Default for Proband {
    
    fn default() -> Self {
      Self {    
        id: String::default(),
        gene_list: Vec::default(),
        disease_list: Vec::default(),
        observed_hpos: Vec::default(),
        severities: BTreeMap::default(),
        }
    }
}