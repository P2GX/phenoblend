use ontolius::TermId;

use serde::{Serialize};
use crate::model::simple_term::SimpleOntologyTerm;
use crate::util::serde_util::serialize_term_id_vec;

#[derive(Serialize, Debug, Clone)]
pub struct Proband {
    pub id: String,
    pub gene_list: Vec<String>,
    pub disease_list: Vec<SimpleOntologyTerm>,
    #[serde(serialize_with = "serialize_term_id_vec")]
    pub observed_hpos: Vec<TermId>
}




impl Default for Proband {
    
    fn default() -> Self {
      Self {    
        id: String::default(),
        gene_list: Vec::default(),
        disease_list: Vec::default(),
        observed_hpos: Vec::default(),
        }
    }
}