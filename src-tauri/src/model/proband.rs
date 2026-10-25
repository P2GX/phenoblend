use std::collections::BTreeMap;

use ontolius::TermId;

use phenopackets::schema::v2::Phenopacket;
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

impl Proband {
    /// Extracts the proband's data from a GA4GH phenopacket, without touching any application
    /// state.
    ///
    /// * `observed_hpos`: the HPO terms of the phenotypic features that are not excluded, in
    ///   file order. Duplicates are kept; later steps collapse them.
    /// * `severities`: the severity of a feature, keyed by its HPO id, when the severity is
    ///   one of the HPO severity terms.
    /// * `disease_list`: the diseases named in the phenopacket.
    /// * `gene_list`: the gene symbols of the genomic interpretations, sorted and without
    ///   duplicates (two variants in one gene give one entry).
    ///
    /// # Errors
    /// Returns an error if the phenopacket has no subject, if a phenotypic feature has no
    /// type or an id that is not a valid term id, or if a disease has no term.
    pub fn from_phenopacket(ppkt: &Phenopacket) -> Result<Self, String> {
        let id = ppkt
            .subject
            .as_ref()
            .map(|subject| subject.id.clone())
            .ok_or_else(|| format!("phenopacket '{}' has no subject", ppkt.id))?;

        let mut observed_hpos: Vec<TermId> = Vec::new();
        let mut severities: BTreeMap<String, Severity> = BTreeMap::new();
        for feature in &ppkt.phenotypic_features {
            if feature.excluded {
                continue;
            }
            let class = feature
                .r#type
                .as_ref()
                .ok_or_else(|| format!("phenopacket '{}' has a phenotypic feature without a type", ppkt.id))?;
            let term_id: TermId = class
                .id
                .parse()
                .map_err(|_| format!("Failed to parse TermId from string '{}'", class.id))?;
            if let Some(severity_class) = &feature.severity {
                match Severity::from_hpo_id(&severity_class.id) {
                    Some(severity) => {
                        severities.insert(class.id.clone(), severity);
                    }
                    None => log::warn!(
                        "Unrecognized severity term {} for feature {}",
                        severity_class.id,
                        class.id
                    ),
                }
            }
            observed_hpos.push(term_id);
        }

        let disease_list = ppkt
            .diseases
            .iter()
            .map(|disease| {
                let term = disease
                    .term
                    .as_ref()
                    .ok_or_else(|| format!("phenopacket '{}' has a disease without a term", ppkt.id))?;
                SimpleOntologyTerm::new(term.id.as_str(), term.label.as_str())
            })
            .collect::<Result<Vec<_>, String>>()?;

        let mut gene_list: Vec<String> = ppkt
            .interpretations
            .iter()
            .filter_map(ga4ghphetools::get_gene_symbol_from_interpretation)
            .collect();
        gene_list.sort();
        gene_list.dedup();

        Ok(Self { id, gene_list, disease_list, observed_hpos, severities })
    }
}