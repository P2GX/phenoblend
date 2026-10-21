use serde::Deserialize;

use crate::blend::dto::Severity;

#[derive(Deserialize, Debug)]
#[serde(rename_all = "camelCase")]
pub struct ObservedHpoInput {
    pub term_id: String,
    #[serde(default)]
    pub modifier_ids: Vec<String>,
}


/// First severity term among `ids`, 
/// By design, ng-hpo-uikit only adds one severity term at most to the modifiers
/// TODO: consider refactor
pub fn severity_from_ids<'a>(
    modifiers: &Vec<String>
) -> Option<Severity> {
    for modifier_tid in modifiers.iter() {
        if let Some(sev) = Severity::from_hpo_id(modifier_tid) {
            return Some(sev);
        }  
    }
    None
}
