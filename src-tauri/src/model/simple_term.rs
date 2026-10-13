
use ontolius::TermId;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use serde::de::Error as _; // brings `custom` into scope without taking the name `Error`

#[derive(Deserialize,Serialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct SimpleOntologyTerm {
     #[serde(serialize_with = "serialize_term_id", deserialize_with="deserialize_term_id")]
    pub term_id: TermId,
    pub term_label: String
}

// A short helper function to cleanly serialize TermId objects via their string representation
fn serialize_term_id<S>(term_id: &TermId, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.serialize_str(&term_id.to_string())
}

fn deserialize_term_id<'de, D>(deserializer: D) -> Result<TermId, D::Error>
where
    D: Deserializer<'de>, {
    let s = String::deserialize(deserializer)?;   // read the JSON string first
    s.parse::<TermId>()  // the last expression is the return value
        .map_err(|e| D::Error::custom(format!("Failed to parse TermId from '{s}': {e}")))
}


impl SimpleOntologyTerm {
    pub fn new(tid: impl Into<String>, label: impl Into<String>) -> Result<Self, String> {
        let s = tid.into();
        let trmid: TermId = s.parse().map_err(|_| format!("Failed to parse TermId from string '{}'.", s))?;
        Ok(Self {
            term_id: trmid,
            term_label: label.into()
        })
    }

    pub fn from_term_id(term_id: &TermId,  label: impl Into<String>) -> Self {
        Self {
            term_id: term_id.clone(),
            term_label: label.into()
        }
    }


    pub fn term_id_as_string(&self) -> String {
        self.term_id.to_string()
    }


 
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_parses_valid_term_id() {
        let t = SimpleOntologyTerm::new("HP:0001250", "Seizure").unwrap();
        assert_eq!(t.term_label, "Seizure");
        assert_eq!(t.term_id_as_string(), "HP:0001250");
    }

    #[test]
    fn new_accepts_owned_strings() {
        let t = SimpleOntologyTerm::new(String::from("HP:0000118"), String::from("Phenotypic abnormality")).unwrap();
        assert_eq!(t.term_id_as_string(), "HP:0000118");
        assert_eq!(t.term_label, "Phenotypic abnormality");
    }

    #[test]
    fn new_rejects_invalid_term_ids() {
        for bad in ["", "not a term id", "HP0001250"] {
            let err = SimpleOntologyTerm::new(bad, "label").unwrap_err();
            assert!(
                err.contains(bad),
                "error message should echo the bad input '{bad}', got: {err}"
            );
        }
    }

    #[test]
    fn from_term_id_copies_id_and_label() {
        let tid: TermId = "HP:0001250".parse().unwrap();
        let t = SimpleOntologyTerm::from_term_id(&tid, "Seizure");
        assert_eq!(t.term_id, tid);
        assert_eq!(t.term_label, "Seizure");
        // original is still usable (function borrows, doesn't consume)
        assert_eq!(tid.to_string(), "HP:0001250");
    }

    #[test]
    fn new_and_from_term_id_agree() {
        let tid: TermId = "HP:0001250".parse().unwrap();
        let a = SimpleOntologyTerm::new("HP:0001250", "Seizure").unwrap();
        let b = SimpleOntologyTerm::from_term_id(&tid, "Seizure");
        assert_eq!(a.term_id, b.term_id);
        assert_eq!(a.term_label, b.term_label);
    }

    #[test]
    fn serializes_term_id_as_curie_string() {
        let t = SimpleOntologyTerm::new("HP:0001250", "Seizure").unwrap();
        let v = serde_json::to_value(&t).unwrap();
        assert_eq!(v["termId"], "HP:0001250");
        assert_eq!(v["termLabel"], "Seizure");
        // exactly two fields, so the frontend contract doesn't silently change
        assert_eq!(v.as_object().unwrap().len(), 2);
    }

    #[test]
    fn serializes_vec_of_terms() {
        let terms = vec![
            SimpleOntologyTerm::new("HP:0001250", "Seizure").unwrap(),
            SimpleOntologyTerm::new("HP:0000118", "Phenotypic abnormality").unwrap(),
        ];
        let json = serde_json::to_string(&terms).unwrap();
        assert_eq!(
            json,
            r#"[{"termId":"HP:0001250","termLabel":"Seizure"},{"termId":"HP:0000118","termLabel":"Phenotypic abnormality"}]"#
        );
    }
}