use std::collections::{HashMap, HashSet};

use ontolius::TermId;
use serde::ser::SerializeMap;
use serde::{Deserialize, Serializer};
use serde::de::{Error, MapAccess};
use serde::de::{Deserializer, SeqAccess, Visitor};
use std::fmt;



// Custom deserializer helper function:
pub fn parse_term_id<'de, D>(deserializer: D) -> Result<TermId, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    let tid: TermId = s.parse().map_err(|e| {
        D::Error::custom(format!("Failed to parse TermId from string '{}': {}", s, e))
    })?;
    Ok(tid) // Or TermId::new(s) depending on your implementation
}


pub fn serialize_term_id<S>(term_id: &TermId, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.serialize_str(&term_id.to_string())
}

pub fn serialize_term_id_set<S>(term_ids: &HashSet<TermId>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    use serde::ser::SerializeSeq;
    let mut seq = serializer.serialize_seq(Some(term_ids.len()))?;
    for id in term_ids {
        seq.serialize_element(&id.to_string())?;
    }
    seq.end()
}



pub fn deserialize_term_id_set<'de, D>(deserializer: D) -> Result<HashSet<TermId>, D::Error>
where
    D: Deserializer<'de>,
{
    // Define a Visitor struct to handle the sequence parsing sequence
    struct TermIdSetVisitor;

    impl<'de> Visitor<'de> for TermIdSetVisitor {
        type Value = HashSet<TermId>;

        fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
            formatter.write_str("a sequence of HPO TermId strings")
        }

        fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
        where
            A: SeqAccess<'de>,
        {
            let mut set = HashSet::with_capacity(seq.size_hint().unwrap_or(0));

            // Iterate over every element in the JSON array
            while let Some(s) = seq.next_element::<String>()? {
                // Parse the string into a TermId using your implementation's FromStr / parse
                let tid: TermId = s.parse().map_err(|e| {
                    serde::de::Error::custom(format!(
                        "Failed to parse TermId from set element '{}': {}", 
                        s, e
                    ))
                })?;
                set.insert(tid);
            }

            Ok(set)
        }
    }

    // Pass our visitor to the deserializer
    deserializer.deserialize_seq(TermIdSetVisitor)
}


// ---------- HashMap<TermId, f64> ----------
// Analogous, but as a JSON *object* {"HP:0009737": 0.78, ...} rather than
// a sequence — this is what was missing, and without it `term_frequencies`
// either fails to compile or panics at runtime with "key must be a string".

pub fn serialize_term_id_freq_map<S>(
    map: &HashMap<TermId, f64>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let mut m = serializer.serialize_map(Some(map.len()))?;
    for (id, freq) in map {
        m.serialize_entry(&id.to_string(), freq)?;
    }
    m.end()
}

pub fn deserialize_term_id_freq_map<'de, D>(
    deserializer: D,
) -> Result<HashMap<TermId, f64>, D::Error>
where
    D: Deserializer<'de>,
{
    struct TermIdFreqMapVisitor;

    impl<'de> Visitor<'de> for TermIdFreqMapVisitor {
        type Value = HashMap<TermId, f64>;

        fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
            formatter.write_str("a map of HPO TermId strings to frequency floats")
        }

        fn visit_map<A>(self, mut access: A) -> Result<Self::Value, A::Error>
        where
            A: MapAccess<'de>,
        {
            let mut map = HashMap::with_capacity(access.size_hint().unwrap_or(0));
            while let Some((key, value)) = access.next_entry::<String, f64>()? {
                let tid: TermId = key.parse().map_err(|e| {
                    serde::de::Error::custom(format!(
                        "Failed to parse TermId key '{}': {}",
                        key, e
                    ))
                })?;
                map.insert(tid, value);
            }
            Ok(map)
        }
    }

    deserializer.deserialize_map(TermIdFreqMapVisitor)
}

#[cfg(test)]
mod tests {
    use crate::hpoa::disease_model::SimpleDiseaseModel;

use super::*;

    fn tid(s: &str) -> TermId {
        s.parse().unwrap()
    }

    fn sample_model() -> SimpleDiseaseModel {
        SimpleDiseaseModel {
            omim_disease_id: tid("OMIM:162200"),
            omim_disease_name: "Neurofibromatosis type 1".to_string(),
            observed_hpo_ids: HashSet::from([tid("HP:0009737")]), // Lisch nodules
            term_frequencies: HashMap::from([(tid("HP:0009737"), 74.0 / 95.0)]),
            excluded_hpo_ids: HashSet::from([tid("HP:0001250")]),
        }
    }

    #[test]
    fn round_trips_through_json() {
        let model = sample_model();
        let json = serde_json::to_string(&model).expect("serialize");
        let back: SimpleDiseaseModel = serde_json::from_str(&json).expect("deserialize");

        assert_eq!(back.omim_disease_id, model.omim_disease_id);
        assert_eq!(back.observed_hpo_ids, model.observed_hpo_ids);
        assert_eq!(back.excluded_hpo_ids, model.excluded_hpo_ids);
        assert_eq!(
            back.term_frequencies[&tid("HP:0009737")],
            model.term_frequencies[&tid("HP:0009737")]
        );
    }

    #[test]
    fn term_frequencies_serializes_as_object_not_array() {
        let model = sample_model();
        let json = serde_json::to_value(&model).expect("serialize");
        let freqs = &json["termFrequencies"];

        // Must be a JSON object keyed by term id string, e.g. {"HP:0009737": 0.778...}
        assert!(freqs.is_object(), "expected object, got {freqs:?}");
        assert!(freqs.get("HP:0009737").is_some());
    }

    #[test]
    fn rejects_malformed_term_id_in_map_key() {
        // the frequency CURIE is not valid
        let bad_json = r#"{
            "omimDiseaseId": "OMIM:162200",
            "omimDiseaseName": "NF1",
            "observedHpoIds": ["HP:0009737"],
            "termFrequencies": {"NOT A VALID TERM ID": 0.5},
            "excludedHpoIds": []
        }"#;
        let result: Result<SimpleDiseaseModel, _> = serde_json::from_str(bad_json);
        assert!(result.is_err());
    }
}