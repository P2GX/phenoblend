use phenopackets::schema::v2::Phenopacket;
use serde::{Deserialize, Serialize};
use serde_json::Value;

// If you have explicit domain models or use an HPO/Phenopacket crate:
#[derive(Serialize, Deserialize, Debug)]
pub struct PhenopacketPayload {
    id: String,
    // Add other fields according to GA4GH spec, or use fallback values
    subject: Option<Value>,
    phenotypic_features: Option<Vec<Value>>,
}


/// Inserts `survivalTimeInDays: 0` into `subject.vitalStatus` when the field is absent.
///
/// Anything that is not a JSON object at that path is left untouched, so malformed
/// input is rejected later by schema validation instead of causing a panic here.
fn default_missing_survival_time(json: &mut Value) {
    if let Some(vital_status) = json
        .pointer_mut("/subject/vitalStatus")
        .and_then(Value::as_object_mut)
    {
        vital_status
            .entry("survivalTimeInDays")
            .or_insert_with(|| Value::from(0));
    }
}

/// Parses a GA4GH Phenopacket from a JSON string.
///
/// The strict Phenopacket deserializer requires `subject.vitalStatus.survivalTimeInDays`
/// whenever a `vitalStatus` block is present, but many real-world files omit it.
/// Before deserializing, this function therefore inserts a default of `0` days if that
/// field is missing. Note that `0` is also the protobuf default for an unset integer,
/// so the inserted value is indistinguishable from "not specified".
///
/// # Errors
///
/// Returns an error message if:
///
/// * `ppkt_string` is not syntactically valid JSON, or
/// * the JSON does not conform to the Phenopacket schema (for example, a field has
///   the wrong type or `vitalStatus` is not an object).
///
/// # Examples
///
/// ```ignore
/// let json = r#"{ "id": "example", "subject": { "id": "patient1" } }"#;
/// let ppkt = load_from_json(json)?;
/// assert_eq!(ppkt.id, "example");
/// ```
pub fn load_from_json(ppkt_string: &str) -> Result<Phenopacket, String> {
    let mut json: Value = serde_json::from_str(ppkt_string)
        .map_err(|e| format!("Invalid JSON syntax structure: {e}"))?;
    default_missing_survival_time(&mut json);
    serde_json::from_value(json)
        .map_err(|e| format!("Phenopacket schema validation error: {e}"))
}