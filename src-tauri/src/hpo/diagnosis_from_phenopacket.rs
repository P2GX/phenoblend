use std::collections::BTreeSet;

use phenopackets::schema::v2::Phenopacket;

use crate::model::simple_term::SimpleOntologyTerm;

/// A diagnosed disease together with the gene that the genomic interpretation points to.
#[derive(Debug, Clone)]
pub struct Diagnosis {
    pub gene_symbol: String,
    pub disease: SimpleOntologyTerm,
}

/// The (gene, disease) pairs of the phenopacket's interpretations, sorted by gene symbol and
/// then disease id, without duplicates.
///
/// Every interpretation must give both a gene symbol (as read by
/// `ga4ghphetools::get_gene_symbol_from_interpretation`, which returns one symbol per
/// interpretation) and a disease (the `disease` of its `diagnosis`).
///
/// # Errors
/// Returns an error naming the phenopacket and the interpretation if an interpretation has no
/// gene symbol or no disease, or if a disease id or label cannot form a valid
/// [`SimpleOntologyTerm`].
pub fn diagnoses_from_phenopacket(ppkt: &Phenopacket) -> Result<Vec<Diagnosis>, String> {
    let mut diagnoses: Vec<Diagnosis> = Vec::new();
    let mut seen: BTreeSet<(String, String)> = BTreeSet::new();
    for interpretation in &ppkt.interpretations {
        let symbol = ga4ghphetools::get_gene_symbol_from_interpretation(interpretation)
            .ok_or_else(|| {
                format!(
                    "interpretation '{}' of phenopacket '{}' has no gene symbol",
                    interpretation.id, ppkt.id
                )
            })?;
        let class = interpretation
            .diagnosis
            .as_ref()
            .and_then(|d| d.disease.as_ref())
            .ok_or_else(|| {
                format!(
                    "interpretation '{}' of phenopacket '{}' has no disease",
                    interpretation.id, ppkt.id
                )
            })?;
        if !seen.insert((symbol.clone(), class.id.clone())) {
            continue; // the same gene and disease listed twice
        }
        diagnoses.push(Diagnosis {
            gene_symbol: symbol,
            disease: SimpleOntologyTerm::new(class.id.as_str(), class.label.as_str())?,
        });
    }
    diagnoses.sort_by(|a, b| {
        (a.gene_symbol.as_str(), a.disease.term_id_as_string())
            .cmp(&(b.gene_symbol.as_str(), b.disease.term_id_as_string()))
    });
    Ok(diagnoses)
}

/// The demand of the multiple-diagnosis analysis: at least two distinct diseases, and at
/// least two distinct genes (the test compares the terms of different genes).
///
/// # Errors
/// Returns a message saying what is missing, suitable for the "reason skipped" column of a
/// batch result.
pub fn require_multiple_diagnoses(diagnoses: &[Diagnosis]) -> Result<(), String> {
    let diseases: BTreeSet<String> =
        diagnoses.iter().map(|d| d.disease.term_id_as_string()).collect();
    let genes: BTreeSet<&str> = diagnoses.iter().map(|d| d.gene_symbol.as_str()).collect();
    if diseases.len() < 2 {
        return Err(format!(
            "need at least two diseases with a gene in the interpretations, found {}",
            diseases.len()
        ));
    }
    if genes.len() < 2 {
        return Err(format!("need at least two distinct genes, found {genes:?}"));
    }
    Ok(())
}

// =====================================================================================
// Tests. The phenopacket is written as JSON, so no nested protobuf structs are built by hand.
// If it fails to parse, adjust the JSON to what your `load_from_json` accepts (for example
// add the `metaData` block that your real files have).
// =====================================================================================

#[cfg(test)]
mod diagnosis_tests {
    use crate::hpo::ppkt_loader::load_from_json;

use super::*;

    fn interpretation(id: &str, disease: Option<(&str, &str)>, gene: Option<&str>) -> String {
        let diagnosis = match (disease, gene) {
            (None, _) => String::new(),
            (Some((did, label)), gene) => {
                let genomic = match gene {
                    Some(symbol) => format!(
                        r#""genomicInterpretations": [{{
                            "subjectOrBiosampleId": "s1",
                            "interpretationStatus": "CAUSATIVE",
                            "variantInterpretation": {{
                                "acmgPathogenicityClassification": "PATHOGENIC",
                                "variationDescriptor": {{
                                    "id": "var-{symbol}",
                                    "geneContext": {{ "valueId": "HGNC:1", "symbol": "{symbol}" }}
                                }}
                            }}
                        }}]"#
                    ),
                    None => r#""genomicInterpretations": []"#.to_string(),
                };
                format!(r#""diagnosis": {{ "disease": {{ "id": "{did}", "label": "{label}" }}, {genomic} }}"#)
            }
        };
        format!(r#"{{ "id": "{id}", "progressStatus": "SOLVED"{} }}"#, if diagnosis.is_empty() { String::new() } else { format!(", {diagnosis}") })
    }

    fn ppkt(interpretations: &[String]) -> Phenopacket {
        let json = format!(
            r#"{{ "id": "ppkt-1", "subject": {{ "id": "s1" }}, "interpretations": [{}] }}"#,
            interpretations.join(",")
        );
        load_from_json(&json).expect("the test phenopacket must parse")
    }

    fn pairs(diagnoses: &[Diagnosis]) -> Vec<(String, String)> {
        diagnoses
            .iter()
            .map(|d| (d.gene_symbol.clone(), d.disease.term_id_as_string()))
            .collect()
    }

    const NF1: (&str, &str) = ("OMIM:162200", "Neurofibromatosis, type 1");
    const NOONAN: (&str, &str) = ("OMIM:163950", "Noonan syndrome 1");

    #[test]
    fn each_interpretation_gives_a_gene_and_its_disease_sorted_by_gene() {
        let p = ppkt(&[
            interpretation("i1", Some(NOONAN), Some("PTPN11")),
            interpretation("i2", Some(NF1), Some("NF1")),
        ]);
        let diagnoses = diagnoses_from_phenopacket(&p).unwrap();
        assert_eq!(
            pairs(&diagnoses),
            [
                ("NF1".to_string(), "OMIM:162200".to_string()),
                ("PTPN11".to_string(), "OMIM:163950".to_string()),
            ]
        );
        assert!(require_multiple_diagnoses(&diagnoses).is_ok());
    }

    #[test]
    fn the_same_gene_and_disease_listed_twice_counts_once() {
        let p = ppkt(&[
            interpretation("i1", Some(NF1), Some("NF1")),
            interpretation("i2", Some(NF1), Some("NF1")),
        ]);
        assert_eq!(diagnoses_from_phenopacket(&p).unwrap().len(), 1);
    }

    #[test]
    fn an_interpretation_without_a_gene_is_an_error_naming_it() {
        let p = ppkt(&[
            interpretation("i1", Some(NOONAN), Some("PTPN11")),
            interpretation("i2", Some(NF1), None),
        ]);
        let err = diagnoses_from_phenopacket(&p).unwrap_err();
        assert!(err.contains("i2") && err.contains("ppkt-1") && err.contains("no gene symbol"), "{err}");
    }

    #[test]
    fn an_interpretation_without_a_disease_is_an_error_naming_it() {
        let p = ppkt(&[
            interpretation("i1", Some(NOONAN), Some("PTPN11")),
            interpretation("i2", None, None),
        ]);
        let err = diagnoses_from_phenopacket(&p).unwrap_err();
        assert!(err.contains("i2") && err.contains("ppkt-1"), "{err}");
    }

    #[test]
    fn one_diagnosis_is_not_enough() {
        let p = ppkt(&[interpretation("i1", Some(NF1), Some("NF1"))]);
        let err = require_multiple_diagnoses(&diagnoses_from_phenopacket(&p).unwrap()).unwrap_err();
        assert!(err.contains("at least two diseases"), "{err}");
    }

    #[test]
    fn two_diseases_of_the_same_gene_are_not_enough() {
        let p = ppkt(&[
            interpretation("i1", Some(NF1), Some("NF1")),
            interpretation("i2", Some(NOONAN), Some("NF1")),
        ]);
        let err = require_multiple_diagnoses(&diagnoses_from_phenopacket(&p).unwrap()).unwrap_err();
        assert!(err.contains("two distinct genes"), "{err}");
    }

    #[test]
    fn no_interpretations_gives_no_diagnoses() {
        let p = ppkt(&[]);
        let diagnoses = diagnoses_from_phenopacket(&p).unwrap();
        assert!(diagnoses.is_empty());
        assert!(require_multiple_diagnoses(&diagnoses).is_err());
    }
}
