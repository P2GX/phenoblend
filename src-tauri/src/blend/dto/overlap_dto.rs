//! Overlap
//! Test the overlap of HPO terms of a phenopacket with the distribution of HPO terms in the disease model.
//! 
//! 
use std::{collections::{BTreeMap, HashMap}, sync::Arc};

use serde::{Deserialize,Serialize};
use ontolius::{TermId, ontology::{OntologyTerms, csr::FullCsrOntology}, term::MinimalTerm};
use crate::{blend::dto::constants::ROOT_PHENOTYPE_ID, model::simple_term::SimpleOntologyTerm};


/// How a gene-disease entity relates to one proband term.
#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "camelCase")]
pub enum MatchType {
    /// The disease is annotated to exactly the proband term.
    Exact,
    /// The disease is annotated to one or more descendants of the proband term
    /// (true path rule, e.g. proband "Cataract", disease "Nuclear cataract").
    Descendant,
    /// No perfect match; `score` is a partial-similarity score. 
    Partial,
    /// Only the root term is shared. Score=0.0 means match at root (no match!)
    NoMatch,
}

impl MatchType {
    /// True if the disease annotation accounts for the proband term
    /// (`Exact` or `Descendant`). `Partial` and `NoMatch` do not count.
    pub fn is_match(self) -> bool {
        matches!(self, MatchType::Exact | MatchType::Descendant)
    }

     pub fn is_partial_match(self) -> bool {
        matches!(self, MatchType::Partial)
    }
}


/// The outcome of matching one proband term against one gene-disease entity.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct GeneMatch {
    pub match_type: MatchType,
    /// 1.0 for `Exact` and `Descendant`; in `[0.0, 1.0]` for `Partial`.
    pub score: f64,
    /// The single term responsible for the match: the identical disease term for `Exact`;
    /// the closest matching descendant for `Descendant` (ties broken by term ID);
    /// the best shared term for `Partial`; `HP:0000118` for `None` (i.e., if there is no specific match, take the term id for 
    /// the root Phenotypic abnormality, HP:0000118).
    pub matched_term: SimpleOntologyTerm,
}


impl GeneMatch {
    pub fn exact(term: SimpleOntologyTerm) -> Self {
        Self { match_type: MatchType::Exact, score: 1.0, matched_term: term }
    }

    pub fn descendant(term: SimpleOntologyTerm) -> Self {
        Self { match_type: MatchType::Descendant, score: 1.0, matched_term: term }
    }

    /// Builds a partial match, or a `None` match when the shared term is the root.
    ///
    /// # Errors
    /// Returns an error if `score` is not finite or lies outside `[0, 1]`, or if the
    /// term is the root but the score is not 0 (a contradiction).
    pub fn partial(score: f64, term: SimpleOntologyTerm) -> Result<Self, String> {
        if !score.is_finite() || !(0.0..=1.0).contains(&score) {
            return Err(format!("partial match score out of range: {score}"));
        }
        if term.term_id_as_string() == ROOT_PHENOTYPE_ID {
            if score != 0.0 {
                return Err(format!("root term {ROOT_PHENOTYPE_ID} with non-zero score {score}"));
            }
            return Ok(Self { match_type: MatchType::NoMatch, score: 0.0, matched_term: term });
        }
        Ok(Self { match_type: MatchType::Partial, score, matched_term: term })
    }

    pub fn is_match(&self) -> bool {
        self.match_type.is_match()
    }

    pub fn is_partial_match(&self) -> bool {
        self.match_type.is_partial_match()
    }
}



/// Severity of an observed feature (children of HP:0012824, "Severity").
#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "camelCase")]
pub enum Severity { Borderline, Mild, Moderate, Severe, Profound }

impl Severity {
    pub fn from_hpo_id(id: &str) -> Option<Self> {
        match id {
            "HP:0012827" => Some(Self::Borderline),
            "HP:0012825" => Some(Self::Mild),
            "HP:0012826" => Some(Self::Moderate),
            "HP:0012828" => Some(Self::Severe),
            "HP:0012829" => Some(Self::Profound),
            _ => None,
        }
    }
}


/// All match results for one observed proband term, across all gene entities.
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct TermMatches {
    pub term: SimpleOntologyTerm,
    /// Optionally indicate the degree of severity reported for the finding.
    pub severity: Option<Severity>,
    /// Gene symbol -> match
    pub matches: BTreeMap<String, GeneMatch>,
}

impl TermMatches {
    pub fn new(term: SimpleOntologyTerm) -> Self {
        Self { term, severity: None, matches: BTreeMap::new() }
    }

    /// Sets the severity recorded for this observed term (`None` leaves it unset).
    #[must_use]
    pub fn with_severity(mut self, severity: Option<Severity>) -> Self {
        self.severity = severity;
        self
    }

    pub fn severe(term: SimpleOntologyTerm) -> Self {
        Self { term, severity: Some(Severity::Severe), matches: BTreeMap::new() }
    }

     pub fn profound(term: SimpleOntologyTerm) -> Self {
        Self { term, severity: Some(Severity::Profound), matches: BTreeMap::new() }
    }
    
    pub fn moderate(term: SimpleOntologyTerm) -> Self {
        Self { term, severity: Some(Severity::Moderate), matches: BTreeMap::new() }
    }

    pub fn mild(term: SimpleOntologyTerm) -> Self {
        Self { term, severity: Some(Severity::Mild), matches: BTreeMap::new() }
    }

       pub fn borderline(term: SimpleOntologyTerm) -> Self {
        Self { term, severity: Some(Severity::Borderline), matches: BTreeMap::new() }
    }


    pub fn insert(&mut self, gene_symbol: &str, m: GeneMatch) {
        self.matches.insert(gene_symbol.to_string(), m);
    }

    /// The wire format your overlap graphic already reads (`hpoId`, `hpoName`, `scores`),
    /// so the frontend contract doesn't change.
    pub fn to_overlap_item(&self) -> OverlapPlotItem {
        OverlapPlotItem {
            hpo_id: self.term.term_id_as_string(),
            hpo_name: self.term.term_label.clone(),
            severity: self.severity,
            scores: self.matches.iter().map(|(g, m)| (g.clone(), m.score)).collect(),
        }
    }
}

/// Represents the scores for one HPO term in the presence matrix
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct OverlapPlotItem {
    pub hpo_id: String,
    pub hpo_name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub severity: Option<Severity>,
    /// A map matching gene symbols to their calculated matrix match scores [0.0 - 1.0]
    pub scores: HashMap<String, f64>,
}

impl OverlapPlotItem {
    pub fn new(hpo_id: &TermId, hpo: Arc<FullCsrOntology>) -> Self {
        let term = hpo.term_by_id(hpo_id);
        let label = match term {
            Some(t) => t.name().to_string(),
            None => "n/a".to_string() // should never happen
        }; 
        Self {
            hpo_id: hpo_id.to_string(),
            hpo_name: label,
            severity: None,
            scores: HashMap::default()
        }
    }

    pub fn add_perfect_match(&mut self, entity: &str) {
        self.scores.insert(entity.to_string(), 1.0);
    }

    pub fn add_match(&mut self, entity: &str, score: f64) {
        self.scores.insert(entity.to_string(), score);
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "camelCase")]
pub struct OverlapPlotPayload {
    /// Ordered array of gene symbols or diseases following column sorting block rules
    pub entities: Vec<String>,
    /// Ordered rows following phenotype grouping tiers
    pub columns: Vec<OverlapPlotItem>,
}


impl OverlapPlotPayload {
    pub fn n_columns(&self) -> usize {
        self.columns.len()
    }

   
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct OverlapData {
     /// Ordered array of gene symbols following column sorting block rules (the genes represent one or more diseases associated with the gene)
    pub entities: Vec<String>,
    /// Ordered rows following phenotype grouping tiers
    pub rows: Vec<TermMatches>,
     /// `n_matched_terms[i]` is the number of rows (proband terms) that have a
    /// match for `entities[i]`. Same length and order as `entities`.
    pub n_matched_terms: Vec<usize>,
    /// `n_partially_matched_terms[i]` is the number of rows (proband terms) that have a
    /// partial match for `entities[i]`.
    pub n_partially_matched_terms: Vec<usize>
}




impl OverlapData {
    pub fn to_payload(&self) -> OverlapPlotPayload {
        let overlap_list: Vec<OverlapPlotItem> = self.rows.iter()
            .map(|entry| entry.to_overlap_item())
            .collect();
        OverlapPlotPayload {
            entities: self.entities.clone(),
            columns: overlap_list
        }
    }

     /// Derives the overlap graphic's payload from the stored match records.
    ///
    /// `entities` is the sorted, de-duplicated list of gene symbols; `columns` has one
    /// entry per record, in the same order.
    pub fn from_term_matches<I>(gene_symbols: I, term_matches: &[TermMatches]) -> Self
    where
        I: IntoIterator<Item = String>,
    {
        let mut entities: Vec<String> = gene_symbols.into_iter().collect();
        entities.sort();
        entities.dedup();
        let mut n_matched_terms = vec![0usize; entities.len()];
        let mut n_partially_matched_terms = vec![0usize; entities.len()];
        for tm in term_matches {
            for (count, gene) in n_matched_terms.iter_mut().zip(&entities) {
                if tm.matches.get(gene).is_some_and(GeneMatch::is_match) {
                    *count += 1;
                } 
            }
            for (count, gene) in n_partially_matched_terms.iter_mut().zip(&entities) {
                if tm.matches.get(gene).is_some_and(GeneMatch::is_partial_match) {
                    *count += 1;
                }
            }
        }
        Self {
            entities,
            rows: term_matches.to_vec(),
            n_matched_terms,
            n_partially_matched_terms
        }
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use rstest::{fixture, rstest};
    use serde_json::json;

    fn term(id: &str, label: &str) -> SimpleOntologyTerm {
        SimpleOntologyTerm::new(id, label).unwrap()
    }

    /// The root term is the only one for which a partial match becomes `NoMatch`.
    fn root() -> SimpleOntologyTerm {
        term(ROOT_PHENOTYPE_ID, "Phenotypic abnormality")
    }

    /// One observed term with one match per `(gene, match)` pair.
    fn row(term: SimpleOntologyTerm, matches: &[(&str, &GeneMatch)]) -> TermMatches {
        let mut tm = TermMatches::new(term);
        for (gene, m) in matches {
            tm.insert(gene, (*m).clone());
        }
        tm
    }

    fn genes(names: &[&str]) -> Vec<String> {
        names.iter().map(|g| g.to_string()).collect()
    }

    /// A match of the requested type, so a case can say what it means and not how it is built.
    fn match_of(match_type: MatchType) -> GeneMatch {
        let m = match match_type {
            MatchType::Exact => GeneMatch::exact(term("HP:0002650", "Scoliosis")),
            MatchType::Descendant => GeneMatch::descendant(term("HP:0002650", "Scoliosis")),
            MatchType::Partial => GeneMatch::partial(0.5, term("HP:0000316", "Hypertelorism")).unwrap(),
            MatchType::NoMatch => GeneMatch::partial(0.0, root()).unwrap(),
        };
        assert_eq!(m.match_type, match_type, "test setup built the wrong match type");
        m
    }

    /// Invariants that must hold for every `OverlapData`, whatever its contents.
    fn assert_consistent(od: &OverlapData) {
        assert_eq!(od.n_matched_terms.len(), od.entities.len());
        assert_eq!(od.n_partially_matched_terms.len(), od.entities.len());
        for (i, entity) in od.entities.iter().enumerate() {
            // A term cannot be both a match and a partial match for the same gene.
            assert!(od.n_matched_terms[i] + od.n_partially_matched_terms[i] <= od.rows.len(),
                "{entity}: more matches than rows");
        }
        assert!(od.entities.windows(2).all(|w| w[0] < w[1]), "entities must be sorted and unique");
    }

    // ---------------------------------------------------------------------------------
    // Shared data. Each gene gets a different pattern, so a bug that swapped the genes
    // (or mixed up the two counters) cannot pass by accident.
    //
    //                   NF1          PTPN11
    //   scoliosis       exact        exact
    //   cleft palate    descendant   exact
    //   hypertelorism   partial      descendant
    //   macrocephaly    partial      (no entry)
    //
    //   => NF1: 2 matched, 2 partial      PTPN11: 3 matched, 0 partial
    // ---------------------------------------------------------------------------------

    #[fixture]
    fn four_term_matches() -> Vec<TermMatches> {
        let scoliosis = term("HP:0002650", "Scoliosis");
        let palate = term("HP:0000175", "Cleft palate");
        let hypertelorism = term("HP:0000316", "Hypertelorism");
        let macrocephaly = term("HP:0000256", "Macrocephaly");
        vec![
            row(scoliosis.clone(), &[
                ("NF1", &GeneMatch::exact(scoliosis.clone())),
                ("PTPN11", &GeneMatch::exact(scoliosis)),
            ]),
            row(palate.clone(), &[
                ("NF1", &GeneMatch::descendant(palate.clone())),
                ("PTPN11", &GeneMatch::exact(palate)),
            ]),
            row(hypertelorism.clone(), &[
                ("NF1", &GeneMatch::partial(0.25, hypertelorism.clone()).unwrap()),
                ("PTPN11", &GeneMatch::descendant(hypertelorism)),
            ]),
            row(macrocephaly.clone(), &[
                ("NF1", &GeneMatch::partial(0.5, macrocephaly).unwrap()),
            ]),
        ]
    }

    #[fixture]
    fn overlap_data(four_term_matches: Vec<TermMatches>) -> OverlapData {
        OverlapData::from_term_matches(genes(&["NF1", "PTPN11"]), &four_term_matches)
    }

    // ---------------------------------------------------------------------------------
    // OverlapData::from_term_matches
    // ---------------------------------------------------------------------------------

    #[rstest]
    fn counts_matches_and_partial_matches_per_gene(overlap_data: OverlapData) {
        assert_eq!(overlap_data.entities, ["NF1", "PTPN11"]);
        assert_eq!(overlap_data.n_matched_terms, [2, 3]);
        assert_eq!(overlap_data.n_partially_matched_terms, [2, 0]);
        assert_eq!(overlap_data.rows.len(), 4);
        assert_consistent(&overlap_data);
    }

    /// One row, one gene, one match of each type: only the right counter moves.
    #[rstest]
    #[case(MatchType::Exact, 1, 0)]
    #[case(MatchType::Descendant, 1, 0)]
    #[case(MatchType::Partial, 0, 1)]
    #[case(MatchType::NoMatch, 0, 0)]
    fn each_match_type_moves_the_right_counter(
        #[case] match_type: MatchType,
        #[case] expected_matched: usize,
        #[case] expected_partial: usize,
    ) {
        let rows = [row(term("HP:0002650", "Scoliosis"), &[("NF1", &match_of(match_type))])];
        let od = OverlapData::from_term_matches(genes(&["NF1"]), &rows);
        assert_eq!(od.n_matched_terms, [expected_matched]);
        assert_eq!(od.n_partially_matched_terms, [expected_partial]);
        assert_consistent(&od);
    }

    #[rstest]
    #[case(&["NF1", "PTPN11"], &["NF1", "PTPN11"])]
    #[case(&["PTPN11", "NF1"], &["NF1", "PTPN11"])] // sorted
    #[case(&["NF1", "PTPN11", "NF1"], &["NF1", "PTPN11"])] // de-duplicated
    #[case(&[], &[])]
    fn entities_are_sorted_and_unique(#[case] input: &[&str], #[case] expected: &[&str]) {
        let od = OverlapData::from_term_matches(genes(input), &[]);
        assert_eq!(od.entities, expected);
        assert_consistent(&od);
    }

    #[test]
    fn gene_missing_from_a_row_counts_as_unmatched() {
        // PTPN11 is an entity, but this row has no entry for it.
        let rows = [row(term("HP:0002650", "Scoliosis"),
                        &[("NF1", &GeneMatch::exact(term("HP:0002650", "Scoliosis")))])];
        let od = OverlapData::from_term_matches(genes(&["NF1", "PTPN11"]), &rows);
        assert_eq!(od.n_matched_terms, [1, 0]);
        assert_eq!(od.n_partially_matched_terms, [0, 0]);
    }

    #[test]
    fn match_for_an_unlisted_gene_is_ignored() {
        // The row mentions OTHER, which is not among the entities, so it must not be counted.
        let rows = [row(term("HP:0002650", "Scoliosis"),
                        &[("OTHER", &GeneMatch::exact(term("HP:0002650", "Scoliosis")))])];
        let od = OverlapData::from_term_matches(genes(&["NF1"]), &rows);
        assert_eq!(od.n_matched_terms, [0]);
    }

    #[test]
    fn no_rows_gives_zero_counts() {
        let od = OverlapData::from_term_matches(genes(&["NF1", "PTPN11"]), &[]);
        assert_eq!(od.n_matched_terms, [0, 0]);
        assert_eq!(od.n_partially_matched_terms, [0, 0]);
        assert!(od.rows.is_empty());
        assert_consistent(&od);
    }

    // ---------------------------------------------------------------------------------
    // GeneMatch::partial: every rule of the constructor, as a table
    // ---------------------------------------------------------------------------------

    /// `expected` is the resulting type, or `None` when construction must fail.
    #[rstest]
    #[case(0.5, false, Some(MatchType::Partial))]
    #[case(0.0, false, Some(MatchType::Partial))] // zero on a non-root term is still partial
    #[case(1.0, false, Some(MatchType::Partial))] // upper bound is inclusive
    #[case(0.0, true, Some(MatchType::NoMatch))] // only the root turns zero into NoMatch
    #[case(0.3, true, None)] // root with a non-zero score contradicts itself
    #[case(-0.1, false, None)]
    #[case(1.1, false, None)]
    #[case(f64::NAN, false, None)]
    #[case(f64::INFINITY, false, None)]
    fn partial_constructor_rules(
        #[case] score: f64,
        #[case] use_root: bool,
        #[case] expected: Option<MatchType>,
    ) {
        let t = if use_root { root() } else { term("HP:0000316", "Hypertelorism") };
        let result = GeneMatch::partial(score, t);
        assert_eq!(result.as_ref().ok().map(|m| m.match_type), expected, "{result:?}");
    }

    // ---------------------------------------------------------------------------------
    // Frontend contract: key names and value shapes the Angular code depends on
    // ---------------------------------------------------------------------------------

    #[rstest]
    fn overlap_data_json_keys(overlap_data: OverlapData) {
        let v = serde_json::to_value(&overlap_data).unwrap();
        let mut keys: Vec<&str> = v.as_object().unwrap().keys().map(String::as_str).collect();
        keys.sort_unstable();
        assert_eq!(keys, ["entities", "nMatchedTerms", "nPartiallyMatchedTerms", "rows"]);
        assert_eq!(v["nMatchedTerms"], json!([2, 3]));
    }

    #[test]
    fn gene_match_json_shape() {
        let m = GeneMatch::exact(term("HP:0002650", "Scoliosis"));
        assert_eq!(
            serde_json::to_value(&m).unwrap(),
            json!({
                "matchType": "exact",
                "score": 1.0,
                "matchedTerm": { "termId": "HP:0002650", "termLabel": "Scoliosis" }
            })
        );
    }

    // ---------------------------------------------------------------------------------
    // Payload for the overlap graphic
    // ---------------------------------------------------------------------------------

    #[rstest]
    fn payload_keeps_row_order_and_scores(overlap_data: OverlapData) {
        let payload = overlap_data.to_payload();
        let ids: Vec<&str> = payload.columns.iter().map(|c| c.hpo_id.as_str()).collect();
        assert_eq!(ids, ["HP:0002650", "HP:0000175", "HP:0000316", "HP:0000256"]);
        assert_eq!(payload.n_columns(), 4);
        assert_eq!(payload.entities, overlap_data.entities);
        let hypertelorism = &payload.columns[2];
        assert_eq!(hypertelorism.scores["NF1"], 0.25);
        assert_eq!(hypertelorism.scores["PTPN11"], 1.0);
        assert!(!payload.columns[3].scores.contains_key("PTPN11")); // absent, not zero
    }

    // ---------------------------------------------------------------------------------
    // Severity
    // ---------------------------------------------------------------------------------

    #[test]
    fn with_severity_sets_and_clears() {
        let tm = TermMatches::new(term("HP:0001639", "Hypertrophic cardiomyopathy"))
            .with_severity(Some(Severity::Severe));
        assert_eq!(tm.severity, Some(Severity::Severe));
        assert_eq!(tm.with_severity(None).severity, None);
    }

    #[rstest]
    #[case("HP:0012827", Some(Severity::Borderline))]
    #[case("HP:0012825", Some(Severity::Mild))]
    #[case("HP:0012826", Some(Severity::Moderate))]
    #[case("HP:0012828", Some(Severity::Severe))]
    #[case("HP:0012829", Some(Severity::Profound))]
    #[case("HP:0000001", None)]
    fn severity_from_hpo_id(#[case] id: &str, #[case] expected: Option<Severity>) {
        assert_eq!(Severity::from_hpo_id(id), expected);
    }
}
