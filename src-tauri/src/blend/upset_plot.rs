//! Data behind the phenotypic UpSet plot.
//!
//! The plot compares the HPO annotations of two or more disease models (A, B, ...) with the
//! phenotypes observed in the proband:
//!
//! 1. Each **horizontal bar** is one disease: all of its explicit HPO annotations. The part the
//!    proband also has is green ("observed"), the rest purple.
//! 2. Each **vertical bar** is an exclusive combination of diseases (A only, B only, A+B, ...):
//!    the annotations that belong to exactly those diseases. Only non-empty combinations exist.
//! 3. Within each vertical bar, green is the part of those annotations the proband has.
//!
//! A disease annotation counts as "the proband has it" when the proband has that term or a
//! descendant of it (true path rule). The bars always show the disease's own explicit
//! annotations; ancestors are never added to them.
//!
//! The pipeline is `build_upset_data` (full detail, for display and debugging) followed by
//! `UpsetPlotPayload::from_upset_data` (the counts the graphic reads).

use std::collections::{BTreeMap, HashMap, HashSet};

use ontolius::{
    ontology::{csr::FullCsrOntology, TaxonomyWalk},
    TermId,
};
use serde::{Deserialize, Serialize};

use crate::{
    blend::{
        disease_gene_entity::GeneDiseaseEntity,
        dto::{UpsetPlotPayload, ROOT_PHENOTYPE_TERM_ID, SIMPLE_ROOT_TERM},
    },
    hpo::hpo_util::ontolius_to_simple,
    model::simple_term::SimpleOntologyTerm,
};

/// A disease annotation that the proband has, with the proband term that accounts for it.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct MatchedTerm {
    /// The disease's own (explicit) annotation, as shown in the bar.
    pub term: SimpleOntologyTerm,
    /// The proband term that accounts for it: the same term if the proband has it, otherwise
    /// a descendant of it (the smallest term ID if there are several). Observed terms are
    /// non-redundant, so an exact match is never listed together with descendants.
    pub matched_by: SimpleOntologyTerm,
}

/// A group of disease annotations split by whether the proband has them.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct TermSplit {
    /// Green: annotations the proband has (exactly, or through a descendant).
    pub proband_overlapped: Vec<MatchedTerm>,
    /// Purple: annotations the proband does not have.
    pub remaining_annotated: Vec<SimpleOntologyTerm>,
}

impl TermSplit {
    /// Height of the green part.
    pub fn n_observed(&self) -> u32 {
        self.proband_overlapped.len() as u32
    }

    /// Height of the whole bar (purple + green).
    pub fn n_annotated(&self) -> u32 {
        (self.proband_overlapped.len() + self.remaining_annotated.len()) as u32
    }
}

/// One **horizontal** bar: all explicit annotations of one disease/gene entity, including
/// terms it shares with the others.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetGene {
    pub gene_symbol: String,
    #[serde(flatten)]
    pub terms: TermSplit,
}

/// One **vertical** bar: the annotations that belong to exactly the entities in
/// `gene_symbol_list` and to no other entity. Each annotation appears in at most one
/// intersection.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetIntersection {
    pub gene_symbol_list: Vec<String>,
    #[serde(flatten)]
    pub terms: TermSplit,
}

/// Everything behind one UpSet plot.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetData {
    /// Horizontal bars, sorted by gene symbol.
    pub genes: Vec<UpsetGene>,
    /// Vertical bars, sorted by number of entities and then alphabetically.
    pub intersections: Vec<UpsetIntersection>,
    /// Proband terms that no disease annotation accounts for (no disease term equals the
    /// term or is one of its ancestors). They are not drawn in the plot.
    pub proband_only: Vec<SimpleOntologyTerm>,
    /// Number of distinct observed terms of the proband.
    pub n_proband_terms: usize,
}

/// Terms implied by every phenotype annotation (the HPO root terms); sharing them says nothing.
fn is_uninformative(t: &TermId) -> bool {
    *t == *ROOT_PHENOTYPE_TERM_ID || t.to_string() == "HP:0000001"
}

/// Term with label; falls back to the root term if the ontology does not know it.
pub(crate) fn simple_term(hpo: &FullCsrOntology, tid: &TermId) -> SimpleOntologyTerm {
    ontolius_to_simple(hpo, tid).unwrap_or_else(|e| {
        log::error!("{e}; substituting the root term");
        SIMPLE_ROOT_TERM.clone()
    })
}

/// One disease annotation and where it belongs.
struct UpsetRow {
    term: SimpleOntologyTerm,
    /// Entities annotated to exactly this term, sorted.
    genes: Vec<String>,
    /// The proband term that accounts for it; `None` if the proband does not have it.
    matched_by: Option<SimpleOntologyTerm>,
}

fn split<'a>(rows: impl Iterator<Item = &'a UpsetRow>) -> TermSplit {
    let mut s = TermSplit::default();
    for r in rows {
        match &r.matched_by {
            Some(m) => s.proband_overlapped.push(MatchedTerm {
                term: r.term.clone(),
                matched_by: m.clone(),
            }),
            None => s.remaining_annotated.push(r.term.clone()),
        }
    }
    s // rows are sorted by term ID, so both lists are too
}

/// Builds the UpSet data.
///
/// `observed` are the proband's terms (non-redundant). `expand` returns a term together with
/// its ancestors; `label` makes the displayable term. Both are parameters so the logic can be
/// tested without an ontology.
pub fn upset_data(
    observed: &[TermId],
    entities: &[GeneDiseaseEntity],
    expand: impl Fn(&TermId) -> Vec<TermId>,
    label: impl Fn(&TermId) -> SimpleOntologyTerm,
) -> UpsetData {
    let mut proband: Vec<TermId> = observed.to_vec();
    proband.sort_by_key(|t| t.to_string());
    proband.dedup();

    // term (a proband term or any of its ancestors) -> the proband terms it accounts for;
    // the lists are sorted because `proband` is
    let mut accounted_by: HashMap<TermId, Vec<TermId>> = HashMap::new();
    for p in &proband {
        let mut implied: HashSet<TermId> = expand(p).into_iter().collect();
        implied.insert(p.clone());
        for a in implied {
            accounted_by.entry(a).or_default().push(p.clone());
        }
    }

    // one entry per distinct disease annotation, with the entities that carry it
    let mut by_term: HashMap<&TermId, Vec<&str>> = HashMap::new();
    for e in entities {
        for t in e.disease_hpo_ids.iter().filter(|t| !is_uninformative(t)) {
            by_term.entry(t).or_default().push(e.gene_symbol.as_str());
        }
    }

    let mut explained: HashSet<TermId> = HashSet::new(); // proband terms matching some annotation
    let mut rows: Vec<UpsetRow> = Vec::new();
    for (t, mut genes) in by_term {
        genes.sort();
        genes.dedup();
        let ps = accounted_by.get(t);
        explained.extend(ps.into_iter().flatten().cloned());
        let matched_by = ps.and_then(|v| v.first()).map(|p| label(p));
        rows.push(UpsetRow {
            term: label(t),
            genes: genes.into_iter().map(String::from).collect(),
            matched_by,
        });
    }
    rows.sort_by_key(|r| r.term.term_id_as_string());

    // horizontal bars
    let mut symbols: Vec<String> = entities.iter().map(|e| e.gene_symbol.clone()).collect();
    symbols.sort();
    symbols.dedup();
    let genes: Vec<UpsetGene> = symbols
        .iter()
        .map(|g| UpsetGene {
            gene_symbol: g.clone(),
            terms: split(rows.iter().filter(|r| r.genes.contains(g))),
        })
        .collect();

    // vertical bars: group the annotations by the exact set of entities that carry them
    let mut regions: BTreeMap<Vec<String>, Vec<&UpsetRow>> = BTreeMap::new();
    for r in &rows {
        regions.entry(r.genes.clone()).or_default().push(r);
    }
    let mut intersections: Vec<UpsetIntersection> = regions
        .into_iter()
        .map(|(genes, rs)| UpsetIntersection {
            gene_symbol_list: genes,
            terms: split(rs.into_iter()),
        })
        .collect();
    intersections.sort_by(|a, b| {
        a.gene_symbol_list
            .len()
            .cmp(&b.gene_symbol_list.len())
            .then_with(|| a.gene_symbol_list.cmp(&b.gene_symbol_list))
    });

    let proband_only: Vec<SimpleOntologyTerm> = proband
        .iter()
        .filter(|p| !explained.contains(*p))
        .map(|p| label(p))
        .collect();

    UpsetData {
        genes,
        intersections,
        proband_only,
        n_proband_terms: proband.len(),
    }
}

/// Verifies that no observed term is an ancestor of another observed term.
///
/// # Errors
/// Names the first redundant term and the descendant that makes it redundant.
pub fn check_no_redundant_terms(
    observed: &[TermId],
    expand: &impl Fn(&TermId) -> Vec<TermId>,
) -> Result<(), String> {
    for q in observed {
        let ancestors = expand(q);
        if let Some(p) = observed.iter().find(|p| *p != q && ancestors.contains(p)) {
            return Err(format!(
                "observed term {p} is redundant: its descendant {q} is also observed"
            ));
        }
    }
    Ok(())
}

/// Builds and validates the UpSet data using the real ontology.
///
/// # Errors
/// Fails if the observed terms are redundant or if the resulting bars are inconsistent.
pub fn build_upset_data(
    entities: &[GeneDiseaseEntity],
    observed: &[TermId],
    hpo: &FullCsrOntology,
) -> Result<UpsetData, String> {
    let expand = |t: &TermId| -> Vec<TermId> { hpo.iter_term_and_ancestor_ids(t).cloned().collect() };
    check_no_redundant_terms(observed, &expand)?;
    let data = upset_data(observed, entities, &expand, |t| simple_term(hpo, t));
    data.validate()?;
    Ok(data)
}



impl UpsetPlotPayload {
    /// The counts-only view the graphic reads: `*_annotated` is the full bar height
    /// (purple + green), `*_observed` the green part.
    pub fn from_upset_data(d: &UpsetData) -> Self {
        Self {
            genes: d.genes.iter().map(|g| g.gene_symbol.clone()).collect(),
            gene_annotated: d.genes.iter().map(|g| g.terms.n_annotated()).collect(),
            gene_observed: d.genes.iter().map(|g| g.terms.n_observed()).collect(),
            combinations: d.intersections.iter().map(|i| i.gene_symbol_list.clone()).collect(),
            combination_annotated: d.intersections.iter().map(|i| i.terms.n_annotated()).collect(),
            combination_observed: d.intersections.iter().map(|i| i.terms.n_observed()).collect(),
        }
    }

  
}

impl UpsetData {
    /// Each entity's totals must equal the sum over the intersections containing it, and no
    /// annotation may sit in two intersections.
    pub fn validate(&self) -> Result<(), String> {
        for g in &self.genes {
            let (mut annotated, mut observed) = (0, 0);
            for i in self
                .intersections
                .iter()
                .filter(|i| i.gene_symbol_list.contains(&g.gene_symbol))
            {
                annotated += i.terms.n_annotated();
                observed += i.terms.n_observed();
            }
            if annotated != g.terms.n_annotated() || observed != g.terms.n_observed() {
                return Err(format!(
                    "{}: horizontal bar ({}/{} observed/annotated) != sum of vertical bars ({observed}/{annotated})",
                    g.gene_symbol,
                    g.terms.n_observed(),
                    g.terms.n_annotated(),
                ));
            }
        }
        let mut seen = HashSet::new();
        for i in &self.intersections {
            let ids = i
                .terms
                .proband_overlapped
                .iter()
                .map(|m| m.term.term_id_as_string())
                .chain(i.terms.remaining_annotated.iter().map(|t| t.term_id_as_string()));
            for id in ids {
                if !seen.insert(id.clone()) {
                    return Err(format!("{id} appears in more than one intersection"));
                }
            }
        }
        Ok(())
    }

    /// Tab-separated table, one line per annotation, to check each step by eye or in a
    /// spreadsheet. Patient-only terms come last.
    pub fn to_tsv(&self) -> String {
        let mut out = String::from("term_id\tlabel\tintersection\tobserved\tmatched_by\n");
        for i in &self.intersections {
            let set = i.gene_symbol_list.join("+");
            for m in &i.terms.proband_overlapped {
                out.push_str(&format!(
                    "{}\t{}\t{set}\ttrue\t{}\n",
                    m.term.term_id_as_string(),
                    m.term.term_label,
                    m.matched_by.term_id_as_string()
                ));
            }
            for t in &i.terms.remaining_annotated {
                out.push_str(&format!(
                    "{}\t{}\t{set}\tfalse\t\n",
                    t.term_id_as_string(),
                    t.term_label
                ));
            }
        }
        for t in &self.proband_only {
            out.push_str(&format!(
                "{}\t{}\t(proband only)\ttrue\t\n",
                t.term_id_as_string(),
                t.term_label
            ));
        }
        out
    }

    /// The payload the graphic reads.
    ///
    /// # Errors
    /// See [`build_upset_data`].
    pub fn build_upset_payload(
        &self
    ) -> Result<UpsetPlotPayload, String> {
        Ok(UpsetPlotPayload::from_upset_data(self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tid(s: &str) -> TermId {
        s.parse().unwrap()
    }

    // Zonular cataract -> Cataract -> Abnormality of the eye -> root -> All
    // Seizure -> Abnormality of the nervous system -> root -> All
    fn fake_expand(t: &TermId) -> Vec<TermId> {
        let chain: &[&str] = match t.to_string().as_str() {
            "HP:0010920" => &["HP:0010920", "HP:0000518", "HP:0000478", "HP:0000118", "HP:0000001"],
            "HP:0000518" => &["HP:0000518", "HP:0000478", "HP:0000118", "HP:0000001"],
            "HP:0001250" => &["HP:0001250", "HP:0000707", "HP:0000118", "HP:0000001"],
            _ => return vec![t.clone()],
        };
        chain.iter().map(|s| tid(s)).collect()
    }

    fn label(t: &TermId) -> SimpleOntologyTerm {
        SimpleOntologyTerm::from_term_id(t, t.to_string())
    }

    fn entity(symbol: &str, ids: &[&str]) -> GeneDiseaseEntity {
        GeneDiseaseEntity {
            gene_symbol: symbol.to_string(),
            disease_hpo_ids: ids.iter().map(|s| tid(s)).collect(),
            ncbi_gene_id: "0".to_string(),
            disease_list: vec![],
        }
    }

    fn data(observed: &[&str], ents: &[GeneDiseaseEntity]) -> UpsetData {
        let obs: Vec<TermId> = observed.iter().map(|s| tid(s)).collect();
        let d = upset_data(&obs, ents, fake_expand, label);
        d.validate().unwrap();
        d
    }

    fn combo(g: &[&str]) -> Vec<String> {
        g.iter().map(|s| s.to_string()).collect()
    }

    fn ids(v: &[SimpleOntologyTerm]) -> Vec<String> {
        v.iter().map(|t| t.term_id_as_string()).collect()
    }

    fn pairs(v: &[MatchedTerm]) -> Vec<String> {
        v.iter()
            .map(|m| format!("{}<-{}", m.term.term_id_as_string(), m.matched_by.term_id_as_string()))
            .collect()
    }

    #[test]
    fn two_diseases_a_b_and_shared() {
        let ents = [
            entity("A", &["HP:0000010", "HP:0000020", "HP:0000030"]),
            entity("B", &["HP:0000030", "HP:0000040"]),
        ];
        let d = data(&["HP:0000010", "HP:0000030", "HP:0000060"], &ents);
        let p = UpsetPlotPayload::from_upset_data(&d);

        assert_eq!(p.genes, ["A", "B"]);
        assert_eq!(p.gene_annotated, [3, 2]);
        assert_eq!(p.gene_observed, [2, 1]);
        assert_eq!(p.combinations, [combo(&["A"]), combo(&["B"]), combo(&["A", "B"])]);
        assert_eq!(p.combination_annotated, [2, 1, 1]); // A only: 10,20; B only: 40; both: 30
        assert_eq!(p.combination_observed, [1, 0, 1]);
        assert_eq!(ids(&d.proband_only), ["HP:0000060"]);
        assert_eq!(d.n_proband_terms, 3);
    }

    #[test]
    fn three_diseases_only_non_empty_regions_appear() {
        let ents = [
            entity("A", &["HP:0000010", "HP:0000020", "HP:0000030"]),
            entity("B", &["HP:0000030", "HP:0000040"]),
            entity("C", &["HP:0000030", "HP:0000050"]),
        ];
        let d = data(&["HP:0000010", "HP:0000030", "HP:0000060"], &ents);
        let p = UpsetPlotPayload::from_upset_data(&d);
        assert_eq!(
            p.combinations,
            [combo(&["A"]), combo(&["B"]), combo(&["C"]), combo(&["A", "B", "C"])]
        );
        assert_eq!(p.combination_annotated, [2, 1, 1, 1]);
        assert_eq!(p.combination_observed, [1, 0, 0, 1]);
        assert_eq!(p.gene_annotated, [3, 2, 2]);
        assert_eq!(p.gene_observed, [2, 1, 1]);
    }

    #[test]
    fn patient_descendant_counts_as_a_match_on_the_disease_term() {
        // A and B both annotated to Cataract; the patient has Zonular cataract
        let ents = [
            entity("A", &["HP:0000518", "HP:0000020"]),
            entity("B", &["HP:0000518", "HP:0001250"]),
        ];
        let d = data(&["HP:0010920"], &ents);

        // the term shown is the disease's own (Cataract), matched by the patient's descendant
        assert_eq!(pairs(&d.genes[0].terms.proband_overlapped), ["HP:0000518<-HP:0010920"]);
        assert_eq!(ids(&d.genes[0].terms.remaining_annotated), ["HP:0000020"]);
        assert_eq!(ids(&d.genes[1].terms.remaining_annotated), ["HP:0001250"]);

        let sets: Vec<_> = d.intersections.iter().map(|i| i.gene_symbol_list.join("+")).collect();
        assert_eq!(sets, ["A", "B", "A+B"]);
        assert!(d.intersections[0].terms.proband_overlapped.is_empty()); // A only: all purple
        assert_eq!(
            pairs(&d.intersections[2].terms.proband_overlapped),
            ["HP:0000518<-HP:0010920"]
        );
        assert!(d.proband_only.is_empty()); // Zonular cataract is accounted for by Cataract
    }

    #[test]
    fn the_patient_having_the_general_term_does_not_make_a_specific_disease_term_green() {
        let d = data(&["HP:0000518"], &[entity("A", &["HP:0010920"])]);
        assert!(d.genes[0].terms.proband_overlapped.is_empty());
        assert_eq!(ids(&d.proband_only), ["HP:0000518"]);
    }

    #[test]
    fn one_patient_term_can_turn_several_disease_terms_green() {
        let d = data(&["HP:0010920"], &[entity("A", &["HP:0000518", "HP:0000478"])]);
        assert_eq!(
            pairs(&d.genes[0].terms.proband_overlapped),
            ["HP:0000478<-HP:0010920", "HP:0000518<-HP:0010920"]
        );
        // so green may exceed the number of observed terms
        assert!(d.genes[0].terms.n_observed() as usize > d.n_proband_terms);
    }

    #[test]
    fn root_terms_are_never_shown() {
        let ents = [entity("A", &["HP:0000118", "HP:0000518"]), entity("B", &["HP:0000118"])];
        let d = data(&["HP:0010920"], &ents);
        let shown: Vec<String> = d
            .intersections
            .iter()
            .flat_map(|i| {
                i.terms
                    .proband_overlapped
                    .iter()
                    .map(|m| m.term.term_id_as_string())
                    .chain(i.terms.remaining_annotated.iter().map(|t| t.term_id_as_string()))
            })
            .collect();
        assert_eq!(shown, ["HP:0000518"]);
        assert_eq!(d.genes.len(), 2); // B still has a bar (empty)
        assert_eq!(d.genes[1].terms.n_annotated(), 0);
    }

    #[test]
    fn validate_catches_inconsistent_bars() {
        let ents = [entity("A", &["HP:0000010", "HP:0000030"]), entity("B", &["HP:0000030"])];
        let good = data(&["HP:0000010"], &ents);

        let mut missing = good.clone();
        missing.intersections[0].terms.proband_overlapped.pop(); // drop A's green term from [A]
        assert!(missing.validate().is_err());

        let mut duplicated = good.clone();
        let t = duplicated.intersections[0].terms.proband_overlapped[0].term.clone();
        duplicated.intersections[1].terms.remaining_annotated.push(t);
        assert!(duplicated.validate().is_err());
    }

    #[test]
    fn output_does_not_depend_on_entity_order() {
        let a = entity("A", &["HP:0000518", "HP:0001250"]);
        let b = entity("B", &["HP:0001250"]);
        let fwd = data(&["HP:0010920"], &[a.clone(), b.clone()]);
        let rev = data(&["HP:0010920"], &[b, a]);
        assert_eq!(serde_json::to_value(&fwd).unwrap(), serde_json::to_value(&rev).unwrap());
    }

    #[test]
    fn payload_json_keys_match_the_frontend_contract() {
        let ents = [entity("A", &["HP:0000010"])];
        let p = UpsetPlotPayload::from_upset_data(&data(&["HP:0000010"], &ents));
        let v = serde_json::to_value(&p).unwrap();
        let mut keys: Vec<&str> = v.as_object().unwrap().keys().map(String::as_str).collect();
        keys.sort();
        assert_eq!(
            keys,
            [
                "combinationAnnotated",
                "combinationObserved",
                "combinations",
                "geneAnnotated",
                "geneObserved",
                "genes"
            ]
        );
    }

    #[test]
    fn upset_data_json_is_flat_and_camel_case() {
        let d = data(&["HP:0010920"], &[entity("A", &["HP:0000518"])]);
        let v = serde_json::to_value(&d).unwrap();
        assert_eq!(v["nProbandTerms"], 1);
        let g = &v["genes"][0];
        assert_eq!(g["geneSymbol"], "A");
        assert_eq!(g["probandOverlapped"][0]["term"]["termId"], "HP:0000518");
        assert_eq!(g["probandOverlapped"][0]["matchedBy"]["termId"], "HP:0010920");
        assert!(g.get("terms").is_none()); // flattened
    }

    #[test]
    fn no_entities_gives_an_empty_payload() {
        let d = data(&["HP:0010920"], &[]);
        let p = UpsetPlotPayload::from_upset_data(&d);
        assert!(p.genes.is_empty() && p.combinations.is_empty());
        assert_eq!(ids(&d.proband_only), ["HP:0010920"]);
    }

    #[test]
    fn redundant_observed_terms_are_detected() {
        let ok = [tid("HP:0010920"), tid("HP:0001250")];
        assert!(check_no_redundant_terms(&ok, &fake_expand).is_ok());
        let bad = [tid("HP:0000518"), tid("HP:0010920")]; // Cataract + its descendant
        let err = check_no_redundant_terms(&bad, &fake_expand).unwrap_err();
        assert!(err.contains("HP:0000518") && err.contains("HP:0010920"));
    }

    #[test]
    fn tsv_has_one_line_per_annotation_plus_patient_only_terms() {
        let ents = [entity("A", &["HP:0000518", "HP:0000020"])];
        let d = data(&["HP:0010920", "HP:0001250"], &ents);
        let tsv = d.to_tsv();
        assert!(tsv.starts_with("term_id\tlabel\tintersection\tobserved\tmatched_by\n"));
        assert_eq!(tsv.lines().count(), 1 + 2 + 1); // header, two annotations, one patient-only term
        assert!(tsv.contains("HP:0000518\tHP:0000518\tA\ttrue\tHP:0010920"));
        assert!(tsv.contains("(proband only)"));
    }
}