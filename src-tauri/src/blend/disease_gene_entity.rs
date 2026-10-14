use std::{collections::{HashMap, HashSet},  sync::Arc};

use ontolius::{TermId, ontology::{TaxonomyQuery, OntologyTerms, csr::FullCsrOntology}, term::MinimalTerm};
use serde::{Serialize, Serializer};
use crate::{blend::dto::{GeneMatch, OverlapPlotPayload, ROOT_PHENOTYPE_TERM_ID, SIMPLE_ROOT_TERM, SpreadPlotPayload, TermMatches, UpsetPlotPayload}, hpo::hpo_util::ontolius_to_simple, hpoa::disease_model::GeneDiseaseAssociation, model::{proband::Proband, simple_term::SimpleOntologyTerm}};
use log::{trace};

fn serialize_term_id_set<S>(ids: &HashSet<TermId>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let mut strs: Vec<String> = ids.iter().map(|t| t.to_string()).collect();
    strs.sort(); // HashSet iteration order is random; sort for stable output
    strs.serialize(serializer)
}

/// Picks the broadest annotated ancestor and computes the partial-match score from it.
///
/// `ancestors` are the disease terms that are ancestors of the query term, each paired
/// with its disease count. Ancestors with a count of 0 are ignored. The broadest one
/// (highest count) is chosen; ties are broken by the smaller term ID so the result does
/// not depend on iteration order. Returns `None` if `pheno_count` is 0 or no ancestor
/// has a positive count.
fn best_ancestor_score(
    pheno_count: usize,
    ancestors: impl IntoIterator<Item = (TermId, usize)>,
) -> Option<(f64, TermId)> {
    if pheno_count == 0 {
        return None;
    }
    let (term, max_count) = ancestors
        .into_iter()
        .filter(|(_, count)| *count > 0)
        .max_by(|(ta, ca), (tb, cb)| {
            ca.cmp(cb).then_with(|| tb.to_string().cmp(&ta.to_string()))
        })?;
    Some((pheno_count as f64 / max_count as f64, term))
}

/// Scores a *partial* match between a query phenotype and the HPO annotations
/// of a gene's disease.
///
/// A partial match means the disease is not annotated with the query term itself but
/// with one or more of its ancestors (broader terms). The score measures how much more
/// specific the query term is than the broadest such ancestor, as a ratio of disease
/// counts:
///
/// ```text
/// score = disease_counts[query_hpo] / max(disease_counts[a])
/// ```
///
/// where `a` ranges over the terms in `gde.disease_hpo_ids` that are ancestors of
/// `query_hpo`, excluding the root `HP:0000118`. Taking the maximum means the broadest
/// (most common) ancestor determines the score, which is the most conservative value
/// when several ancestors match. That ancestor is also returned, so callers can record
/// which annotation produced the score. If two ancestors tie on count, the one with the
/// smaller term ID is returned.
///
/// Interpretation: a score near 1 means the ancestor is about as specific as the query
/// term, so the match is nearly as good as an exact one. A score near 0 means the
/// ancestor is far more general. For example, if the query term is annotated to 5
/// diseases and the broadest matching ancestor is annotated to 50, the score is `0.1`.
///
/// # Arguments
///
/// * `query_hpo` - The phenotype term being scored.
/// * `gde` - The gene-disease entity whose `disease_hpo_ids` are searched for
///   ancestors of `query_hpo`.
/// * `hpo` - The ontology, used for ancestor queries.
/// * `disease_counts` - Number of diseases annotated with each HPO term (with the true
///   path rule applied). Terms missing from the map count as `0`.
///
/// # Returns
///
/// A pair `(score, term)`:
///
/// * `score` in `(0.0, 1.0]` and the broadest matching ancestor, when a partial match
///   exists and `disease_counts` is propagated (see Notes);
/// * `(0.0, HP:0000118)` when there is no partial match: the query term is not annotated
///   to any disease, or no disease term other than the root is an ancestor of it.
///
/// # Notes
///
/// * Only ancestors count. A disease term identical to `query_hpo` is not considered
///   here, so exact and descendant matches must be handled by the caller first.
/// * The root `HP:0000118` is never returned as a match, even if the disease is annotated
///   to it, because sharing only the root means no meaningful match.
/// * The `(0.0, 1.0]` bound assumes `disease_counts` is propagated up the hierarchy, so
///   an ancestor's count is never smaller than that of its descendants. With
///   direct-annotation counts, or counts that are out of date, the ratio can exceed 1.
pub fn partial_match_score(
    query_hpo: &TermId,
    gde: &GeneDiseaseEntity,
    hpo: &FullCsrOntology,
    disease_counts: &HashMap<TermId, usize>,
) -> (f64, TermId) {
    let pheno_count = disease_counts.get(query_hpo).copied().unwrap_or(0);
    let ancestors = gde
        .disease_hpo_ids
        .iter()
        .filter(|d| **d != *ROOT_PHENOTYPE_TERM_ID && hpo.is_descendant_of(query_hpo, *d))
        .map(|d| (d.clone(), disease_counts.get(d).copied().unwrap_or(0)));

    let (score, term) = best_ancestor_score(pheno_count, ancestors)
        .unwrap_or_else(|| (0.0, ROOT_PHENOTYPE_TERM_ID.clone()));
    trace!("partial_match_score query={query_hpo} score={score} via={term}");
    (score, term)
}

/// Create a SimpleTerm from a TermId. This should never fail because we are retrieving terms that we know are in the Ontology,
/// but if for some reason it does fail, log it and return the Root.
pub fn simple_term(hpo: &FullCsrOntology, tid: &TermId) -> SimpleOntologyTerm {
    ontolius_to_simple(hpo, tid).unwrap_or_else(|e| {
        log::error!("{e}; substituting root term {}", *ROOT_PHENOTYPE_TERM_ID);
        SIMPLE_ROOT_TERM.clone()
    })
}

/// Returns `Some` for an exact or descendant match, `None` when only a partial match is possible.
fn classify_perfect_match(
    ppkt_tid: &TermId,
    gde: &GeneDiseaseEntity,
    hpo: &FullCsrOntology,
) -> Option<GeneMatch> {
    if gde.disease_hpo_ids.contains(ppkt_tid) {
        return Some(GeneMatch::exact(simple_term(hpo, ppkt_tid)));
    }
    let descendants: Vec<&TermId> = gde
        .disease_hpo_ids
        .iter()
        .filter(|id| hpo.is_descendant_of(*id, ppkt_tid))
        .collect();
    descendants
        .iter()
        .copied()
        // closest: no other matching term lies between this one and the proband term
        .filter(|d| !descendants.iter().any(|o| *o != *d && hpo.is_descendant_of(*d, *o)))
        .min_by_key(|t| t.to_string())
        .map(|t| GeneMatch::descendant(simple_term(hpo, t)))
}

/// Computes, for every distinct observed proband term, how it matches each gene entity.
///
/// Columns follow the order of `proband.observed_hpos`; duplicate terms are skipped.
///
/// # Errors
/// Returns an error if a partial-match score is invalid (see [`GeneMatch::partial`]),
/// naming the proband term and gene.
pub fn compute_term_matches(
    proband: &Proband,
    gde_list: &[GeneDiseaseEntity],
    disease_counts: &HashMap<TermId, usize>,
    hpo: &FullCsrOntology,
) -> Result<Vec<TermMatches>, String> {
    let mut seen = HashSet::new();
    let mut all = Vec::new();
    for tid in proband.observed_hpos.iter().filter(|t| seen.insert((*t).clone())) {
        let mut tm = TermMatches::new(crate::hpo::hpo_util::ontolius_to_simple(hpo, tid)?);
        for gde in gde_list {
            let m = match classify_perfect_match(tid, gde, hpo) {
                Some(m) => m,
                None => {
                    let (score, via) = partial_match_score(tid, gde, hpo, disease_counts);
                    GeneMatch::partial(score, crate::hpo::hpo_util::ontolius_to_simple(hpo, &via)?)
                        .map_err(|e| format!("{tid} / {}: {e}", gde.gene_symbol))?
                }
            };
            tm.insert(&gde.gene_symbol, m);
        }
        all.push(tm);
    }
    Ok(all)
}

/// This is the structure we use for the phenoblend analysis. We have one such entity for each gene, and the entities can have one or more diseases associated with the gene
/// disease_hpo_ids then has all the observed HPO identifiers for all of the diseases
/// individual_explicit_ids has all of the explicit annotations from the phenopacket/case we are analyzing.
#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct GeneDiseaseEntity {
    pub ncbi_gene_id: String,
    pub gene_symbol: String,
    pub disease_list: Vec<SimpleOntologyTerm>,
    #[serde(serialize_with = "serialize_term_id_set")]
    pub(crate) disease_hpo_ids: HashSet<TermId> 
}

impl GeneDiseaseEntity {
    pub fn new(gda_list: &Vec<GeneDiseaseAssociation>) -> Result<Self, String> {
        let first_gda = gda_list.first().ok_or("gda_list cannot be empty")?;
        let ncbi_gene_id = first_gda.ncbi_gene_id.clone();
        let gene_symbol = first_gda.gene_symbol.clone();
        
        let mut disease_list = Vec::with_capacity(gda_list.len());
        let mut disease_hpo_ids = HashSet::new();

        for d in gda_list {
            if d.ncbi_gene_id != ncbi_gene_id || d.gene_symbol != gene_symbol {
                return Err("gda_list contains mismatched Gene id/symbol".to_string());
            }
            let model = d.disease_model.clone()
                .ok_or_else(|| format!("Could not find model for gene {}", gene_symbol))?;
            disease_hpo_ids.extend(model.observed_hpo_ids);
            let simple_disease_term = SimpleOntologyTerm::new(
                model.omim_disease_id.to_string(),
                model.omim_disease_name.clone()
            )?;
            disease_list.push(simple_disease_term);
        }
        Ok(Self {
            ncbi_gene_id: ncbi_gene_id,
            gene_symbol: gene_symbol,
            disease_list: disease_list,
            disease_hpo_ids: disease_hpo_ids
        })
    }


    fn debug_print_term(tid: &TermId,
        hpo: Arc<FullCsrOntology>, message: &str) {
            match hpo.term_by_id(tid) {
                Some(t) => { trace!("{}: {} ({})", message, t.name(), tid);},
                None => { trace!("{}: Could not find label for {}", message, tid); }
            }
        }


    /// Returns `true` if the gene-disease entity is implicitly or explicitly annotated
    /// to `ppkt_tid`.
    ///
    /// By the true path rule, an annotation to a descendant also counts as an annotation
    /// to the ancestor. For example, a disease annotated to "Nuclear cataract" is also
    /// annotated to "Cataract", so a proband with "Cataract" is a perfect match.
     pub fn is_perfect_match(
        ppkt_tid: &TermId,
        gde: &GeneDiseaseEntity,
        hpo: Arc<FullCsrOntology>) -> bool {
            let matched = if gde.disease_hpo_ids.contains(ppkt_tid) {
                Some(ppkt_tid)
            } else {
                gde.disease_hpo_ids
                    .iter()
                    .find(|id| hpo.is_descendant_of(*id, ppkt_tid))
            };
           matched.is_some()
    }



    pub fn get_overlap_matrix_payload(
        proband: Proband, 
        gde_list: &Vec<GeneDiseaseEntity>, 
        disease_counts: &HashMap<TermId, usize>,
        hpo: Arc<FullCsrOntology>
    ) -> Result<OverlapPlotPayload, String> {
        let tm = compute_term_matches(&proband, gde_list, disease_counts, hpo.as_ref())?;
            Ok(OverlapPlotPayload::from_term_matches(
                gde_list.iter().map(|g| g.gene_symbol.clone()),
                &tm,
            ))
    }

    pub fn get_upset_plot_payloadOLD(
        proband: Proband, 
        gda_map: &HashMap<String, Vec<GeneDiseaseAssociation>>,
        hpo: Arc<FullCsrOntology>
    ) -> Result<UpsetPlotPayload, String> {

       /* let mut gd_entry_list: Vec<GeneDiseaseEntity> = Vec::new();
        for (symbol, gda_list) in gda_map.into_iter() {
            let gd_entity =   GeneDiseaseEntity::new(gda_list)?;
            gd_entry_list.push(gd_entity);
        }
        Ok(crate::blend::upset_plot::build_upset_payload(&gd_entry_list, &proband, hpo.as_ref())?)
        */ 
        Err("Needs refactor".to_ascii_lowercase())
    }


    pub fn get_spread_plot_payload(
        proband: Proband, 
        gda_map: &HashMap<String, Vec<GeneDiseaseAssociation>>,
        hpo: Arc<FullCsrOntology>
    ) -> Result<SpreadPlotPayload, String> {
        let mut gd_entry_list: Vec<GeneDiseaseEntity> = Vec::new();
        for (symbol, gda_list) in gda_map.into_iter() {
            let gd_entity =   GeneDiseaseEntity::new(gda_list)?;
            gd_entry_list.push(gd_entity);
        }
       crate::blend::spread_plot::get_spread_plot_payload(proband, gda_map, hpo)
    }


   


}


#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::sync::Arc;
    use std::collections::HashMap;
    use ontolius::TermId;
    use crate::tests::hpo_fixture::hpo;
    

    #[rstest::fixture]
    pub fn brca_disease() -> SimpleOntologyTerm {
        SimpleOntologyTerm::new("OMIM:604370", "Breast-ovarian cancer, familial, 1").unwrap()
    }

    #[rstest]
    #[case(0, 0, 0.0)]     // Zero proband count -> should be 0.0
    #[case(10, 0, 0.0)]    // Zero descendant count -> should be 0.0
    #[case(5, 10, 0.5)]    // Normal calculation: 5 / 10 = 0.5
    fn test_partial_match_score_edge_cases(
        hpo: Arc<FullCsrOntology>,
        #[case] pheno_count: usize,
        #[case] max_geno_count: usize,
        #[case] expected_score: f64,
    ) {
        // Mock or construct your inputs here. 
        // Since `partial_match_score` relies on `disease_counts` and `hpo`, 
        // you can populate a mock disease_counts map to mimic these counts.
        
        let mut disease_counts = HashMap::new();
        let disease_hpo: TermId = "HP:0000518".parse().unwrap(); // Cataract  HP:0000518,  general ancestor term
        let query_hpo : TermId = "HP:0010920".parse().unwrap(); // Zonular cataract HP:0010920
        assert!(
            hpo.is_descendant_of(&query_hpo, &disease_hpo),
            "fixture assumption: HP:0010920 is a descendant of HP:0000518"
        );

        if pheno_count > 0 {
            disease_counts.insert(query_hpo.clone(), pheno_count);
        }
        if max_geno_count > 0 {
            disease_counts.insert(disease_hpo.clone(), max_geno_count);
        }
        let bca = SimpleOntologyTerm::new("OMIM:604370", "Breast-ovarian cancer, familial, 1").unwrap();
        let hpo_set = HashSet::from([
            (disease_hpo),
        ]);
        let gde = GeneDiseaseEntity {
            gene_symbol:"BRCA1".to_string(),
            disease_hpo_ids: hpo_set, 
            ncbi_gene_id: "672".to_string(), 
            disease_list: vec![bca]
        };


        let score = partial_match_score(
            &query_hpo,
            &gde,
            hpo.as_ref(),
            &disease_counts,
        );

        assert!((score.0 - expected_score).abs() < f64::EPSILON);
        assert!(score.0.is_finite(), "Score must never be INF or NaN!");
    }

    fn tid(s: &str) -> TermId { s.parse().unwrap() }

    #[test]
    fn broadest_ancestor_determines_score() {
        let anc = vec![(tid("HP:0000002"), 50), (tid("HP:0000003"), 20)];
        let (score, term) = best_ancestor_score(5, anc).unwrap();
        assert!((score - 0.1).abs() < 1e-12); // 5 / 50
        assert_eq!(term, tid("HP:0000002"));
    }

    #[test]
    fn ties_are_broken_by_smaller_term_id_regardless_of_order() {
        let a = vec![(tid("HP:0000002"), 10), (tid("HP:0000001"), 10)];
        let b = vec![(tid("HP:0000001"), 10), (tid("HP:0000002"), 10)];
        assert_eq!(best_ancestor_score(5, a).unwrap().1, tid("HP:0000001"));
        assert_eq!(best_ancestor_score(5, b).unwrap().1, tid("HP:0000001"));
    }

    #[test]
    fn no_score_without_query_count_or_usable_ancestor() {
        assert!(best_ancestor_score(0, vec![(tid("HP:0000002"), 10)]).is_none());
        assert!(best_ancestor_score(5, Vec::new()).is_none());
        assert!(best_ancestor_score(5, vec![(tid("HP:0000002"), 0)]).is_none());
    }

    #[test]
    fn score_of_one_when_counts_equal() {
        let (score, _) = best_ancestor_score(7, vec![(tid("HP:0000002"), 7)]).unwrap();
        assert_eq!(score, 1.0);
    }
}