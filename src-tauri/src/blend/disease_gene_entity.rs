use std::{collections::{HashMap, HashSet},  sync::Arc};

use ontolius::{TermId, ontology::{TaxonomyQuery, OntologyTerms, csr::FullCsrOntology}, term::MinimalTerm};
use serde::{Serialize, Serializer};
use crate::{blend::dto::{OverlapPlotItem, OverlapPlotPayload, SpreadPlotPayload, UpsetPlotPayload}, hpoa::disease_model::GeneDiseaseAssociation, model::{proband::Proband, simple_term::SimpleOntologyTerm}};
use log::{trace};

fn serialize_term_id_set<S>(ids: &HashSet<TermId>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let mut strs: Vec<String> = ids.iter().map(|t| t.to_string()).collect();
    strs.sort(); // HashSet iteration order is random; sort for stable output
    strs.serialize(serializer)
}

/// This is the structure we use for the phenoblend analysis. We have one such entity for each gene, and the entities can have one or more diseases associated with the gene
/// disease_hpo_ids then has all the observed HPO identifiers for all of the diseases
/// individual_explicit_ids has all of the explicit annotations from the phenopacket/case we are analyzing.
#[derive(Clone, Debug, Serialize)]
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


    /// We consider a term from the phenopacket (tid) and ask if there is an 
    /// exact match in a GeneDisease entity. If the GDE is annotated to the same
    /// term, obviously there is an exact match. By the true path rule, if the 
    /// GDE is annotated to a descendent of the term, then there is also an exact 
    /// match. For instance, if the proband is annotated to "Cataract" and the GDE
    /// is annotated to "Nuclear cataract", then the disease is implicitly also annotated
    /// to "Cataract" and we have an exact match.
     pub fn is_perfect_match(
        ppkt_tid: &TermId,
        gde: &GeneDiseaseEntity,
        hpo: Arc<FullCsrOntology>) -> bool {
            // 1. exact match
            if gde.disease_hpo_ids.contains(ppkt_tid) {
                 trace!(
                        "is_perfect_match term={} disease={} perfect_match={}",
                        ppkt_tid,
                        gde.gene_symbol,
                        ppkt_tid
                    );
                return true;
            }
            // 2. match to a descendent
            for disease_hpo_id in gde.disease_hpo_ids.iter() {
                if hpo.is_descendant_of(disease_hpo_id, ppkt_tid) {
                    trace!(
                        "is_perfect_match term={} disease={} perfect_match={}",
                        ppkt_tid,
                        gde.gene_symbol,
                        disease_hpo_id
                    );
                    return true;
                }
            }
            return false;
    }


    pub fn get_overlap_matrix_payload(
        proband: Proband, 
        gde_list: &Vec<GeneDiseaseEntity>, 
        disease_counts: &HashMap<TermId, usize>,
        hpo: Arc<FullCsrOntology>
    ) -> OverlapPlotPayload {
        let mut term_to_item_d: HashMap<TermId, OverlapPlotItem> = HashMap::new();
        let mut gene_entities: Vec<String> = gde_list
            .iter()
            .map(|gde| gde.gene_symbol.clone())
            .collect();
        gene_entities.sort();
        trace!("get_overlap_matrix_payload");
        trace!("gene_entities {:?}", gene_entities);

        for tid in proband.observed_hpos.iter() {
            for gde in gde_list {
                let pm_item = term_to_item_d
                    .entry(tid.clone())
                    .or_insert_with(|| OverlapPlotItem::new(tid, Arc::clone(&hpo)));
                if Self::is_perfect_match(tid, gde, Arc::clone(&hpo)) {
                    pm_item.add_perfect_match(&gde.gene_symbol);
                } else {
                    let match_score = Self::partial_match_score(tid, gde, Arc::clone(&hpo), disease_counts);
                    pm_item.add_match(&gde.gene_symbol, match_score);
                }
            }
        }
        let rows: Vec<OverlapPlotItem> = term_to_item_d.into_values().collect();
        OverlapPlotPayload { entities: gene_entities, columns: rows }
    }

    pub fn get_upset_plot_payload(
        proband: Proband, 
        gda_map: &HashMap<String, Vec<GeneDiseaseAssociation>>,
        hpo: Arc<FullCsrOntology>
    ) -> Result<UpsetPlotPayload, String> {

        let mut gd_entry_list: Vec<GeneDiseaseEntity> = Vec::new();
        for (symbol, gda_list) in gda_map.into_iter() {
            let gd_entity =   GeneDiseaseEntity::new(gda_list)?;
            gd_entry_list.push(gd_entity);
        }
        Ok(crate::blend::upset_plot::build_upset_payload(&gd_entry_list, &proband, hpo.clone()))
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


    /// Scores a *partial* match between a query phenotype and the HPO annotations
    /// of a gene's disease.
    ///
    /// A partial match means the gene's disease is not annotated with the query term
    /// itself but with one of its ancestors (a broader term). The score measures how
    /// much more specific the query term is than that ancestor, as a ratio of disease
    /// counts:
    ///
    /// ```text
    /// score = disease_counts[query_hpo] / max(disease_counts[a])
    /// ```
    ///
    /// where `a` ranges over the terms in `gde.disease_hpo_ids` that are ancestors of
    /// `query_hpo`. Taking the maximum means the broadest (most common) ancestor
    /// determines the score, which gives the most conservative value when several
    /// ancestors match.
    ///
    /// Interpretation: a score near 1 means the ancestor is about as specific as the
    /// query term, so the match is nearly as good as an exact one. A score near 0 means
    /// the ancestor is far more general. For example, if the query term is annotated to
    /// 5 diseases and the broadest matching ancestor is annotated to 50, the score is
    /// `0.1`.
    ///
    /// # Arguments
    ///
    /// * `query_hpo` - The phenotype term being scored.
    /// * `gde` - The gene-disease entity whose `disease_hpo_ids` are searched for
    ///   ancestors of `query_hpo`.
    /// * `hpo` - The ontology, used for ancestor queries.
    /// * `disease_counts` - Number of diseases annotated with each HPO term. Terms
    ///   missing from the map count as `0`.
    ///
    /// # Returns
    ///
    /// A score in `(0.0, 1.0]` if the counts are propagated (see Notes), or `0.0` if
    /// the query term is not annotated to any disease or if none of the gene's
    /// disease terms is an ancestor of the query term.
    ///
    /// # Notes
    ///
    /// * Only ancestors count. A disease term identical to `query_hpo` is not
    ///   considered here, so exact matches must be handled by the caller.
    /// * The `(0.0, 1.0]` bound assumes `disease_counts` is propagated up the
    ///   hierarchy, so an ancestor's count is never smaller than that of its
    ///   descendants. With direct-annotation counts the ratio can exceed 1.
    pub fn partial_match_score(
        query_hpo:&TermId, 
        gde: &GeneDiseaseEntity,
        hpo: Arc<FullCsrOntology>,
        disease_counts: &HashMap<TermId, usize>) -> f64 {
            let pheno_hpo_count = *disease_counts.get(query_hpo).unwrap_or(&0);
            //trace!("partial_match_score-pheno_hpo_count={}", pheno_hpo_count);
            let mut not_anc = 0;
            if pheno_hpo_count == 0 {
                return 0.0 as f64;
            }
            // Find the maximum disease count among the more specific descendant terms
            let mut max_geno_hpo_count = 0;
            for disease_hpo_id in &gde.disease_hpo_ids {
                if hpo.is_descendant_of(query_hpo, disease_hpo_id ) {
                    let anc_count = *disease_counts.get(disease_hpo_id).unwrap_or(&0);
                    if anc_count > max_geno_hpo_count {
                        max_geno_hpo_count = anc_count;
                    }
                    trace!(
                    "{} descendant of ? {} (Y/N:{}) count={} pheno_hpo_count={}",
                    disease_hpo_id,
                    query_hpo,
                    hpo.is_descendant_of(query_hpo, disease_hpo_id),
                    disease_counts.get(disease_hpo_id).unwrap_or(&0),
                    pheno_hpo_count
                );
                } else {
                    not_anc += 1;
                }
                
            }
           
            let score = if max_geno_hpo_count == 0 {
                0.0
            } else {
                pheno_hpo_count  as f64 /  max_geno_hpo_count as f64
            };
            trace!("score={} - Not ancestor={} - max_geno_hpo_count={} - ",score, not_anc, max_geno_hpo_count);
            score
            
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
        let disease_hpo: TermId = "HP:0000118".parse().unwrap(); // Phenotypic abnormality, meant just to be a more general ancestor term
        let query_hpo : TermId = "HP:0001234".parse().unwrap(); // Hitchhiker thumb

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


        let score = GeneDiseaseEntity::partial_match_score(
            &query_hpo,
            &gde,
            hpo,
            &disease_counts,
        );

        assert!((score - expected_score).abs() < f64::EPSILON);
        assert!(score.is_finite(), "Score must never be INF or NaN!");
    }
}