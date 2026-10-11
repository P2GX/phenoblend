use std::{collections::HashSet, sync::Arc};
use itertools::Itertools;
use ontolius::{TermId, ontology::{TaxonomyWalk, csr::FullCsrOntology}};

use crate::{blend::{disease_gene_entity::GeneDiseaseEntity, dto::UpsetPlotPayload}, model::proband::Proband};

/// Expands a set of TermIds to include all of their ancestors (inclusive of the originals).
fn propagate_ancestors(terms: &HashSet<TermId>, hpo: &FullCsrOntology) -> HashSet<TermId> {
    let mut expanded = terms.clone();
    for term_id in terms {
        expanded.extend(hpo.iter_term_and_ancestor_ids(term_id).cloned());
    }
    expanded
}

/// Precomputes an ancestor-propagated term closure for each entity, once.
fn build_expanded_entity_terms(
    entities: &[GeneDiseaseEntity],
    hpo: &FullCsrOntology,
) -> Vec<HashSet<TermId>> {
    entities
        .iter()
        .map(|e| propagate_ancestors(&e.disease_hpo_ids, hpo))
        .collect()
}

/// Marginal statistics for a single entity: raw (non-propagated) annotated
/// and observed counts.
fn gene_marginal(
    entity: &GeneDiseaseEntity,
    observed_phenotypes: &HashSet<TermId>,
) -> (u32, u32) {
    let annotated = entity.disease_hpo_ids.len() as u32;
    let observed = entity.disease_hpo_ids
        .intersection(observed_phenotypes)
        .count() as u32;
    (annotated, observed)
}

/// Per-gene marginal statistics across all entities: raw (non-propagated)
/// annotated and observed counts.
fn build_gene_marginals(
    entities: &[GeneDiseaseEntity],
    observed_phenotypes: &HashSet<TermId>,
) -> (Vec<u32>, Vec<u32>) {
    entities
        .iter()
        .map(|entity| gene_marginal(entity, observed_phenotypes))
        .unzip()
}

/// Computes the exclusive intersection of expanded term sets for one specific
/// combination of entity indices — terms shared by every included entity,
/// minus any term also covered by an excluded entity.
fn exclusive_intersection(
    combo_indices: &[usize],
    expanded_terms: &[HashSet<TermId>],
) -> HashSet<TermId> {
    let n = expanded_terms.len();
    let included: &[usize] = combo_indices;

    let mut shared = expanded_terms[included[0]].clone();
    for &idx in &included[1..] {
        shared = shared.intersection(&expanded_terms[idx]).cloned().collect();
    }

    for idx in 0..n {
        if !included.contains(&idx) {
            shared = shared.difference(&expanded_terms[idx]).cloned().collect();
        }
    }

    shared
}

/// Builds every exclusive combination (of every size, 1..=n) across all entities,
/// dropping combinations whose exclusive intersection is empty.
fn build_combinations(
    entities: &[GeneDiseaseEntity],
    expanded_terms: &[HashSet<TermId>],
    combo_observed_set: &HashSet<TermId>,
) -> (Vec<Vec<String>>, Vec<u32>, Vec<u32>) {
    let n_genes = entities.len();
    let mut combinations = Vec::new();
    let mut combination_annotated = Vec::new();
    let mut combination_observed = Vec::new();

    for size in 1..=n_genes {
        for combo_indices in (0..n_genes).combinations(size) {
            let shared_terms = exclusive_intersection(&combo_indices, expanded_terms);
            if shared_terms.is_empty() {
                continue;
            }

            let combo_symbols: Vec<String> = combo_indices
                .iter()
                .map(|&idx| entities[idx].gene_symbol.clone())
                .collect();

            combinations.push(combo_symbols);
            combination_annotated.push(shared_terms.len() as u32);
            combination_observed.push(shared_terms.intersection(combo_observed_set).count() as u32);
        }
    }

    (combinations, combination_annotated, combination_observed)
}

pub fn build_upset_payload(
    entities: &[GeneDiseaseEntity],
    proband: &Proband,
    hpo: Arc<FullCsrOntology>,
) -> UpsetPlotPayload {
    let gene_symbols: Vec<String> = entities.iter().map(|e| e.gene_symbol.clone()).collect();
    let observed_phenotypes: HashSet<TermId> = proband.observed_hpos.iter().cloned().collect();
    let combo_observed_set = propagate_ancestors(&observed_phenotypes, &hpo);

    let expanded_terms = build_expanded_entity_terms(entities, &hpo);
    let (gene_annotated, gene_observed) = build_gene_marginals(entities, &observed_phenotypes);
    let (combinations, combination_annotated, combination_observed) =
        build_combinations(entities, &expanded_terms, &combo_observed_set);

    UpsetPlotPayload {
        genes: gene_symbols,
        combinations,
        combination_annotated,
        combination_observed,
        gene_annotated,
        gene_observed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    fn term(id: &str) -> TermId {
        let tid: TermId = id.parse().unwrap();
        tid
    }

    fn set(ids: &[&str]) -> HashSet<TermId> {
        ids.iter().map(|s| term(s)).collect()
    }

    fn entity(gene_symbol: &str, hpo_ids: &[&str]) -> GeneDiseaseEntity {
        GeneDiseaseEntity {
            gene_symbol: gene_symbol.to_string(),
            disease_hpo_ids: set(hpo_ids),
            ncbi_gene_id: String::default(),
            disease_list: vec![],
        }
    }

    #[rstest]
    fn no_overlap_with_observed_yields_zero_observed_count() {
        let e = entity("GENE1", &["HP:0001", "HP:0002"]);
        let observed = set(&["HP:0099"]);
        let (annotated, observed_count) = gene_marginal(&e, &observed);
        assert_eq!(annotated, 2);
        assert_eq!(observed_count, 0);
    }

    #[rstest]
    fn partial_overlap_with_observed_counts_only_shared_terms() {
        let e = entity("GENE1", &["HP:0001", "HP:0002", "HP:0003"]);
        let observed = set(&["HP:0001", "HP:0003", "HP:0099"]);
        let (annotated, observed_count) = gene_marginal(&e, &observed);
        assert_eq!(annotated, 3);
        assert_eq!(observed_count, 2);
    }

    #[rstest]
    fn empty_disease_hpo_ids_yields_zero_for_both() {
        let e = entity("GENE1", &[]);
        let observed = set(&["HP:0001"]);
        let (annotated, observed_count) = gene_marginal(&e, &observed);
        assert_eq!(annotated, 0);
        assert_eq!(observed_count, 0);
    }
}