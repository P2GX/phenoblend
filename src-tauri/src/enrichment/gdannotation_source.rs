use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use ontolius::TermId;
use ontolius::ontology::csr::FullCsrOntology;
use ontolius::ontology::TaxonomyWalk;
use super::annotation_source::{AnnotationSource, DescendantsProvider};
use crate::hpoa::disease_model::{GeneDiseaseAssociation};




pub struct OntoliusDescendants {
    hpo: Arc<FullCsrOntology>,
}

impl OntoliusDescendants {
    pub fn new(hpo: Arc<FullCsrOntology>) -> Self {
        Self { hpo }
    }
}

impl DescendantsProvider for OntoliusDescendants {
    fn term_and_descendants(&self, term_id: &TermId) -> HashSet<TermId> {
        let mut set: HashSet<TermId> = self.hpo.iter_descendant_ids(term_id).cloned().collect();
        set.insert(term_id.clone()); // iter_descendant_ids excludes the term itself
        set
    }
}




pub struct GeneDiseaseAnnotationSource<'a, D: DescendantsProvider> {
    gene_to_disease: &'a HashMap<String, Vec<GeneDiseaseAssociation>>,
    descendants: &'a D,
}
impl<'a, D: DescendantsProvider> GeneDiseaseAnnotationSource<'a, D> {
    pub fn new(
        gene_to_disease: &'a HashMap<String, Vec<GeneDiseaseAssociation>>,
        descendants: &'a D,
    ) -> Self {
        Self { gene_to_disease, descendants }
    }

    fn associations_for_gene(&self, gene: &str) -> &[GeneDiseaseAssociation] {
        self.gene_to_disease.get(gene).map(Vec::as_slice).unwrap_or(&[])
    }

    fn disease_ids_for_gene(&self, gene: &str) -> HashSet<TermId> {
        self.associations_for_gene(gene)
            .iter()
            .map(|assoc| assoc.disease_id.clone())
            .collect()
    }
}

impl<'a, D: DescendantsProvider> AnnotationSource for GeneDiseaseAnnotationSource<'a, D> {
    fn direct_terms_for_gene(&self, gene: &str) -> HashSet<TermId> {
        self.associations_for_gene(gene)
            .iter()
            .filter_map(|assoc| assoc.disease_model.as_ref())
            .flat_map(|model| model.observed_hpo_ids.iter().cloned())
            .collect()
    }

    fn descendant_closure(&self, terms: &HashSet<TermId>) -> HashSet<TermId> {
        terms
            .iter()
            .flat_map(|t| self.descendants.term_and_descendants(t))
            .collect()
    }

    fn diseases_for_gene(&self, gene: &str) -> HashSet<TermId> {
        self.disease_ids_for_gene(gene)
    }

    fn term_frequency_in_disease(&self, disease_id: &TermId, term_id: &TermId) -> Option<f64> {
        // Scans across all genes' associations for the matching disease_id,
        // since `gene_to_disease` is keyed by gene, not disease. Fine for a
        // single duo-enrichment call's small gene set; if this needs to run
        // over a large batch later, build a flattened disease_id -> model
        // index once up front instead of scanning per lookup.
        self.gene_to_disease
            .values()
            .flatten()
            .find(|assoc| &assoc.disease_id == disease_id)
            .and_then(|assoc| assoc.disease_model.as_ref())
            .and_then(|model| model.term_frequencies.get(term_id).copied())
    }
}