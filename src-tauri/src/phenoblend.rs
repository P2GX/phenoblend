
use std::collections::{BTreeMap, HashSet};
use std::{collections::HashMap, sync::Arc};

use fenominal::{AutoCompleter, Fenominal, FenominalSentence, OntologyMatch};
use ga4ghphetools::{dto::hpo_term_dto::HpoTermDuplet, load_ontology};
use ga4ghphetools::HierarchyMapItem;
use ontolius::ontology::{MetadataAware, OntologyTerms};
use ontolius::{TermId,  ontology::csr::FullCsrOntology};
use phenopackets::schema::v2::Phenopacket;


use crate::blend::disease_gene_entity::GeneDiseaseEntity;
use crate::blend::dto::{ObservedHpoInput, OverlapData, PhenoBlendAnalysis, Severity, SpreadPlotPayload, UpsetPlotData, severity_from_ids};
use crate::blend::spread_plot::build_spread_payload;
use crate::model::Notebook;
use crate::model::status::InitializationStatusDto;
use crate::{blend::dto::OverlapPlotPayload, hpoa::disease_model::SimpleDiseaseModel, model::{proband::Proband, simple_term::SimpleOntologyTerm}};
use crate::hpoa::disease_model::GeneDiseaseAssociation;
use crate::util::settings::PhenoblendSettings;

pub struct PhenoblendSingleton {
    pub(crate) settings: PhenoblendSettings,
    pub(crate) individual: Proband,
    pub(crate) hpo: Option<Arc<FullCsrOntology>>,
    pub(crate) omim_disease_models:  Option<HashMap<TermId, SimpleDiseaseModel>>,
    pub(crate) gene_to_disease_d: Option<HashMap<String, Vec<GeneDiseaseAssociation>>>,
    /// Counts of explicit annotations of each HPO term to diseases (ranges from 0 to ca. 9000 and used to calculate similarity)
    pub(crate) hpo_term_to_disease_count_d: HashMap<TermId, usize>,
    /// Autocompletion functionality from Fenominal
    pub(crate) autocompleter: Option<AutoCompleter>,
    /// Object that records input data and analysis parameters and results for Q/C or downstream use.
    pub(crate) notebook: Option<Notebook>
}


impl PhenoblendSingleton {
    pub fn new() -> Self {
        PhenoblendSingleton::default()
    }


    pub fn set_hpo(&mut self, ontology: Arc<FullCsrOntology>, hpo_json_path: &str) {
        self.hpo = Some(ontology.clone());
        self.autocompleter = Some(AutoCompleter::new(ontology.clone()));
        let _ = self.settings.set_hp_json_path(hpo_json_path);
    }

    pub fn get_hpo(&self) -> Option<Arc<FullCsrOntology>> {
        match &self.hpo {
            Some(hpo) => Some(hpo.clone()),
            None => None,
        }
    }



    pub fn set_hpoa_d(&mut self, hpoa_d: HashMap<TermId, SimpleDiseaseModel>, hpoa_path: &str) {
        self.omim_disease_models = Some(hpoa_d);
        let _ = self.calculate_disease_counts();
        let _ = self.settings.set_hpoa_path(hpoa_path);
    }
    

    pub fn set_gene_to_disease(&mut self, g2d: HashMap<String, Vec<GeneDiseaseAssociation>>, g2d_path: &str) {
        self.gene_to_disease_d = Some(g2d);
        let _ = self.settings.set_g2d_path(g2d_path);
    }

    
    /// Input a GA4GH Phenopacket.
    pub fn ingest_ppkt(&mut self, ppkt: Phenopacket) -> Result<(), String> {
        let proband = Proband::from_phenopacket(&ppkt)?;
        let hpo = self
            .hpo
            .clone()
            .ok_or_else(|| "HPO ontology not loaded".to_string())?;
        let simple_observed =
            crate::hpo::hpo_util::ontolius_to_simple_list(hpo.as_ref(), &proband.observed_hpos)?;
        self.individual = proband;
        self.init_notebook_from_ppkt(&ppkt.id, &simple_observed);
        Ok(())
    }


    fn init_notebook_from_ppkt(
        &mut self, 
        id: impl Into<String>, 
        observed_hpo: &Vec<SimpleOntologyTerm>) {
            let notebook = Notebook::from_ppkt(observed_hpo.clone(), id.into());
            self.notebook = Some(notebook);
    }

    fn init_notebook_from_ner(
        &mut self,
         observed_hpo: &Vec<SimpleOntologyTerm>) {
            let notebook = Notebook::from_ner(observed_hpo.clone());
            self.notebook = Some(notebook);
    }


    pub fn add_observed_hpos_from_ner(
        &mut self,
        observed: Vec<ObservedHpoInput>
    ) -> Result<(), String>{
        let mut observed_tid: Vec<TermId> = Vec::new();
         let mut severities: BTreeMap<String, Severity> = BTreeMap::new();
        for obs in &observed {
            let term_id_str = &obs.term_id;
            let tid: TermId = term_id_str.parse().map_err(|_| format!("Could not parse {term_id_str}"))?;
            observed_tid.push(tid);
            if let Some(sev) = severity_from_ids(&obs.modifier_ids){
                severities.insert(term_id_str.clone(), sev);
            }
           
        }
        let proband = Proband { 
            id: String::default(), 
            gene_list: vec![], 
            disease_list: vec![], 
            observed_hpos: observed_tid.clone(),
            severities: severities.clone(),
         };
        self.individual = proband;
        let hpo = self.hpo
            .clone()
            .ok_or_else(|| "HPO ontology not loaded".to_string())?;
       
        let simple_observed = crate::hpo::hpo_util::ontolius_to_simple_list(hpo.as_ref(), &observed_tid)?;
        self.init_notebook_from_ner(&simple_observed);
        Ok(())
    }

    pub fn calculate_overlap_matrix(
        &mut self, 
        annotations: HashMap<String, Vec<GeneDiseaseAssociation>>
    ) -> Result<OverlapPlotPayload, String> {
        let hpo = self.hpo.as_ref()
            .ok_or_else(|| "Missing required resource: HPO Ontology".to_string())?;
        
        let proband = self.individual.clone();
        let pm = crate::blend::overlap_payload::calculate_overlap_matrix(
            hpo.clone(), 
            &annotations, 
            &self.hpo_term_to_disease_count_d, 
            proband)?;
        Ok(crate::blend::overlap_payload::sort_presence_payload(pm))
    }

    pub fn get_upset_plot_data(
        &self,
    )  -> Result<UpsetPlotData, String> {
        let nb = self.notebook.as_ref().ok_or("Notebook not initialized")?;
        let ana = nb.analysis.as_ref().ok_or("Analysis not initialized")?;
        ana.upset_payload()
    }

    pub fn get_overlap_plot_payload(
        &self
    ) -> Result<OverlapPlotPayload, String> {
        let nb = self.notebook.as_ref().ok_or("Notebook not initialized")?;
        let ana = nb.analysis.as_ref().ok_or("Analysis not initialized")?;
        ana.overlap_payload()
    }


    pub fn get_spread_plot_payload(
        &self,
    ) -> Result<SpreadPlotPayload, String> {
        let nb = self.notebook.as_ref().ok_or("Notebook not initialized")?;
        let ana = nb.analysis.as_ref().ok_or("Analysis not initialized")?;
         Ok(ana.spread.clone())
  }



        pub fn calculate_disease_counts(&mut self) -> Result<(), String>{
            let mut disease_counts: HashMap<TermId, usize> = HashMap::new();
            if let Some(models) = &self.omim_disease_models {
                for model in models.values() {
                    for hpo_id in &model.observed_hpo_ids {
                        *disease_counts.entry(hpo_id.clone()).or_insert(0) += 1;
                    }
                }
            }
            self.hpo_term_to_disease_count_d = disease_counts;
            Ok(())
        }
        


    pub fn disease_count(&self, term_id: &TermId) -> usize {
        self.hpo_term_to_disease_count_d
            .get(term_id)
            .copied()
            .unwrap_or(0)
    }

     /// Provide Strings with TermId - Label that will be used for autocompletion
    /// fenominal functionality
    pub fn search_hpo(&self, query: &str, limit: usize) -> Vec<OntologyMatch> {
        self.autocompleter
            .as_ref()
            .map(|ac| ac.search_hpo(query, limit))
            .unwrap_or_default()
    }


    pub fn get_hpo_parent_and_children_terms(&self, term_id: &str) -> Result<HierarchyMapItem, String> {
        match &self.hpo {
            Some(hpo) => {
                let hm = ga4ghphetools::get_hpo_parent_and_children_terms(term_id, hpo.clone());
                Ok(hm)
            },
            None => Err("Could not retrieve parent/child hierarchy".to_string())
        }
    }

    pub fn mine_clinical_text(
        &self,
        text: &str
     ) -> Result<Vec<FenominalSentence>, String> {
        let hpo = self.hpo.as_ref().ok_or_else(|| "HPO not initialized".to_string())?;
        let fenominal = Fenominal::new_hpo(hpo.clone());
        fenominal.mine_sentences(text).map_err(|e|e.to_string())
    }

    pub fn get_modifiers(&self) -> Result<Vec<HpoTermDuplet>, String> {
        let hpo = self.hpo.as_ref().ok_or_else(|| "HPO not initialized".to_string())?;
        ga4ghphetools::get_modifiers(hpo.clone())
    }

    pub fn perform_hpo_autocomplete(&self, query: String) -> Result<Vec<OntologyMatch>, String> {
        let autocompleter = self.autocompleter.as_ref().ok_or_else(|| "Autocomplete not initialized".to_string())?;
        let n_term_limit = 20;
        Ok(autocompleter.search_hpo(&query, n_term_limit))
    }

    /// Returns all GeneDiseaseAssociation entries for genes whose symbol
    /// starts with `query` (case-insensitive). Since a single gene symbol
    /// can map to multiple diseases, and a query can match multiple gene
    /// symbols, results are flattened across both dimensions.
    pub fn autocomplete_gene_symbol(&self, query: &str, limit: usize) -> Result<Vec<GeneDiseaseAssociation>, String> {
        let map = self.gene_to_disease_d
            .as_ref()
            .ok_or_else(|| "Gene-to-disease map not initialized".to_string())?;
        let disease_map = self.omim_disease_models
            .as_ref()
            .ok_or_else(|| "disease models map n ot initialized!".to_string())?;

        let trimmed = query.trim();
        if trimmed.is_empty() {
            return Ok(vec![]);
        }

        let query_lower = trimmed.to_lowercase();

        let mut matches: Vec<GeneDiseaseAssociation> = map
            .iter()
            .filter(|(gene_symbol, _)| gene_symbol.to_lowercase().starts_with(&query_lower))
            .flat_map(|(_, associations)| associations.iter().cloned())
            .collect();

        for m in &mut matches {
            match disease_map.get(&m.disease_id) {
                Some(simple_disease) => {
                    m.disease_model = Some(simple_disease.clone());
                }
                None => {
                    m.disease_model = Some(SimpleDiseaseModel::fallback_from_id(m.disease_id.clone()));
                }
            }
        }

        matches.sort_by(|a, b| {
            a.gene_symbol.cmp(&b.gene_symbol)
                .then_with(|| a.disease_id.to_string().cmp(&b.disease_id.to_string()))
        });

        matches.truncate(limit);

        Ok(matches)
    }

   pub fn get_observed_hpo_count(&self) -> usize {
        self.individual.observed_hpos.len()
    }

   /// What the frontend shows about the loaded resources, read from the current state.
    pub fn initialization_status(&self) -> InitializationStatusDto {
        InitializationStatusDto {
            hpo_loaded: self.hpo.is_some(),
            hpo_terms: self.hpo.as_ref().map_or(0, |o| o.len()),
            hpo_version: self.hpo.as_ref().map_or("n/a".to_string(), |o| o.version().to_string()),
            hpoa_loaded: self.omim_disease_models.is_some(),
            hpoa_diseases: self.omim_disease_models.as_ref().map_or(0, |omim_d|omim_d.len()),
            g2d_loaded: self.gene_to_disease_d.is_some(),
            n_genes: self.gene_to_disease_d.as_ref().map_or(0, |gd_d|gd_d.len()),
        }
    }



   

     pub(crate) fn notebook_ref(&self) -> Result<&Notebook, String> {
        self.notebook
            .as_ref()
            .ok_or_else(|| "No notebook available yet".to_string())
    }

    pub(crate) fn notebook_mut(&mut self) -> Result<&mut Notebook, String> {
        self.notebook
            .as_mut()
            .ok_or_else(|| "No notebook available yet".to_string())
    }


    /// Runs the full analysis for the notebook's current observed terms against `entities`.
    ///
    /// Read-only: returns the result without touching `self`, so a failure leaves the
    /// notebook as it was. Requires the ontology and the disease counts to be loaded.
    fn run_analysis(&self, entities: &[GeneDiseaseEntity]) -> Result<PhenoBlendAnalysis, String> {
        let hpo = self.get_hpo().ok_or("Could not get ontology".to_string())?; // Arc<FullCsrOntology>; adapt if get_hpo still returns Option
        let observed: Vec<TermId> = self
            .notebook_ref()?
            .observed_hpo
            .iter()
            .map(|t| t.term_id.clone())
            .collect();

        let term_matches =
            crate::blend::disease_gene_entity::compute_term_matches(&self.individual, entities, &self.hpo_term_to_disease_count_d, hpo.as_ref())?;
    
        let overlap = OverlapData::from_term_matches(
            entities.iter().map(|g| g.gene_symbol.clone()),
            &term_matches,
        );
       
        let mut ids_by_gene: BTreeMap<String, BTreeMap<String, &TermId>> = BTreeMap::new();
        for g in entities {
            let ids = ids_by_gene.entry(g.gene_symbol.clone()).or_default();
            for t in &g.disease_hpo_ids {
                ids.insert(t.to_string(), t);
            }
        }
        let gene_annotations = ids_by_gene
            .into_iter()
            .map(|(gene, ids)| {
                let terms = ids
                    .into_values()
                    .map(|t| crate::blend::disease_gene_entity::simple_term(hpo.as_ref(), t))
                    .collect();
                (gene, terms)
            })
            .collect();
        let upset_input = crate::blend::upset_plot::plot_input(entities, &observed, hpo.as_ref());
        let upset = crate::blend::upset_plot::plot_data(&upset_input)?;
        let mut genotypes: Vec<String> = entities.iter().map(|g| g.gene_symbol.clone()).collect();
        genotypes.sort();
        genotypes.dedup();
        let spread = build_spread_payload(&observed, &genotypes, entities, hpo.as_ref())?;
        let enrichment = crate::blend::shared_enrichment::all_enrichments(
            &upset_input,
            crate::blend::shared_enrichment::DEFAULT_N_SIM,
            &self.individual.id,
        )?;
        Ok(PhenoBlendAnalysis { term_matches, gene_annotations, overlap, upset, spread, enrichment })
    }

    /// Builds the gene-disease entities, runs the analysis, and stores both in the notebook.
    ///
    /// The disease counts are recalculated first so they always reflect the loaded
    /// disease models. If any step fails, the notebook is left unchanged.
    ///
    /// # Errors
    /// Fails if the ontology or disease models are not loaded, if there is no notebook,
    /// or if an entity or match cannot be built.
    pub(crate) fn initialize_notebook(
        &mut self,
        gda_list: HashMap<String, Vec<GeneDiseaseAssociation>>,
    ) -> Result<Notebook, String> {
        let entities = crate::blend::overlap_payload::build_entity_list(&gda_list)?;
        let analysis = self.run_analysis(&entities)?;   // may fail: nothing changed yet
        let nb = self.notebook_mut()?;
        nb.gene_disease_list = entities;
        nb.analysis = Some(analysis);
        Ok(nb.clone())
    }

}






impl Default for PhenoblendSingleton {
    fn default() -> Self {
        let settings = PhenoblendSettings::load_settings();
        let mut singleton = Self { 
            settings: settings.clone(),
            individual: Proband::default(),
            hpo: None,
            omim_disease_models: None,
            gene_to_disease_d: None,
            hpo_term_to_disease_count_d: HashMap::new(),
            autocompleter: None,
            notebook: None
        };
       if let Some(ontology) = settings.get_hp_json_path().ok().and_then(|path| load_ontology(&path).ok()) {
            singleton.hpo = Some(ontology.clone());
            let ac = AutoCompleter::new(ontology.clone());
            singleton.autocompleter = Some(ac);
        } else {
            println!("Did not get ontology");
            println!("Oath: {:?}", settings.get_hp_json_path());
        }
        if let Some(omim_map) = settings.get_hpoa_path().ok().and_then(|path| crate::hpoa::hpoa_ingest::load_hpoa_d(&path).ok()) {
            singleton.omim_disease_models = Some(omim_map);
            let _ = singleton.calculate_disease_counts();
        }
        if let Some(g2d) = settings.get_g2d_path().ok().and_then(|path| crate::hpoa::gene_to_disease::load_gene_disease_associations(&path).ok()) {
            singleton.gene_to_disease_d = Some(g2d);
        }
        singleton
    }
}








#[cfg(test)]
mod tests {
    use rstest::rstest;

use super::*;
    use crate::tests::hpo_fixture::a13;
    use crate::tests::hpo_fixture::a13_ppkt;
    use phenopackets::schema::v2::core::{Disease, Individual, OntologyClass, PhenotypicFeature};


    #[rstest]
    fn a13_severity_reaches_proband(a13: Phenopacket) {
        let mut pblend = PhenoblendSingleton::new();
        let result = pblend.ingest_ppkt(a13);
        assert!(result.is_ok());
        let p = pblend.individual.clone();
       
        assert_eq!(p.observed_hpos.len(), 5);
        assert_eq!(p.severities.len(), 1);
        assert_eq!(p.severities.get("HP:0001639"), Some(&Severity::Severe));
        assert_eq!(p.gene_list, ["MYH7", "PRKAG2"]);
    }

  fn class(id: &str, label: &str) -> OntologyClass {
        OntologyClass { id: id.to_string(), label: label.to_string() }
    }

    fn feature(id: &str, excluded: bool) -> PhenotypicFeature {
        PhenotypicFeature {
            r#type: Some(class(id, "label")),
            excluded,
            ..Default::default()
        }
    }

    fn ppkt_with(features: Vec<PhenotypicFeature>) -> Phenopacket {
        Phenopacket {
            id: "ppkt-1".to_string(),
            subject: Some(Individual { id: "patient-1".to_string(), ..Default::default() }),
            phenotypic_features: features,
            ..Default::default()
        }
    }

    #[test]
    fn excluded_features_are_skipped() {
        let ppkt = ppkt_with(vec![feature("HP:0001250", false), feature("HP:0004322", true)]);
        let proband = Proband::from_phenopacket(&ppkt).unwrap();
        let ids: Vec<String> = proband.observed_hpos.iter().map(|t| t.to_string()).collect();
        assert_eq!(ids, ["HP:0001250"]);
        assert_eq!(proband.id, "patient-1");
    }

    #[test]
    fn a_recognized_severity_is_recorded_and_an_unknown_one_ignored() {
        let mut mild = feature("HP:0001250", false);
        mild.severity = Some(class("HP:0012825", "Mild"));
        let mut odd = feature("HP:0004322", false);
        odd.severity = Some(class("HP:0000001", "Not a severity"));
        let proband = Proband::from_phenopacket(&ppkt_with(vec![mild, odd])).unwrap();
        assert_eq!(proband.severities.get("HP:0001250"), Some(&Severity::Mild));
        assert!(!proband.severities.contains_key("HP:0004322"));
        assert_eq!(proband.observed_hpos.len(), 2); // the odd severity does not drop the feature
    }

    #[test]
    fn the_diseases_of_the_phenopacket_are_listed() {
        let mut ppkt = ppkt_with(vec![]);
        ppkt.diseases = vec![Disease {
            term: Some(class("OMIM:162200", "Neurofibromatosis, type 1")),
            ..Default::default()
        }];
        let proband = Proband::from_phenopacket(&ppkt).unwrap();
        assert_eq!(proband.disease_list.len(), 1);
        assert_eq!(proband.disease_list[0].term_id_as_string(), "OMIM:162200");
    }

    #[test]
    fn a_missing_subject_is_an_error_naming_the_phenopacket() {
        let mut ppkt = ppkt_with(vec![]);
        ppkt.subject = None;
        let err = Proband::from_phenopacket(&ppkt).unwrap_err();
        assert!(err.contains("ppkt-1"), "{err}");
    }

    #[test]
    fn a_feature_without_a_type_is_an_error() {
        let mut broken = feature("HP:0001250", false);
        broken.r#type = None;
        let err = Proband::from_phenopacket(&ppkt_with(vec![broken])).unwrap_err();
        assert!(err.contains("without a type"), "{err}");
    }

    #[test]
    fn an_invalid_term_id_is_an_error() {
        let err = Proband::from_phenopacket(&ppkt_with(vec![feature("not-an-id", false)])).unwrap_err();
        assert!(err.contains("not-an-id"), "{err}");
    }

    #[test]
    fn the_a13_fixture_gives_its_non_excluded_features_and_sorted_unique_genes() {
        let ppkt = a13_ppkt();
        let expected = ppkt.phenotypic_features.iter().filter(|f| !f.excluded).count();
        let proband = Proband::from_phenopacket(&ppkt).unwrap();
        assert_eq!(proband.observed_hpos.len(), expected);
        let mut sorted = proband.gene_list.clone();
        sorted.sort();
        sorted.dedup();
        assert_eq!(proband.gene_list, sorted);
    }
}
