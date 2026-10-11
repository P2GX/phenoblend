
use std::collections::HashSet;
use std::{collections::HashMap, sync::Arc};

use fenominal::{AutoCompleter, Fenominal, FenominalSentence, OntologyMatch};
use ga4ghphetools::{dto::hpo_term_dto::HpoTermDuplet, load_ontology};
use ga4ghphetools::HierarchyMapItem;
use ontolius::{TermId,  ontology::csr::FullCsrOntology};
use phenopackets::schema::v2::Phenopacket;
use rand::SeedableRng;
use rand::rngs::StdRng;


use crate::blend::dto::{SpreadPlotPayload, UpsetPlotPayload};
use crate::enrichment::dto::DuoEnrichmentSummary;
use crate::enrichment::gdannotation_source::{GeneDiseaseAnnotationSource, OntoliusDescendants};
use crate::enrichment::melded_ppkt::MeldedPpkt;
use crate::enrichment::test::{EnrichmentTest, Method};
use crate::excel::payload_to_excel::export_upset_plot_to_xlsx;
use crate::model::notebook::Notebook;
use crate::{blend::dto::OverlapPlotPayload, hpoa::disease_model::SimpleDiseaseModel, model::{proband::Proband, simple_term::SimpleOntologyTerm}};
use crate::hpoa::disease_model::GeneDiseaseAssociation;
use crate::util::settings::PhenoblendSettings;

pub struct PhenoblendSingleton {
    pub(crate) settings: PhenoblendSettings,
    pub(crate) individual: Proband,
    pub(crate) hpo: Option<Arc<FullCsrOntology>>,
    pub(crate) omim_disease_models:  Option<HashMap<TermId, SimpleDiseaseModel>>,
    pub(crate) gene_to_disease_d: Option<HashMap<String, Vec<GeneDiseaseAssociation>>>,
    pub(crate) disease_count_d: HashMap<TermId, usize>,
    pub(crate) autocompleter: Option<AutoCompleter>,
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
        let proband_id = match ppkt.subject {
            Some(sjt) => sjt.id.clone(),
            None => { return Err(format!("Could not extract subject from phenopacket: {:?}", ppkt));},
        };
        let mut observed_hpo: Vec<TermId> = Vec::new();
        let mut disease_list: Vec<SimpleOntologyTerm> = Vec::new();
        let mut gene_symbol_list: Vec<String> = Vec::new();
        for pf in ppkt.phenotypic_features.iter() {
            if pf.excluded {
                continue;
            }
            let hpo_id = match &pf.r#type {
                Some(oclass) => oclass.id.clone(),
                None => { return Err(format!("Could not extract id from Ontology Class: {:?}", &pf.r#type));},
            };
            let hpo_tid: TermId = hpo_id.parse::<TermId>().map_err(|_| format!("Failed to parse TermId from string '{}'", hpo_id))?;
            observed_hpo.push(hpo_tid);
        }
        for disease in ppkt.diseases {
            let dterm = disease.term.as_ref()
                .ok_or_else(|| format!("Could not extract term from disease: {:?}", disease))?;
            let sterm = SimpleOntologyTerm::new(dterm.id.as_str(), dterm.label.as_str())?;
            disease_list.push(sterm);
        }
        for interpretation in ppkt.interpretations {
            if let Some(symbol) = ga4ghphetools::get_gene_symbol_from_interpretation(&interpretation) {
                gene_symbol_list.push(symbol);
            }
        }
        let proband = Proband { 
            id: proband_id, 
            gene_list: gene_symbol_list, 
            disease_list: disease_list, 
            observed_hpos: observed_hpo.clone()
        };
        self.individual = proband;
        let hpo = self.hpo
            .clone()
            .ok_or_else(|| "HPO ontology not loaded".to_string())?;
       
        let simple_observed = crate::hpo::hpo_util::ontolius_to_simple_list(hpo.clone(), &observed_hpo)?;
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
        observed: Vec<String>
    ) -> Result<(), String>{
        let mut observed_tid: Vec<TermId> = Vec::new();
        for obs in observed {
            let tid: TermId = obs.parse().map_err(|e| format!("Could not parse {obs}"))?;
            observed_tid.push(tid);
        }
        let proband = Proband { 
            id: String::default(), 
            gene_list: vec![], 
            disease_list: vec![], 
            observed_hpos: observed_tid.clone()
         };
        self.individual = proband;
        let hpo = self.hpo
            .clone()
            .ok_or_else(|| "HPO ontology not loaded".to_string())?;
       
        let simple_observed = crate::hpo::hpo_util::ontolius_to_simple_list(hpo.clone(), &observed_tid)?;
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
        let pm = crate::blend::overlap_matrix::calculate_overlap_matrix(
            hpo.clone(), 
            &annotations, 
            &self.disease_count_d, 
            proband)?;
        Ok(crate::blend::overlap_matrix::sort_presence_payload(pm))
    }

    pub fn get_upset_plot_payload(
        &self,
        annotations: HashMap<String, Vec<GeneDiseaseAssociation>>
    )  -> Result<UpsetPlotPayload, String> {
         let hpo = self.hpo.as_ref()
            .ok_or_else(|| "Missing required resource: HPO Ontology".to_string())?;
        
        let proband = self.individual.clone();
        let upset = crate::blend::disease_gene_entity::GeneDiseaseEntity::get_upset_plot_payload(
            proband, 
            &annotations, 
            hpo.clone())?;

        Ok(upset)
    }


    pub fn get_spread_plot_payload(
        &self,
        annotations: HashMap<String, Vec<GeneDiseaseAssociation>>
    ) -> Result<SpreadPlotPayload, String> {
        let hpo = self.hpo.as_ref()
            .ok_or_else(|| "Missing required resource: HPO Ontology".to_string())?;
        
        let proband = self.individual.clone();
        let spread = crate::blend::disease_gene_entity::GeneDiseaseEntity::get_spread_plot_payload(
            proband, 
            &annotations, 
            hpo.clone())?;
            Ok(spread)
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
            self.disease_count_d = disease_counts;
            Ok(())
        }
        


    pub fn disease_count(&self, term_id: &TermId) -> usize {
        self.disease_count_d
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
                    m.disease_model = Some(SimpleDiseaseModel::from_id(m.disease_id.clone()));
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

    pub fn output_excel(
        &mut self, 
        data_type: String, 
        annotations: HashMap<String, Vec<GeneDiseaseAssociation>>, 
        path_str: String) 
    -> Result<(), String> {
        if data_type == "overlap" {
            let payload = self.calculate_overlap_matrix(annotations)?;
            crate::excel::payload_to_excel::export_overlap_plot_to_xlsx(&payload, &path_str)
                .map_err(|e| format!("Could not output overlap Excel file: {e}"))?;
            Ok(())
        } else if data_type == "upset" {
            let payload = self.get_upset_plot_payload(annotations)?;
            export_upset_plot_to_xlsx(&payload, &path_str)
                .map_err(|e| format!("Could not output upset Excel file: {e}"))?;
            Ok(())
        } else {
            Err(format!("{data_type} excel export not supported yet"))
        }
    }



     pub fn get_duo_enrichment_payload(
        &self,
        annotations: HashMap<String, Vec<GeneDiseaseAssociation>>,
        n_sim: usize,
    ) -> Result<Vec<DuoEnrichmentSummary>, String> {
        let hpo = self.hpo.as_ref().ok_or("Ontology not loaded")?;

        let genes: Vec<String> = annotations.keys().cloned().collect();
        let pairs = gene_pairs(&genes);
        if pairs.is_empty() {
            return Err(format!(
                "Need at least 2 genes to test for shared enrichment, got {}: {:?}",
                genes.len(),
                genes
            ));
        }

        let observed: HashSet<TermId> = self.individual.observed_hpos.iter().cloned().collect();
        let descendants = OntoliusDescendants::new(hpo.clone());
        // annotations passed in from the frontend now supplies both the gene
        // list AND the per-gene disease models directly, so the adapter
        // reads from `annotations` rather than `self.gene_to_disease_d` /
        // `self.omim_disease_models` -- see the AnnotationSource impl below.
        let source = GeneDiseaseAnnotationSource::new(&annotations, &descendants);

        let mut results = Vec::with_capacity(pairs.len());
        for (gene_a, gene_b) in pairs {
            let melded = MeldedPpkt::new(gene_a.clone(), gene_b.clone(), &observed, &source);
            let mut rng = StdRng::seed_from_u64(seed_for_pair(&self.individual.id, &gene_a, &gene_b));

            let unweighted = EnrichmentTest::new(Method::Unweighted, n_sim).run(&melded, &source, &mut rng);
            let frequency = EnrichmentTest::new(Method::Frequency, n_sim).run(&melded, &source, &mut rng);

            results.push(DuoEnrichmentSummary {
                gene_a,
                gene_b,
                n_union: melded.n_union(),
                n_shared: melded.n_shared(),
                n_obs_union: melded.n_obs_union(),
                n_obs_shared: melded.n_obs_shared(),
                has_overlap: melded.has_overlap(),
                unweighted,
                frequency,
            });
        }
        Ok(results)
    }


}


    /// All unordered 2-combinations of `genes`, deduplicated. A proband with
/// exactly 2 genes yields exactly 1 pair (the common case); 3+ genes yields
/// every candidate pair, so the enrichment test doesn't need to know in
/// advance how many genes were solved for a given proband.
pub fn gene_pairs(genes: &[String]) -> Vec<(String, String)> {
    let mut pairs = Vec::new();
    for i in 0..genes.len() {
        for j in (i + 1)..genes.len() {
            pairs.push((genes[i].clone(), genes[j].clone()));
        }
    }
    pairs
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
            disease_count_d: HashMap::new(),
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

/// Deterministic per-(proband, gene pair) seed, so re-running the analysis
/// on the same proband and pair reproduces the same Monte Carlo draw, but
/// different pairs within one proband don't share an RNG stream.
fn seed_for_pair(proband_id: &str, gene_a: &str, gene_b: &str) -> u64 {
    use std::hash::{Hash, Hasher};
    use std::collections::hash_map::DefaultHasher;
    let mut hasher = DefaultHasher::new();
    proband_id.hash(&mut hasher);
    gene_a.hash(&mut hasher);
    gene_b.hash(&mut hasher);
    hasher.finish()
}