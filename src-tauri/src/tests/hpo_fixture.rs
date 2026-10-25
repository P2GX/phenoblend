use flate2::bufread::GzDecoder;
use ontolius::{TermId, io::OntologyLoaderBuilder, ontology::csr::FullCsrOntology};
use phenopackets::schema::v2::Phenopacket;
use rstest::fixture;
use std::{
    collections::HashMap, fs::File, io::BufReader, sync::{Arc, LazyLock},
};

use crate::{blend::disease_gene_entity::GeneDiseaseEntity, model::simple_term::SimpleOntologyTerm};

/// A thread-safe, lazily-loaded singleton instance of the FullCsrOntology (HPO).
/// 
/// This ensures the heavy HPO graph is parsed from disk exactly once across the 
/// entire test suite, saving substantial I/O and parsing time during test runs.
/// It uses `CARGO_MANIFEST_DIR` to locate the resource file reliably regardless 
/// of which subdirectory `cargo test` is executed from.
pub static HPO: LazyLock<Arc<FullCsrOntology>> = LazyLock::new(|| {
    let path = format!("{}/resources/hp.v2026-09-01.json.gz", env!("CARGO_MANIFEST_DIR"));
    let file = File::open(&path).unwrap_or_else(|_| {
        panic!("Failed to open HPO ontology resource file at expected path: {}", path)
    });
    
    let reader = GzDecoder::new(BufReader::new(file));
    let loader = OntologyLoaderBuilder::new().obographs_parser().build();
    let hpo = loader.load_from_read(reader).unwrap_or_else(|e| {
        panic!("Failed to parse HPO obographs ontology: {:?}", e)
    });
    
    Arc::new(hpo)
});

/// An `rstest` fixture that provides a cheap clone of the shared HPO singleton `Arc`.
/// 
/// Injecting this into your test functions lets `rstest` automatically supply 
/// the preloaded ontology without reloading it for every test.
#[fixture]
pub fn hpo() -> Arc<FullCsrOntology> {
    Arc::clone(&HPO)
}


const A13_JSON: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"), "/resources/data/PMID_27082122_A13.json"
));

pub fn a13_ppkt() -> Phenopacket {
    serde_json::from_str(A13_JSON).expect("fixture must parse as a Phenopacket")
}

#[rstest::fixture]
pub fn a13() -> Phenopacket { 
    a13_ppkt() 
}


// Real HPO terms used by the tests. `fixture_relationships_hold` (in the upset_plot tests)
// checks the ontology facts they rely on.
pub const ROOT: &str = "HP:0000118"; // Phenotypic abnormality
pub const SEIZURE: &str = "HP:0001250";
pub const FOCAL: &str = "HP:0007359"; // Focal-onset seizure, a seizure type
pub const TONIC: &str = "HP:0010818"; // Generalized tonic seizure, another seizure type
pub const SHORT_STATURE: &str = "HP:0004322"; // unrelated to seizures
pub const SCOLIOSIS: &str = "HP:0002650"; // unrelated to seizures and to short stature

pub fn tid(s: &str) -> TermId {
    s.parse().expect("valid CURIE")
}

pub fn tids(ts: &[&str]) -> Vec<TermId> {
    ts.iter().map(|t| tid(t)).collect()
}

/// One gene with the given directly annotated terms. The gene id and the disease are
/// placeholders: the tests only use the gene symbol and the HPO term ids.
pub fn entity(gene: &str, ts: &[&str]) -> GeneDiseaseEntity {
    GeneDiseaseEntity {
        ncbi_gene_id: format!("test-{gene}"),
        gene_symbol: gene.to_string(),
        disease_list: vec![
            SimpleOntologyTerm::new("OMIM:000000", format!("Test disease of {gene}")).unwrap(),
        ],
        disease_hpo_ids: tids(ts).into_iter().collect(),
        term_frequencies: HashMap::new()
    }
}