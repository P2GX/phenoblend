use flate2::bufread::GzDecoder;
use ontolius::{ontology::csr::FullCsrOntology, io::OntologyLoaderBuilder};
use rstest::fixture;
use std::{
    fs::File,
    io::BufReader,
    sync::{Arc, LazyLock},
};

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