use clap::Parser;
use ontolius::{io::OntologyLoaderBuilder, ontology::csr::FullCsrOntology};
use phenoblend_lib::write_batch;
use phenoblend_lib::load_hpoa_d;
use std::{path::PathBuf, sync::Arc};

/// Shared-term enrichment for a folder of GA4GH phenopackets.
#[derive(Parser)]
#[command(name = "phenoblend-batch", version)]
struct Args {
    /// Folder with the phenopackets (.json files)
    #[arg(long)] phenopackets: PathBuf,
    /// The HPO file (hp.json)
    #[arg(long)] hpo: String,
    /// The HPOA file
    #[arg(long)] hpoa: String,
    /// Output table (TSV)
    #[arg(long)] out: PathBuf,
    /// Simulated individuals per test
    #[arg(long, default_value_t = 10_000)] n_sim: usize,
    /// Also run the shuffled-pair control with this seed
    #[arg(long)] shuffle_seed: Option<u64>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let hpo = load_hpo(&args.hpo)?;
    let models = load_hpoa_d(&args.hpoa)?;

    let phenopackets = args.phenopackets.display().to_string();
    let n_sim = args.n_sim.to_string();
    let shuffle_seed = args.shuffle_seed.map_or("none".to_string(), |s| s.to_string());
    let provenance = [
        ("app_version", env!("CARGO_PKG_VERSION")),
        ("hpo_file", args.hpo.as_str()),
        ("hpoa_file", args.hpoa.as_str()),
        ("phenopackets", phenopackets.as_str()),
        ("n_sim", n_sim.as_str()),
        ("shuffle_seed", shuffle_seed.as_str()),
    ];

    let summary = write_batch(
        &args.phenopackets, &args.out, hpo.as_ref(), &models,
        args.n_sim, args.shuffle_seed, &provenance,
    )?;
    eprintln!(
        "{} files, {} analyzed, {} skipped, {} errors",
        summary.n_files, summary.n_files_analyzed, summary.n_skipped_rows, summary.n_error_rows
    );
    Ok(())
}

/// Load HPO JSON
fn load_hpo(json_path: &str) -> Result<Arc<FullCsrOntology>, Box<dyn std::error::Error>> {
    let loader = OntologyLoaderBuilder::new().obographs_parser().build();
    let hpo: FullCsrOntology = loader.load_from_path(json_path)?;
    Ok(Arc::new(hpo))
}