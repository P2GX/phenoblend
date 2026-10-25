use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use ontolius::{TermId, ontology::csr::FullCsrOntology};
use phenopackets::schema::v2::Phenopacket;
use serde::Serialize;

use crate::{blend::dto::shared_term_enrichment::SharedTermEnrichment, hpoa::disease_model::{GeneDiseaseAssociation, SimpleDiseaseModel}, model::proband::Proband, phenoblend::PhenoblendSingleton};
use crate::blend::batch_tsv::{RowContext, error_row, provenance_lines, result_row, shuffle_assignment, skipped_row, tsv_header};
use crate::blend::shared_enrichment::all_enrichments;
use crate::blend::upset_plot::plot_input;
use crate::blend::overlap_payload::build_entity_list;
use crate::hpo::diagnosis_from_phenopacket::{diagnoses_from_phenopacket, require_multiple_diagnoses, Diagnosis};
use crate::hpo::ppkt_loader::load_from_json;

/// The result of analyzing one phenopacket: the individual and the enrichment results for
/// their diagnosed genes.
#[derive(Debug, Clone)]
pub struct PhenopacketAnalysis {
    /// Identifier of the phenopacket.
    pub phenopacket_id: String,
    /// Identifier of the subject; it also seeds the simulations.
    pub individual_id: String,
    /// Number of distinct observed (not excluded) HPO terms of the individual.
    pub n_observed_terms: usize,
    /// One result for every pair of genes and, for three or four genes, one for all of them
    /// (see [`all_enrichments`]).
    pub results: Vec<SharedTermEnrichment>,
}

/// What is read from a phenopacket before any statistics: who the individual is, which terms
/// they have, and which (gene, disease) pairs were diagnosed. Kept apart so that the shuffled
/// control can combine the terms of one individual with the diagnoses of another.
#[derive(Debug, Clone)]
pub struct PreparedPhenopacket {
    pub phenopacket_id: String,
    pub individual_id: String,
    pub observed: Vec<TermId>,
    pub diagnoses: Vec<Diagnosis>,
}

/// Reads a phenopacket and checks that it can be analyzed (two diseases, two genes).
///
/// # Errors
/// Returns `Skipped` for a diagnosis without a gene, and `Failed` for anything that cannot be
/// read, as [`Proband::from_phenopacket`], [`diagnoses_from_phenopacket`] and
/// [`require_multiple_diagnoses`] do.
pub fn prepare_phenopacket(ppkt: &Phenopacket) -> Result<PreparedPhenopacket, PrepareProblem> {
    if let Some(reason) = non_genetic_diagnosis(ppkt) {
        return Err(PrepareProblem::Skipped(reason));
    }
    let proband = Proband::from_phenopacket(ppkt)?;
    let diagnoses = diagnoses_from_phenopacket(ppkt)?;
    require_multiple_diagnoses(&diagnoses)?;
    Ok(PreparedPhenopacket {
        phenopacket_id: ppkt.id.clone(),
        individual_id: proband.id,
        observed: proband.observed_hpos,
        diagnoses,
    })
}
/// The enrichment results for the terms of `terms_of` against the diagnoses of `diagnoses_of`.
/// The two are the same phenopacket for the real analysis, and different ones for the shuffled
/// control. The simulations are seeded from the id of `terms_of` and the genes tested.
///
/// # Errors
/// Returns an error if a disease is not in `disease_models`, or as the functions it calls.
pub fn analyze_prepared(
    terms_of: &PreparedPhenopacket,
    diagnoses_of: &PreparedPhenopacket,
    hpo: &FullCsrOntology,
    disease_models: &HashMap<TermId, SimpleDiseaseModel>,
    n_sim: usize,
) -> Result<PhenopacketAnalysis, String> {
    let mut annotation_map: HashMap<String, Vec<GeneDiseaseAssociation>> = HashMap::new();
    for diagnosis in &diagnoses_of.diagnoses {
        let model = disease_models.get(&diagnosis.disease.term_id).ok_or_else(|| {
            format!(
                "disease {} of gene {} (phenopacket '{}') was not found in the disease models",
                diagnosis.disease.term_id_as_string(),
                diagnosis.gene_symbol,
                diagnoses_of.phenopacket_id
            )
        })?;
        annotation_map
            .entry(diagnosis.gene_symbol.clone())
            .or_default()
            .push(GeneDiseaseAssociation {
                // the gene id is not needed here, but must be the same for one gene
                ncbi_gene_id: String::new(),
                gene_symbol: diagnosis.gene_symbol.clone(),
                association_type: "PHENOPACKET_INTERPRETATION".to_string(),
                disease_id: diagnosis.disease.term_id.clone(),
                disease_model: Some(model.clone()),
                source: format!("phenopacket {}", diagnoses_of.phenopacket_id),
            });
    }

    let entities = build_entity_list(&annotation_map)?;
    let input = plot_input(&entities, &terms_of.observed, hpo);
    let results = all_enrichments(&input, n_sim, &terms_of.individual_id)?;

    Ok(PhenopacketAnalysis {
        phenopacket_id: terms_of.phenopacket_id.clone(),
        individual_id: terms_of.individual_id.clone(),
        n_observed_terms: input.proband_closures.len(),
        results,
    })
}

/// Analyzes one phenopacket: the shared-term enrichment of the genes named in its
/// interpretations, using the annotations of the diseases named there.
///
/// Every interpretation must give a gene and a disease, and there must be at least two
/// diseases and two distinct genes. Each disease must be in `disease_models`. The result is
/// reproducible: the simulations are seeded from the subject's id and the genes.
///
/// # Errors
/// Returns an error describing the first problem: a phenopacket that cannot be read, an
/// interpretation without a gene or disease, fewer than two diagnoses, or a disease that is
/// not in `disease_models`.
pub fn analyze_phenopacket(
    ppkt: &Phenopacket,
    hpo: &FullCsrOntology,
    disease_models: &HashMap<TermId, SimpleDiseaseModel>,
    n_sim: usize,
) -> Result<PhenopacketAnalysis, String> {
    let prepared = prepare_phenopacket(ppkt)?;
    analyze_prepared(&prepared, &prepared, hpo, disease_models, n_sim)
}

impl PhenoblendSingleton {
    /// [`analyze_phenopacket`] with the ontology and the disease models of the application.
    /// Does not change the individual or the notebook.
    ///
    /// # Errors
    /// Returns an error if the ontology or the disease models are not loaded, or as
    /// [`analyze_phenopacket`].
    pub fn analyze_phenopacket(
        &self,
        ppkt: &Phenopacket,
        n_sim: usize,
    ) -> Result<PhenopacketAnalysis, String> {
        let hpo = self.hpo.as_ref().ok_or("HPO ontology not loaded")?;
        let models = self.omim_disease_models.as_ref().ok_or("Disease models not loaded")?;
        analyze_phenopacket(ppkt, hpo.as_ref(), models, n_sim)
    }
}

// -------------------------------------------------------------------------------------
// The folder
// -------------------------------------------------------------------------------------

/// What a batch run did.
/// What a batch run did.
#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct BatchSummary {
    /// Number of `.json` files found.
    pub n_files: usize,
    /// Number of files that gave results in the real analysis.
    pub n_files_analyzed: usize,
    /// Number of rows with `status = skipped` (left out on purpose, for example a non-genetic
    /// diagnosis).
    pub n_skipped_rows: usize,
    /// Number of rows with `status = error` (real analysis and control).
    pub n_error_rows: usize,
    /// Number of rows in the table, all kinds included.
    pub n_rows: usize,
}


/// The `.json` files of a directory, in name order (so a run is reproducible).
fn json_files(directory: &Path) -> Result<Vec<PathBuf>, String> {
    let mut files: Vec<PathBuf> = fs::read_dir(directory)
        .map_err(|e| format!("cannot read directory {}: {e}", directory.display()))?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.extension().is_some_and(|ext| ext.eq_ignore_ascii_case("json")))
        .collect();
    files.sort();
    Ok(files)
}

fn read_and_prepare(path: &Path) -> Result<PreparedPhenopacket, PrepareProblem> {
    let text = fs::read_to_string(path)
        .map_err(|e| PrepareProblem::Failed(format!("cannot read file: {e}")))?;
    let ppkt = load_from_json(&text)?;
    prepare_phenopacket(&ppkt)
}

fn push_result_rows(lines: &mut Vec<String>, analysis: &PhenopacketAnalysis, mode: &str) {
    let context = RowContext {
        phenopacket_id: &analysis.phenopacket_id,
        individual_id: &analysis.individual_id,
        mode,
        n_observed_terms: analysis.n_observed_terms,
    };
    for result in &analysis.results {
        lines.push(result_row(&context, result));
    }
}

/// Analyzes every `.json` phenopacket of `directory` and returns the table as text (comment
/// lines with `provenance`, a header, one row per result or error), with a summary.
///
/// With `shuffle_seed`, a control is added (`mode = shuffled`): each individual that was
/// analyzed keeps their own terms but is paired with the genes and diseases of another
/// individual, chosen by a seeded shuffle in which nobody keeps their own.
///
/// A file that cannot be read or analyzed gives one error row and does not stop the run. A
/// file that cannot be read or prepared is identified in its row by the file name; one that
/// fails later is identified by its phenopacket id.
///
/// # Errors
/// Returns an error only if the directory cannot be read.
pub fn run_batch(
    directory: &Path,
    hpo: &FullCsrOntology,
    disease_models: &HashMap<TermId, SimpleDiseaseModel>,
    n_sim: usize,
    shuffle_seed: Option<u64>,
    provenance: &[(&str, &str)],
) -> Result<(String, BatchSummary), String> {
    let files = json_files(directory)?;
    let mut lines: Vec<String> = Vec::new();

    // the real analysis
    let mut analyzed: Vec<PreparedPhenopacket> = Vec::new();
    for path in &files {
        let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
        match read_and_prepare(path) {
            Err(PrepareProblem::Skipped(message)) => {
                lines.push(skipped_row(&name, "observed", &message));
            }
            Err(PrepareProblem::Failed(message)) => {
                lines.push(error_row(&name, "observed", &message));
            }
            Ok(prepared) => match analyze_prepared(&prepared, &prepared, hpo, disease_models, n_sim) {
                Ok(analysis) => {
                    push_result_rows(&mut lines, &analysis, "observed");
                    analyzed.push(prepared);
                }
                Err(message) => {
                    lines.push(error_row(&prepared.phenopacket_id, "observed", &message));
                }
            },
        }
    }
    let n_files_analyzed = analyzed.len();

    // the control
    if let Some(seed) = shuffle_seed {
        match shuffle_assignment(analyzed.len(), seed) {
            Err(message) => lines.push(error_row("(control)", "shuffled", &message)),
            Ok(assignment) => {
                for (i, terms_of) in analyzed.iter().enumerate() {
                    let diagnoses_of = &analyzed[assignment[i]];
                    match analyze_prepared(terms_of, diagnoses_of, hpo, disease_models, n_sim) {
                        Ok(analysis) => push_result_rows(&mut lines, &analysis, "shuffled"),
                        Err(message) => {
                            lines.push(error_row(&terms_of.phenopacket_id, "shuffled", &message));
                        }
                    }
                }
            }
        }
    }

    let count_status = |status: &str| {
        lines.iter().filter(|line| line.split('\t').nth(2) == Some(status)).count()
    };
    let summary = BatchSummary {
        n_files: files.len(),
        n_files_analyzed,
        n_skipped_rows: count_status("skipped"),
        n_error_rows: count_status("error"),
        n_rows: lines.len(),
    };
    let mut text = provenance_lines(provenance);
    text.push_str(&tsv_header());
    text.push('\n');
    for line in &lines {
        text.push_str(line);
        text.push('\n');
    }
    Ok((text, summary))
}



/// [`run_batch`], writing the table to `output`.
///
/// # Errors
/// Returns an error if the directory cannot be read or the file cannot be written.
pub fn write_batch(
    directory: &Path,
    output: &Path,
    hpo: &FullCsrOntology,
    disease_models: &HashMap<TermId, SimpleDiseaseModel>,
    n_sim: usize,
    shuffle_seed: Option<u64>,
    provenance: &[(&str, &str)],
) -> Result<BatchSummary, String> {
    let (text, summary) = run_batch(directory, hpo, disease_models, n_sim, shuffle_seed, provenance)?;
    fs::write(output, text).map_err(|e| format!("cannot write {}: {e}", output.display()))?;
    Ok(summary)
}

/// Why a phenopacket was not analyzed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PrepareProblem {
    /// Left out on purpose; the message says why.
    Skipped(String),
    /// Could not be read, or is inconsistent.
    Failed(String),
}

impl From<String> for PrepareProblem {
    fn from(message: String) -> Self {
        Self::Failed(message)
    }
}

impl From<PrepareProblem> for String {
    fn from(problem: PrepareProblem) -> Self {
        match problem {
            PrepareProblem::Skipped(m) | PrepareProblem::Failed(m) => m,
        }
    }
}

/// The message for a phenopacket with an interpretation that gives no gene, presumably a
/// non-genetic diagnosis.
fn non_genetic_diagnosis(ppkt: &Phenopacket) -> Option<String> {
    ppkt.interpretations.iter().find_map(|interpretation| {
        if ga4ghphetools::get_gene_symbol_from_interpretation(interpretation).is_some() {
            return None;
        }
        let disease = interpretation
            .diagnosis
            .as_ref()
            .and_then(|d| d.disease.as_ref())
            .map(|class| format!("{} ({})", class.label, class.id))
            .unwrap_or_else(|| "no disease".to_string());
        Some(format!(
            "skipped phenopacket '{}': interpretation '{}' has no gene symbol (disease: {}); \
             this is presumably a non-genetic diagnosis, and the analysis needs a gene for every diagnosis",
            ppkt.id, interpretation.id, disease
        ))
    })
}