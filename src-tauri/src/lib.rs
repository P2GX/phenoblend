pub mod phenoblend;

mod enrichment;
mod excel;
mod hpoa;
mod hpo;
mod model;
mod blend;
mod util;
#[cfg(test)]
pub mod tests;


use ontolius::ontology::MetadataAware;
use serde::{self,Serialize};

use fenominal::OntologyMatch;
use fenominal::FenominalSentence;
use ga4ghphetools::dto::hpo_term_dto::HpoTermDuplet;
use ga4ghphetools::HierarchyMapItem;
use ontolius::ontology::OntologyTerms;
use tauri_plugin_dialog::DialogExt;
use tauri::{AppHandle, Emitter, WindowEvent};
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use ga4ghphetools::{pick_file_and_process, load_ontology, OntologyLoadEvent};


use crate::blend::dto::UpsetPlotPayload;
use crate::blend::dto::SpreadPlotPayload;
use crate::enrichment::dto::DuoEnrichmentSummary;
use crate::model::Notebook;
use crate::{blend::dto::OverlapPlotPayload, phenoblend::PhenoblendSingleton};
use crate::model::status::InitializationStatusDto;
use crate::model::emit_notebook;
use crate::hpoa::disease_model::GeneDiseaseAssociation;
use tauri_plugin_log::{Target, TargetKind};


struct AppState {
    phenoblend: Mutex<PhenoblendSingleton>,
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    let app_state = Arc::new(AppState {
        phenoblend: Mutex::new(PhenoblendSingleton::new()),
    });
    let log_level = log::LevelFilter::Trace;

    tauri::Builder::default()
        .manage(app_state)
        .plugin(tauri_plugin_fs::init())
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_opener::init())
        .plugin(tauri_plugin_shell::init())     
        .plugin(
            tauri_plugin_log::Builder::new()
                .level(log_level)
                //.level_for("tao", log::LevelFilter::Warn)
               // .level_for("winit", log::LevelFilter::Warn)
                //.level_for("tauri", log::LevelFilter::Warn)
                .targets([
                    Target::new(TargetKind::Stdout),
                   // Target::new(TargetKind::LogDir { file_name: Some("phenoblendtk.log".to_string()) }), 
                   // Target::new(TargetKind::Webview), // forwards logs to the browser devtools console
                ])
                .build(),
        )
        .invoke_handler(tauri::generate_handler![
            add_observed_hpos_from_ner,
            analyze_duo_enrichment,
            autocomplete_gene_symbol,
            check_initialization_status,
            download_excel_summary,
            export_svg_to_pdf,
            get_hpo_autocomplete,
            get_hpo_modifiers,
            get_hpo_parent_and_children_terms,
            get_observed_hpo_count,
            get_overlap_plot,
            get_spread_plot_payload,
            get_upset_plot_payload,
            ingest_phenopacket,
            initialize_notebook,
            load_hpo,
            load_hpoas,
            load_gene_disease_associations,
            mine_clinical_text,
            perform_hpo_autocomplete,
            save_svg_file
        ])
        .setup(|_app| {
            Ok(())
        })
        .on_window_event(|window, event| {
            if let WindowEvent::CloseRequested { api, .. } = event {
                api.prevent_close();
                window.emit("close-requested", ()).unwrap_or_default();
            }
        })
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}




#[tauri::command]
fn ingest_phenopacket(
    app: AppHandle,
    state: tauri::State<'_, Arc<AppState>>,
    ppkt_string: String
) -> Result<(), String> {
    let result: Result<Notebook, String> = (|| {
        let mut singleton = state.phenoblend.lock().map_err(|e|e.to_string())?;
        let phenopacket = crate::hpo::ppkt_loader::load_from_json(&ppkt_string)?;
        singleton.ingest_ppkt(phenopacket)?;
        Ok(singleton.notebook_ref()?.clone())
    })();
    eprintln!("result={:?}", result);
    match result {
        Ok(snapshot) => {
            if let Err(e) = emit_notebook(&app, &snapshot) {
                log::error!("Failed to emit notebook: {e}");
            }
            Ok(())
        },
        Err(e) => {Err(e)},
    }
}


/// Load the Human Phenotype Ontology (HPO)
#[tauri::command]
async fn load_hpo(
    app: AppHandle,
    state: tauri::State<'_, Arc<AppState>>,
) -> Result<(), String> {
    let state_handle = state.inner().clone();
    let _ = app.emit("hpo-load-event", OntologyLoadEvent::loading());
    pick_file_and_process(app, "hpo-load-event", move |hpo_json_path, app_handle| async move {
        match load_ontology(&hpo_json_path) {
            Ok(ontology) => {
                let mut singleton = state_handle.phenoblend.lock().unwrap();
                let n_terms = ontology.len();
                let version = ontology.version().to_string();
                singleton.set_hpo(ontology, &hpo_json_path);
                let _ = app_handle.emit(
                    "hpo-load-event", 
                    OntologyLoadEvent::success(version, n_terms)
                );
            },
            Err(e) => { 
                let _ = app_handle.emit("hpo-load-event", OntologyLoadEvent::error(e.to_string()));
            }
        }
    });

    Ok(())
}


#[tauri::command]
async fn add_observed_hpos_from_ner(
    state: tauri::State<'_, Arc<AppState>>,
    observed: Vec<String>
) -> Result<(), String>{
    let state_handle = state.inner().clone();
    let mut singleton = state_handle.phenoblend.lock().map_err(|e| e.to_string())?;
    singleton.add_observed_hpos_from_ner(observed)?;
    Ok(())
}

#[tauri::command]
async fn load_hpoas(
    app: AppHandle,
    state: tauri::State<'_, Arc<AppState>>,
) -> Result<(), String> {
 let state_handle = state.inner().clone();
    let _ = app.emit("hpoa-load-event", OntologyLoadEvent::loading());
    pick_file_and_process(app, "hpoa-load-event", move |hpoa_path, app_handle| async move {
        match crate::hpoa::hpoa_ingest::load_hpoa_d(&hpoa_path) {
            Ok(disease_model_map) => {
                let mut singleton = state_handle.phenoblend.lock().unwrap();
                let n_terms = disease_model_map.len();
                singleton.set_hpoa_d(disease_model_map, &hpoa_path);
                let _ = app_handle.emit(
                    "hpoa-load-event", 
                    OntologyLoadEvent::success("HPOAs loaded".to_string(), n_terms)
                );
            },
            Err(e) => { 
                let _ = app_handle.emit("hpoa-load-event", OntologyLoadEvent::error(e.to_string()));
            }
        }
    });

    Ok(())
}


#[tauri::command]
async fn load_gene_disease_associations(
    app: AppHandle,
    state: tauri::State<'_, Arc<AppState>>,
) -> Result<(), String> {
 let state_handle = state.inner().clone();
    let _ = app.emit("g2d-load-event", OntologyLoadEvent::loading());
    pick_file_and_process(app, "g2d-load-event", move |g2d_path, app_handle| async move {
        match crate::hpoa::gene_to_disease::load_gene_disease_associations(&g2d_path) {
            Ok(g2d) => {
                let mut singleton = state_handle.phenoblend.lock().unwrap();
                let n_terms = g2d.len();
                singleton.set_gene_to_disease(g2d, &g2d_path);
                let _ = app_handle.emit(
                    "g2d-load-event", 
                    OntologyLoadEvent::success("gene to disease loaded".to_string(), n_terms)
                );
                println!("emitting g2d");
            },
            Err(e) => { 
                let _ = app_handle.emit("g2d-load-event", OntologyLoadEvent::error(e.to_string()));
            }
        }
    });

    Ok(())
}

#[tauri::command]
fn get_overlap_plot(
    state: tauri::State<'_, Arc<AppState>>,
) -> Result<OverlapPlotPayload, String> {
    
    let singleton = state.phenoblend.lock().map_err(|e| e.to_string())?;
    singleton.get_overlap_plot_payload()
}

#[tauri::command]
fn get_upset_plot_payload(
    state: tauri::State<'_, Arc<AppState>>,
) -> Result<UpsetPlotPayload, String> {
    let state_handle = state.inner().clone();
    let singleton = state_handle.phenoblend.lock().map_err(|e| e.to_string())?;
    singleton.get_upset_plot_payload()
}

#[tauri::command]
fn get_spread_plot_payload( state: tauri::State<'_, Arc<AppState>>,
) -> Result<SpreadPlotPayload, String> {
    let state_handle = state.inner().clone();
    let singleton = state_handle.phenoblend.lock().map_err(|e| e.to_string())?;
    singleton.get_spread_plot_payload()
  }
  

#[tauri::command]
fn get_observed_hpo_count( state: tauri::State<'_, Arc<AppState>>,
) -> usize {
    let state_handle = state.inner().clone();
    // We use this to show the number of HPOs in. the GUI. If this fails, just return zero, but we never expect this to happen
    let singleton = match state_handle.phenoblend.lock() {
        Ok(guard) => guard,
        Err(_) => return 0,
    };
    singleton.get_observed_hpo_count()
}


/// This function supplies the autocompletion candidates for angular for the HPO
/// The JavaScript ensures that query is at least 3 letters
#[tauri::command]
fn get_hpo_autocomplete(
    state: tauri::State<'_, Arc<AppState>>,
    query: String
) -> Vec<OntologyMatch> {
    let singleton = match state.phenoblend.lock() {
        Ok(s) => s,
        Err(_) => return vec![],
    };

    // If query is too short, don't even bother searching
    if query.len() < 3 {
        return vec![];
    }
    singleton.search_hpo(&query, 20)
}

#[tauri::command]
async fn mine_clinical_text(
    state: tauri::State<'_, Arc<AppState>>,
    text: String,
) -> Result<Vec<FenominalSentence>, String> {
    let singleton = match state.phenoblend.lock() {
        Ok(s) => s,
        Err(_) => return Err("Failed to acquire application state lock".to_string()),
    };
    singleton
        .mine_clinical_text(&text)
}


#[tauri::command]
async fn check_initialization_status(
    state: tauri::State<'_, Arc<AppState>>,
) -> Result<InitializationStatusDto, String> {
    let singleton = state.phenoblend.lock()
        .map_err(|_| "Failed to lock state".to_string())?;

    Ok(InitializationStatusDto {
        hpo_loaded: singleton.hpo.is_some(),
        hpo_terms: singleton.hpo.as_ref().map(|h| h.len()).unwrap_or(0), // adjust based on your ontology structure length method
        hpoa_loaded: singleton.omim_disease_models.is_some(),
        hpoa_diseases: singleton.omim_disease_models.as_ref().map(|m| m.len()).unwrap_or(0),
        g2d_loaded: singleton.gene_to_disease_d.is_some(),
        n_genes: singleton.gene_to_disease_d.as_ref().map(|g| g.len()).unwrap_or(0),
    })
}



#[tauri::command]
async fn get_hpo_parent_and_children_terms(
    state: tauri::State<'_, Arc<AppState>>,
    term_id: &str,
) -> Result<HierarchyMapItem, String> {
    let singleton = state.phenoblend.lock()
        .map_err(|_| "Failed to lock state".to_string())?;
    singleton.get_hpo_parent_and_children_terms(term_id)
}

/// format matching the TypeScript `HpoTermMinimal` interface in ng-hpo-uikit.
#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct HpoTermMinimalDto {
    pub term_id: String,
    pub label: String,
}

impl From<HpoTermDuplet> for HpoTermMinimalDto {
    fn from(d: HpoTermDuplet) -> Self {
        HpoTermMinimalDto {
            term_id: d.hpo_id,
            label: d.hpo_label,
        }
    }
}

#[tauri::command]
async fn get_hpo_modifiers(
    state: tauri::State<'_, Arc<AppState>>
) -> Result<Vec<HpoTermMinimalDto>, String> {
    let singleton = state.phenoblend.lock()
        .map_err(|_| "Failed to lock state".to_string())?;
    let duplets = singleton.get_modifiers()?;

    Ok(duplets.into_iter().map(HpoTermMinimalDto::from).collect())
}




#[tauri::command]
async fn perform_hpo_autocomplete(state: tauri::State<'_, Arc<AppState>>, query: String) -> Result<Vec<OntologyMatch>, String> {
    let singleton = state.phenoblend.lock()
        .map_err(|_| "Failed to lock state".to_string())?;
    singleton.perform_hpo_autocomplete(query)
}


#[tauri::command]
async fn autocomplete_gene_symbol(
    state: tauri::State<'_, Arc<AppState>>,
    query: &str,
) -> Result<Vec<GeneDiseaseAssociation>, String> {
    let singleton = state.phenoblend.lock()
        .map_err(|_| "Failed to lock state".to_string())?;
    let limit = 20;
    singleton.autocomplete_gene_symbol(query, limit)
}


#[tauri::command]
async fn save_svg_file(
    app: tauri::AppHandle,
    svg_content: String,
    default_filename: String,
) -> Result<bool, String> {
    let file_path = app
        .dialog()
        .file()
        .set_file_name(&default_filename)
        .add_filter("SVG Image", &["svg"])
        .blocking_save_file();

    match file_path {
        Some(path) => {
            let path = path.as_path().ok_or("Invalid path returned from dialog")?;
            std::fs::write(path, svg_content).map_err(|e| e.to_string())?;
            Ok(true)
        }
        None => Ok(false), // user cancelled — not an error
    }
}




#[tauri::command]
async fn export_svg_to_pdf(
    app: tauri::AppHandle,
    svg_content: String,
    default_filename: String,
) -> Result<bool, String> {
    // 1. Parse the SVG into a usvg tree
    let mut font_db = svg2pdf::usvg::fontdb::Database::new();
    font_db.load_system_fonts();
    let font_db = Arc::new(font_db);

    let opt = svg2pdf::usvg::Options {
        fontdb: font_db.clone(),
        ..Default::default()
    };

    let tree = svg2pdf::usvg::Tree::from_str(&svg_content, &opt)
        .map_err(|e| format!("SVG parse error: {e}"))?;

    // 2. Convert to PDF bytes
    let pdf_bytes = svg2pdf::to_pdf(
        &tree,
        svg2pdf::ConversionOptions::default(),
        svg2pdf::PageOptions::default(),
    )
    .map_err(|e| format!("PDF conversion error: {e}"))?;

    // 3. Show native save dialog
    let file_path = app
        .dialog()
        .file()
        .set_file_name(&default_filename)
        .add_filter("PDF Document", &["pdf"])
        .blocking_save_file();

    match file_path {
        Some(path) => {
            let path = path.as_path().ok_or("Invalid path returned from dialog")?;
            std::fs::write(path, pdf_bytes).map_err(|e| e.to_string())?;
            Ok(true)
        }
        None => Ok(false), // user cancelled — not an error
    }
}


#[tauri::command]
async fn download_excel_summary(
    app: tauri::AppHandle,
    state: tauri::State<'_, Arc<AppState>>,
    data_type: String,
    annotations: HashMap<String, Vec<GeneDiseaseAssociation>>) -> Result<(), String>
{
    let mut singleton = state.phenoblend.lock()
        .map_err(|_| "Failed to lock state".to_string())?;
    let fname = format!("{}-annotations.xlsx", data_type);
    let file_path = app
        .dialog()
        .file()
        .add_filter("Excel Workbook", &["xlsx"])
        .set_file_name(fname)
        .blocking_save_file();

    let Some(path) = file_path else {
        return Err("Export cancelled".to_string());
    };

    let path_str = path.to_string();
    singleton.output_excel(data_type, annotations, path_str)?;

    Ok(())

}




const DEFAULT_N_SIM: usize = 100_000;

#[tauri::command]
fn analyze_duo_enrichment(
    state: tauri::State<'_, Arc<AppState>>,
    annotations: HashMap<String, Vec<GeneDiseaseAssociation>>,
    n_sim: Option<usize>,
) -> Result<Vec<DuoEnrichmentSummary>, String> {
    let state_handle = state.inner().clone();
    let singleton = state_handle.phenoblend.lock().map_err(|e| e.to_string())?;
    singleton.get_duo_enrichment_payload(annotations, n_sim.unwrap_or(DEFAULT_N_SIM))
}

#[tauri::command]
fn initialize_notebook(
    app: tauri::AppHandle,
    state: tauri::State<'_, Arc<AppState>>,
    gda_list: HashMap<String, Vec<GeneDiseaseAssociation>>) -> Result<(), String> {
    let mut singleton = state.phenoblend.lock().map_err(|e| e.to_string())?;
    let notebook = singleton.initialize_notebook(gda_list)?;
    if let Err(e) = emit_notebook(&app, &notebook) {
        log::error!("Failed to emit notebook: {e}");
    }
    Ok(())
}