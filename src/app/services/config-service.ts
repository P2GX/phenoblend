import { Injectable } from '@angular/core';
import { invoke } from "@tauri-apps/api/core";

import { FenominalSentence, HierarchyMapItem, OntologyMatch, HpoTermMinimal } from 'ng-hpo-uikit';
import { InitializationStatusDto } from '../models/status_dto'; 
import { DuoEnrichmentSummary, GeneDiseaseAnnotations, GeneDiseaseAssociation } from '../models/interfaces';
import { PresenceMatrixPayload, SpreadPlotPayload, UpsetPlotPayload } from 'ngx-phenoprofile';

@Injectable({
  providedIn: 'root'
})
export class ConfigService {
 

  async loadHPO(): Promise<void> {
    return await invoke("load_hpo");
  }

  async loadHpoas(): Promise<void> {
    return await invoke("load_hpoas");
  }

  async loadGeneToDisease(): Promise<void> {
    return await invoke("load_gene_disease_associations");
  }

  async ingestPhenopacket(ppkt: string): Promise<void> {
    return await invoke("ingest_phenopacket", {'ppktString': ppkt});
  }

  async getOverlapPlotData(annotationMap: Record<string, GeneDiseaseAssociation[]>): Promise<PresenceMatrixPayload> {
    return await invoke("get_overlap_plot", {annotations: annotationMap});
  }

  async getUpsetPlotPayload(annotationMap: Record<string, GeneDiseaseAssociation[]>): Promise<UpsetPlotPayload> {
    return await invoke<UpsetPlotPayload>("get_upset_plot_payload", {annotations: annotationMap});
  }

  async getSpreadPlotPayload(annotationMap: Record<string, GeneDiseaseAssociation[]>): Promise<SpreadPlotPayload> {
    return await invoke<SpreadPlotPayload>("get_spread_plot_payload", {annotations: annotationMap});
  }

  async getAutocompleteHpo(value: string): Promise<OntologyMatch[]> {
    return invoke<OntologyMatch[]>('get_hpo_autocomplete', { query: value });
  }

  async mineClinicalText(text: string): Promise<FenominalSentence[]> {
    console.log('[mineClinicalText] invoking mine_clinical_text with text:', text);
    try {
      const result = await invoke<FenominalSentence[]>('mine_clinical_text', { text });
      console.log('[mineClinicalText] received result:', result);
      return result;
    } catch (error) {
      console.error('[mineClinicalText] invoke failed:', error);
      throw error;
    }
  }

  async checkInitializationStatus(): Promise<InitializationStatusDto> {
    return invoke<InitializationStatusDto>('check_initialization_status');
  }

  async getHpoParentAndChildrenTerms(termId: string): Promise<HierarchyMapItem> {
    return invoke<HierarchyMapItem>('get_hpo_parent_and_children_terms', { termId });
  }

  async getHpoModifiers(): Promise<HpoTermMinimal[]> {
    return invoke<HpoTermMinimal[]>('get_hpo_modifiers');
  }

  async performHpoAutocomplete(query: string): Promise<OntologyMatch[]> {
    return invoke<OntologyMatch[]>('perform_hpo_autocomplete', { query });
  }

  async  autocompleteGeneSymbol(query: string): Promise<GeneDiseaseAssociation[]> {
    return invoke<GeneDiseaseAssociation[]>('autocomplete_gene_symbol', { query });
  }

  async getObservedHpoCount(): Promise<number> {
    return invoke<number>('get_observed_hpo_count');
  }

  async saveSvgFile(svgContent: string, defaultFilename: string): Promise<boolean> {
    return invoke('save_svg_file', { svgContent, defaultFilename });
  }

  async exportSvgToPdf(svgContent: string, defaultFilename: string): Promise<boolean> {
    return invoke('export_svg_to_pdf', { svgContent, defaultFilename });
  }

  async addObservedHposFromNER(observed: string[]): Promise<null> {
    return invoke<null>('add_observed_hpos_from_ner', {observed});
  }

  async downloadExcelSummary(
    dataType: string,
    annotations: Record<string, GeneDiseaseAssociation[]>,
  ): Promise<void> {
    await invoke<void>('download_excel_summary', {
      dataType,
      annotations,
    });
  }

  /**
 * Runs the duo-enrichment test for every pairwise combination of genes
 * present in `annotations`. Observed HPO terms and the ontology are read
 * from the backend singleton; only the candidate gene/disease annotations
 * and (optionally) the Monte Carlo sample size need to be supplied.
 *
 * @param annotations gene symbol -> its GeneDiseaseAssociation list
 * @param nSim        Monte Carlo draws for the FREQUENCY null (backend
 *                    defaults to 100_000 if omitted; consider passing a
 *                    smaller value, e.g. 5_000-10_000, for snappier
 *                    interactive use)
 */
async analyzeDuoEnrichment(
  annotations: GeneDiseaseAnnotations,
  nSim?: number,
): Promise<DuoEnrichmentSummary[]> {
  return invoke<DuoEnrichmentSummary[]>('analyze_duo_enrichment', {
    annotations,
    nSim,
  });
}
  
}