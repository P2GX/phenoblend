export interface SimpleDiseaseModel {
  omimDiseaseId: string;
  omimDiseaseName: string;
  observedHpoIds: string[];
  excludedHpoIds: string[];
}

export interface GeneDiseaseAssociation {
  ncbiGeneId: string;
  geneSymbol: string;
  associationType: string;
  diseaseId: string;
  diseaseModel: SimpleDiseaseModel | null;
  source: string;
}



/** Mirrors Rust's Method enum (serialized as lowercase strings). */
export type EnrichmentMethod = 'unweighted' | 'frequency';

/** Mirrors Rust's EnrichmentResult. */
export interface EnrichmentResult {
  /** pmf[k] = P(duo count == k) under the null, for k in 0..=m_eff */
  pmf: number[];
  pValue: number;
  /** true if the draw count would exhaust the urn (null model is degenerate) */
  degenerate: boolean;
}

/** Mirrors Rust's DuoEnrichmentSummary. */
export interface DuoEnrichmentSummary {
  geneA: string;
  geneB: string;
  nUnion: number;
  nShared: number;
  nObsUnion: number;
  nObsShared: number;
  hasOverlap: boolean;
  /** null only if EnrichmentTest::run returned None (shouldn't happen for UNWEIGHTED) */
  unweighted: EnrichmentResult | null;
  /** null if no term in the urn had a positive frequency weight */
  frequency: EnrichmentResult | null;
}

export interface ObservedHpoInput {
  termId: string;
  /** All modifier term IDs; the backend picks out the severity terms. */
  modifierIds: string[];
}

export type GeneDiseaseAnnotations = Record<string, GeneDiseaseAssociation[]>;