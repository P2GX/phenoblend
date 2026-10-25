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





export interface ObservedHpoInput {
  termId: string;
  /** All modifier term IDs; the backend picks out the severity terms. */
  modifierIds: string[];
}

export type GeneDiseaseAnnotations = Record<string, GeneDiseaseAssociation[]>;