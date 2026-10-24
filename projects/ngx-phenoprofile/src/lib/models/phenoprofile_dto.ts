
export type Severity = 'borderline' | 'mild' | 'moderate' | 'severe' | 'profound';


export interface OverlapPlotItem {
  hpoId: string;
  hpoName: string;
  severity?: Severity;
  // Keyed by Gene Symbol, value ranges from [0.0 to 1.0]
  scores: { [geneSymbol: string]: number };
}

export interface OverlapPlotPayload {
  entities: string[]; // Already sorted by your column ranking rules
  columns: OverlapPlotItem[]; // Already sorted by your tier blocks
}

export interface UpsetPlotPayload {
  genes: string[];
  combinations: string[][];
  combinationAnnotated: number[];
  combinationObserved: number[];
  geneAnnotated: number[];
  geneObserved: number[];
}



export interface UpsetRowData { 
  geneSymbol: string; 
  nProbandTerms: number; 
  nOverlappingTerms: number; 
}


export interface SimpleOntologyTerm { termId: string; termLabel: string; }
export interface ProbandTermMatch { term: SimpleOntologyTerm; annotatedBy: SimpleOntologyTerm[]; }
export interface MatchedTerm { term: SimpleOntologyTerm; matchedBy: SimpleOntologyTerm[]; }

export interface UpsetRowData {
  geneSymbol: string; nProbandTerms: number; nOverlappingTerms: number;
  overlappingTerms: ProbandTermMatch[]; nonOverlappingTerms: SimpleOntologyTerm[];
}

export interface UpsetColumnData { 
  rowIndices: number[]; 
  nAnnotatedTerms: number; 
  nOverlappingTerms: number; 
  matchedTerms: MatchedTerm[];
  unmatchedTerms: SimpleOntologyTerm[];
}

export interface ProbandSegment { 
  rowIndices: number[]; 
  nTerms: number; 
  terms: SimpleOntologyTerm[];
}

export interface UpsetSummary {
  nProbandTerms: number; nProbandExplained: number; nProbandUnexplained: number;
  maxColumn: number; maxRow: number;
}
export interface UpsetPlotData { 
  rows: UpsetRowData[]; 
  columns: UpsetColumnData[]; 
  probandBreakdown: ProbandSegment[];
  summary: UpsetSummary; 
}


export interface SpreadPlotCategory {
  id: string;
  name: string;
  alias?: string;
  ppktValue: number;         // Patient/Phenopacket fraction
  geneValues: number[];      // Array of values for each gene combo/series
}

export interface SpreadPlotPayload {
  seriesLabels: string[];    // e.g., ["Ppkt", "MBNL1", "DMPK", "MBNL1+DMPK"]
  categories: SpreadPlotCategory[];
}