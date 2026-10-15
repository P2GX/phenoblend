// notebook.model.ts
export interface SimpleOntologyTerm {
  termId: string;
  termLabel: string;
}

export interface GeneDiseaseEntity {
  // TODO: mirror Rust struct; 
}

export interface SimpleOntologyTerm {
    termId: string;
     termLabel: string;
}

export type MatchType = 'exact' | 'descendant' | 'partial' | 'noMatch';

export interface GeneMatch {
    matchType: MatchType,
    score: number;
    matchedTerm: SimpleOntologyTerm;
}


export interface TermMatches {
    term: SimpleOntologyTerm,
    matches: Record<string, GeneMatch>;
}


export interface OverlapData {
    entities: string[];
     columns: TermMatches[];
}

export interface MatchedTerm {
    term: SimpleOntologyTerm,
    matchedBy: SimpleOntologyTerm,
}

/// A group of disease annotations split by whether the proband has them.
export interface TermSplit {
  /** Green: annotations the proband has (exactly, or through a descendant). */
  probandOverlapped: MatchedTerm[];
  /** Purple: annotations the proband does not have. */
  remainingAnnotated: SimpleOntologyTerm[];
}

/// One **horizontal** bar: all explicit annotations of one disease/gene entity, including
/// terms it shares with the others.
export interface UpsetGene extends TermSplit {
  geneSymbol: string;
}



export interface UpsetIntersection {
     geneSymbolList: string[];
     probandOverlappedTids: string[];
    /** HPO Terms annotation to the gene but not reported in the proband (explicit only) */
    remainingAnnotatedTids: string[];
}

export interface UpsetData {
  /** Horizontal bars, sorted by gene symbol. */
  genes: UpsetGene[];
  /** Vertical bars, sorted by number of entities and then alphabetically. */
  intersections: UpsetIntersection[];
  /** 
   * Proband terms that no disease annotation accounts for (no disease term equals the 
   * term or is one of its ancestors). They are not drawn in the plot. 
   */
  probandOnly: SimpleOntologyTerm[];
  /** Number of distinct observed terms of the proband. */
  nProbandTerms: number;
}

export interface SpreadPlotCategory {
  id: string;
  name: string;
  /** Optional field mapping to `alias?: string` */
  alias?: string;
  /** Patient/Phenopacket fraction (`ppktValue`) */
  ppktValue: number;
  /** Array of values for each gene combo/series (`geneValues`) */
  geneValues: number[];
}

export interface SpreadPlotPayload {
  /** e.g., ["Ppkt", "MBNL1", "DMPK", "MBNL1+DMPK"] */
  seriesLabels: string[];
  categories: SpreadPlotCategory[];
}

export interface PhenoBlendAnalysis {
    /// Source of truth: how each observed proband term matches each gene entity.
    termMatches: TermMatches[];
    /** All annotated terms per gene entity (needed later for the UpSet "annotated" bars). */
  geneAnnotations: Record<string, SimpleOntologyTerm[]>;
  overlap: OverlapData;
  upset: UpsetData;
  spread: SpreadPlotPayload;
}

export interface Notebook {
  id: string;
  timestamp: number;
  observedHpo: SimpleOntologyTerm[];
  dataInputHistory: string;
  geneDiseaseList: GeneDiseaseEntity[];
  analysis?: PhenoBlendAnalysis
}


export function initialNotebook(): Notebook {
    let n: Notebook = {
        id: "",
        timestamp: 0,
        observedHpo: [],
        dataInputHistory: "not initialized",
        geneDiseaseList: []
    };
    return n;
}