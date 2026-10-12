// notebook.model.ts
export interface SimpleOntologyTerm {
  termId: string;
  termLabel: string;
}

export interface GeneDiseaseEntity {
  // TODO: mirror your Rust struct; I haven't seen it
}

export interface Notebook {
  id: string;
  timestamp: number;
  observedHpo: SimpleOntologyTerm[];
  dataInputHistory: string;
  geneDiseaseList: GeneDiseaseEntity[];
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