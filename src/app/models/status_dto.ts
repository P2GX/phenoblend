export interface StatusDto {
    hpoLoaded: boolean;
    hpoVersion: string;
    nHpoTerms: number;
    ptTemplatePath: string;
    ptTemplateLoaded: boolean;
    biocuratorOrcid: string;
    hpoJsonPath: string;
    hasError: boolean;
    errorMessage: string;
}



export function defaultStatusDto(): StatusDto {
    return {
        hpoLoaded: false,
        hpoVersion: "not initialized",
        nHpoTerms: 0,
        ptTemplatePath: "not initialized",
        ptTemplateLoaded: false,
        biocuratorOrcid: "not initialized",
        hpoJsonPath: "not initialized",
        hasError: false,
        errorMessage: ""
    };
}




export interface InitializationStatusDto {
  hpo_loaded: boolean;
  hpo_terms: number;
  hpo_version: string;
  hpoa_loaded: boolean;
  hpoa_diseases: number;
  g2d_loaded: boolean;
  n_genes: number;
}