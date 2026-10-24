
use std::collections::{BTreeMap, BTreeSet};
use ontolius::{TermId, ontology::{TaxonomyWalk, csr::FullCsrOntology}};
use serde::{Serialize, Deserialize};


use crate::model::simple_term::SimpleOntologyTerm;


pub type TermKey = String; // CURIE, e.g. "HP:0001250"
pub type GeneKey = String; // gene symbol, e.g. "NF1"



/// HPO annotations of one gene, pooled across all diseases associated with it.
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct GeneTerms {
    /// Terms directly annotated to the gene's disease(s).
    pub annotated: BTreeSet<TermKey>,
    /// `annotated` plus all ancestors (true-path closure); always a superset of `annotated`.
    pub closure: BTreeSet<TermKey>,
}


/// The true-path closure of `t`: `t` itself plus all of its ancestors.
pub fn closure_of(hpo: &FullCsrOntology, t: &TermId) -> BTreeSet<TermKey> {
    let closure: BTreeSet<TermKey> = hpo
        .iter_term_and_ancestor_ids(t)
        .map(|a| a.to_string())
        .collect();
    closure
}

impl GeneTerms {
    pub fn add_term(&mut self, hpo: &FullCsrOntology, t: &TermId) {
        self.annotated.insert(t.to_string());
        self.closure.extend(closure_of(hpo, t));
    }
}


/// One part of the proband's terms, grouped by the exact set of genes that overlap them.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ProbandSegment {
    /// Indices into `rows` of the overlapping genes, ascending. Empty for terms that no gene
    /// overlaps.
    pub row_indices: Vec<usize>,
    /// Number of proband terms whose overlapping genes are exactly `row_indices`.
    pub n_terms: usize,
    /// The proband terms, sorted by term id
    pub terms: Vec<SimpleOntologyTerm>,
}



/// Term closures of the proband and of every gene: the input to [`plot_data`].
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetPlotInput {
    /// Key: an HPO term observed in the proband; value: that term plus all of its ancestors.
    /// Duplicate observed terms collapse into one entry.
    pub proband_closures: BTreeMap<TermKey, BTreeSet<TermKey>>,
    /// Key: gene symbol; value: HPO annotations of all diseases associated with the gene.
    pub gene_annotations: BTreeMap<GeneKey, GeneTerms>,
    /// Display data for every observed and every annotated term. Ancestors that are neither
    /// are never shown, so they are not included.
    pub terms: BTreeMap<TermKey, SimpleOntologyTerm>,
}






/// A disease annotation that the proband has, with the proband term that accounts for it.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct MatchedTerm {
    /// The disease's own (explicit) annotation, as shown in the bar.
    pub term: SimpleOntologyTerm,
    /// Observed proband terms equal to `term` or one of its descendants. Never empty.
    /// If the disease has “Seizure” and the individual has “Focal seizure” and “Tonic seizure”, 
    /// both are observed and both fall under the one disease term.
    pub matched_by: Vec<SimpleOntologyTerm>,
}





#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetPlotData {
    /// Horizontal bars (one per gene), sorted by gene symbol.
    ///
    /// The full length of every bar is the number of distinct HPO terms of the proband;
    /// the highlighted part is the number of those terms that the gene accounts for
    /// (see the module documentation for the matching rule).
    ///
    /// Example: the proband has 20 HPO terms, 4 of which overlap gene A and 7 gene B.
    /// Both bars are 20 units long, with highlights of 4 units for A and 7 units for B.
    pub rows: Vec<UpsetRowData>,

    /// Vertical bars, one per gene and one per combination of 2 to 4 genes, sorted by
    /// `n_overlapping_terms` descending, then `n_annotated_terms` descending, then
    /// `row_indices`.
    ///
    /// Combinations are *inclusive*: the A+B column also counts terms that gene C has.
    /// Only genes that overlap at least one proband term are combined, and combinations
    /// without any overlap are omitted.
    ///
    /// Example: gene A has 40 annotations, of which 4 overlap the proband, so its column is
    /// 40 units high with a highlight of 4. For A+B, the annotations common to both genes
    /// (see [`UpsetColumnData::n_annotated_terms`]) set the height.
    pub columns: Vec<UpsetColumnData>,
    /// The proband's terms partitioned by exactly which genes overlap them, for the bar
    /// beneath the plot. The counts sum to `summary.n_proband_terms`. The first segment is
    /// always the terms no gene explains (empty `row_indices`); the rest are ordered by the
    /// number of genes, then by `row_indices`. With at most four genes every combination is
    /// listed, even with a count of zero; with more genes only combinations that occur are.
    pub proband_breakdown: Vec<ProbandSegment>,

    pub summary: UpsetSummary,
}





#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ProbandTermMatch {
    pub term: SimpleOntologyTerm,
    /// Gene annotations equal to `term` or one of its ancestors. Never empty.
    pub annotated_by: Vec<SimpleOntologyTerm>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetRowData {
    pub gene_symbol: String,
    /// Number of distinct HPO terms observed in the proband (an individual represented by a
    /// phenopacket or a list of observed HPO terms). Identical in every row; it is the full
    /// bar length.
    pub n_proband_terms: usize,
    /// How many of those proband terms overlap the gene's annotations (pooled across all of
    /// the gene's diseases): the term is annotated to the gene or is a descendant of an
    /// annotated term (true-path rule).
    ///
    /// Formally, the number of proband terms `p` with
    /// `gene.annotated ∩ proband_closures[p] ≠ ∅`.
    pub n_overlapping_terms: usize,
    pub overlapping_terms: Vec<ProbandTermMatch>,
    pub non_overlapping_terms: Vec<SimpleOntologyTerm>,
}

/// One column of the pheno-upset plot: the phenotypic profile of one gene (e.g. NF1) or of a
/// combination of genes (e.g. NF1+PTPN11).
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetColumnData {
    /// Filled dots: indices into `rows`, ascending.
    pub row_indices: Vec<usize>,
    /// Number of distinct HPO terms shared by all genes of the column under the true-path
    /// rule. A term counts if every gene is annotated with it or with one of its descendants,
    /// and it is explicitly annotated to at least one of the genes. The more general term
    /// represents the shared annotation: if gene A has HP:1, gene B has HP:2 and HP:1 is a
    /// descendant of HP:2, then HP:2 is counted once, but ancestors of HP:2 that no gene
    /// annotates are not. For a single gene this is the number of its annotated terms.
    pub n_annotated_terms: usize,
    /// How many of the `n_annotated_terms` shared terms the proband matches, that is, that
    /// equal or are an ancestor of an observed proband term. Each shared term counts once,
    /// however many proband terms fall under it, so this never exceeds `n_annotated_terms`.
    pub n_overlapping_terms: usize,
     /// Shared terms the proband matches, each with the observed terms that match it.
    pub matched_terms: Vec<MatchedTerm>,
    /// Shared terms the proband does not match.
    pub unmatched_terms: Vec<SimpleOntologyTerm>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpsetSummary {
    /// Number of distinct HPO terms observed in the proband.
    pub n_proband_terms: usize,
    /// Distinct proband terms that overlap at least one gene's annotations.
    pub n_proband_explained: usize,
    /// Distinct proband terms that overlap no gene (`n_proband_terms - n_proband_explained`).
    pub n_proband_unexplained: usize,
    /// Axis maximum for the columns: the largest `n_annotated_terms`.
    pub max_column: usize,
    /// Axis maximum for the rows: every bar has length `n_proband_terms`.
    pub max_row: usize,
}
