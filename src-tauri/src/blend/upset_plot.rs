//! Data for the pheno-upset plot and for the notebook page that documents it: how a
//! proband's HPO terms relate to the HPO annotations of individual genes and of combinations
//! of genes (multiple genetic diagnoses, MGD).
//!
//! The analysis has two steps:
//!
//! 1. [`plot_input`] expands every term to its true-path closure (the term plus all of its
//!    ancestors), for the proband and for each gene, and collects labels for display.
//! 2. [`plot_data`] derives [`UpsetPlotData`]. Every count in it is the length of a term list
//!    stored next to it, so the plot and the notebook page cannot disagree.
//!
//! # Matching rule
//!
//! A proband term `p` *overlaps* a set of annotations `S` if `p` is in `S` or is a descendant
//! of a term in `S` (true-path rule). With `proband_closures[p]` being `p` plus its ancestors,
//! this is `S ∩ proband_closures[p] ≠ ∅`. The reverse direction (the gene is annotated with a
//! term more specific than the proband's) is deliberately *not* a match.
use std::collections::{BTreeMap, BTreeSet};

use ontolius::{TermId, ontology::{TaxonomyWalk, csr::FullCsrOntology}};
use serde::{Deserialize, Serialize};

use crate::{blend::{disease_gene_entity::{GeneDiseaseEntity, simple_term}, dto::{GeneKey, MatchedTerm, ProbandSegment, ProbandTermMatch, TermKey, UpsetColumnData, UpsetPlotData, UpsetRowData, UpsetSummary}}, model::simple_term::SimpleOntologyTerm};


/// Most genes combined in one column. Multiple genetic diagnoses with more than four genes
/// are not expected, and the number of combinations grows quickly.
const MAX_GENES_PER_COLUMN: usize = 4;

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

/// HPO annotations of one gene, pooled across all diseases associated with it.
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct GeneTerms {
    /// Terms directly annotated to the gene's disease(s).
    pub annotated: BTreeSet<TermKey>,
    /// `annotated` plus all ancestors (true-path closure); always a superset of `annotated`.
    pub closure: BTreeSet<TermKey>,
}

impl GeneTerms {
    /// Records `t` as annotated and adds it, with its ancestors, to the closure.
    pub fn add_term(&mut self, hpo: &FullCsrOntology, t: &TermId) {
        self.add_term_with_closure(t.to_string(), closure_of(hpo, t));
    }

    /// Ontology-free core of [`GeneTerms::add_term`].
    fn add_term_with_closure(&mut self, term: TermKey, closure: BTreeSet<TermKey>) {
        self.annotated.insert(term);
        self.closure.extend(closure);
    }
}

/// The true-path closure of `t`: `t` itself plus all of its ancestors.
fn closure_of(hpo: &FullCsrOntology, t: &TermId) -> BTreeSet<TermKey> {
    let mut closure: BTreeSet<TermKey> = hpo
        .iter_term_and_ancestor_ids(t)
        .map(|a| a.to_string())
        .collect();
    closure.insert(t.to_string()); // harmless if the iterator already includes the term
    closure
}

/// Step 1: expands the proband's observed terms and the genes' annotations to closures.
///
/// Diseases that share a gene are pooled into one [`GeneTerms`].
pub fn plot_input(
    entities: &[GeneDiseaseEntity],
    observed: &[TermId],
    hpo: &FullCsrOntology,
) -> UpsetPlotInput {
    let mut terms: BTreeMap<TermKey, SimpleOntologyTerm> = BTreeMap::new();

    let mut proband_closures: BTreeMap<TermKey, BTreeSet<TermKey>> = BTreeMap::new();
    for p in observed {
        terms.entry(p.to_string()).or_insert_with(|| simple_term(hpo, p));
        proband_closures.insert(p.to_string(), closure_of(hpo, p));
    }

    let mut gene_annotations: BTreeMap<GeneKey, GeneTerms> = BTreeMap::new();
    for e in entities {
        let gt = gene_annotations.entry(e.gene_symbol.clone()).or_default();
        for t in &e.disease_hpo_ids {
            terms.entry(t.to_string()).or_insert_with(|| simple_term(hpo, t));
            gt.add_term(hpo, t);
        }
    }
    UpsetPlotInput { proband_closures, gene_annotations, terms }
}






type Closures = BTreeMap<TermKey, BTreeSet<TermKey>>;

/// Read-only data shared by the builders below.
struct Ctx<'a> {
    proband: &'a Closures,
    terms: &'a BTreeMap<TermKey, SimpleOntologyTerm>,
    /// Union of all proband closures: every term that is an observed proband term or one of
    /// its ancestors.
    proband_closure_union: BTreeSet<TermKey>,
}

impl Ctx<'_> {
    /// The display data for `key`.
    ///
    /// # Errors
    /// Fails if no label was collected for the term. [`plot_input`] collects one for every
    /// observed term and gene annotation, so this indicates an input assembled by hand or a
    /// bug in the caller; the error names the term.
    fn term(&self, key: &TermKey) -> Result<SimpleOntologyTerm, String> {
        self.terms
            .get(key)
            .cloned()
            .ok_or_else(|| format!("no display data (label) was collected for term {key}"))
    }

    fn terms_of<'k>(
        &self,
        keys: impl IntoIterator<Item = &'k TermKey>,
    ) -> Result<Vec<SimpleOntologyTerm>, String> {
        keys.into_iter().map(|k| self.term(k)).collect()
    }
}

/// The overlap of one gene with the proband's terms.
fn row_for(gene: &GeneKey, gt: &GeneTerms, ctx: &Ctx) -> Result<UpsetRowData, String> {
    let mut overlapping_terms = Vec::new();
    let mut non_overlapping_terms = Vec::new();
    for (p, closure) in ctx.proband {
        let annotated_by = ctx.terms_of(gt.annotated.intersection(closure))?;
        if annotated_by.is_empty() {
            non_overlapping_terms.push(ctx.term(p)?);
        } else {
            overlapping_terms.push(ProbandTermMatch { term: ctx.term(p)?, annotated_by });
        }
    }
    Ok(UpsetRowData {
        gene_symbol: gene.clone(),
        n_proband_terms: ctx.proband.len(),
        n_overlapping_terms: overlapping_terms.len(),
        overlapping_terms,
        non_overlapping_terms,
    })
}

/// Builds the column for the genes at `indices`, which must be non-empty.
fn column_for(
    indices: &[usize],
    genes: &[&GeneTerms],
    ctx: &Ctx,
) -> Result<UpsetColumnData, String> {
    let mut common = genes[indices[0]].closure.clone();
    let mut annotated_by_any: BTreeSet<TermKey> = BTreeSet::new();
    for &i in indices {
        common = common.intersection(&genes[i].closure).cloned().collect();
        annotated_by_any.extend(genes[i].annotated.iter().cloned());
    }
    common.retain(|t| annotated_by_any.contains(t));

    let mut matched_terms = Vec::new();
    let mut unmatched_terms = Vec::new();
    for t in &common {
        if ctx.proband_closure_union.contains(t) {
            let matched_by = ctx.terms_of(
                ctx.proband.iter().filter(|(_, c)| c.contains(t)).map(|(p, _)| p),
            )?;
            matched_terms.push(MatchedTerm { term: ctx.term(t)?, matched_by });
        } else {
            unmatched_terms.push(ctx.term(t)?);
        }
    }
    Ok(UpsetColumnData {
        row_indices: indices.to_vec(),
        n_annotated_terms: matched_terms.len() + unmatched_terms.len(),
        n_overlapping_terms: matched_terms.len(),
        matched_terms,
        unmatched_terms,
    })
}

/// All `k`-element subsets of `items`, each in ascending order.
fn combinations(items: &[usize], k: usize) -> Vec<Vec<usize>> {
    fn rec(items: &[usize], k: usize, start: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if cur.len() == k {
            out.push(cur.clone());
            return;
        }
        for i in start..items.len() {
            cur.push(items[i]);
            rec(items, k, i + 1, cur, out);
            cur.pop();
        }
    }
    let mut out = Vec::new();
    rec(items, k, 0, &mut Vec::new(), &mut out);
    out
}

/// Partitions the proband's terms by the exact set of genes that overlap each of them.
fn proband_breakdown(genes: &[&GeneTerms], ctx: &Ctx) -> Result<Vec<ProbandSegment>, String> {
    let mut by_signature: BTreeMap<Vec<usize>, Vec<TermKey>> = BTreeMap::new();
    for (p, closure) in ctx.proband {
        let signature: Vec<usize> = (0..genes.len())
            .filter(|&i| !genes[i].annotated.is_disjoint(closure))
            .collect();
        by_signature.entry(signature).or_default().push(p.clone());
    }

    // Ordered by number of genes, then lexicographically.
    let mut signatures: BTreeSet<(usize, Vec<usize>)> = BTreeSet::new();
    if genes.len() <= MAX_GENES_PER_COLUMN {
        let all: Vec<usize> = (0..genes.len()).collect();
        for k in 0..=genes.len() {
            signatures.extend(combinations(&all, k).into_iter().map(|s| (k, s)));
        }
    } else {
        signatures.insert((0, Vec::new()));
        signatures.extend(by_signature.keys().map(|s| (s.len(), s.clone())));
    }

    signatures
        .into_iter()
        .map(|(_, row_indices)| {
            let terms = ctx.terms_of(by_signature.get(&row_indices).into_iter().flatten())?;
            Ok(ProbandSegment { row_indices, n_terms: terms.len(), terms })
        })
        .collect()
}

/// Step 2: derives the rows, columns, breakdown and summary from the closures.
///
/// # Errors
/// Returns an error naming the term if a displayed term (an observed term or a gene
/// annotation) has no entry in `input.terms`.
pub fn plot_data(input: &UpsetPlotInput) -> Result<UpsetPlotData, String> {
    let ctx = Ctx {
        proband: &input.proband_closures,
        terms: &input.terms,
        proband_closure_union: input.proband_closures.values().flatten().cloned().collect(),
    };
    let n_proband_terms = ctx.proband.len();

    // BTreeMap iterates in key order, so rows are sorted by gene symbol.
    let rows: Vec<UpsetRowData> = input
        .gene_annotations
        .iter()
        .map(|(gene, gt)| row_for(gene, gt, &ctx))
        .collect::<Result<_, _>>()?;
    let genes: Vec<&GeneTerms> = input.gene_annotations.values().collect();

    let mut columns: Vec<UpsetColumnData> = (0..genes.len())
        .map(|i| column_for(&[i], &genes, &ctx))
        .collect::<Result<_, _>>()?;

    // Combine only genes that overlap the proband, to keep the number of combinations small.
    let candidates: Vec<usize> = (0..rows.len())
        .filter(|&i| rows[i].n_overlapping_terms > 0)
        .collect();
    for k in 2..=MAX_GENES_PER_COLUMN.min(candidates.len()) {
        for combo in combinations(&candidates, k) {
            let col = column_for(&combo, &genes, &ctx)?;
            if col.n_overlapping_terms > 0 {
                columns.push(col);
            }
        }
    }
    columns.sort_by(|a, b| {
        b.n_overlapping_terms
            .cmp(&a.n_overlapping_terms)
            .then(b.n_annotated_terms.cmp(&a.n_annotated_terms))
            .then(a.row_indices.cmp(&b.row_indices))
    });

    let proband_breakdown = proband_breakdown(&genes, &ctx)?;
    let n_proband_unexplained = proband_breakdown
        .iter()
        .find(|s| s.row_indices.is_empty())
        .map_or(0, |s| s.n_terms);

    let summary = UpsetSummary {
        n_proband_terms,
        n_proband_explained: n_proband_terms - n_proband_unexplained,
        n_proband_unexplained,
        max_column: columns.iter().map(|c| c.n_annotated_terms).max().unwrap_or(0),
        max_row: n_proband_terms,
    };

    Ok(UpsetPlotData { rows, columns, proband_breakdown, summary })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::hpo_fixture::hpo;
    use rstest::rstest;
    use std::{collections::HashSet, sync::Arc};

    // ---------------------------------------------------------------------------------
    // The real HPO is loaded once by the fixture. These tests use only the terms below;
    // `fixture_relationships_hold` checks the ontology facts they depend on, so a wrong
    // id fails there with a clear message and not somewhere downstream.
    // ---------------------------------------------------------------------------------

    const ROOT: &str = "HP:0000118"; // Phenotypic abnormality
    const SEIZURE: &str = "HP:0001250";
    const FOCAL: &str = "HP:0007359"; // Focal-onset seizure, a seizure type
    const TONIC: &str = "HP:0010818"; // Generalized tonic seizure, another seizure type
    const SHORT_STATURE: &str = "HP:0004322"; // unrelated to seizures
    const SCOLIOSIS: &str = "HP:0002650"; // unrelated to seizures and to short stature

    // ---------------------------------------------------------------------------------
    // Helpers
    // ---------------------------------------------------------------------------------

    fn tid(s: &str) -> TermId {
        s.parse().expect("valid CURIE")
    }

    fn tids(ts: &[&str]) -> Vec<TermId> {
        ts.iter().map(|t| tid(t)).collect()
    }

    /// One gene with the given annotated terms.
    ///
    /// Adapt this one function to however `GeneDiseaseEntity` is constructed in your code;
    /// the tests only need the gene symbol and the disease term ids.
    fn entity(gene: &str, ts: &[&str]) -> GeneDiseaseEntity {
        let m: HashSet<TermId> = tids(ts).into_iter().collect();
        GeneDiseaseEntity {
            gene_symbol: gene.to_string(),
            disease_hpo_ids: m,
            ncbi_gene_id: "123".to_string(),
            disease_list: Vec::default(),
            
        }
    }

    fn ids_of(ts: &[SimpleOntologyTerm]) -> Vec<String> {
        ts.iter().map(|t| t.term_id_as_string()).collect()
    }

    fn set(ts: &[&str]) -> BTreeSet<TermKey> {
        ts.iter().map(|t| t.to_string()).collect()
    }

    /// Runs both steps.
    fn run(hpo: &FullCsrOntology, entities: &[GeneDiseaseEntity], observed: &[&str]) -> UpsetPlotData {
        plot_data(&plot_input(entities, &tids(observed), hpo)).expect("valid input")
    }

    fn row<'a>(data: &'a UpsetPlotData, gene: &str) -> &'a UpsetRowData {
        data.rows.iter().find(|r| r.gene_symbol == gene).expect("row")
    }

    /// The column whose dots are exactly the given genes.
    fn col<'a>(data: &'a UpsetPlotData, genes: &[&str]) -> Option<&'a UpsetColumnData> {
        let idx: Vec<usize> = genes
            .iter()
            .map(|g| data.rows.iter().position(|r| r.gene_symbol == *g).expect("gene"))
            .collect();
        data.columns.iter().find(|c| c.row_indices == idx)
    }

    // ---------------------------------------------------------------------------------
    // The ontology facts the other tests rely on
    // ---------------------------------------------------------------------------------

    #[rstest]
    fn fixture_relationships_hold(hpo: Arc<FullCsrOntology>) {
        let closure = |t: &str| closure_of(&hpo, &tid(t));
        for t in [SEIZURE, FOCAL, TONIC, SHORT_STATURE, SCOLIOSIS] {
            assert!(closure(t).contains(ROOT), "{t} should be below {ROOT}");
        }
        assert!(closure(FOCAL).contains(SEIZURE), "FOCAL must be a kind of SEIZURE");
        assert!(closure(TONIC).contains(SEIZURE), "TONIC must be a kind of SEIZURE");
        assert!(!closure(SEIZURE).contains(FOCAL));
        assert!(!closure(SHORT_STATURE).contains(SEIZURE));
        assert!(!closure(SCOLIOSIS).contains(SHORT_STATURE));
        assert!(!closure(SHORT_STATURE).contains(SCOLIOSIS));
    }

    // ---------------------------------------------------------------------------------
    // closure_of and plot_input
    // ---------------------------------------------------------------------------------

    #[rstest]
    #[case(ROOT)]
    #[case(SEIZURE)]
    #[case(FOCAL)]
    #[case(SHORT_STATURE)]
    fn closure_contains_the_term_itself(hpo: Arc<FullCsrOntology>, #[case] term: &str) {
        assert!(closure_of(&hpo, &tid(term)).contains(term));
    }

    /// Regression: gene annotations used to get no label, which broke `plot_data`.
    #[rstest]
    fn plot_input_labels_observed_and_annotated_terms(hpo: Arc<FullCsrOntology>) {
        let input = plot_input(&[entity("A", &[FOCAL])], &tids(&[SHORT_STATURE]), &hpo);
        let keys: BTreeSet<&str> = input.terms.keys().map(String::as_str).collect();
        assert_eq!(keys, BTreeSet::from([FOCAL, SHORT_STATURE])); // ancestors are not labelled
        assert!(!input.terms[FOCAL].term_label.is_empty());
    }

    #[rstest]
    fn plot_input_collapses_duplicate_observed_terms(hpo: Arc<FullCsrOntology>) {
        let input = plot_input(&[], &tids(&[FOCAL, FOCAL, SHORT_STATURE]), &hpo);
        assert_eq!(input.proband_closures.len(), 2);
    }

    #[rstest]
    fn plot_input_pools_diseases_of_the_same_gene(hpo: Arc<FullCsrOntology>) {
        let entities = [entity("A", &[FOCAL]), entity("A", &[SHORT_STATURE]), entity("B", &[TONIC])];
        let input = plot_input(&entities, &[], &hpo);
        assert_eq!(input.gene_annotations["A"].annotated, set(&[FOCAL, SHORT_STATURE]));
        assert_eq!(input.gene_annotations["B"].annotated, set(&[TONIC]));
        assert_eq!(input.gene_annotations.len(), 2);
    }

    // ---------------------------------------------------------------------------------
    // Rows (left bars): proband-side counts and lists
    // ---------------------------------------------------------------------------------

    // Proband: focal seizure, tonic seizure, short stature.
    #[rstest]
    #[case(&[SEIZURE], 2)] // general gene term: both seizure types are descendants
    #[case(&[FOCAL], 1)] // exact match
    #[case(&[FOCAL, TONIC], 2)]
    #[case(&[SHORT_STATURE], 1)]
    #[case(&[ROOT], 3)] // everything is below the root phenotype term
    fn row_overlap_counts_descendants_of_annotations(
        hpo: Arc<FullCsrOntology>,
        #[case] annotated: &[&str],
        #[case] expected: usize,
    ) {
        let data = run(&hpo, &[entity("A", annotated)], &[FOCAL, TONIC, SHORT_STATURE]);
        assert_eq!(row(&data, "A").n_proband_terms, 3);
        assert_eq!(row(&data, "A").n_overlapping_terms, expected);
    }

    #[rstest]
    fn gene_more_specific_than_proband_is_not_a_match(hpo: Arc<FullCsrOntology>) {
        // The gene has the specific term; the proband only the general one.
        let data = run(&hpo, &[entity("A", &[FOCAL])], &[SEIZURE]);
        assert_eq!(row(&data, "A").n_overlapping_terms, 0);
        assert_eq!(data.summary.n_proband_explained, 0);
    }

    #[rstest]
    fn row_lists_overlapping_terms_with_annotating_terms(hpo: Arc<FullCsrOntology>) {
        let data = run(
            &hpo,
            &[entity("A", &[SEIZURE, SHORT_STATURE])],
            &[FOCAL, TONIC, SCOLIOSIS],
        );
        let r = row(&data, "A");
        let overlapping: Vec<(String, Vec<String>)> = r
            .overlapping_terms
            .iter()
            .map(|m| (m.term.term_id_as_string(), ids_of(&m.annotated_by)))
            .collect();
        assert_eq!(
            overlapping,
            [
                (FOCAL.to_string(), vec![SEIZURE.to_string()]),
                (TONIC.to_string(), vec![SEIZURE.to_string()]),
            ]
        );
        assert_eq!(ids_of(&r.non_overlapping_terms), [SCOLIOSIS]);
        assert_eq!(r.n_overlapping_terms, 2);
    }

    // ---------------------------------------------------------------------------------
    // Columns (top bars): gene-side counts and lists
    // ---------------------------------------------------------------------------------

    #[rstest]
    fn single_gene_column_lists_matched_and_unmatched_terms(hpo: Arc<FullCsrOntology>) {
        let data = run(&hpo, &[entity("A", &[SEIZURE, SHORT_STATURE])], &[FOCAL, TONIC]);
        let c = col(&data, &["A"]).expect("column");
        assert_eq!((c.n_annotated_terms, c.n_overlapping_terms), (2, 1));
        assert_eq!(c.matched_terms.len(), 1);
        assert_eq!(c.matched_terms[0].term.term_id_as_string(), SEIZURE);
        assert_eq!(ids_of(&c.matched_terms[0].matched_by), [FOCAL, TONIC]); // counted once
        assert_eq!(ids_of(&c.unmatched_terms), [SHORT_STATURE]);
    }

    #[rstest]
    fn column_does_not_match_terms_more_specific_than_the_proband(hpo: Arc<FullCsrOntology>) {
        let data = run(&hpo, &[entity("A", &[FOCAL])], &[SEIZURE]);
        let c = col(&data, &["A"]).expect("column");
        assert_eq!((c.n_annotated_terms, c.n_overlapping_terms), (1, 0));
        assert_eq!(ids_of(&c.unmatched_terms), [FOCAL]);
    }

    #[rstest]
    fn pair_column_counts_the_more_general_term_once(hpo: Arc<FullCsrOntology>) {
        // A has the specific term, B the general one: the shared annotation is the general one.
        let entities = [entity("A", &[FOCAL]), entity("B", &[SEIZURE])];
        let data = run(&hpo, &entities, &[FOCAL, TONIC, SHORT_STATURE]);
        let c = col(&data, &["A", "B"]).expect("pair column");
        assert_eq!((c.n_annotated_terms, c.n_overlapping_terms), (1, 1));
        assert_eq!(c.matched_terms[0].term.term_id_as_string(), SEIZURE);
    }

    #[rstest]
    fn pair_column_ignores_unannotated_shared_ancestors(hpo: Arc<FullCsrOntology>) {
        // Two different seizure types share many ancestors, but no gene annotates them.
        let entities = [entity("A", &[FOCAL]), entity("B", &[TONIC])];
        let data = run(&hpo, &entities, &[FOCAL, TONIC]);
        assert!(col(&data, &["A", "B"]).is_none());
        assert!(col(&data, &["A"]).is_some());
        assert!(col(&data, &["B"]).is_some());
    }

    #[rstest]
    fn pair_column_is_inclusive_of_other_genes(hpo: Arc<FullCsrOntology>) {
        let entities = [entity("A", &[FOCAL]), entity("B", &[FOCAL]), entity("C", &[FOCAL])];
        let data = run(&hpo, &entities, &[FOCAL]);
        assert_eq!(col(&data, &["A", "B"]).expect("A+B").n_annotated_terms, 1);
        assert_eq!(col(&data, &["A", "B", "C"]).expect("A+B+C").n_annotated_terms, 1);
    }

    #[rstest]
    fn column_overlap_never_exceeds_annotated_terms(hpo: Arc<FullCsrOntology>) {
        // One shared annotation (SEIZURE) with two proband terms below it.
        let entities = [entity("A", &[SEIZURE]), entity("B", &[SEIZURE])];
        let data = run(&hpo, &entities, &[FOCAL, TONIC]);
        let c = col(&data, &["A", "B"]).expect("pair column");
        assert_eq!((c.n_annotated_terms, c.n_overlapping_terms), (1, 1));
        assert!(data.columns.iter().all(|c| c.n_overlapping_terms <= c.n_annotated_terms));
    }

    #[rstest]
    fn genes_without_overlap_are_not_combined(hpo: Arc<FullCsrOntology>) {
        let entities = [entity("A", &[FOCAL]), entity("B", &[FOCAL]), entity("C", &[SHORT_STATURE])];
        let data = run(&hpo, &entities, &[FOCAL]);
        assert!(col(&data, &["A", "B"]).is_some());
        assert!(col(&data, &["A", "C"]).is_none());
        assert!(col(&data, &["C"]).is_some()); // single-gene columns are always present
    }

    #[rstest]
    fn at_most_four_genes_are_combined(hpo: Arc<FullCsrOntology>) {
        let entities: Vec<_> = ["A", "B", "C", "D", "E"].iter().map(|g| entity(g, &[FOCAL])).collect();
        let data = run(&hpo, &entities, &[FOCAL]);
        // C(5,1) + C(5,2) + C(5,3) + C(5,4) = 30; the 5-gene combination is not built.
        assert_eq!(data.columns.len(), 30);
        assert!(data.columns.iter().all(|c| c.row_indices.len() <= 4));
    }

    // ---------------------------------------------------------------------------------
    // Breakdown (bottom bar) and summary
    // ---------------------------------------------------------------------------------

    #[rstest]
    fn breakdown_partitions_terms_by_exact_gene_set(hpo: Arc<FullCsrOntology>) {
        // FOCAL is overlapped by A (via SEIZURE) and B; TONIC only by A; SHORT_STATURE by nobody.
        let entities = [entity("A", &[SEIZURE]), entity("B", &[FOCAL])];
        let data = run(&hpo, &entities, &[FOCAL, TONIC, SHORT_STATURE]);
        let seen: Vec<(Vec<usize>, Vec<String>)> = data
            .proband_breakdown
            .iter()
            .map(|s| (s.row_indices.clone(), ids_of(&s.terms)))
            .collect();
        assert_eq!(
            seen,
            [
                (vec![], vec![SHORT_STATURE.to_string()]),
                (vec![0], vec![TONIC.to_string()]),
                (vec![1], vec![]),
                (vec![0, 1], vec![FOCAL.to_string()]),
            ]
        );
        assert!(data.proband_breakdown.iter().all(|s| s.n_terms == s.terms.len()));
    }

    #[rstest]
    fn summary_counts_explained_and_unexplained_terms(hpo: Arc<FullCsrOntology>) {
        let entities = [entity("A", &[FOCAL]), entity("B", &[SHORT_STATURE])];
        let data = run(&hpo, &entities, &[FOCAL, TONIC, SHORT_STATURE]);
        let s = &data.summary;
        assert_eq!(s.n_proband_terms, 3);
        assert_eq!(s.n_proband_explained, 2); // FOCAL via A, SHORT_STATURE via B; TONIC by nobody
        assert_eq!(s.n_proband_unexplained, 1);
        assert_eq!(s.max_row, 3);
        assert_eq!(s.max_column, data.columns.iter().map(|c| c.n_annotated_terms).max().unwrap());
    }

    #[rstest]
    fn no_genes_gives_empty_plot(hpo: Arc<FullCsrOntology>) {
        let data = run(&hpo, &[], &[FOCAL, TONIC]);
        assert!(data.rows.is_empty() && data.columns.is_empty());
        assert_eq!(data.summary.n_proband_unexplained, 2);
    }

    #[rstest]
    fn no_proband_terms_gives_zero_overlaps(hpo: Arc<FullCsrOntology>) {
        let data = run(&hpo, &[entity("A", &[FOCAL])], &[]);
        assert_eq!(row(&data, "A").n_proband_terms, 0);
        assert_eq!(row(&data, "A").n_overlapping_terms, 0);
    }

    // ---------------------------------------------------------------------------------
    // Errors, not panics
    // ---------------------------------------------------------------------------------

    /// A missing label is an error that names the term.
    #[rstest]
    #[case(SEIZURE)] // a gene annotation without a label
    #[case(SHORT_STATURE)] // an observed term without a label
    fn plot_data_reports_the_term_with_a_missing_label(
        hpo: Arc<FullCsrOntology>,
        #[case] unlabelled: &str,
    ) {
        let mut input = plot_input(&[entity("A", &[SEIZURE])], &tids(&[SHORT_STATURE]), &hpo);
        input.terms.remove(unlabelled);
        let err = plot_data(&input).unwrap_err();
        assert!(err.contains(unlabelled), "{err}");
    }

    #[rstest]
    fn every_count_equals_the_length_of_its_list(hpo: Arc<FullCsrOntology>) {
        let entities = [
            entity("A", &[SEIZURE, SHORT_STATURE]),
            entity("B", &[FOCAL]),
            entity("C", &[TONIC, SCOLIOSIS]),
        ];
        let data = run(&hpo, &entities, &[FOCAL, TONIC, SHORT_STATURE]);
        for r in &data.rows {
            assert_eq!(r.n_overlapping_terms, r.overlapping_terms.len());
            assert_eq!(r.n_proband_terms, r.overlapping_terms.len() + r.non_overlapping_terms.len());
        }
        for c in &data.columns {
            assert_eq!(c.n_overlapping_terms, c.matched_terms.len());
            assert_eq!(c.n_annotated_terms, c.matched_terms.len() + c.unmatched_terms.len());
            assert!(c.matched_terms.iter().all(|m| !m.matched_by.is_empty()));
        }
        let s = &data.summary;
        assert_eq!(s.n_proband_explained + s.n_proband_unexplained, s.n_proband_terms);
    }
}
