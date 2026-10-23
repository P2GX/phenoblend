//! Data for the pheno-upset plot: how a proband's HPO terms relate to the HPO annotations
//! of individual genes and of combinations of genes (multiple genetic diagnoses, MGD).
//!
//! The analysis has two steps:
//!
//! 1. [`plot_input`] expands every term to its true-path closure (the term plus all of its
//!    ancestors), for the proband and for each gene.
//! 2. [`plot_data`] derives the rows, columns and summary that the frontend draws.
//!
//! # Matching rule
//!
//! A proband term `p` *overlaps* a set of annotations `S` if `p` is in `S` or is a descendant
//! of a term in `S` (true-path rule). With `proband_closures[p]` being `p` plus its ancestors,
//! this is `S ∩ proband_closures[p] ≠ ∅`. The reverse direction (the gene is annotated with a
//! term more specific than the proband's) is deliberately *not* a match.
use std::collections::{BTreeMap, BTreeSet};

use ontolius::{TermId, ontology:: csr::FullCsrOntology};


use crate::blend::{disease_gene_entity::GeneDiseaseEntity, dto::{GeneKey, GeneTerms, ProbandSegment, ROOT_PHENOTYPE_TERM_ID, TermKey, UpsetColumnData, UpsetPlotData, UpsetPlotInput, UpsetRowData, UpsetSummary, closure_of}};

/// Most genes combined in one column. Multiple genetic diagnoses with more than four genes
/// are not expected, and the number of combinations grows quickly.
const MAX_GENES_PER_COLUMN: usize = 4;

    

pub fn plot_input(
    entities: &[GeneDiseaseEntity],
    observed: &[TermId],
    hpo: &FullCsrOntology,
) -> UpsetPlotInput {
    let proband_closures: BTreeMap<TermKey, BTreeSet<TermKey>> = observed
        .iter()
        .map(|p| (p.to_string(), closure_of(hpo, p)))
        .collect();
    let mut gene_annotations: BTreeMap<GeneKey, GeneTerms> = BTreeMap::new();
    for e in entities {
        let gt = gene_annotations.entry(e.gene_symbol.clone()).or_default();
        for t in &e.disease_hpo_ids {
            gt.add_term(hpo, t);
        }
    }
   UpsetPlotInput { proband_closures, gene_annotations }
}



type Closures = BTreeMap<TermKey, BTreeSet<TermKey>>;

/// Number of proband terms p for which `terms ∩ proband_closures[p] ≠ ∅`,
/// i.e. p is one of `terms` or a descendant of one of them (true-path rule).
fn count_overlapping(proband: &Closures, terms: &BTreeSet<TermKey>) -> usize {
    proband.values().filter(|c| !terms.is_disjoint(c)).count()
}

/// Terms shared by all the given genes under the true-path rule, restricted to terms
/// that are explicitly annotated to at least one of them.
fn column_for(
    indices: &[usize], 
    genes: &[&GeneTerms], 
    proband_term_closure: &BTreeSet<TermKey>) -> UpsetColumnData {
    let mut common: BTreeSet<TermKey> = genes[indices[0]].closure.clone();
    let mut annotated_by_any: BTreeSet<TermKey> = BTreeSet::new();
    for &i in indices {
        common = common.intersection(&genes[i].closure).cloned().collect();
        annotated_by_any.extend(genes[i].annotated.iter().cloned());
    }
    common.retain(|t| annotated_by_any.contains(t));

    UpsetColumnData {
        row_indices: indices.to_vec(),
        n_annotated_terms: common.len(),
        n_overlapping_terms: common.intersection(proband_term_closure).count(), // TODO CHECK
    }
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
fn proband_breakdown(proband: &Closures, genes: &[&GeneTerms]) -> Vec<ProbandSegment> {
    let mut counts: BTreeMap<Vec<usize>, usize> = BTreeMap::new();
    for closure in proband.values() {
        let signature: Vec<usize> = (0..genes.len())
            .filter(|&i| !genes[i].annotated.is_disjoint(closure))
            .collect();
        *counts.entry(signature).or_default() += 1;
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
        signatures.extend(counts.keys().map(|s| (s.len(), s.clone())));
    }

    signatures
        .into_iter()
        .map(|(_, row_indices)| {
            let n_terms = counts.get(&row_indices).copied().unwrap_or(0);
            ProbandSegment { row_indices, n_terms }
        })
        .collect()
}

/// Step 2: derives the rows, columns and summary of the plot from the closures.
pub fn plot_data(input: &UpsetPlotInput) -> UpsetPlotData {
    let proband = &input.proband_closures;
    let n_proband_terms = proband.len();

    // BTreeMap iterates in key order, so rows are sorted by gene symbol.
    let rows: Vec<UpsetRowData> = input
        .gene_annotations
        .iter()
        .map(|(gene, gt)| UpsetRowData {
            gene_symbol: gene.clone(),
            n_proband_terms,
            n_overlapping_terms: count_overlapping(proband, &gt.annotated),
        })
        .collect();
    let genes: Vec<&GeneTerms> = input.gene_annotations.values().collect();

    let proband_term_closure: BTreeSet<TermKey> = proband.values().flatten().cloned().collect();
    let mut columns: Vec<UpsetColumnData> = (0..genes.len())
        .map(|i| column_for(&[i], &genes, &proband_term_closure))
        .collect();

    // Combine only genes that overlap the proband, to keep the number of combinations small.
    let candidates: Vec<usize> = (0..rows.len())
        .filter(|&i| rows[i].n_overlapping_terms > 0)
        .collect();
    for k in 2..=MAX_GENES_PER_COLUMN.min(candidates.len()) {
        for combo in combinations(&candidates, k) {
            let col = column_for(&combo, &genes, &proband_term_closure);
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

    let n_proband_explained = proband
        .values()
        .filter(|c| genes.iter().any(|g| !g.annotated.is_disjoint(c)))
        .count();

    let summary = UpsetSummary {
        n_proband_terms,
        n_proband_explained,
        n_proband_unexplained: n_proband_terms - n_proband_explained,
        max_column: columns.iter().map(|c| c.n_annotated_terms).max().unwrap_or(0),
        max_row: n_proband_terms,
    };
    let proband_breakdown = proband_breakdown(proband, &genes);

    UpsetPlotData { rows, columns, proband_breakdown, summary }
}




#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use crate::tests::hpo_fixture::hpo;


}