use std::collections::{HashMap, BTreeSet};
use std::sync::Arc;
use itertools::Itertools;

use ontolius::TermId;
use ontolius::ontology::{TaxonomyQuery, TaxonomyWalk, OntologyTerms};
use ontolius::ontology::csr::FullCsrOntology;
use ontolius::term::MinimalTerm;

use crate::blend::disease_gene_entity::GeneDiseaseEntity;
use crate::blend::dto::{SpreadPlotCategory, SpreadPlotPayload};
use crate::hpoa::disease_model::GeneDiseaseAssociation;
use crate::model::proband::Proband;

/// Converts raw gene->GDA mappings into clean entities, alongside a parallel
/// list of gene symbols in the same (HashMap-iteration) order.
fn build_entities(
    gda_map: &HashMap<String, Vec<GeneDiseaseAssociation>>,
) -> Result<(Vec<GeneDiseaseEntity>, Vec<String>), String> {
    let mut gd_entry_list = Vec::new();
    let mut genotypes = Vec::new();

    for (symbol, gda_list) in gda_map.iter() {
        let gd_entity = GeneDiseaseEntity::new(gda_list)?;
        gd_entry_list.push(gd_entity);
        genotypes.push(symbol.clone());
    }

    Ok((gd_entry_list, genotypes))
}

/// Fetches the top-level HPO organ-system categories (direct children of HP:0000118).
fn top_level_categories(hpo: &FullCsrOntology) -> Result<Vec<TermId>, String> {
    let phenotype_root: TermId = "HP:0000118".parse().map_err(|_| "Could not map root term".to_string())?;
    Ok(hpo.iter_child_ids(&phenotype_root).cloned().collect())
}

/// Computes every non-empty subset (of every size) of the given genotype symbols.
fn genotype_power_set(genotypes: &[String]) -> Vec<Vec<String>> {
    let mut power_set = Vec::new();
    for size in 1..=genotypes.len() {
        for combo in genotypes.iter().cloned().combinations(size) {
            power_set.push(combo);
        }
    }
    power_set
}

/// Computes, for each top-level category, the fraction of the proband's observed
/// phenotypes that fall under that category (self or descendant).
fn ppkt_category_fractions(
    observed_hpos: &[TermId],
    categories: &[TermId],
    hpo: &FullCsrOntology,
) -> HashMap<TermId, f64> {
    let mut counts: HashMap<TermId, f64> = categories.iter().map(|c| (c.clone(), 0.0)).collect();

    for pheno in observed_hpos {
        for cat in categories {
            if cat == pheno || hpo.is_ancestor_of(cat, pheno) {
                *counts.entry(cat.clone()).or_insert(0.0) += 1.0;
            }
        }
    }

    let total = observed_hpos.len() as f64;
    for val in counts.values_mut() {
        *val /= total;
    }

    counts
}

/// Computes the category fractions for a single genotype combination: the
/// fraction of the combination's union of HPO terms falling under each category.
/// Returns NaN for every category if the combination's term union is empty.
fn combo_category_fractions(
    combo: &[String],
    entity_lookup: &HashMap<String, &GeneDiseaseEntity>,
    categories: &[TermId],
    hpo: &FullCsrOntology,
) -> HashMap<TermId, f64> {
    let mut union_hpos: BTreeSet<TermId> = BTreeSet::new();
    for gene in combo {
        if let Some(entity) = entity_lookup.get(gene) {
            for term in &entity.disease_hpo_ids {
                union_hpos.insert(term.clone());
            }
        }
    }

    let mut cat_scores: HashMap<TermId, f64> = categories.iter().map(|c| (c.clone(), 0.0)).collect();
    let denominator = union_hpos.len() as f64;

    if denominator > 0.0 {
        for union_hpo in &union_hpos {
            for cat in categories {
                if hpo.is_ancestor_of(cat, union_hpo) || cat == union_hpo {
                    *cat_scores.entry(cat.clone()).or_insert(0.0) += 1.0;
                }
            }
        }
        for val in cat_scores.values_mut() {
            *val /= denominator;
        }
    } else {
        for val in cat_scores.values_mut() {
            *val = f64::NAN;
        }
    }

    cat_scores
}

/// Computes category fractions for every combination in the power set, keyed
/// by the combination's canonical "+"-joined label.
fn build_combo_column_data(
    genes_power_set: &[Vec<String>],
    entity_lookup: &HashMap<String, &GeneDiseaseEntity>,
    categories: &[TermId],
    hpo: &FullCsrOntology,
) -> (HashMap<String, HashMap<TermId, f64>>, Vec<String>) {
    let mut combo_column_data = HashMap::new();
    let mut labels_power_set = Vec::new();

    for combo in genes_power_set {
        let label = combo.join("+");
        labels_power_set.push(label.clone());

        let cat_scores = combo_category_fractions(combo, entity_lookup, categories, hpo);
        combo_column_data.insert(label, cat_scores);
    }

    (combo_column_data, labels_power_set)
}

/// Scores each individual gene by the dot product of its category fractions
/// with the patient's category fractions, then returns gene symbols ranked
/// descending by that similarity score.
fn rank_genes_by_similarity(
    genotypes: &[String],
    combo_column_data: &HashMap<String, HashMap<TermId, f64>>,
    ppkt_counts: &HashMap<TermId, f64>,
    categories: &[TermId],
) -> Vec<String> {
    let mut single_gene_scores: Vec<(String, f64)> = Vec::new();

    for gene in genotypes {
        let mut score = 0.0;
        if let Some(cat_scores) = combo_column_data.get(gene) {
            for cat in categories {
                let ppkt_val = ppkt_counts.get(cat).cloned().unwrap_or(0.0);
                let gene_val = cat_scores.get(cat).cloned().unwrap_or(0.0);
                if !gene_val.is_nan() {
                    score += gene_val * ppkt_val;
                }
            }
        }
        single_gene_scores.push((gene.clone(), score));
    }

    single_gene_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    single_gene_scores.into_iter().map(|(g, _)| g).collect()
}

/// Orders every power-set combination label according to the single-gene
/// similarity ranking, with each combo's internal gene order kept canonical
/// (matching original input index positions).
fn order_combination_labels(
    ranked_genes: &[String],
    gene_pos_map: &HashMap<String, usize>,
) -> Vec<String> {
    let mut ordered_labels = Vec::new();

    for size in 1..=ranked_genes.len() {
        for mut combo in ranked_genes.iter().cloned().combinations(size) {
            combo.sort_by_key(|g| gene_pos_map.get(g).cloned().unwrap_or(0));
            ordered_labels.push(combo.join("+"));
        }
    }

    ordered_labels
}

/// Sorts categories descending by their patient-profile ("Ppkt") fraction.
fn sort_categories_by_ppkt_value(
    categories: &[TermId],
    ppkt_counts: &HashMap<TermId, f64>,
) -> Vec<TermId> {
    let mut sorted = categories.to_vec();
    sorted.sort_by(|a, b| {
        let a_val = ppkt_counts.get(a).unwrap_or(&0.0);
        let b_val = ppkt_counts.get(b).unwrap_or(&0.0);
        b_val.partial_cmp(a_val).unwrap_or(std::cmp::Ordering::Equal)
    });
    sorted
}

/// Builds the final per-category payload rows, in the given sorted category order.
fn build_final_categories(
    sorted_categories: &[TermId],
    ordered_labels: &[String],
    ppkt_counts: &HashMap<TermId, f64>,
    combo_column_data: &HashMap<String, HashMap<TermId, f64>>,
    hpo: &FullCsrOntology,
) -> Vec<SpreadPlotCategory> {
    let mut final_categories = Vec::with_capacity(sorted_categories.len());

    for cat_id in sorted_categories {
        let ppkt_val = ppkt_counts.get(cat_id).cloned().unwrap_or(0.0);

        let gene_values: Vec<f64> = ordered_labels
            .iter()
            .map(|label| {
                combo_column_data
                    .get(label)
                    .and_then(|m| m.get(cat_id))
                    .cloned()
                    .unwrap_or(f64::NAN)
            })
            .collect();

        let hpo_name = hpo.term_by_id(cat_id)
            .map(|term| term.name().to_string())
            .unwrap_or_else(|| cat_id.to_string());

        // Skip organs with no observations in our data
        let epsilon = 1e-8;
        let max_gene_val = gene_values.iter().copied().reduce(f64::max).unwrap_or(0.0);
        if ppkt_val < epsilon && max_gene_val < epsilon {
            continue;
        }

        final_categories.push(SpreadPlotCategory {
            id: cat_id.to_string(),
            name: hpo_name,
            alias: None,
            ppkt_value: ppkt_val,
            gene_values,
        });
    }

    final_categories
}





pub fn get_spread_plot_payload(
    proband: Proband,
    genotypes: &[String],
    gd_entities: &[GeneDiseaseEntity],
    hpo: &FullCsrOntology
) -> Result<SpreadPlotPayload, String> {
    //let (gd_entry_list, genotypes) = build_entities(gda_map)?;

    if proband.observed_hpos.is_empty() || genotypes.is_empty() {
        return Err("Phenotypes and genotypes must not be empty.".to_string());
    }

    let categories = top_level_categories(&hpo)?;
    let genes_power_set = genotype_power_set(&genotypes);

    let gene_pos_map: HashMap<String, usize> = genotypes.iter()
        .enumerate()
        .map(|(i, g)| (g.clone(), i))
        .collect();

    let ppkt_counts = ppkt_category_fractions(&proband.observed_hpos, &categories, &hpo);

    let entity_lookup: HashMap<String, &GeneDiseaseEntity> = gd_entities.iter()
        .map(|gde| (gde.gene_symbol.clone(), gde))
        .collect();

    let (combo_column_data, _labels_power_set) =
        build_combo_column_data(&genes_power_set, &entity_lookup, &categories, &hpo);

    let ranked_genes = rank_genes_by_similarity(&genotypes, &combo_column_data, &ppkt_counts, &categories);
    let ordered_labels = order_combination_labels(&ranked_genes, &gene_pos_map);
    let sorted_categories = sort_categories_by_ppkt_value(&categories, &ppkt_counts);

    let final_categories = build_final_categories(
        &sorted_categories,
        &ordered_labels,
        &ppkt_counts,
        &combo_column_data,
        &hpo,
    );

    let mut series_labels = vec!["Ppkt".to_string()];
    series_labels.extend(ordered_labels);

    Ok(SpreadPlotPayload {
        series_labels,
        categories: final_categories,
    })
}



#[cfg(test)]
mod tests {
    use super::*;
    use rstest::{fixture, rstest};
    use std::collections::HashMap;

    // ---------------------------------------------------------------------
    // genotype_power_set
    // ---------------------------------------------------------------------

    #[fixture]
    fn three_genes() -> Vec<String> {
        vec!["GENE_A".to_string(), "GENE_B".to_string(), "GENE_C".to_string()]
    }

    #[rstest]
    fn power_set_of_empty_input_is_empty(
    ) {
        let result = genotype_power_set(&[]);
        assert!(result.is_empty());
    }

    #[rstest]
    fn power_set_of_single_gene_yields_one_combo() {
        let genes = vec!["GENE_A".to_string()];
        let result = genotype_power_set(&genes);
        assert_eq!(result, vec![vec!["GENE_A".to_string()]]);
    }

    #[rstest]
    fn power_set_has_correct_total_count(three_genes: Vec<String>) {
        // For n=3: C(3,1) + C(3,2) + C(3,3) = 3 + 3 + 1 = 7
        let result = genotype_power_set(&three_genes);
        assert_eq!(result.len(), 7);
    }

    #[rstest]
    fn power_set_contains_full_combination(three_genes: Vec<String>) {
        let result = genotype_power_set(&three_genes);
        assert!(result.contains(&three_genes));
    }

    #[rstest]
    fn power_set_contains_each_single_gene_alone(three_genes: Vec<String>) {
        let result = genotype_power_set(&three_genes);
        for gene in &three_genes {
            assert!(result.contains(&vec![gene.clone()]));
        }
    }

    #[rstest]
    fn power_set_contains_no_duplicate_combinations(three_genes: Vec<String>) {
        let result = genotype_power_set(&three_genes);
        let unique: std::collections::HashSet<_> = result.iter().cloned().collect();
        assert_eq!(unique.len(), result.len());
    }

    // ---------------------------------------------------------------------
    // order_combination_labels
    // ---------------------------------------------------------------------

    #[fixture]
    fn gene_pos_map() -> HashMap<String, usize> {
        // Original input order: GENE_A=0, GENE_B=1, GENE_C=2
        [
            ("GENE_A".to_string(), 0),
            ("GENE_B".to_string(), 1),
            ("GENE_C".to_string(), 2),
        ]
        .into_iter()
        .collect()
    }

    #[rstest]
    fn single_gene_ranking_yields_single_labels(gene_pos_map: HashMap<String, usize>) {
        let ranked = vec!["GENE_A".to_string()];
        let result = order_combination_labels(&ranked, &gene_pos_map);
        assert_eq!(result, vec!["GENE_A".to_string()]);
    }

    #[rstest]
    fn combo_labels_use_canonical_original_order_not_rank_order(
        gene_pos_map: HashMap<String, usize>,
    ) {
        // Ranked order puts GENE_C first (e.g. it scored highest), but the
        // combo label should still reflect original input order (A, B, C),
        // not the ranking order.
        let ranked = vec!["GENE_C".to_string(), "GENE_A".to_string(), "GENE_B".to_string()];
        let result = order_combination_labels(&ranked, &gene_pos_map);

        // The full 3-gene combination must be labeled in original order.
        assert!(result.contains(&"GENE_A+GENE_B+GENE_C".to_string()));
        // It must NOT appear in ranked order.
        assert!(!result.contains(&"GENE_C+GENE_A+GENE_B".to_string()));
    }

    #[rstest]
    fn total_label_count_matches_power_set_size(gene_pos_map: HashMap<String, usize>) {
        let ranked = vec!["GENE_A".to_string(), "GENE_B".to_string(), "GENE_C".to_string()];
        let result = order_combination_labels(&ranked, &gene_pos_map);
        // Same combinatorial count as genotype_power_set: 7 for n=3
        assert_eq!(result.len(), 7);
    }

    #[rstest]
    fn unknown_gene_defaults_to_position_zero(gene_pos_map: HashMap<String, usize>) {
        // A gene missing from gene_pos_map should sort as if at position 0,
        // per the `.unwrap_or(0)` fallback in the implementation.
        let ranked = vec!["GENE_B".to_string(), "UNKNOWN_GENE".to_string()];
        let result = order_combination_labels(&ranked, &gene_pos_map);

        // UNKNOWN_GENE (falls back to pos 0) should sort before GENE_B (pos 1)
        assert!(result.contains(&"UNKNOWN_GENE+GENE_B".to_string()));
    }

    // ---------------------------------------------------------------------
    // build_entities
    // ---------------------------------------------------------------------
    // NOTE: placeholder — I don't know GeneDiseaseAssociation's fields or
    // GeneDiseaseEntity::new's exact validation behavior (when does it
    // return Err?). Fill in `sample_gda(...)` once you share the real shape.

    // #[fixture]
    // fn sample_gda_map() -> HashMap<String, Vec<GeneDiseaseAssociation>> {
    //     let mut map = HashMap::new();
    //     map.insert("GENE_A".to_string(), vec![/* GeneDiseaseAssociation { ... } */]);
    //     map
    // }
    //
    // #[rstest]
    // fn build_entities_preserves_gene_symbols(sample_gda_map: HashMap<String, Vec<GeneDiseaseAssociation>>) {
    //     let (entities, genotypes) = build_entities(&sample_gda_map).unwrap();
    //     assert_eq!(entities.len(), genotypes.len());
    //     assert!(genotypes.contains(&"GENE_A".to_string()));
    // }
    //
    // #[rstest]
    // fn build_entities_propagates_entity_construction_error() {
    //     // Construct a GeneDiseaseAssociation that should cause
    //     // GeneDiseaseEntity::new(...) to return Err, and assert build_entities
    //     // surfaces that Err rather than panicking or silently dropping it.
    // }
}