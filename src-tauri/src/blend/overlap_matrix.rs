use std::collections::HashMap;
use std::sync::Arc;

use ontolius::{TermId, ontology::csr::FullCsrOntology};

use crate::blend::disease_gene_entity::GeneDiseaseEntity;
use crate::blend::dto::OverlapPlotItem;
use crate::blend::dto::OverlapPlotPayload;
use crate::hpoa::disease_model::GeneDiseaseAssociation;
use crate::model::proband::Proband;



/// Builds one [`GeneDiseaseEntity`] per gene from the gene-symbol-keyed association map.
///
/// Each value in `annotation_map` (all associations for a single gene) is passed to
/// [`GeneDiseaseEntity::new`]. The map keys are used only to order the output and to
/// label errors; the entities take their gene information from the associations.
///
/// The result is sorted by map key (gene symbol), so repeated calls on the same input
/// give the same order. This matters because the list ends up in the notebook that is
/// re-emitted to the frontend, and `HashMap` iteration order would otherwise vary
/// between runs.
///
/// # Errors
///
/// Returns the first error produced by [`GeneDiseaseEntity::new`], prefixed with the
/// gene symbol it occurred for. No partial list is returned.
pub fn build_entity_list(
    annotation_map: &HashMap<String, Vec<GeneDiseaseAssociation>>,
) -> Result<Vec<GeneDiseaseEntity>, String> {
    // transform the k+v pairs of the map into entries of a Vec
    let mut genes: Vec<_> = annotation_map.iter().collect();
    // sort by first entry, i.e., gene symbols
    genes.sort_unstable_by(|a, b| a.0.cmp(b.0));
    // GeneDiseaseEntity::new produces a Result<GeneDiseaseEntity,String>
    // given the return type, collect creates Result<Vec<GeneDiseaseEntity>, String>
    genes
        .into_iter()
        .map(|(gene, gdas)| {
            GeneDiseaseEntity::new(gdas).map_err(|e| format!("gene '{gene}': {e}"))
        })
        .collect()
}

pub fn calculate_overlap_matrix(
    hpo: Arc<FullCsrOntology>,
    annotation_map: &HashMap<String, Vec<GeneDiseaseAssociation>>,
    disease_counts: &HashMap<TermId, usize>,
    proband: Proband,
) -> Result<OverlapPlotPayload, String> {
    let gd_entry_list = build_entity_list(annotation_map)?;

    // Now we have one gene disease entry for each gene or gene combination. This entry contains
    // HPOs for all of the gene-associated diseases that the user chose in the GUI
    let payload = GeneDiseaseEntity::get_overlap_matrix_payload(
        proband,
        &gd_entry_list,
        disease_counts,
        hpo.clone(),
    )?;

    Ok(payload)
}

/// Computes the total score (summed across all rows) for each entity/column.
fn compute_entity_sums(payload: &OverlapPlotPayload) -> HashMap<String, f64> {
    payload.entities.iter()
        .map(|entity| {
            let sum: f64 = payload.columns.iter()
                .map(|item| item.scores.get(entity).copied().unwrap_or(0.0))
                .sum();
            (entity.clone(), sum)
        })
        .collect()
}

/// Sorts entities (columns) descending by total score. Ties preserve original
/// input order for caller stability.
fn sort_entities_by_score(
    entities: &[String],
    entity_sums: &HashMap<String, f64>,
) -> Vec<String> {
    let mut sorted = entities.to_vec();
    sorted.sort_by(|a, b| {
        let sum_a = entity_sums.get(a).unwrap_or(&0.0);
        let sum_b = entity_sums.get(b).unwrap_or(&0.0);
        sum_b.total_cmp(sum_a) // descending
    });
    sorted
}

/// Sort key + payload for a row with at least one full (score == 1.0) match.
type FullMatchKey = (isize, Vec<usize>, usize, OverlapPlotItem);
/// Sort key + payload for a row with a best partial (0.0 < score < 1.0) match.
type PartialMatchKey = (usize, f64, usize, OverlapPlotItem);
/// Sort key + payload for a row with no matches at all.
type ZeroMatchKey = (usize, OverlapPlotItem);

enum RowBucket {
    Full(FullMatchKey),
    Partial(PartialMatchKey),
    Zero(ZeroMatchKey),
}

/// Classifies a single row into the full/partial/zero-match bucket, computing
/// whatever sort key that bucket needs relative to the newly sorted entity order.
fn classify_row(
    item: OverlapPlotItem,
    sorted_entities: &[String],
    original_positions: &HashMap<String, usize>,
) -> RowBucket {
    let input_idx = *original_positions.get(&item.hpo_id).unwrap_or(&0);

    let full_positions: Vec<usize> = sorted_entities.iter()
        .enumerate()
        .filter(|(_, entity)| {
            let score = item.scores.get(*entity).copied().unwrap_or(0.0);
            (score - 1.0).abs() < f64::EPSILON
        })
        .map(|(i, _)| i)
        .collect();

    if !full_positions.is_empty() {
        let count_key = -(full_positions.len() as isize);
        return RowBucket::Full((count_key, full_positions, input_idx, item));
    }

    let mut best_score = 0.0;
    let mut primary_pos = None;
    for (i, entity) in sorted_entities.iter().enumerate() {
        let score = item.scores.get(entity).copied().unwrap_or(0.0);
        if score > best_score {
            best_score = score;
            primary_pos = Some(i);
        }
    }

    match primary_pos {
        Some(pos) => RowBucket::Partial((pos, best_score, input_idx, item)),
        None => RowBucket::Zero((input_idx, item)),
    }
}

/// Splits all rows into their three sort buckets.
fn classify_rows(
    columns: Vec<OverlapPlotItem>,
    sorted_entities: &[String],
    original_positions: &HashMap<String, usize>,
) -> (Vec<FullMatchKey>, Vec<PartialMatchKey>, Vec<ZeroMatchKey>) {
    let mut full_keys = Vec::new();
    let mut partial_keys = Vec::new();
    let mut zero_keys = Vec::new();

    for item in columns {
        match classify_row(item, sorted_entities, original_positions) {
            RowBucket::Full(k) => full_keys.push(k),
            RowBucket::Partial(k) => partial_keys.push(k),
            RowBucket::Zero(k) => zero_keys.push(k),
        }
    }

    (full_keys, partial_keys, zero_keys)
}

/// Sorts full-match rows: by descending match count, then lexicographic
/// match-position vector, then stable original position.
fn sort_full_matches(mut keys: Vec<FullMatchKey>) -> Vec<FullMatchKey> {
    keys.sort_by(|a, b| {
        match a.0.cmp(&b.0) {
            std::cmp::Ordering::Equal => match a.1.cmp(&b.1) {
                std::cmp::Ordering::Equal => a.2.cmp(&b.2),
                other => other,
            },
            other => other,
        }
    });
    keys
}

/// Sorts partial-match rows: by ascending leftmost primary position, then
/// descending best score, then stable original position.
fn sort_partial_matches(mut keys: Vec<PartialMatchKey>) -> Vec<PartialMatchKey> {
    keys.sort_by(|a, b| {
        match a.0.cmp(&b.0) {
            std::cmp::Ordering::Equal => match b.1.total_cmp(&a.1) {
                std::cmp::Ordering::Equal => a.2.cmp(&b.2),
                other => other,
            },
            other => other,
        }
    });
    keys
}

/// Sorts zero-match rows: stable original position order.
fn sort_zero_matches(mut keys: Vec<ZeroMatchKey>) -> Vec<ZeroMatchKey> {
    keys.sort_by_key(|k| k.0);
    keys
}

/// Sorts the columns (gene entities) and rows (PresenceMatrixItem) following the
/// explicit hierarchical block rules established in the Python pipeline.
pub fn sort_presence_payload(payload: OverlapPlotPayload) -> OverlapPlotPayload {
    if payload.entities.is_empty() || payload.columns.is_empty() {
        return payload;
    }
    let n_columns = payload.n_columns();

    let entity_sums = compute_entity_sums(&payload);
    let sorted_entities = sort_entities_by_score(&payload.entities, &entity_sums);

    let original_positions: HashMap<String, usize> = payload.columns.iter()
        .enumerate()
        .map(|(i, item)| (item.hpo_id.clone(), i))
        .collect();

    let (full_keys, partial_keys, zero_keys) =
        classify_rows(payload.columns, &sorted_entities, &original_positions);

    let full_keys = sort_full_matches(full_keys);
    let partial_keys = sort_partial_matches(partial_keys);
    let zero_keys = sort_zero_matches(zero_keys);

    let mut sorted_columns = Vec::with_capacity(n_columns);
    sorted_columns.extend(full_keys.into_iter().map(|k| k.3));
    sorted_columns.extend(partial_keys.into_iter().map(|k| k.3));
    sorted_columns.extend(zero_keys.into_iter().map(|k| k.1));

    OverlapPlotPayload {
        entities: sorted_entities,
        columns: sorted_columns,
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use rstest::{fixture, rstest};
    use std::collections::HashMap;

    // -----------------------------------------------------------------
    // Fixtures
    // -----------------------------------------------------------------

    fn item(hpo_id: &str, name: &str, scores: &[(&str, f64)]) -> OverlapPlotItem {
        OverlapPlotItem {
            hpo_id: hpo_id.to_string(),
            hpo_name: name.to_string(),
            scores: scores.iter().map(|(k, v)| (k.to_string(), *v)).collect(),
        }
    }

    #[fixture]
    fn three_entities() -> Vec<String> {
        vec!["GENE_A".to_string(), "GENE_B".to_string(), "GENE_C".to_string()]
    }

    // -----------------------------------------------------------------
    // compute_entity_sums
    // -----------------------------------------------------------------

    #[rstest]
    fn entity_sums_add_scores_across_all_rows(three_entities: Vec<String>) {
        let payload = OverlapPlotPayload {
            entities: three_entities,
            columns: vec![
                item("HP:0001", "Term 1", &[("GENE_A", 1.0), ("GENE_B", 0.5)]),
                item("HP:0002", "Term 2", &[("GENE_A", 0.0), ("GENE_B", 0.5)]),
            ],
        };

        let sums = compute_entity_sums(&payload);
        assert_eq!(sums.get("GENE_A"), Some(&1.0));
        assert_eq!(sums.get("GENE_B"), Some(&1.0));
        assert_eq!(sums.get("GENE_C"), Some(&0.0)); // no rows reference it at all
    }

    #[rstest]
    fn entity_sums_missing_score_defaults_to_zero() {
        let payload = OverlapPlotPayload {
            entities: vec!["GENE_A".to_string()],
            columns: vec![item("HP:0001", "Term 1", &[])], // no scores recorded at all
        };
        let sums = compute_entity_sums(&payload);
        assert_eq!(sums.get("GENE_A"), Some(&0.0));
    }

    // -----------------------------------------------------------------
    // sort_entities_by_score
    // -----------------------------------------------------------------

    #[rstest]
    fn entities_sort_descending_by_score(three_entities: Vec<String>) {
        let sums: HashMap<String, f64> = [
            ("GENE_A".to_string(), 0.4),
            ("GENE_B".to_string(), 0.9),
            ("GENE_C".to_string(), 0.1),
        ].into_iter().collect();

        let sorted = sort_entities_by_score(&three_entities, &sums);
        assert_eq!(sorted, vec!["GENE_B", "GENE_A", "GENE_C"]);
    }

    #[rstest]
    fn tied_scores_preserve_original_order(three_entities: Vec<String>) {
        let sums: HashMap<String, f64> = [
            ("GENE_A".to_string(), 0.5),
            ("GENE_B".to_string(), 0.5),
            ("GENE_C".to_string(), 0.5),
        ].into_iter().collect();

        let sorted = sort_entities_by_score(&three_entities, &sums);
        // All tied: stable sort should preserve input order exactly.
        assert_eq!(sorted, vec!["GENE_A", "GENE_B", "GENE_C"]);
    }

    #[rstest]
    fn entity_missing_from_sums_defaults_to_zero(three_entities: Vec<String>) {
        let sums: HashMap<String, f64> = [
            ("GENE_A".to_string(), 0.9),
            // GENE_B, GENE_C intentionally absent
        ].into_iter().collect();

        let sorted = sort_entities_by_score(&three_entities, &sums);
        assert_eq!(sorted[0], "GENE_A");
    }

    // -----------------------------------------------------------------
    // classify_row
    // -----------------------------------------------------------------

    #[rstest]
    fn row_with_one_full_match_classifies_as_full(three_entities: Vec<String>) {
        let original_positions: HashMap<String, usize> =
            [("HP:0001".to_string(), 0)].into_iter().collect();

        let row = item("HP:0001", "Term 1",  &[("GENE_A", 1.0), ("GENE_B", 0.3)]);
        match classify_row(row, &three_entities, &original_positions) {
            RowBucket::Full((count_key, positions, idx, _)) => {
                assert_eq!(count_key, -1); // one full match => count_key = -1
                assert_eq!(positions, vec![0]); // GENE_A is at sorted-index 0
                assert_eq!(idx, 0);
            }
            _ => panic!("expected Full bucket"),
        }
    }

    #[rstest]
    fn row_with_two_full_matches_records_both_positions(three_entities: Vec<String>) {
        let original_positions: HashMap<String, usize> =
            [("HP:0001".to_string(), 2)].into_iter().collect();

        let row = item("HP:0001", "Term 1", &[("GENE_A", 1.0), ("GENE_C", 1.0), ("GENE_B", 0.0)]);
        match classify_row(row, &three_entities, &original_positions) {
            RowBucket::Full((count_key, positions, idx, _)) => {
                assert_eq!(count_key, -2);
                assert_eq!(positions, vec![0, 2]); // GENE_A=0, GENE_C=2 in sorted_entities order
                assert_eq!(idx, 2);
            }
            _ => panic!("expected Full bucket"),
        }
    }

    #[rstest]
    fn row_with_no_full_match_but_partial_classifies_as_partial(three_entities: Vec<String>) {
        let original_positions: HashMap<String, usize> =
            [("HP:0001".to_string(), 1)].into_iter().collect();

        let row = item("HP:0001", "Term 1", &[("GENE_A", 0.4), ("GENE_B", 0.7), ("GENE_C", 0.2)]);
        match classify_row(row, &three_entities, &original_positions) {
            RowBucket::Partial((pos, best_score, idx, _)) => {
                assert_eq!(pos, 1); // GENE_B is the best score, at sorted-index 1
                assert_eq!(best_score, 0.7);
                assert_eq!(idx, 1);
            }
            _ => panic!("expected Partial bucket"),
        }
    }

    #[rstest]
    fn row_with_all_zero_scores_classifies_as_zero(three_entities: Vec<String>) {
        let original_positions: HashMap<String, usize> =
            [("HP:0001".to_string(), 3)].into_iter().collect();

        let row = item("HP:0001", "Term 1", &[("GENE_A", 0.0), ("GENE_B", 0.0), ("GENE_C", 0.0)]);
        match classify_row(row, &three_entities, &original_positions) {
            RowBucket::Zero((idx, _)) => assert_eq!(idx, 3),
            _ => panic!("expected Zero bucket"),
        }
    }

    #[rstest]
    fn row_missing_from_original_positions_defaults_to_zero_index(three_entities: Vec<String>) {
        let original_positions: HashMap<String, usize> = HashMap::new(); // empty
        let row = item("HP:9999", "Term 999", &[("GENE_A", 1.0)]);
        match classify_row(row, &three_entities, &original_positions) {
            RowBucket::Full((_, _, idx, _)) => assert_eq!(idx, 0),
            _ => panic!("expected Full bucket"),
        }
    }

    // -----------------------------------------------------------------
    // sort_full_matches
    // -----------------------------------------------------------------

    #[rstest]
    fn full_matches_sort_by_descending_match_count() {
        let keys = vec![
            (-1isize, vec![0usize], 0usize, item("HP:0001", "Term 1", &[])),
            (-2isize, vec![0usize, 1usize], 1usize, item("HP:0002", "Term 2", &[])),
        ];
        let sorted = sort_full_matches(keys);
        assert_eq!(sorted[0].3.hpo_id, "HP:0002"); // 2 matches sorts before 1 match
        assert_eq!(sorted[1].3.hpo_id, "HP:0001");
    }

    #[rstest]
    fn full_matches_tie_break_by_position_vector_then_original_index() {
        let keys = vec![
            (-1isize, vec![1usize], 0usize, item("HP:0001", "Term 1", &[])), // match at position 1
            (-1isize, vec![0usize], 1usize, item("HP:0002", "Term 2", &[])), // match at position 0
        ];
        let sorted = sort_full_matches(keys);
        // Same count (-1), so compares position vectors: [0] < [1]
        assert_eq!(sorted[0].3.hpo_id, "HP:0002");
        assert_eq!(sorted[1].3.hpo_id, "HP:0001");
    }

    // -----------------------------------------------------------------
    // sort_partial_matches
    // -----------------------------------------------------------------

    #[rstest]
    fn partial_matches_sort_by_ascending_primary_position() {
        let keys = vec![
            (2usize, 0.5f64, 0usize, item("HP:0001", "Term 1", &[])),
            (0usize, 0.3f64, 1usize, item("HP:0002", "Term 2", &[])),
        ];
        let sorted = sort_partial_matches(keys);
        assert_eq!(sorted[0].3.hpo_id, "HP:0002"); // position 0 comes first
        assert_eq!(sorted[1].3.hpo_id, "HP:0001");
    }

    #[rstest]
    fn partial_matches_tie_break_by_descending_score() {
        let keys = vec![
            (0usize, 0.3f64, 0usize, item("HP:0001", "Term 1", &[])),
            (0usize, 0.9f64, 1usize, item("HP:0002", "Term 2", &[])),
        ];
        let sorted = sort_partial_matches(keys);
        // Same position (0), so higher score (0.9) sorts first
        assert_eq!(sorted[0].3.hpo_id, "HP:0002");
        assert_eq!(sorted[1].3.hpo_id, "HP:0001");
    }

    #[rstest]
    fn partial_matches_final_tiebreak_by_original_index() {
        let keys = vec![
            (0usize, 0.5f64, 5usize, item("HP:0001", "Term 1", &[])),
            (0usize, 0.5f64, 2usize, item("HP:0002", "Term 2", &[])),
        ];
        let sorted = sort_partial_matches(keys);
        // Same position and score: lower original index (2) sorts first
        assert_eq!(sorted[0].3.hpo_id, "HP:0002");
        assert_eq!(sorted[1].3.hpo_id, "HP:0001");
    }

    // -----------------------------------------------------------------
    // sort_zero_matches
    // -----------------------------------------------------------------

    #[rstest]
    fn zero_matches_sort_by_original_index() {
        let keys = vec![
            (3usize, item("HP:0001", "Term 1", &[])),
            (1usize, item("HP:0002", "Term 2", &[])),
        ];
        let sorted = sort_zero_matches(keys);
        assert_eq!(sorted[0].1.hpo_id, "HP:0002");
        assert_eq!(sorted[1].1.hpo_id, "HP:0001");
    }
}