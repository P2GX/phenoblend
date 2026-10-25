// New file `blend/batch_tsv.rs` (add `pub mod batch_tsv;` to blend's mod.rs).
// Uses `serde_json` and `rstest` (already in your crate).

//! The tab-separated output of a batch run, and the assignment for the shuffled-pair control.
//!
//! One row per phenopacket and gene set. A phenopacket that could not be analyzed gets one row
//! with `status = error` and the reason, so a single bad file never stops a run and the
//! failures can be counted afterwards.
//!
//! The columns are chosen so that the cohort statistics can be computed outside the
//! application (for example in R or Python) from the file alone: besides the p-values, each
//! row has `P(X > k)` and `P(X = k)`, from which a randomized p-value
//! `p_greater + U · p_equal` (exactly uniform under the null) can be made with any `U`.
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

use crate::blend::dto::shared_term_enrichment::SharedTermEnrichment;

/// The columns of the output, in order.
pub const TSV_COLUMNS: [&str; 27] = [
    "phenopacket_id",
    "individual_id",
    "status",
    "message",
    "mode",
    "genes",
    "n_genes",
    "n_observed_terms",
    "n_urn",
    "n_shared",
    "n_drawn",
    "n_observed",
    "expected_shared",
    "p_value",
    "p_mid",
    "p_greater",
    "p_equal",
    "degenerate",
    "sim_n",
    "sim_seed",
    "sim_p_value",
    "sim_expected_shared",
    "sim_n_greater",
    "sim_n_equal",
    "sim_mean_weight_shared",
    "sim_mean_weight_unshared",
    "sim_n_downweighted_terms",
];

/// The header line of the table.
pub fn tsv_header() -> String {
    TSV_COLUMNS.join("\t")
}

/// Comment lines (`# key: value`) that record how the file was made, for the top of the file.
/// R and Python read the table with `comment = "#"` / `comment="#"`.
pub fn provenance_lines(items: &[(&str, &str)]) -> String {
    items
        .iter()
        .map(|(key, value)| format!("# {}: {}\n", clean(key), clean(value)))
        .collect()
}

/// Tabs and line breaks would break the table, so they become spaces.
fn clean(text: &str) -> String {
    text.replace(['\t', '\n', '\r'], " ")
}

/// What identifies a row.
pub struct RowContext<'a> {
    pub phenopacket_id: &'a str,
    pub individual_id: &'a str,
    /// `observed` for the real analysis, `shuffled` for the control.
    pub mode: &'a str,
    /// Number of distinct observed terms of the individual.
    pub n_observed_terms: usize,
}

/// The row for one enrichment result.
pub fn result_row(ctx: &RowContext, r: &SharedTermEnrichment) -> String {
    let k = r.n_observed;
    let p_equal = r.pmf[k];
    let p_greater: f64 = r.pmf[k + 1..].iter().sum();

    let simulation = match &r.simulated {
        Some(s) => {
            let n = s.n_sim as f64;
            let n_equal = (s.pmf[k] * n).round() as usize;
            let n_greater = (s.pmf[k + 1..].iter().sum::<f64>() * n).round() as usize;
            [
                s.n_sim.to_string(),
                s.seed.to_string(),
                s.p_value.to_string(),
                s.expected_shared.to_string(),
                n_greater.to_string(),
                n_equal.to_string(),
                s.mean_weight_shared.map(|w| w.to_string()).unwrap_or_default(),
                s.mean_weight_unshared.map(|w| w.to_string()).unwrap_or_default(),
                s.n_downweighted_terms.to_string(),
            ]
        }
        None => Default::default(),
    };

    let mut cells: Vec<String> = vec![
        clean(ctx.phenopacket_id),
        clean(ctx.individual_id),
        "ok".to_string(),
        String::new(),
        clean(ctx.mode),
        clean(&r.genes.join("+")),
        r.genes.len().to_string(),
        ctx.n_observed_terms.to_string(),
        r.n_urn.to_string(),
        r.n_shared.to_string(),
        r.n_drawn.to_string(),
        r.n_observed.to_string(),
        r.expected_shared.to_string(),
        r.p_value.to_string(),
        r.p_mid.to_string(),
        p_greater.to_string(),
        p_equal.to_string(),
        r.degenerate.to_string(),
    ];
    cells.extend(simulation);
    cells.join("\t")
}

/// The row for a phenopacket that could not be analyzed.
pub fn error_row(phenopacket_id: &str, mode: &str, message: &str) -> String {
    let mut cells = vec![String::new(); TSV_COLUMNS.len()];
    cells[0] = clean(phenopacket_id);
    cells[2] = "error".to_string();
    cells[3] = clean(message);
    cells[4] = clean(mode);
    cells.join("\t")
}

/// For the shuffled-pair control: `assignment[i]` is the index of the phenopacket whose gene
/// set individual `i` is paired with. It is a permutation of `0..n` without fixed points (no
/// individual keeps their own genes), and it is the same for the same `seed`.
///
/// The indices are shuffled (Fisher–Yates) and then each is sent to the next one in the
/// shuffled order, which makes a single cycle.
///
/// # Errors
/// Returns an error if `n < 2`, because then no individual can be given other genes.
pub fn shuffle_assignment(n: usize, seed: u64) -> Result<Vec<usize>, String> {
    if n < 2 {
        return Err(format!("need at least two phenopackets to shuffle gene sets, got {n}"));
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let mut order: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        let j = rng.random_range(0..=i);
        order.swap(i, j);
    }
    let mut assignment = vec![0usize; n];
    for position in 0..n {
        assignment[order[position]] = order[(position + 1) % n];
    }
    Ok(assignment)
}

/// A row without results: the phenopacket was not analyzed, for the stated reason.
fn problem_row(phenopacket_id: &str, mode: &str, status: &str, message: &str) -> String {
    let mut cells = vec![String::new(); TSV_COLUMNS.len()];
    cells[0] = clean(phenopacket_id);
    cells[2] = status.to_string();
    cells[3] = clean(message);
    cells[4] = clean(mode);
    cells.join("\t")
}

/// A phenopacket left out on purpose, for example a non-genetic diagnosis.
pub fn skipped_row(phenopacket_id: &str, mode: &str, message: &str) -> String {
    problem_row(phenopacket_id, mode, "skipped", message)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    /// A result with simulation, built from JSON so the test does not depend on how the
    /// application builds its terms. N=105, K=19, m=20, k=2; the exact pmf has P(X=2) = 0.25 and
    /// P(X>2) = 0.375 + 0.125 = 0.5, so P(X>=2) = 0.75; the simulation has 1000 individuals,
    /// 300 at 2 and 400 above.
    fn sample() -> SharedTermEnrichment {
        let mut pmf = vec![0.0; 21];
        pmf[0] = 0.125;
        pmf[1] = 0.125;
        pmf[2] = 0.25;
        pmf[3] = 0.375;
        pmf[4] = 0.125;
        let mut sim_pmf = vec![0.0; 21];
        sim_pmf[1] = 0.3;
        sim_pmf[2] = 0.3;
        sim_pmf[3] = 0.4;
        let json = serde_json::json!({
            "genes": ["NF1", "PTPN11"],
            "nUrn": 105, "nShared": 19, "nDrawn": 20, "nObserved": 2,
            "pmf": pmf, "pValue": 0.75, "pMid": 0.625,
            "expectedShared": 3.6, "degenerate": false,
            "simulated": {
                "nSim": 1000, "seed": "1470397748979433974", "pmf": sim_pmf,
                "pValue": 0.7, "expectedShared": 2.1, "nDownweightedTerms": 40,
                "meanWeightShared": 0.31, "meanWeightUnshared": null, "urn": []
            }
        });
        serde_json::from_value(json).expect("the sample result must deserialize")
    }

    fn ctx() -> RowContext<'static> {
        RowContext {
            phenopacket_id: "ppkt-1",
            individual_id: "patient-1",
            mode: "observed",
            n_observed_terms: 23,
        }
    }

    fn cells(row: &str) -> Vec<&str> {
        row.split('\t').collect()
    }

    fn cell<'a>(row: &'a [&'a str], column: &str) -> &'a str {
        row[TSV_COLUMNS.iter().position(|c| *c == column).expect("known column")]
    }

    #[test]
    fn the_header_and_every_row_have_the_same_number_of_columns() {
        assert_eq!(tsv_header().split('\t').count(), TSV_COLUMNS.len());
        assert_eq!(cells(&result_row(&ctx(), &sample())).len(), TSV_COLUMNS.len());
        assert_eq!(cells(&error_row("p", "observed", "boom")).len(), TSV_COLUMNS.len());
    }

    #[test]
    fn column_names_are_unique() {
        let mut names: Vec<&str> = TSV_COLUMNS.to_vec();
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), TSV_COLUMNS.len());
    }

    #[test]
    fn a_result_row_holds_the_values_of_the_result() {
        let row = result_row(&ctx(), &sample());
        let c = cells(&row);
        assert_eq!(cell(&c, "phenopacket_id"), "ppkt-1");
        assert_eq!(cell(&c, "individual_id"), "patient-1");
        assert_eq!(cell(&c, "status"), "ok");
        assert_eq!(cell(&c, "message"), "");
        assert_eq!(cell(&c, "mode"), "observed");
        assert_eq!(cell(&c, "genes"), "NF1+PTPN11");
        assert_eq!(cell(&c, "n_genes"), "2");
        assert_eq!(cell(&c, "n_observed_terms"), "23");
        assert_eq!(cell(&c, "n_urn"), "105");
        assert_eq!(cell(&c, "n_shared"), "19");
        assert_eq!(cell(&c, "n_drawn"), "20");
        assert_eq!(cell(&c, "n_observed"), "2");
        assert_eq!(cell(&c, "p_value"), "0.75");
        assert_eq!(cell(&c, "degenerate"), "false");
    }

    #[test]
    fn p_greater_and_p_equal_split_the_p_value() {
        // P(X >= 2) = P(X = 2) + P(X > 2) = 0.25 + 0.5 = 0.75
        let row = result_row(&ctx(), &sample());
        let c = cells(&row);
        let p_equal: f64 = cell(&c, "p_equal").parse().unwrap();
        let p_greater: f64 = cell(&c, "p_greater").parse().unwrap();
        assert!((p_equal - 0.25).abs() < 1e-12);
        assert!((p_greater - 0.5).abs() < 1e-12);
        let p_value: f64 = cell(&c, "p_value").parse().unwrap();
        assert!((p_equal + p_greater - p_value).abs() < 1e-12);
    }

    #[test]
    fn the_simulation_columns_hold_counts_recovered_from_the_pmf() {
        let row = result_row(&ctx(), &sample());
        let c = cells(&row);
        assert_eq!(cell(&c, "sim_n"), "1000");
        assert_eq!(cell(&c, "sim_seed"), "1470397748979433974"); // every digit
        assert_eq!(cell(&c, "sim_n_equal"), "300");
        assert_eq!(cell(&c, "sim_n_greater"), "400");
        assert_eq!(cell(&c, "sim_mean_weight_shared"), "0.31");
        assert_eq!(cell(&c, "sim_mean_weight_unshared"), ""); // none
        assert_eq!(cell(&c, "sim_n_downweighted_terms"), "40");
    }

    #[test]
    fn without_a_simulation_its_columns_are_empty() {
        let mut result = sample();
        result.simulated = None;
        let row = result_row(&ctx(), &result);
        let c = cells(&row);
        assert_eq!(c.len(), TSV_COLUMNS.len());
        assert!(c[TSV_COLUMNS.iter().position(|n| *n == "sim_n").unwrap()..].iter().all(|x| x.is_empty()));
    }

    #[test]
    fn an_error_row_has_the_reason_and_no_numbers() {
        let row = error_row("ppkt-9", "observed", "disease OMIM:1 not found");
        let c = cells(&row);
        assert_eq!(cell(&c, "phenopacket_id"), "ppkt-9");
        assert_eq!(cell(&c, "status"), "error");
        assert_eq!(cell(&c, "message"), "disease OMIM:1 not found");
        assert_eq!(cell(&c, "n_urn"), "");
        assert_eq!(cell(&c, "p_value"), "");
    }

    #[test]
    fn tabs_and_line_breaks_in_text_cannot_break_the_table() {
        let row = error_row("a\tb", "observed", "line one\nline two\ttab");
        assert_eq!(row.matches('\n').count(), 0);
        assert_eq!(cells(&row).len(), TSV_COLUMNS.len());
        assert!(row.contains("line one line two tab"));
    }

    #[test]
    fn provenance_is_written_as_comment_lines() {
        let text = provenance_lines(&[("hpo_version", "2026-09-01"), ("note", "a\nb")]);
        assert_eq!(text, "# hpo_version: 2026-09-01\n# note: a b\n");
    }

    // ---------------------------------------------------------------------------------
    // The shuffled-pair control
    // ---------------------------------------------------------------------------------

    #[rstest]
    #[case(2)]
    #[case(3)]
    #[case(10)]
    #[case(250)]
    fn the_assignment_is_a_permutation_without_fixed_points(#[case] n: usize) {
        for seed in [1u64, 2, 99, 2026] {
            let assignment = shuffle_assignment(n, seed).unwrap();
            assert_eq!(assignment.len(), n);
            let mut sorted = assignment.clone();
            sorted.sort_unstable();
            assert_eq!(sorted, (0..n).collect::<Vec<_>>(), "must use every individual once");
            assert!(assignment.iter().enumerate().all(|(i, a)| i != *a), "no one keeps their own");
        }
    }

    #[test]
    fn the_assignment_is_reproducible_and_depends_on_the_seed() {
        assert_eq!(shuffle_assignment(50, 7).unwrap(), shuffle_assignment(50, 7).unwrap());
        assert_ne!(shuffle_assignment(50, 7).unwrap(), shuffle_assignment(50, 8).unwrap());
    }

    #[test]
    fn two_individuals_swap_their_gene_sets() {
        assert_eq!(shuffle_assignment(2, 5).unwrap(), [1, 0]);
    }

    #[rstest]
    #[case(0)]
    #[case(1)]
    fn fewer_than_two_individuals_cannot_be_shuffled(#[case] n: usize) {
        let err = shuffle_assignment(n, 1).unwrap_err();
        assert!(err.contains("at least two"), "{err}");
    }
}
