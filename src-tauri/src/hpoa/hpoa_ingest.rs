use std::collections::{HashMap, HashSet};

use log::warn;

use oboannotation::hpo::{FrequencyData, HpoAnnotations};
use oboannotation::io::AnnotationLoader;


use ontolius::TermId;

use crate::util::errors::PhenoblendError;
use crate::hpoa::disease_model::SimpleDiseaseModel;

use std::sync::LazyLock;
static EXCLUDED_TERM_ID: LazyLock<TermId> = LazyLock::new(|| {
    let tid: TermId  = "HP:0040285".parse().unwrap();
    tid
});

static OMIM_PREFIX: &str = "OMIM";

/// Key used to track the max frequency seen for a given disease/term pair
/// across all HPOA lines (multiple publications can annotate the same
/// disease/term with different frequencies).
type DiseaseTermKey = (TermId, TermId);


/// Convert an HPOA `Frequency` value into a numeric frequency in [0.0, 1.0].
///
/// - Missing frequency defaults to 1.0 (treated as "always present" —
///   NOTE: this conflates "unspecified" with "obligate"; revisit if that
///   distinction matters downstream).
/// - `Ratio` is computed as a true floating point division. Both operands
///   MUST be cast to f64 individually before dividing — dividing as
///   integers first (`(numerator / denominator) as f64`) truncates any
///   ratio < 1 down to 0.0, which silently reclassifies real annotations
///   (e.g. "74/95") as excluded. This was the root cause of terms like
///   Lisch nodules / NF1 disappearing.
fn frequency_to_f64(freq: Option<oboannotation::hpo::FrequencyData>) -> f64 {
    freq.map_or(1.0, |f| match f {
        oboannotation::hpo::FrequencyData::TermId(term_id) if term_id == *EXCLUDED_TERM_ID => 0.0,
        oboannotation::hpo::FrequencyData::TermId(_) => 1.0,
        oboannotation::hpo::FrequencyData::Ratio{ n:numerator, m:denominator } => {
            if denominator == 0 {
                0.0
            } else {
                numerator as f64 / denominator as f64
            }
        }
        oboannotation::hpo::FrequencyData::Percentage(perc) => 0.01*perc,
    })
}

struct ParsedRow {
    disease_id: TermId,
    disease_name: String,
    hpo_term_id: TermId,
    freq: f64,
}


/// Folds all rows into:
/// - `max_freq`: the maximum frequency ever observed for each (disease, term) pair
/// - `disease_id_to_name`: disease id -> display name
///
/// Keeping the max, rather than the first or last row seen, means a term
/// annotated as rare in one publication and common in another is reported
/// at its most-supported (highest) frequency, and a term with any positive
/// evidence is never lost to a later exclusion row for the same pair.
fn aggregate_rows(
    rows: impl IntoIterator<Item = ParsedRow>,
) -> (HashMap<DiseaseTermKey, f64>, HashMap<TermId, String>) {
    let mut max_freq: HashMap<DiseaseTermKey, f64> = HashMap::new();
    let mut disease_id_to_name: HashMap<TermId, String> = HashMap::new();

    for row in rows {
        let key = (row.disease_id.clone(), row.hpo_term_id.clone());
        max_freq
            .entry(key)
            .and_modify(|f| *f = f.max(row.freq))
            .or_insert(row.freq);

        disease_id_to_name
            .entry(row.disease_id.clone())
            .or_insert(row.disease_name);
    }

    (max_freq, disease_id_to_name)
}

/// Splits the flat (disease, term) -> max_freq map into, per disease:
/// - observed_hpo_ids: terms with freq > 0.0
/// - term_frequencies: freq value for each observed term
/// - excluded_hpo_ids: terms whose max freq across all publications is 0.0
///   (i.e. never positively observed anywhere)
fn build_per_disease_maps(
    max_freq: HashMap<DiseaseTermKey, f64>,
) -> HashMap<TermId, (HashSet<TermId>, HashMap<TermId, f64>, HashSet<TermId>)> {
    let mut per_disease: HashMap<TermId, (HashSet<TermId>, HashMap<TermId, f64>, HashSet<TermId>)> =
        HashMap::new();

    for ((disease_id, term_id), freq) in max_freq {
        let entry = per_disease.entry(disease_id).or_insert_with(|| {
            (HashSet::new(), HashMap::new(), HashSet::new())
        });
        if freq > 0.0 {
            entry.0.insert(term_id.clone());
            entry.1.insert(term_id, freq);
        } else {
            entry.2.insert(term_id);
        }
    }

    per_disease
}

fn build_models(
    per_disease: HashMap<TermId, (HashSet<TermId>, HashMap<TermId, f64>, HashSet<TermId>)>,
    mut disease_id_to_name: HashMap<TermId, String>,
) -> HashMap<TermId, SimpleDiseaseModel> {
    let mut simple_model_d = HashMap::new();

    for (disease_id, (observed_hpo_ids, term_frequencies, excluded_hpo_ids)) in per_disease {
        let Some(disease_name) = disease_id_to_name.remove(&disease_id) else {
            warn!("Skipping {disease_id}: missing disease name");
            continue;
        };

        let model = SimpleDiseaseModel::new(
            disease_id.clone(),
            disease_name,
            observed_hpo_ids,
            term_frequencies,
            excluded_hpo_ids,
        );
        simple_model_d.insert(disease_id, model);
    }
    simple_model_d
}


fn get_term_f(term_id: &TermId) -> Result<f64, String> {
    super::frequency_term::FrequencyTerm::try_from(term_id).map(super::frequency_term::FrequencyTerm::mean_frequency)
}

fn get_freq_float(freqd: &FrequencyData) -> Result<f64, String>  {
    match freqd {
        FrequencyData::TermId(term_id) => get_term_f(term_id),
        FrequencyData::Ratio { n, m } if *m > 0 => Ok(f64::from(*n) / f64::from(*m)),
        FrequencyData::Ratio { .. } => Ok(0.0), // zero denominator: no usable ratio
        FrequencyData::Percentage(perc) => Ok(f64::from(*perc) / 100.0),
    }
}

pub fn load_hpoa_d(fpath_hpoa: &str) -> Result<HashMap<TermId, SimpleDiseaseModel>, PhenoblendError> {
    let data: oboannotation::hpo::HpoAnnotations = HpoAnnotations::load_from_path(fpath_hpoa)
        .map_err(|e| PhenoblendError::hpoa_load_error(format!("{fpath_hpoa}: {e}")))?;
    let rows: Vec<ParsedRow> = data
        .annotations()
        .into_iter()
        .filter(|line| line.disease_id.prefix().to_string() == OMIM_PREFIX)
        .map(|line| {
            let freq: f64 = match &line.frequency {
                Some(freq) => get_freq_float(freq.data()).map_err(|e|PhenoblendError::parse_error(e))?,
                None => 1.0, // by default, if no freq data available it is assumed to be 100%.
            };
            Ok(ParsedRow {
                disease_id: line.disease_id.clone(),
                disease_name: line.disease_name.clone(),
                hpo_term_id: line.phenotype_term_id.clone(),
                freq,
            })
        })
        .collect::<Result<_, PhenoblendError>>()?;
    let (max_freq, disease_id_to_name) = aggregate_rows(rows);
    let per_disease = build_per_disease_maps(max_freq);
    Ok(build_models(per_disease, disease_id_to_name))
}


#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    fn tid(s: &str) -> TermId {
        s.parse().unwrap()
    }

    // --- frequency_to_f64 ---

    #[rstest]
    #[case::none_defaults_to_obligate(None, 1.0)]
    #[case::excluded_term_id_is_zero(
        Some(oboannotation::hpo::FrequencyData::TermId(tid("HP:0040285"))),
        0.0
    )]
    #[case::other_term_id_is_present(
        Some(oboannotation::hpo::FrequencyData::TermId(tid("HP:0040283"))), // e.g. "Occasional"
        1.0
    )]
    #[case::explicit_frequency_passthrough(Some(oboannotation::hpo::FrequencyData::Percentage(42.0)), 0.42)]
    fn test_frequency_to_f64_simple_cases(#[case] input: Option<oboannotation::hpo::FrequencyData>, #[case] expected: f64) {
        assert_eq!(frequency_to_f64(input), expected);
    }

    #[rstest]
    // This is the regression case: a ratio less than 1 must NOT truncate to 0.
    #[case::proper_fraction_lisch_nodule_like(74, 95, 74.0 / 95.0)]
    #[case::proper_fraction_small(1, 2, 0.5)]
    #[case::exact_all_present(5, 5, 1.0)]
    #[case::zero_numerator_true_exclusion(0, 10, 0.0)]
    fn test_frequency_to_f64_ratio(
        #[case] numerator: u32,
        #[case] denominator: u32,
        #[case] expected: f64,
    ) {
        let result = frequency_to_f64(Some(oboannotation::hpo::FrequencyData::Ratio { n:numerator, m:denominator }));
        assert!(
            (result - expected).abs() < 1e-9,
            "expected {expected}, got {result}"
        );
        // Explicitly guard against the integer-truncation regression:
        // any proper fraction must be > 0.0 unless numerator is 0.
        if numerator > 0 {
            assert!(result > 0.0, "ratio {numerator}/{denominator} truncated to 0.0");
        }
    }

    #[rstest]
    fn test_frequency_to_f64_ratio_zero_denominator_no_panic() {
        let result = frequency_to_f64(Some(oboannotation::hpo::FrequencyData::Ratio { n: 3, m: 0 }));
        assert_eq!(result, 0.0);
    }



 
}