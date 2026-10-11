use ontolius::TermId;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum FrequencyTerm {
    #[serde(rename = "HP:0040280")] Obligate,      // 100%
    #[serde(rename = "HP:0040281")] VeryFrequent,  // 80-99%
    #[serde(rename = "HP:0040282")] Frequent,      // 30-79%
    #[serde(rename = "HP:0040283")] Occasional,    // 5-29%
    #[serde(rename = "HP:0040284")] VeryRare,      // 1-4%
    #[serde(rename = "HP:0040285")] Excluded,      // 0%
}

impl FrequencyTerm {
    pub const ALL: [FrequencyTerm; 6] = [
        Self::Obligate, Self::VeryFrequent, Self::Frequent,
        Self::Occasional, Self::VeryRare, Self::Excluded,
    ];

    pub const fn curie(self) -> &'static str {
        match self {
            Self::Obligate => "HP:0040280",
            Self::VeryFrequent => "HP:0040281",
            Self::Frequent => "HP:0040282",
            Self::Occasional => "HP:0040283",
            Self::VeryRare => "HP:0040284",
            Self::Excluded => "HP:0040285",
        }
    }

    /// Inclusive frequency range covered by the term, as fractions.
    pub const fn range(self) -> (f64, f64) {
        match self {
            Self::Obligate => (1.0, 1.0),
            Self::VeryFrequent => (0.80, 0.99),
            Self::Frequent => (0.30, 0.79),
            Self::Occasional => (0.05, 0.29),
            Self::VeryRare => (0.01, 0.04),
            Self::Excluded => (0.0, 0.0),
        }
    }

    /// Midpoint of the term's frequency range, as a fraction in `[0, 1]`.
    ///
    /// For example, `Frequent` (30-79%) gives `0.545`. `Obligate` and `Excluded`
    /// are the degenerate ranges 100% and 0%, so they return exactly `1.0` and `0.0`.
    pub fn mean_frequency(self) -> f64 {
        let (lo, hi) = self.range();
        (lo + hi) / 2.0
    }
}

impl TryFrom<&TermId> for FrequencyTerm {
    type Error = String;

    fn try_from(tid: &TermId) -> Result<Self, Self::Error> {
        let s = tid.to_string();
        Self::ALL
            .into_iter()
            .find(|t| t.curie() == s)
            .ok_or_else(|| format!("'{s}' is not an HPO frequency term"))
    }
}




#[cfg(test)]
mod tests {

    use super::*;
    use rstest::rstest;


    #[rstest]
    fn only_frequency_terms_accepted() {
        for t in FrequencyTerm::ALL {
            let tid: TermId = t.curie().parse().unwrap();
            assert_eq!(FrequencyTerm::try_from(&tid).unwrap(), t);
        }
        let other: TermId = "HP:0001250".parse().unwrap(); // not a valid Frequency Term in HPO
        assert!(FrequencyTerm::try_from(&other).is_err());
        assert!(serde_json::from_str::<FrequencyTerm>("\"HP:0000118\"").is_err());
        assert_eq!(
            serde_json::to_string(&FrequencyTerm::Frequent).unwrap(),
            "\"HP:0040282\""
        );
    }

    #[rstest]
    fn mean_frequency_is_midpoint_of_range() {
        let expected = [
            (FrequencyTerm::Obligate, 1.0),
            (FrequencyTerm::VeryFrequent, 0.895),
            (FrequencyTerm::Frequent, 0.545),
            (FrequencyTerm::Occasional, 0.17),
            (FrequencyTerm::VeryRare, 0.025),
            (FrequencyTerm::Excluded, 0.0),
        ];
        for (term, mean) in expected {
            assert!(
                (term.mean_frequency() - mean).abs() < 1e-9,
                "{term:?}: got {}, expected {mean}",
                term.mean_frequency()
            );
        }
    }

    #[rstest]
    fn means_decrease_monotonically_with_rarity() {
        let means: Vec<f64> = FrequencyTerm::ALL.iter().map(|t| t.mean_frequency()).collect();
        assert!(means.windows(2).all(|w| w[0] > w[1]));
    }

}