use std::sync::LazyLock;

use crate::model::simple_term::SimpleOntologyTerm;
use ontolius::TermId;

/// Root of the HPO phenotype branch; sharing only this term means no meaningful match.
pub const ROOT_PHENOTYPE_ID: &str = "HP:0000118";
pub static ROOT_PHENOTYPE_TERM_ID: std::sync::LazyLock<TermId> = std::sync::LazyLock::new(|| {
    let tid: TermId  = ROOT_PHENOTYPE_ID.parse().unwrap();
    tid
});
pub static SIMPLE_ROOT_TERM: LazyLock<SimpleOntologyTerm> = LazyLock::new(|| {
    SimpleOntologyTerm::new(ROOT_PHENOTYPE_ID, "Phenotypic abnormality").unwrap()
});
