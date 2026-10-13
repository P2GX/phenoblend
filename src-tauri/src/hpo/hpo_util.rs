use ontolius::{TermId, ontology::{OntologyTerms, csr::FullCsrOntology}, term::MinimalTerm};
use crate::model::simple_term::SimpleOntologyTerm;





pub fn ontolius_to_simple(hpo: &FullCsrOntology, term_id: &TermId) -> Result<SimpleOntologyTerm, String>
{
        match hpo.term_by_id(term_id) {
            Some(term) => Ok(SimpleOntologyTerm::from_term_id(term_id, term.name())),
            None => Err(format!("Could not create Term from {}", term_id)),
        }
    
}

/// Fail-fast: returns the first error, otherwise all terms in input order.
pub fn ontolius_to_simple_list(
    hpo: &FullCsrOntology,
    term_ids: &[TermId],
) -> Result<Vec<SimpleOntologyTerm>, String>
{
    term_ids
        .iter()
        .map(|id| ontolius_to_simple(hpo, id))
        .collect()
}