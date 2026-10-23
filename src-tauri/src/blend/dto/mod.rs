//! dto
//! Module with structures for representing and analysing the melded phenotype data
//! Three main analysis types are supported.

mod constants;
mod overlap_dto;
mod input;
mod phenoblend_analysis;
mod spread_dto;
mod upset_dto;

pub(crate) use constants::{ROOT_PHENOTYPE_TERM_ID, SIMPLE_ROOT_TERM};
pub use input::{ObservedHpoInput,severity_from_ids};
pub use overlap_dto::{GeneMatch, TermMatches, OverlapPlotItem, OverlapPlotPayload, OverlapData, Severity};
pub use spread_dto::{assign_to_categories, CategoryAssignment, SpreadPlotCategory, SpreadPlotPayload,SeriesSummary, SeriesTerms};
pub use upset_dto::{UpsetPlotInput, ProbandSegment, TermKey, GeneKey, GeneTerms, UpsetPlotData, UpsetRowData, UpsetColumnData, closure_of, UpsetSummary};
pub use phenoblend_analysis::PhenoBlendAnalysis;