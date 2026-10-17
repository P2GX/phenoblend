//! dto
//! Module with structures for representing and analysing the melded phenotype data
//! Three main analysis types are supported.

mod constants;
mod overlap_dto;
mod phenoblend_analysis;
mod spread_dto;
mod upset_dto;

pub(crate) use constants::{ROOT_PHENOTYPE_TERM_ID, SIMPLE_ROOT_TERM};

pub use overlap_dto::{GeneMatch, TermMatches, OverlapPlotItem, OverlapPlotPayload, OverlapData};
pub use spread_dto::{assign_to_categories, CategoryAssignment, SpreadPlotCategory, SpreadPlotPayload,SeriesSummary, SeriesTerms};
pub use upset_dto::{UpsetGene, UpsetIntersection, UpsetPlotPayload};
pub use phenoblend_analysis::PhenoBlendAnalysis;