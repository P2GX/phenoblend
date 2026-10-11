use rust_xlsxwriter::{Workbook, Format, FormatAlign};
use crate::blend::dto::{OverlapPlotPayload, UpsetPlotPayload}; // adjust to your actual module path

pub fn export_upset_plot_to_xlsx(
    payload: &UpsetPlotPayload,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut workbook = Workbook::new();
    let worksheet = workbook.add_worksheet();

    let header_format = Format::new()
        .set_bold()
        .set_background_color("#D9E1F2")
        .set_align(FormatAlign::Center);

    // --- Header row: HPO Label, HPO Id, then one column per entity ---
    worksheet.write_string_with_format(0, 0, "HPO Label", &header_format)?;
    worksheet.write_string_with_format(0, 1, "HPO Id", &header_format)?;

    /* 
    for (i, entity) in payload.entities.iter().enumerate() {
        let col = (2 + i) as u16;
        worksheet.write_string_with_format(0, col, entity, &header_format)?;
    }

    // --- Data rows ---
    for (row_idx, item) in payload.columns.iter().enumerate() {
        let r = (row_idx + 1) as u32;

        worksheet.write_string(r, 0, &item.hpo_name)?;
        worksheet.write_string(r, 1, &item.hpo_id)?;

        for (i, entity) in payload.entities.iter().enumerate() {
            let col = (2 + i) as u16;
            // Missing entity in this row's scores map → 0.0, matching the
            // same fallback convention used elsewhere (e.g. compute_entity_sums).
            let score = item.scores.get(entity).copied().unwrap_or(0.0);
            worksheet.write_number(r, col, score)?;
        }
    }

    // Column widths: label/id get more room, score columns stay narrow
    worksheet.set_column_width(0, 35)?;
    worksheet.set_column_width(1, 12)?;
    for i in 0..payload.entities.len() {
        worksheet.set_column_width((2 + i) as u16, 14)?;
    }
    */
    workbook.save(output_path)?;
    Ok(())
}




pub fn export_overlap_plot_to_xlsx(
    payload: &OverlapPlotPayload,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut workbook = Workbook::new();
    let worksheet = workbook.add_worksheet();

    let header_format = Format::new()
        .set_bold()
        .set_background_color("#D9E1F2")
        .set_align(FormatAlign::Center);

    // --- Header row: HPO Label, HPO Id, then one column per entity ---
    worksheet.write_string_with_format(0, 0, "HPO Label", &header_format)?;
    worksheet.write_string_with_format(0, 1, "HPO Id", &header_format)?;

    for (i, entity) in payload.entities.iter().enumerate() {
        let col = (2 + i) as u16;
        worksheet.write_string_with_format(0, col, entity, &header_format)?;
    }

    // --- Data rows ---
    for (row_idx, item) in payload.columns.iter().enumerate() {
        let r = (row_idx + 1) as u32;

        worksheet.write_string(r, 0, &item.hpo_name)?;
        worksheet.write_string(r, 1, &item.hpo_id)?;

        for (i, entity) in payload.entities.iter().enumerate() {
            let col = (2 + i) as u16;
            // Missing entity in this row's scores map → 0.0, matching the
            // same fallback convention used elsewhere (e.g. compute_entity_sums).
            let score = item.scores.get(entity).copied().unwrap_or(0.0);
            worksheet.write_number(r, col, score)?;
        }
    }

    // Column widths: label/id get more room, score columns stay narrow
    worksheet.set_column_width(0, 35)?;
    worksheet.set_column_width(1, 12)?;
    for i in 0..payload.entities.len() {
        worksheet.set_column_width((2 + i) as u16, 14)?;
    }

    workbook.save(output_path)?;
    Ok(())
}