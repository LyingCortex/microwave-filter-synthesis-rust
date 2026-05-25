use super::Table;

pub fn format_markdown_table(table: &Table) -> String {
    let mut lines = Vec::with_capacity(table.rows.len() + 2);
    lines.push(format!("| {} |", table.headers.join(" | ")));
    lines.push(format!(
        "| {} |",
        table.headers
            .iter()
            .map(|_| "---")
            .collect::<Vec<_>>()
            .join(" | ")
    ));
    for row in &table.rows {
        lines.push(format!("| {} |", row.join(" | ")));
    }
    lines.join("\n")
}
