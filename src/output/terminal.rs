use super::Table;

pub fn format_section_title(title: &str) -> String {
    let inner = format!(" {} ", title);
    let line = "═".repeat(inner.chars().count());
    format!("╔{}╗\n║{}║\n╚{}╝", line, inner, line)
}

pub fn format_box_table(table: &Table) -> String {
    let mut widths = table
        .headers
        .iter()
        .map(|header| multiline_display_width(header))
        .collect::<Vec<_>>();
    for row in &table.rows {
        for (index, value) in row.iter().enumerate() {
            widths[index] = widths[index].max(multiline_display_width(value));
        }
    }

    let top = format_border('┌', '┬', '┐', &widths);
    let mid = format_border('├', '┼', '┤', &widths);
    let bottom = format_border('└', '┴', '┘', &widths);
    let header = format_row(&table.headers, &widths);
    let body = table
        .rows
        .iter()
        .map(|row| format_row(row, &widths))
        .collect::<Vec<_>>();

    std::iter::once(top)
        .chain(std::iter::once(header))
        .chain(std::iter::once(mid))
        .chain(body)
        .chain(std::iter::once(bottom))
        .collect::<Vec<_>>()
        .join("\n")
}

fn format_border(left: char, middle: char, right: char, widths: &[usize]) -> String {
    let cells = widths
        .iter()
        .map(|width| "─".repeat(*width + 2))
        .collect::<Vec<_>>()
        .join(&middle.to_string());
    format!("{left}{cells}{right}")
}

fn format_row(values: &[String], widths: &[usize]) -> String {
    let split_values = values
        .iter()
        .map(|value| value.lines().map(|line| line.to_string()).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let row_height = split_values
        .iter()
        .map(|lines| lines.len())
        .max()
        .unwrap_or(1);

    (0..row_height)
        .map(|line_index| {
            let cells = split_values
                .iter()
                .enumerate()
                .map(|(index, lines)| {
                    let value = lines.get(line_index).map(String::as_str).unwrap_or("");
                    let pad =
                        widths[index] + value.chars().count().saturating_sub(display_width(value));
                    format!(" {:<pad$} ", value, pad = pad)
                })
                .collect::<Vec<_>>()
                .join("│");
            format!("│{cells}│")
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn display_width(text: &str) -> usize {
    text.chars()
        .map(|ch| match ch {
            '\u{221E}' => 2,
            _ => 1,
        })
        .sum()
}

fn multiline_display_width(text: &str) -> usize {
    text.lines().map(display_width).max().unwrap_or(0)
}
