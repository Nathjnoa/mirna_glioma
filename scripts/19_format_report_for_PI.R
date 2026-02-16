#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# 19_format_report_for_PI.R
# Convert survGSEA top12 report to user-friendly Excel format for PI
#
# Features:
# - Organized tabs by database and enrichment direction
# - Color-coded Q-values (gradient: green = significant, red = not significant)
# - Auto-filters on all columns
# - Glossary tab with explanations
# - Summary tab with overview statistics
# -----------------------------------------------------------------------------

library(openxlsx)

# --- Parameters ---
input_csv <- "results_29s/tables/reports/GSEA_survGSEA_top12_annotated_20260215_181921.csv"
output_xlsx <- "results_29s/tables/reports/SurvGSEA_top12_formatted_PI.xlsx"

cat("=== Formatting survGSEA Report for PI ===\n")
cat("Input:", input_csv, "\n")
cat("Output:", output_xlsx, "\n\n")

# --- Read data ---
df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
cat("Loaded", nrow(df), "pathways\n\n")

# --- Parse miRNAs column to extract top 3 ---
format_top_mirnas <- function(mirna_col) {
  # Extract miRNA names (before the #)
  mirnas <- strsplit(as.character(mirna_col), ";")[[1]]
  mirnas <- trimws(mirnas)
  mirnas <- sub(" \\(#.*", "", mirnas)  # Remove (# and everything after)

  # Get top 3
  top3 <- head(mirnas, 3)
  paste(top3, collapse = "; ")
}

df$Top_3_miRNAs <- sapply(df$miRNAs_with_rank, format_top_mirnas)

# --- Add relevance stars based on Q-value ---
df$Relevancia <- ifelse(df$Q_value < 0.001, "⭐⭐⭐",
                 ifelse(df$Q_value < 0.01, "⭐⭐",
                 ifelse(df$Q_value < 0.05, "⭐", "")))

# --- Translate enrichment to Spanish ---
df$Direccion <- ifelse(df$Enrichment == "depleted", "Peor pronóstico", "Mejor pronóstico")

# --- Reorder columns for clarity ---
df_formatted <- df[, c("Category", "Subcategory", "Direccion", "Q_value", "Relevancia",
                       "Observed", "Top_3_miRNAs", "miRNAs_with_rank")]
colnames(df_formatted) <- c("Base de Datos", "Vía Biológica", "Asociación Pronóstico",
                            "Q-value (Significancia)", "Relevancia", "Nº miRNAs",
                            "Top 3 miRNAs", "Todos los miRNAs (con ranking)")

# --- Create workbook ---
wb <- createWorkbook()

# --- Tab 1: Summary ---
addWorksheet(wb, "Resumen")

summary_text <- c(
  "RESUMEN - Análisis de Supervivencia y miRNAs en Glioma",
  "",
  paste("Total de vías analizadas:", nrow(df)),
  paste("Vías asociadas a PEOR pronóstico:", sum(df$Enrichment == "depleted")),
  paste("Vías asociadas a MEJOR pronóstico:", sum(df$Enrichment == "enriched")),
  "",
  "Distribución por base de datos:",
  paste("  - KEGG:", sum(grepl("KEGG", df$Category))),
  paste("  - Reactome:", sum(grepl("Reactome", df$Category))),
  paste("  - GO Biological Process:", sum(grepl("GO.*Biological", df$Category))),
  "",
  "INTERPRETACIÓN:",
  "• Peor pronóstico = miRNAs con niveles DEPLETED (bajos) en estos pathways",
  "• Mejor pronóstico = miRNAs con niveles ENRICHED (altos) en estos pathways",
  "• Q-value < 0.05 = resultado estadísticamente significativo",
  "• ⭐⭐⭐ = altamente significativo (Q < 0.001)",
  "",
  "Navegación:",
  "• Ver pestañas abajo para resultados por categoría",
  "• Usar filtros (flecha en encabezados) para buscar vías específicas",
  "• Ver pestaña 'Glosario' para explicación de columnas"
)

writeData(wb, "Resumen", summary_text, startCol = 1, startRow = 1)

# Format summary
summary_style <- createStyle(fontSize = 12, fontColour = "#000000")
addStyle(wb, "Resumen", summary_style, rows = 1:length(summary_text), cols = 1, gridExpand = TRUE)

# Bold title
title_style <- createStyle(fontSize = 14, fontColour = "#000000", textDecoration = "bold")
addStyle(wb, "Resumen", title_style, rows = 1, cols = 1)

# --- Tab 2-7: Data tables by category ---
categories <- list(
  "Peor Pronóstico - KEGG" = df_formatted[df_formatted$`Asociación Pronóstico` == "Peor pronóstico" &
                                           grepl("KEGG", df_formatted$`Base de Datos`), ],
  "Peor Pronóstico - Reactome" = df_formatted[df_formatted$`Asociación Pronóstico` == "Peor pronóstico" &
                                                grepl("Reactome", df_formatted$`Base de Datos`), ],
  "Peor Pronóstico - GO BP" = df_formatted[df_formatted$`Asociación Pronóstico` == "Peor pronóstico" &
                                             grepl("GO.*Biological", df_formatted$`Base de Datos`), ],
  "Mejor Pronóstico - KEGG" = df_formatted[df_formatted$`Asociación Pronóstico` == "Mejor pronóstico" &
                                            grepl("KEGG", df_formatted$`Base de Datos`), ],
  "Mejor Pronóstico - Reactome" = df_formatted[df_formatted$`Asociación Pronóstico` == "Mejor pronóstico" &
                                                 grepl("Reactome", df_formatted$`Base de Datos`), ],
  "Mejor Pronóstico - GO BP" = df_formatted[df_formatted$`Asociación Pronóstico` == "Mejor pronóstico" &
                                              grepl("GO.*Biological", df_formatted$`Base de Datos`), ]
)

for (sheet_name in names(categories)) {
  df_subset <- categories[[sheet_name]]

  if (nrow(df_subset) == 0) next

  # Order by Q-value (most significant first)
  df_subset <- df_subset[order(df_subset$`Q-value (Significancia)`), ]

  # Create sheet
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df_subset, startRow = 1, startCol = 1, headerStyle = createStyle(textDecoration = "bold", fgFill = "#4F81BD", fontColour = "#FFFFFF"))

  # Add filters
  addFilter(wb, sheet_name, rows = 1, cols = 1:ncol(df_subset))

  # Format Q-value column with conditional formatting (gradient)
  qval_col <- which(colnames(df_subset) == "Q-value (Significancia)")

  # Green for low Q-values (significant), red for high
  conditionalFormatting(wb, sheet_name,
                        cols = qval_col,
                        rows = 2:(nrow(df_subset) + 1),
                        type = "colourScale",
                        style = c("#63BE7B", "#FFEB84", "#F8696B"),
                        rule = c(0, 0.025, 0.05))

  # Bold highly significant rows (Q < 0.001)
  sig_rows <- which(df_subset$`Q-value (Significancia)` < 0.001) + 1  # +1 for header
  if (length(sig_rows) > 0) {
    bold_style <- createStyle(textDecoration = "bold")
    addStyle(wb, sheet_name, bold_style, rows = sig_rows, cols = 1:ncol(df_subset), gridExpand = TRUE)
  }

  # Auto-width columns
  setColWidths(wb, sheet_name, cols = 1:ncol(df_subset), widths = "auto")

  # Freeze first row
  freezePane(wb, sheet_name, firstRow = TRUE)
}

# --- Tab 8: miRNA Lookup Table (for easy searching) ---
addWorksheet(wb, "Búsqueda de miRNAs")

cat("\nCreating miRNA lookup table...\n")

# Parse all miRNAs and create long-format table
mirna_lookup <- data.frame()

for (i in 1:nrow(df_formatted)) {
  pathway <- df_formatted[i, "Vía Biológica"]
  database <- df_formatted[i, "Base de Datos"]
  direction <- df_formatted[i, "Asociación Pronóstico"]
  qvalue <- df_formatted[i, "Q-value (Significancia)"]
  relevance <- df_formatted[i, "Relevancia"]

  # Parse miRNAs from original column
  mirna_string <- df[i, "miRNAs_with_rank"]
  mirnas <- strsplit(as.character(mirna_string), ";")[[1]]
  mirnas <- trimws(mirnas)

  # Extract miRNA name and z-score
  for (mirna in mirnas) {
    # Extract name (before #)
    mirna_name <- sub(" \\(#.*", "", mirna)

    # Extract rank and z-score
    rank_match <- regmatches(mirna, regexpr("#[0-9]+", mirna))
    rank <- ifelse(length(rank_match) > 0, gsub("#", "", rank_match), NA)

    z_match <- regmatches(mirna, regexpr("z=[-0-9.]+", mirna))
    z_score <- ifelse(length(z_match) > 0, gsub("z=", "", z_match), NA)

    # Add to lookup table
    mirna_lookup <- rbind(mirna_lookup, data.frame(
      miRNA = mirna_name,
      Ranking = as.numeric(rank),
      `Z-score` = as.numeric(z_score),
      `Vía Biológica` = pathway,
      `Base de Datos` = database,
      `Asociación Pronóstico` = direction,
      `Q-value` = qvalue,
      Relevancia = relevance,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  }
}

cat("  Total filas en tabla de búsqueda:", nrow(mirna_lookup), "\n")

# Sort by miRNA name for easy lookup
mirna_lookup <- mirna_lookup[order(mirna_lookup$miRNA), ]

# Write to sheet
writeData(wb, "Búsqueda de miRNAs", mirna_lookup, startRow = 1, startCol = 1,
          headerStyle = createStyle(textDecoration = "bold", fgFill = "#4F81BD", fontColour = "#FFFFFF"))

# Add filters
addFilter(wb, "Búsqueda de miRNAs", rows = 1, cols = 1:ncol(mirna_lookup))

# Conditional formatting on Q-value
qval_col_lookup <- which(colnames(mirna_lookup) == "Q-value")
conditionalFormatting(wb, "Búsqueda de miRNAs",
                      cols = qval_col_lookup,
                      rows = 2:(nrow(mirna_lookup) + 1),
                      type = "colourScale",
                      style = c("#63BE7B", "#FFEB84", "#F8696B"),
                      rule = c(0, 0.025, 0.05))

# Auto-width columns
setColWidths(wb, "Búsqueda de miRNAs", cols = 1:ncol(mirna_lookup), widths = "auto")

# Freeze first row
freezePane(wb, "Búsqueda de miRNAs", firstRow = TRUE)

# Add instruction at top
instruction_text <- "USO: Use Ctrl+F para buscar cualquier miRNA. Use los filtros (flechas en encabezados) para filtrar por vía, base de datos o dirección."
writeData(wb, "Búsqueda de miRNAs", instruction_text, startRow = 1, startCol = ncol(mirna_lookup) + 2)
instruction_style <- createStyle(fontSize = 10, fontColour = "#0066CC", textDecoration = "italic")
addStyle(wb, "Búsqueda de miRNAs", instruction_style, rows = 1, cols = ncol(mirna_lookup) + 2)

# --- Tab 9: Glossary ---
addWorksheet(wb, "Glosario")

glossary_text <- c(
  "GLOSARIO - Explicación de Columnas",
  "",
  "Base de Datos:",
  "  Origen de la vía biológica (KEGG, Reactome, o Gene Ontology Biological Process)",
  "",
  "Vía Biológica:",
  "  Nombre del proceso biológico o pathway analizado",
  "",
  "Asociación Pronóstico:",
  "  • Peor pronóstico = miRNAs con niveles bajos (depleted) están asociados a peor sobrevida",
  "  • Mejor pronóstico = miRNAs con niveles altos (enriched) están asociados a mejor sobrevida",
  "",
  "Q-value (Significancia):",
  "  Valor estadístico ajustado por comparaciones múltiples",
  "  • Menor valor = más significativo",
  "  • Q < 0.05 = significativo",
  "  • Q < 0.01 = muy significativo",
  "  • Q < 0.001 = altamente significativo",
  "  • Columna tiene gradiente de color: verde = significativo, rojo = no significativo",
  "",
  "Relevancia:",
  "  Sistema de estrellas basado en Q-value:",
  "  • ⭐⭐⭐ = Q < 0.001 (altamente significativo)",
  "  • ⭐⭐ = Q < 0.01 (muy significativo)",
  "  • ⭐ = Q < 0.05 (significativo)",
  "",
  "Nº miRNAs:",
  "  Cantidad de miRNAs que participan en esta vía biológica",
  "",
  "Top 3 miRNAs:",
  "  Los 3 miRNAs más importantes en este pathway",
  "",
  "Todos los miRNAs (con ranking):",
  "  Lista completa de miRNAs con su posición en el ranking de supervivencia",
  "  Formato: nombre (#posición, z=valor)",
  "  • z positivo = asociado a peor pronóstico",
  "  • z negativo = asociado a mejor pronóstico"
)

writeData(wb, "Glosario", glossary_text, startCol = 1, startRow = 1)

# Format glossary
glossary_style <- createStyle(fontSize = 11, fontColour = "#000000")
addStyle(wb, "Glosario", glossary_style, rows = 1:length(glossary_text), cols = 1, gridExpand = TRUE)

# Bold title and headers
glossary_title_style <- createStyle(fontSize = 14, fontColour = "#000000", textDecoration = "bold")
addStyle(wb, "Glosario", glossary_title_style, rows = 1, cols = 1)

header_rows <- grep("^[A-Z].+:$", glossary_text)
header_style <- createStyle(textDecoration = "bold", fontColour = "#4F81BD")
for (row in header_rows) {
  addStyle(wb, "Glosario", header_style, rows = row, cols = 1)
}

# Adjust column width
setColWidths(wb, "Glosario", cols = 1, widths = 100)

# --- Save workbook ---
saveWorkbook(wb, output_xlsx, overwrite = TRUE)

cat("\n✓ Excel formateado guardado en:", output_xlsx, "\n")
cat("\nPestañas creadas:\n")
cat("  1. Resumen (overview)\n")
cat("  2-4. Peor Pronóstico (KEGG, Reactome, GO BP)\n")
cat("  5-7. Mejor Pronóstico (KEGG, Reactome, GO BP)\n")
cat("  8. Búsqueda de miRNAs (tabla long para Ctrl+F) ⭐ NUEVO\n")
cat("  9. Glosario (explicaciones)\n")
cat("\nCaracterísticas:\n")
cat("  ✓ Ordenado por significancia (Q-value)\n")
cat("  ✓ Gradiente de color en Q-values\n")
cat("  ✓ Filtros activos en todas las columnas\n")
cat("  ✓ Sistema de estrellas de relevancia\n")
cat("  ✓ Top 3 miRNAs extraídos\n")
cat("  ✓ Filas altamente significativas en negritas\n")
cat("  ✓ Cada miRNA en celda separada (pestaña Búsqueda) para Ctrl+F fácil\n")
