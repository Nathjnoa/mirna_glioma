#!/usr/bin/env Rscript
# ============================================================================
# Reduce GO term redundancy in survGSEA results using rrvgo
# - Maps GO term names → GO IDs via GO.db
# - Applies rrvgo semantic similarity reduction (GO BP, GO MF separately)
# - Unmapped terms are kept as-is
# - KEGG and Reactome pass through without modification
# - Regenerates bubble plots from reduced table
# ============================================================================
options(stringsAsFactors = FALSE)

# ---- CLI argument parsing ----
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[[hit + 1]]
}
as_num_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}
safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# ---- CLI defaults ----
in_root   <- get_arg("--in_root", file.path("results", "tables", "SurvivalRank_GSEA"))
out_root  <- get_arg("--out_root", file.path("results", "figures", "SurvivalRank_GSEA"))
tbl_root  <- get_arg("--tbl_root", file.path("results", "tables", "SurvivalRank_GSEA"))
run_tag   <- get_arg("--run_tag", "SurvivalRank_CoxZ_miRPathDB")
plot_label <- get_arg("--label", run_tag)

sim_threshold <- as_num_safe(get_arg("--sim_threshold", "0.9"))
sim_method    <- get_arg("--sim_method", "Rel")  # Rel, Lin, Wang
orgdb         <- get_arg("--orgdb", "org.Hs.eg.db")

n_enriched <- as.integer(get_arg("--n_enriched", "12"))
n_depleted <- as.integer(get_arg("--n_depleted", "12"))
wrap_width <- as.integer(get_arg("--wrap_width", "58"))
preset     <- get_arg("--preset", "double_col")
stamp      <- tolower(get_arg("--stamp", "true")) %in% c("true","t","1","yes","y")
seed       <- as.integer(get_arg("--seed", "42"))

# ---- Resolve paths ----
get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}
script_dir <- get_script_dir()
proj_root  <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
resolve_path <- function(p) {
  if (file.exists(p) || dir.exists(p)) return(p)
  alt <- file.path(proj_root, p)
  if (file.exists(alt) || dir.exists(alt)) return(alt)
  p
}
in_root  <- resolve_path(in_root)
out_root <- resolve_path(out_root)
tbl_root <- resolve_path(tbl_root)

# ---- Logging ----
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
out_tbl_dir <- file.path(tbl_root, run_tag, "reduced_rrvgo")
out_fig_dir <- file.path(out_root, run_tag, "reduced_rrvgo")
dir.create(out_tbl_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

log_fp <- file.path("logs", paste0("survGSEA_reduce_redundancy_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")
on.exit({ sink(type = "message"); sink(type = "output"); close(zz) }, add = TRUE)

cat("============================================\n")
cat("SurvivalRank GSEA - Reduce GO Redundancy\n")
cat("============================================\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("Input root:", in_root, "\n")
cat("Table output:", out_tbl_dir, "\n")
cat("Figure output:", out_fig_dir, "\n")
cat("Run tag:", run_tag, "\n")
cat("Similarity threshold:", sim_threshold, "\n")
cat("Similarity method:", sim_method, "\n")
cat("OrgDb:", orgdb, "\n")
cat("N enriched:", n_enriched, " | N depleted:", n_depleted, "\n")
cat("Preset:", preset, "\n")
cat("============================================\n\n")

suppressPackageStartupMessages({
  library(GO.db)
  library(rrvgo)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(patchwork)
})

# Load orgdb
suppressPackageStartupMessages(library(orgdb, character.only = TRUE))

set.seed(seed)

# ===========================================================================
# PART 1: Load miEAA results
# ===========================================================================
find_latest <- function(root, run_tag, pattern) {
  dir_path <- file.path(root, run_tag)
  if (!dir.exists(dir_path)) return(NA_character_)
  ff <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(ff) == 0) return(NA_character_)
  ff[order(file.info(ff)$mtime, decreasing = TRUE)][1]
}

use_fp <- find_latest(in_root, run_tag, "^miEAA_GSEA_all_miRPathDB_expert_.*\\.tsv$")
if (is.na(use_fp) || !file.exists(use_fp)) stop("No miEAA results found in: ", file.path(in_root, run_tag))
cat("Input file:", use_fp, "\n\n")

df <- read.delim(use_fp, check.names = FALSE)
names(df) <- trimws(names(df))
if (nrow(df) == 0) stop("Empty miEAA table.")

# Standardize columns
df$pval      <- as_num_safe(df[["P-value"]])
df$padj      <- as_num_safe(df[["P-adjusted"]])
df$qval      <- as_num_safe(df[["Q-value"]])
df$observed  <- as_num_safe(df[["Observed"]])
df$term      <- str_squish(as.character(df[["Subcategory"]]))
df$category  <- as.character(df[["Category"]])
df$enrichment <- tolower(as.character(df[["Enrichment"]]))

cat("Total terms loaded:", nrow(df), "\n")
cat("Enriched:", sum(df$enrichment == "enriched"), " | Depleted:", sum(df$enrichment == "depleted"), "\n\n")

# ===========================================================================
# PART 2: Map GO term names to GO IDs
# ===========================================================================
cat("--- Mapping GO term names to GO IDs ---\n")

# Build lookup: GO term name → GO ID from GO.db
go_terms_db <- AnnotationDbi::select(GO.db, keys = keys(GO.db, "GOID"),
                                      columns = c("GOID", "TERM", "ONTOLOGY"),
                                      keytype = "GOID")

# Build normalized name lookup (lowercase, trimmed)
go_terms_db$term_norm <- tolower(str_squish(go_terms_db$TERM))

# Identify GO rows in our data
is_go_bp <- grepl("GO.*Biological.process", df$category, ignore.case = TRUE)
is_go_mf <- grepl("GO.*Molecular.function", df$category, ignore.case = TRUE)
is_go <- is_go_bp | is_go_mf

# Map by exact match on normalized name
df$go_id <- NA_character_
df$go_ont <- NA_character_

if (any(is_go)) {
  terms_to_map <- tolower(str_squish(df$term[is_go]))
  match_idx <- match(terms_to_map, go_terms_db$term_norm)

  mapped <- !is.na(match_idx)
  df$go_id[is_go][mapped]  <- go_terms_db$GOID[match_idx[mapped]]
  df$go_ont[is_go][mapped] <- go_terms_db$ONTOLOGY[match_idx[mapped]]

  n_go <- sum(is_go)
  n_mapped <- sum(mapped)
  n_unmapped <- n_go - n_mapped

  cat("GO terms total:", n_go, "\n")
  cat("Mapped to GO ID:", n_mapped, "(", round(100 * n_mapped / n_go, 1), "%)\n")
  cat("Unmapped:", n_unmapped, "(will be kept without reduction)\n")

  # Export unmapped terms for reference
  if (n_unmapped > 0) {
    unmapped_terms <- unique(df$term[is_go][!mapped])
    unmapped_fp <- file.path(out_tbl_dir, paste0("GO_unmapped_terms_", ts, ".txt"))
    writeLines(unmapped_terms, unmapped_fp)
    cat("Unmapped terms list:", unmapped_fp, "\n")
  }
} else {
  cat("WARN: No GO terms found in the data.\n")
}

# ===========================================================================
# PART 3: Apply rrvgo reduction
# ===========================================================================
cat("\n--- Applying rrvgo semantic similarity reduction ---\n")

# Track which rows to keep (TRUE = keep, FALSE = removed by rrvgo)
df$keep <- TRUE
df$rrvgo_parent <- NA_character_
df$rrvgo_cluster <- NA_integer_

reduce_go <- function(df_subset, ont_label, orgdb_obj) {
  go_ids <- df_subset$go_id
  has_id <- !is.na(go_ids)

  if (sum(has_id) < 3) {
    cat("  Too few mapped terms (", sum(has_id), ") for rrvgo; keeping all.\n")
    return(df_subset)
  }

  go_ids_valid <- go_ids[has_id]
  scores <- setNames(-log10(pmax(df_subset$padj[has_id], .Machine$double.xmin)), go_ids_valid)

  # Remove duplicates (same GO ID appearing multiple times)
  if (any(duplicated(go_ids_valid))) {
    dup_mask <- !duplicated(go_ids_valid)
    go_ids_valid <- go_ids_valid[dup_mask]
    scores <- scores[dup_mask]
  }

  sim_matrix <- tryCatch(
    calculateSimMatrix(go_ids_valid,
                       orgdb = orgdb_obj,
                       ont = ont_label,
                       method = sim_method),
    error = function(e) {
      cat("  ERROR in calculateSimMatrix:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(sim_matrix) || nrow(sim_matrix) < 2) {
    cat("  Sim matrix failed or too small; keeping all terms.\n")
    return(df_subset)
  }

  # Filter scores to only include GO IDs present in sim_matrix
  valid_in_matrix <- intersect(names(scores), rownames(sim_matrix))
  if (length(valid_in_matrix) < 2) {
    cat("  Too few valid terms in sim matrix; keeping all.\n")
    return(df_subset)
  }
  scores <- scores[valid_in_matrix]

  reduced <- tryCatch(
    reduceSimMatrix(sim_matrix,
                    scores = scores,
                    threshold = sim_threshold,
                    orgdb = orgdb_obj),
    error = function(e) {
      cat("  ERROR in reduceSimMatrix:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(reduced) || nrow(reduced) == 0) {
    cat("  Reduction failed; keeping all terms.\n")
    return(df_subset)
  }

  # Mark parent terms (representatives) to keep
  # reduced$parent contains the GO ID of the cluster representative
  # reduced$go contains the GO ID of each term
  go_to_parent <- setNames(reduced$parent, reduced$go)
  go_to_cluster <- setNames(reduced$cluster, reduced$go)

  # Update df_subset
  for (i in seq_len(nrow(df_subset))) {
    gid <- df_subset$go_id[i]
    if (is.na(gid)) next  # unmapped → keep
    if (gid %in% names(go_to_parent)) {
      df_subset$rrvgo_parent[i] <- go_to_parent[[gid]]
      df_subset$rrvgo_cluster[i] <- go_to_cluster[[gid]]
      # Keep only parent (representative) terms
      df_subset$keep[i] <- (gid == go_to_parent[[gid]])
    }
    # If gid not in reduced output (dropped by sim matrix), keep it
  }

  n_before <- sum(has_id)
  n_after <- sum(df_subset$keep[has_id])
  cat("  Mapped terms:", n_before, "→ Reduced:", n_after,
      "(removed", n_before - n_after, "redundant)\n")

  df_subset
}

# Process each combination of ontology × direction
for (ont in c("BP", "MF")) {
  ont_pattern <- if (ont == "BP") "GO.*Biological.process" else "GO.*Molecular.function"
  ont_mask <- grepl(ont_pattern, df$category, ignore.case = TRUE)

  for (dir in c("enriched", "depleted")) {
    mask <- ont_mask & df$enrichment == dir
    if (sum(mask) == 0) next

    cat(sprintf("\nProcessing GO %s - %s (n=%d)\n", ont, dir, sum(mask)))
    df[mask, ] <- reduce_go(df[mask, , drop = FALSE], ont, orgdb)
  }
}

# ===========================================================================
# PART 4: Build reduced table
# ===========================================================================
cat("\n--- Building reduced table ---\n")

df_reduced <- df[df$keep, , drop = FALSE]

cat("Before reduction:", nrow(df), "\n")
cat("After reduction:", nrow(df_reduced), "\n")
cat("Removed:", nrow(df) - nrow(df_reduced), "redundant GO terms\n\n")

# Summary by category × direction
summary_tbl <- df_reduced %>%
  group_by(category, enrichment) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(category, enrichment)

cat("Summary after reduction:\n")
print(as.data.frame(summary_tbl), row.names = FALSE)

# Export reduced table (keep original columns + rrvgo annotations)
export_cols <- c(names(read.delim(use_fp, check.names = FALSE, nrows = 1)),
                 "go_id", "rrvgo_parent", "rrvgo_cluster")
export_cols <- intersect(export_cols, names(df_reduced))

reduced_fp <- file.path(out_tbl_dir, paste0("miEAA_GSEA_reduced_rrvgo_", ts, ".tsv"))
write.table(df_reduced[, export_cols], reduced_fp, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nReduced table:", reduced_fp, "\n")

# Export full table with keep/cluster annotations for reference
full_annotated_fp <- file.path(out_tbl_dir, paste0("miEAA_GSEA_full_rrvgo_annotated_", ts, ".tsv"))
write.table(df, full_annotated_fp, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Full annotated table:", full_annotated_fp, "\n")

# Export summary as Markdown table
md_fp <- file.path(out_tbl_dir, paste0("reduction_summary_", ts, ".md"))
md_lines <- c(
  paste0("# rrvgo Reduction Summary"),
  "",
  paste0("- **Run tag**: ", plot_label),
  paste0("- **Similarity threshold**: ", sim_threshold),
  paste0("- **Method**: ", sim_method),
  paste0("- **Date**: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  paste0("- **Before reduction**: ", nrow(df)),
  paste0("- **After reduction**: ", nrow(df_reduced)),
  paste0("- **Removed**: ", nrow(df) - nrow(df_reduced), " redundant GO terms"),
  "",
  "## Terms per category",
  "",
  "| Category | Enriched | Depleted |",
  "|----|----:|----:|"
)
summary_wide <- summary_tbl %>%
  tidyr::pivot_wider(names_from = enrichment, values_from = n, values_fill = 0)
for (i in seq_len(nrow(summary_wide))) {
  row <- summary_wide[i, ]
  md_lines <- c(md_lines, sprintf("| %s | %d | %d |",
                                   row$category,
                                   if ("enriched" %in% names(row)) row$enriched else 0,
                                   if ("depleted" %in% names(row)) row$depleted else 0))
}
writeLines(md_lines, md_fp)
cat("Summary MD:", md_fp, "\n")

# ===========================================================================
# PART 5: Regenerate bubble plots from reduced table
# ===========================================================================
cat("\n--- Generating bubble plots from reduced data ---\n")

# ---- Presets ----
presets <- list(
  single_col = list(
    width_mm = 85, height_mm = 65, base = 8, title = 9, axis = 8,
    legend = 7, ticks = 7, line = 0.6, point = 1.6, spacing_mm = 2,
    margin_mm = 6, dpi_vector = 300, dpi_raster = 600
  ),
  double_col = list(
    width_mm = 220, height_mm = 160, base = 30, title = 35, axis = 30,
    legend = 28, ticks = 28, line = 0.8, point = 2.2, spacing_mm = 3,
    margin_mm = 8, dpi_vector = 300, dpi_raster = 600
  ),
  presentation = list(
    width_mm = 254, height_mm = 143, base = 18, title = 22, axis = 18,
    legend = 16, ticks = 16, line = 1.5, point = 4, spacing_mm = 6,
    margin_mm = 10, dpi_vector = 300, dpi_raster = 300
  ),
  poster = list(
    width_mm = 508, height_mm = 356, base = 28, title = 34, axis = 28,
    legend = 24, ticks = 24, line = 2.0, point = 6, spacing_mm = 10,
    margin_mm = 15, dpi_vector = 300, dpi_raster = 300
  ),
  per_database = list(
    width_mm = 220, height_mm = 160, base = 22, title = 26, axis = 22,
    legend = 20, ticks = 20, line = 0.8, point = 2.2, spacing_mm = 3,
    margin_mm = 8, dpi_vector = 300, dpi_raster = 600
  ),
  enriched_depleted = list(
    width_mm = 220, height_mm = 160, base = 20, title = 24, axis = 20,
    legend = 18, ticks = 18, line = 0.8, point = 2.2, spacing_mm = 3,
    margin_mm = 8, dpi_vector = 300, dpi_raster = 600
  ),
  multipanel = list(
    width_mm = 220, height_mm = 160, base = 16, title = 20, axis = 16,
    legend = 14, ticks = 14, line = 0.8, point = 2.2, spacing_mm = 3,
    margin_mm = 8, dpi_vector = 300, dpi_raster = 600
  )
)
if (!preset %in% names(presets)) stop("Invalid preset: ", preset)
p <- presets[[preset]]

# ---- Theme ----
theme_pub <- function(preset_params = p) {
  theme_minimal(base_size = preset_params$base, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = preset_params$title, face = "bold", hjust = 0),
      axis.title = element_text(size = preset_params$axis, face = "plain"),
      axis.text.x = element_text(size = preset_params$ticks, color = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = preset_params$axis, color = "black", lineheight = 1.1,
                                  margin = margin(r = 3, unit = "mm")),
      strip.text = element_text(size = preset_params$axis, face = "bold"),
      legend.title = element_text(size = preset_params$legend, face = "bold"),
      legend.text = element_text(size = preset_params$legend),
      legend.position = "right",
      legend.key.size = unit(4, "mm"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.border = element_blank(),
      panel.spacing.x = unit(preset_params$spacing_mm, "mm"),
      panel.spacing.y = unit(preset_params$spacing_mm * 2.4, "mm"),
      plot.margin = margin(preset_params$margin_mm, preset_params$margin_mm, preset_params$margin_mm, preset_params$margin_mm + 6, unit = "mm"),
      strip.background = element_rect(fill = "grey95", color = NA)
    )
}

colors_direction <- c("Enriched" = "#009E73", "Depleted" = "#D55E00")

# ---- Prepare plot data ----
map_db <- function(x) {
  x <- tolower(x)
  db <- rep(NA_character_, length(x))
  db[grepl("go.*biological.process", x)] <- "GO Biological process (miRPathDB)"
  db[grepl("go.*molecular.function", x)] <- "GO Molecular function (miRPathDB)"
  db[grepl("kegg", x)] <- "KEGG (miRPathDB)"
  db[grepl("reactome", x)] <- "Reactome (miRPathDB)"
  db
}

df_plot <- df_reduced
df_plot$db <- map_db(df_plot$category)
df_plot$direction <- ifelse(df_plot$enrichment == "enriched", "Enriched", "Depleted")
df_plot$direction <- factor(df_plot$direction, levels = c("Enriched", "Depleted"))
df_plot$logp <- -log10(pmax(df_plot$padj, .Machine$double.xmin))

mirpath_levels <- c(
  "GO Biological process (miRPathDB)",
  "KEGG (miRPathDB)",
  "Reactome (miRPathDB)"
)
df_plot <- df_plot[df_plot$db %in% mirpath_levels & is.finite(df_plot$padj) & nzchar(df_plot$term), ]
df_plot$db <- factor(df_plot$db, levels = mirpath_levels)

# ---- Save function ----
save_plot <- function(plot, out_base, n_rows = 1, width_mult = 1, height_mult = 1) {
  width_mm  <- p$width_mm * width_mult
  height_mm <- p$height_mm * height_mult * n_rows
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else pdf
  pdf_path <- paste0(out_base, ".pdf")
  tryCatch({
    ggsave(pdf_path, plot, width = width_mm, height = height_mm,
           units = "mm", dpi = p$dpi_vector, limitsize = FALSE, device = pdf_dev)
    cat("  Exported:", pdf_path, sprintf("(%d x %d mm)\n", round(width_mm), round(height_mm)))
  }, error = function(e) warning("Failed PDF: ", conditionMessage(e)))
  svg_path <- paste0(out_base, ".svg")
  if (requireNamespace("svglite", quietly = TRUE)) {
    tryCatch({
      ggsave(svg_path, plot, width = width_mm, height = height_mm,
             units = "mm", dpi = p$dpi_vector, limitsize = FALSE, device = svglite::svglite)
      cat("  Exported:", svg_path, "\n")
    }, error = function(e) warning("Failed SVG: ", conditionMessage(e)))
  }
  png_path <- paste0(out_base, ".png")
  tryCatch({
    ggsave(png_path, plot, width = width_mm, height = height_mm,
           units = "mm", dpi = p$dpi_raster, limitsize = FALSE, device = "png")
    cat("  Exported:", png_path, sprintf("(%d dpi)\n", p$dpi_raster))
  }, error = function(e) warning("Failed PNG: ", conditionMessage(e)))
  invisible(NULL)
}

# ---- Create output subdirectories ----
dir.create(file.path(out_fig_dir, "enriched_depleted"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_fig_dir, "per_database"), showWarnings = FALSE, recursive = TRUE)

# ---- Generate plots ----
plots_list <- list()
plots_per_db <- list()

for (dir_type in c("Enriched", "Depleted")) {
  n_terms <- if (dir_type == "Depleted") n_depleted else n_enriched

  df_dir <- df_plot %>%
    filter(direction == dir_type) %>%
    group_by(db) %>%
    arrange(padj, .by_group = TRUE) %>%
    slice_head(n = n_terms) %>%
    ungroup()

  dir_counts <- table(df_dir$db)
  cat(sprintf("  %s terms (reduced): %s\n", dir_type,
              paste(names(dir_counts), dir_counts, sep = "=", collapse = ", ")))

  if (nrow(df_dir) == 0) {
    cat(sprintf("  SKIP %s: no terms available\n", dir_type))
    next
  }

  # Order terms within each DB
  df_dir <- df_dir %>%
    arrange(db, padj) %>%
    mutate(term_plot = paste(term, db, sep = "|||"))
  df_dir$term_plot <- factor(df_dir$term_plot, levels = rev(unique(df_dir$term_plot)))

  # Use enriched_depleted preset for multipanel enriched/depleted plots
  p_ed <- presets[["enriched_depleted"]]
  size_range_ed <- c(p_ed$point * 0.8, p_ed$point * 2.8)
  dir_color <- colors_direction[[dir_type]]

  # Multipanel plot (all databases)
  p_dir <- ggplot(df_dir, aes(x = logp, y = term_plot)) +
    geom_point(aes(size = observed), alpha = 0.85, color = dir_color) +
    facet_wrap(~db, ncol = 3, drop = FALSE, scales = "free_y",
               labeller = labeller(db = label_wrap_gen(width = 16))) +
    scale_y_discrete(labels = function(x) str_wrap(gsub("\\|\\|\\|.*$", "", x), width = wrap_width),
                     expand = expansion(mult = 0.06, add = 1.2)) +
    scale_size(range = size_range_ed, name = "Observed") +
    labs(
      x = "-log10(P-adjusted)",
      y = NULL,
      title = paste0("miRPathDB ", dir_type, ": ", plot_label)
    ) +
    theme_pub(p_ed)

  # Save multipanel plot in enriched_depleted/
  dir_tag <- tolower(dir_type)
  out_base <- file.path(out_fig_dir, "enriched_depleted",
                        paste0("SurvivalRank_Bubble_miRPathDB_reduced_", dir_tag, "_",
                               safe_name(plot_label), if (stamp) paste0("_", ts) else ""))
  save_plot(p_dir, out_base, n_rows = 1, width_mult = 4.2, height_mult = 2.2)

  plots_list[[dir_type]] <- p_dir

  # ---- Per-database plots ----
  # Use per_database preset for individual database plots
  p_db_preset <- presets[["per_database"]]
  size_range_db <- c(p_db_preset$point * 0.8, p_db_preset$point * 2.8)

  for (db_name in unique(df_dir$db)) {
    df_db <- df_dir %>% filter(db == db_name)
    if (nrow(df_db) == 0) next

    # Re-order for single DB
    df_db$term_plot <- factor(df_db$term, levels = rev(unique(df_db$term)))

    p_db <- ggplot(df_db, aes(x = logp, y = term_plot)) +
      geom_point(aes(size = observed), alpha = 0.85, color = dir_color) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width),
                       expand = expansion(mult = 0.06, add = 1.0)) +
      scale_size(range = size_range_db, name = "Observed") +
      labs(
        x = "-log10(P-adjusted)",
        y = NULL,
        title = paste0(db_name, " - ", dir_type)
      ) +
      theme_pub(p_db_preset)

    # Save per-database plot
    db_safe <- safe_name(db_name)
    out_base_db <- file.path(out_fig_dir, "per_database",
                             paste0("SurvivalRank_", db_safe, "_", dir_tag, "_",
                                    safe_name(plot_label), if (stamp) paste0("_", ts) else ""))
    save_plot(p_db, out_base_db, n_rows = 1, width_mult = 2.0, height_mult = max(1.2, nrow(df_db) * 0.24))

    # Store for potential use
    key <- paste(dir_type, db_name, sep = "_")
    plots_per_db[[key]] <- p_db
  }
}


cat("\n============================================\n")
cat("DONE\n")
cat("Reduced table:", reduced_fp, "\n")
cat("Full annotated:", full_annotated_fp, "\n")
cat("Figures:", out_fig_dir, "\n")
cat("Log:", log_fp, "\n")
