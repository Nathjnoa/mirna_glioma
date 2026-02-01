#!/usr/bin/env Rscript
# ============================================================================
# miEAA GSEA Bubble Plots - Publication-ready figures
# Generates separate bubble plots for enriched vs depleted processes per DB
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

# ---- Input parameters with validation ----
spec_fp <- get_arg("--spec", file.path("config", "de_specs.csv"))
in_root <- get_arg("--in_root", file.path("results", "tables", "MiEAA_GSEA"))
out_root <- get_arg("--out_root", file.path("results", "figures", "MiEAA_GSEA_bubble"))
run_tag <- get_arg("--run_tag", "A_conservative")
n_mirpathdb <- as.integer(get_arg("--n_mirpathdb", "12"))
wrap_width <- as.integer(get_arg("--wrap_width", "40"))
preset <- get_arg("--preset", "double_col")
seed <- as.integer(get_arg("--seed", "42"))

# Input validation
if (!nzchar(run_tag)) stop("run_tag cannot be empty")
if (n_mirpathdb < 1) stop("n_mirpathdb must be >= 1")
if (wrap_width < 10) stop("wrap_width must be >= 10")

get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
}
script_dir <- get_script_dir()
proj_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
resolve_path <- function(p) {
  if (file.exists(p) || dir.exists(p)) return(p)
  alt <- file.path(proj_root, p)
  if (file.exists(alt) || dir.exists(alt)) return(alt)
  p
}
spec_fp <- resolve_path(spec_fp)
in_root <- resolve_path(in_root)
out_root <- resolve_path(out_root)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_root, run_tag), showWarnings = FALSE, recursive = TRUE)
log_fp <- file.path("logs", paste0("mieaa_gsea_bubble_", ts, ".txt"))
log_fp_copy <- file.path(out_root, run_tag, paste0("mieaa_gsea_bubble_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(zz)
  if (file.exists(log_fp)) {
    file.copy(log_fp, log_fp_copy, overwrite = TRUE)
  }
}, add = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(patchwork)  # For multi-panel combined figures
})

# ---- Publication-ready presets (explicit dimensions and styling) ----
presets <- list(
  single_col = list(
    width_mm = 85, height_mm = 65, base = 8, title = 9, axis = 8,
    legend = 7, ticks = 7, line = 0.6, point = 1.6, spacing_mm = 2,
    margin_mm = 6, dpi_vector = 300, dpi_raster = 600
  ),
  double_col = list(
    width_mm = 180, height_mm = 120, base = 8.5, title = 10, axis = 8.5,
    legend = 7.5, ticks = 7.5, line = 0.7, point = 1.8, spacing_mm = 2.5,
    margin_mm = 7, dpi_vector = 300, dpi_raster = 600
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
  )
)
if (!preset %in% names(presets)) stop("preset invalido: ", preset, " (esperado: ", paste(names(presets), collapse = ", "), ")")
p <- presets[[preset]]

# Set seed for reproducibility (affects jitter if used)
set.seed(seed)

# ---- Publication-ready theme (colorblind-safe, explicit sizing) ----
theme_pub <- function() {
  theme_minimal(base_size = p$base, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = p$title, face = "bold", hjust = 0),
      axis.title = element_text(size = p$axis, face = "plain"),
      axis.text.x = element_text(size = p$ticks, color = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = p$axis, color = "black"),  # Increased from ticks to axis for better readability
      strip.text = element_text(size = p$axis, face = "bold"),
      legend.title = element_text(size = p$legend, face = "bold"),
      legend.text = element_text(size = p$legend),
      legend.position = "right",
      legend.key.size = unit(4, "mm"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.border = element_blank(),
      panel.spacing = unit(p$spacing_mm, "mm"),
      plot.margin = margin(p$margin_mm, p$margin_mm, p$margin_mm, p$margin_mm + 6, unit = "mm"),
      strip.background = element_rect(fill = "grey95", color = NA)
    )
}

# ---- Colorblind-safe palette (Okabe-Ito inspired) ----
# Enriched: green (#009E73), Depleted: orange (#D55E00), Other: grey
colors_direction <- c(
  "Enriched" = "#009E73",
  "Depleted" = "#D55E00",
  "Other" = "#999999"
)

find_latest_all <- function(root, aid_fs, run_tag) {
  dir_path <- file.path(root, aid_fs, run_tag)
  if (!dir.exists(dir_path)) return(NA_character_)
  pat <- paste0("^MiEAA_GSEA_all_", aid_fs, "_", run_tag, "_.*\\.tsv$")
  ff <- list.files(dir_path, pattern = pat, full.names = TRUE)
  if (length(ff) == 0) return(NA_character_)
  ff[order(file.info(ff)$mtime, decreasing = TRUE)][1]
}

norm_name <- function(x) {
  x <- tolower(x)
  gsub("[^a-z0-9]+", "", x)
}
find_col <- function(df, patterns, exclude = character()) {
  nn <- norm_name(names(df))
  if (length(exclude) == 1 && is.na(exclude)) exclude <- character()
  for (pat in patterns) {
    idx <- which(grepl(pat, nn))
    if (length(idx) > 0) {
      cand <- names(df)[idx[1]]
      if (!cand %in% exclude) return(cand)
    }
  }
  NA_character_
}

map_db <- function(x) {
  x <- tolower(x)
  db <- rep(NA_character_, length(x))
  db[grepl("mirpathdb.*go.*biological", x) | grepl("go.*biological.*mirpathdb", x) |
        grepl("go.*biologicalprocess", x) | grepl("go_biological_process", x)] <- "GO Biological process (miRPathDB)"
  db[grepl("mirpathdb.*go.*molecular", x) | grepl("go.*molecular.*mirpathdb", x) |
        grepl("go.*molecularfunction", x) | grepl("go_molecular_function", x)] <- "GO Molecular function (miRPathDB)"
  db[grepl("mirpathdb.*kegg", x) | grepl("kegg.*mirpathdb", x)] <- "KEGG (miRPathDB)"
  db[grepl("mirpathdb.*reactome", x) | grepl("reactome.*mirpathdb", x)] <- "Reactome (miRPathDB)"
  db[grepl("mirtarbase", x)] <- "miRTarBase"
  db
}

select_pval <- function(df) {
  padj_col <- find_col(df, c("padjusted", "padj", "adjustedp", "padjust", "p\\.?adj"))
  q_col <- find_col(df, c("qvalue", "qval", "fdr"))
  p_col <- find_col(df, c("pvalue", "pval"))
  if (!is.na(padj_col)) return(list(col = padj_col, source = "P-adjusted"))
  if (!is.na(q_col)) return(list(col = q_col, source = "Q-value"))
  if (!is.na(p_col)) return(list(col = p_col, source = "P-value"))
  stop("No se encontro columna de p/Q-value en la tabla miEAA.")
}

order_terms <- function(df, group_col = NULL) {
  if (!is.null(group_col)) {
    df <- df %>% arrange(.data[[group_col]], desc(logp))
    df$term_plot <- paste(df$term, df[[group_col]], sep = "|||")
    levels <- unique(df$term_plot)
    df$term_plot <- factor(df$term_plot, levels = rev(levels))
  } else {
    df <- df %>% arrange(desc(logp))
    df$term_plot <- factor(df$term, levels = df$term)
    df$term_plot <- factor(df$term_plot, levels = rev(levels(df$term_plot)))
  }
  df
}

# ---- Deterministic export function (PDF + SVG + PNG) ----
save_plot <- function(plot, out_base, n_rows = 1, width_mult = 1, height_mult = 1) {
  width_mm <- p$width_mm * width_mult
  height_mm <- p$height_mm * height_mult * n_rows

  # PDF with font embedding
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else pdf
  pdf_path <- paste0(out_base, ".pdf")
  tryCatch({
    ggsave(pdf_path, plot, width = width_mm, height = height_mm,
           units = "mm", dpi = p$dpi_vector, limitsize = FALSE, device = pdf_dev)
    cat("  Exported:", pdf_path, sprintf("(%d x %d mm)\n", round(width_mm), round(height_mm)))
  }, error = function(e) {
    warning("Failed to create PDF: ", conditionMessage(e))
  })

  # SVG (vector, editable text)
  svg_path <- paste0(out_base, ".svg")
  if (requireNamespace("svglite", quietly = TRUE)) {
    tryCatch({
      ggsave(svg_path, plot, width = width_mm, height = height_mm,
             units = "mm", dpi = p$dpi_vector, limitsize = FALSE, device = svglite::svglite)
      cat("  Exported:", svg_path, "\n")
    }, error = function(e) {
      warning("Failed to create SVG: ", conditionMessage(e))
    })
  } else {
    cat("  SKIP SVG (svglite not installed)\n")
  }

  # PNG (raster, high-res)
  png_path <- paste0(out_base, ".png")
  tryCatch({
    ggsave(png_path, plot, width = width_mm, height = height_mm,
           units = "mm", dpi = p$dpi_raster, limitsize = FALSE, device = "png")
    cat("  Exported:", png_path, sprintf("(%d dpi)\n", p$dpi_raster))
  }, error = function(e) {
    warning("Failed to create PNG: ", conditionMessage(e))
  })

  # Verify at least PDF was created
  if (!file.exists(pdf_path)) {
    stop("CRITICAL: Failed to create PDF export: ", pdf_path)
  }

  invisible(list(pdf = pdf_path, svg = svg_path, png = png_path))
}

# ---- Validate inputs ----
if (!file.exists(spec_fp)) stop("Spec file not found: ", spec_fp)
spec <- read.csv(spec_fp, check.names = FALSE)
if (!is.data.frame(spec) || nrow(spec) == 0) stop("Spec file is empty or invalid")
if (!("analysis_id" %in% names(spec))) stop("Spec file missing 'analysis_id' column")
analysis_ids <- as.character(spec$analysis_id)
analysis_ids <- analysis_ids[nzchar(analysis_ids)]
if (length(analysis_ids) == 0) stop("No valid analysis_ids found in spec")

if (!dir.exists(in_root)) {
  stop("Input directory not found: ", in_root, "\nCheck --in_root parameter")
}

# ---- Structured logging ----
cat("============================================\n")
cat("miEAA GSEA Bubble Plots - Publication Ready\n")
cat("============================================\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("Spec file:", spec_fp, "\n")
cat("Input root:", in_root, "\n")
cat("Output root:", out_root, "\n")
cat("Run tag:", run_tag, "\n")
cat("Top N per DB:", n_mirpathdb, "\n")
cat("Label wrap width:", wrap_width, "\n")
cat("Preset:", preset, "\n")
cat("Seed (reproducibility):", seed, "\n")
cat("Figure dimensions:", sprintf("%d x %d mm", p$width_mm, p$height_mm), "\n")
cat("DPI (vector/raster):", sprintf("%d / %d", p$dpi_vector, p$dpi_raster), "\n")
cat("Analysis IDs:", length(analysis_ids), "\n")
cat("Log file:", log_fp, "\n")
cat("Log copy:", log_fp_copy, "\n")
cat("============================================\n\n")

mirpath_levels <- c(
  "GO Biological process (miRPathDB)",
  "KEGG (miRPathDB)",
  "Reactome (miRPathDB)"
)

for (aid in analysis_ids) {
  aid_fs <- safe_name(aid)
  top50_fp <- find_latest_all(in_root, aid_fs, run_tag)
  if (is.na(top50_fp) || !file.exists(top50_fp)) {
    cat("SKIP", aid, "-> no hay MiEAA_GSEA_all para run_tag=", run_tag, "\n")
    next
  }

  cat("\n--------------------------------------------\n")
  cat("Analysis:", aid, "\n")
  cat("Top50 file:", top50_fp, "\n")

  df <- read.delim(top50_fp, check.names = FALSE)
  names(df) <- trimws(names(df))
  if (!is.data.frame(df) || nrow(df) == 0) {
    cat("WARN", aid, "-> top50 vacio\n")
    next
  }

  category_col <- find_col(df, c("category", "database", "source", "collection", "set"))
  enrich_col <- find_col(df, c("enrichment"))
  obs_col <- find_col(df, c("observed", "hits", "count"))
  term_col <- find_col(df, c("term", "name", "description", "label", "pathway", "subcategory", "subcat"), exclude = category_col)
  if (is.na(term_col)) {
    name_lower <- tolower(names(df))
    sub_idx <- which(name_lower == "subcategory")
    if (length(sub_idx) > 0) term_col <- names(df)[sub_idx[1]]
  }
  if (is.na(term_col)) {
    name_lower <- tolower(names(df))
    drop_cols <- c(category_col, obs_col)
    drop_cols <- drop_cols[!is.na(drop_cols)]
    drop_cols <- unique(c(drop_cols, names(df)[grepl("p\\s*value|p\\s*adjusted|q\\s*value", name_lower)]))
    drop_cols <- unique(c(drop_cols, "miRNAs/precursors", "miRNAs.precursors"))
    cand <- setdiff(names(df), drop_cols)
    if (length(cand) > 0) {
      term_col <- cand[1]
      cat("WARN", aid, "-> term_col inferido:", term_col, "\n")
    }
  }

  if (is.na(term_col)) {
    cat("DEBUG columns:", paste(names(df), collapse = ", "), "\n")
    stop("No se encontro columna de terminos en: ", top50_fp)
  }
  if (is.na(category_col)) {
    cat("WARN", aid, "-> no se encontro columna de categoria; no se puede inferir DB.\n")
  }

  p_sel <- select_pval(df)
  if (!is.na(enrich_col)) {
    dir_raw <- tolower(as.character(df[[enrich_col]]))
    df$direction <- ifelse(grepl("enrich|over", dir_raw), "Enriched",
                           ifelse(grepl("deplet|under", dir_raw), "Depleted", "Other"))
    df$direction <- factor(df$direction, levels = c("Enriched", "Depleted", "Other"))
  } else {
    df$direction <- NA_character_
  }
  df$term <- str_squish(as.character(df[[term_col]]))
  df$pval_use <- as_num_safe(df[[p_sel$col]])
  if (!is.na(obs_col)) {
    df$observed <- as_num_safe(df[[obs_col]])
  } else {
    df$observed <- 1
    cat("WARN", aid, "-> no se encontro 'Observed'; usando size=1.\n")
  }
  df$observed <- ifelse(is.finite(df$observed), df$observed, NA_real_)
  df$db_raw <- if (!is.na(category_col)) as.character(df[[category_col]]) else NA_character_
  df$db <- map_db(df$db_raw)
  df$logp <- -log10(pmax(df$pval_use, .Machine$double.xmin))

  df <- df[is.finite(df$pval_use) & nzchar(df$term), , drop = FALSE]
  if (nrow(df) == 0) {
    cat("WARN", aid, "-> sin filas validas tras filtros.\n")
    next
  }

  # ---- miRPathDB bubble plots: SEPARATE enriched vs depleted ----
  df_mirpath <- df[df$db %in% mirpath_levels, , drop = FALSE]
  df_mirpath$db <- factor(df_mirpath$db, levels = mirpath_levels)

  # Check if direction column is available
  has_direction <- !is.na(enrich_col) && "direction" %in% names(df_mirpath) && any(df_mirpath$direction %in% c("Enriched", "Depleted"))

  if (!has_direction) {
    cat("WARN", aid, "-> No direction column; generating single combined plot.\n")
    # Fallback: single plot with all terms
    df_mirpath <- df_mirpath %>%
      group_by(db) %>%
      arrange(pval_use, .by_group = TRUE) %>%
      slice_head(n = n_mirpathdb) %>%
      ungroup()

    if (nrow(df_mirpath) > 0) {
      df_mirpath <- order_terms(df_mirpath, group_col = "db")
      size_range <- c(p$point * 0.8, p$point * 2.8)
      p_mirpath <- ggplot(df_mirpath, aes(x = logp, y = term_plot)) +
        geom_point(aes(size = observed), alpha = 0.85, color = "#009E73") +
        facet_wrap(~db, ncol = 2, drop = FALSE, scales = "free_y",
                   labeller = labeller(db = label_wrap_gen(width = 16))) +
        scale_y_discrete(labels = function(x) str_wrap(gsub("\\|\\|\\|.*$", "", x), width = wrap_width)) +
        scale_size(range = size_range, name = "Observed") +
        labs(x = paste0("-log10(", p_sel$source, ")"), y = NULL,
             title = paste0("miRPathDB: ", aid)) +
        theme_pub()

      out_base <- file.path(out_root, run_tag, paste0("Bubble_miRPathDB_combined_", aid_fs, "_", run_tag))
      save_plot(p_mirpath, out_base, n_rows = 2, width_mult = 1.25)
    }
  } else {
    # ---- CREATE OUTPUT SUBDIRECTORIES ----
    dir_enriched_depleted <- file.path(out_root, run_tag, "enriched_depleted")
    dir_per_database <- file.path(out_root, run_tag, "per_database")
    dir.create(dir_enriched_depleted, showWarnings = FALSE, recursive = TRUE)
    dir.create(dir_per_database, showWarnings = FALSE, recursive = TRUE)

    # ---- SEPARATE PLOTS: Enriched and Depleted ----
    # Store plots for combined figure
    plots_list <- list()

    for (dir_type in c("Enriched", "Depleted")) {
      df_dir <- df_mirpath %>%
        filter(direction == dir_type) %>%
        group_by(db) %>%
        arrange(pval_use, .by_group = TRUE) %>%
        slice_head(n = n_mirpathdb) %>%
        ungroup()

      dir_counts <- table(df_dir$db)
      cat(sprintf("  %s terms: %s\n", dir_type, paste(names(dir_counts), dir_counts, sep = "=", collapse = ", ")))

      if (nrow(df_dir) == 0) {
        cat(sprintf("  SKIP %s: no terms available\n", dir_type))
        next
      }

      # Order terms within each DB
      df_dir <- order_terms(df_dir, group_col = "db")

      # Size range for bubbles
      size_range <- c(p$point * 0.8, p$point * 2.8)

      # Color: single color per direction (no legend needed)
      dir_color <- colors_direction[[dir_type]]

      # Create plot with 3-panel layout (ncol=3 for GO BP, KEGG, Reactome)
      p_dir <- ggplot(df_dir, aes(x = logp, y = term_plot)) +
        geom_point(aes(size = observed), alpha = 0.85, color = dir_color) +
        facet_wrap(~db, ncol = 3, drop = FALSE, scales = "free_y",
                   labeller = labeller(db = label_wrap_gen(width = 16))) +
        scale_y_discrete(labels = function(x) str_wrap(gsub("\\|\\|\\|.*$", "", x), width = wrap_width)) +
        scale_size(range = size_range, name = "Observed") +
        labs(
          x = paste0("-log10(", p_sel$source, ")"),
          y = NULL,
          title = paste0("miRPathDB ", dir_type, ": ", aid)
        ) +
        theme_pub()

      # Save individual plot to enriched_depleted/
      dir_tag <- tolower(dir_type)
      out_base <- file.path(dir_enriched_depleted, paste0("Bubble_miRPathDB_", dir_tag, "_", aid_fs, "_", run_tag, "_", ts))
      save_plot(p_dir, out_base, n_rows = 1, width_mult = 1.5)
      cat(sprintf("  Saved %s plot to enriched_depleted/\n", dir_type))

      # Store plot for combined figure
      plots_list[[dir_type]] <- p_dir

      # ---- PER-DATABASE PLOTS ----
      for (db_name in mirpath_levels) {
        df_db <- df_dir %>% filter(db == db_name)
        if (nrow(df_db) == 0) next

        # Order terms for single-DB plot
        df_db <- df_db %>% arrange(desc(logp))
        df_db$term_plot <- factor(df_db$term, levels = df_db$term)
        df_db$term_plot <- factor(df_db$term_plot, levels = rev(levels(df_db$term_plot)))

        # Single-panel plot for this database
        p_db <- ggplot(df_db, aes(x = logp, y = term_plot)) +
          geom_point(aes(size = observed), alpha = 0.85, color = dir_color) +
          scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width)) +
          scale_size(range = size_range, name = "Observed") +
          labs(
            x = paste0("-log10(", p_sel$source, ")"),
            y = NULL,
            title = paste0(db_name, " (", dir_type, "): ", aid)
          ) +
          theme_pub()

        # Save to per_database/
        db_tag <- safe_name(gsub("\\(miRPathDB\\)", "", db_name))
        out_base_db <- file.path(dir_per_database, paste0(aid_fs, "_", db_tag, "_", dir_tag, "_", run_tag, "_", ts))
        save_plot(p_db, out_base_db, n_rows = 1, width_mult = 0.7)
        cat(sprintf("    Saved %s %s to per_database/\n", db_name, dir_type))
      }
    }

    # ---- COMBINED PLOT: Enriched + Depleted side-by-side ----
    if (length(plots_list) == 2) {
      cat("  Creating combined figure (Enriched + Depleted)...\n")

      # Combine plots horizontally with patchwork
      p_combined <- plots_list[["Enriched"]] + plots_list[["Depleted"]] +
        plot_layout(guides = "collect") +
        plot_annotation(
          title = paste0("miRPathDB Enrichment Summary: ", aid),
          theme = theme(
            plot.title = element_text(size = p$title + 1, face = "bold", hjust = 0.5)
          )
        )

      # Save combined plot to enriched_depleted/ (wider, 2x width for side-by-side)
      out_base_combined <- file.path(dir_enriched_depleted, paste0("Bubble_miRPathDB_combined_", aid_fs, "_", run_tag, "_", ts))
      save_plot(p_combined, out_base_combined, n_rows = 1, width_mult = 3)
      cat("  Saved combined plot to enriched_depleted/\n")
    } else {
      cat("  SKIP combined plot: need both Enriched and Depleted\n")
    }
  }

  # QC summary
  max_label <- max(nchar(df$term), na.rm = TRUE)
  cat("QC: term_max_chars=", max_label,
      " | p_source=", p_sel$source,
      " | n_total=", nrow(df), "\n")
}
