#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

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

spec_fp <- get_arg("--spec", file.path("config", "de_specs.csv"))
in_root <- get_arg("--in_root", file.path("results", "tables", "MiEAA_GSEA"))
out_root <- get_arg("--out_root", file.path("results", "figures", "MiEAA_GSEA_bubble"))
run_tag <- get_arg("--run_tag", "A_conservative")
n_mirpathdb <- as.integer(get_arg("--n_mirpathdb", "12"))
wrap_width <- as.integer(get_arg("--wrap_width", "40"))
preset <- get_arg("--preset", "double_col")

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
})

presets <- list(
  single_col = list(width_mm = 85, height_mm = 65, base = 8, title = 9, axis = 8,
                    legend = 7, ticks = 7, line = 0.6, point = 1.6, spacing_mm = 2, margin_mm = 6),
  double_col = list(width_mm = 180, height_mm = 120, base = 8.5, title = 10, axis = 8.5,
                    legend = 7.5, ticks = 7.5, line = 0.7, point = 1.8, spacing_mm = 2.5, margin_mm = 7),
  presentation = list(width_mm = 254, height_mm = 143, base = 18, title = 22, axis = 18,
                      legend = 16, ticks = 16, line = 1.5, point = 4, spacing_mm = 6, margin_mm = 10),
  poster = list(width_mm = 508, height_mm = 356, base = 28, title = 34, axis = 28,
                legend = 24, ticks = 24, line = 2.0, point = 6, spacing_mm = 10, margin_mm = 15)
)
if (!preset %in% names(presets)) stop("preset invalido: ", preset)
p <- presets[[preset]]

theme_pub <- function() {
  theme_minimal(base_size = p$base, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = p$title, face = "bold"),
      axis.title = element_text(size = p$axis),
      axis.text = element_text(size = p$ticks, color = "black"),
      strip.text = element_text(size = p$axis, face = "bold"),
      legend.title = element_text(size = p$legend),
      legend.text = element_text(size = p$legend),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.spacing = unit(p$spacing_mm, "mm"),
      plot.margin = margin(p$margin_mm, p$margin_mm, p$margin_mm, p$margin_mm + 6, unit = "mm")
    )
}

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
  q_col <- find_col(df, c("qvalue", "qval", "fdr"))
  padj_col <- find_col(df, c("padjusted", "padj", "adjustedp", "padjust", "p\\.?adj"))
  p_col <- find_col(df, c("pvalue", "pval"))
  if (!is.na(q_col)) return(list(col = q_col, source = "Q-value"))
  if (!is.na(padj_col)) return(list(col = padj_col, source = "P-adjusted"))
  if (!is.na(p_col)) return(list(col = p_col, source = "P-value"))
  stop("No se encontro columna de p/Q-value en la tabla miEAA.")
}

order_terms <- function(df, group_col = NULL) {
  if (!is.null(group_col)) {
    df <- df %>% arrange(.data[[group_col]], desc(logp))
    df$term_plot <- paste(df$term, df[[group_col]], sep = "|||")
    levels <- unique(df$term_plot)
    df$term_plot <- factor(df$term_plot, levels = levels)
  } else {
    df <- df %>% arrange(desc(logp))
    df$term_plot <- factor(df$term, levels = df$term)
  }
  df
}

save_plot <- function(plot, out_base, n_rows = 1, width_mult = 1) {
  width_mm <- p$width_mm * width_mult
  height_mm <- p$height_mm * n_rows
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(paste0(out_base, ".pdf"), plot, width = width_mm, height = height_mm,
         units = "mm", dpi = 300, limitsize = FALSE, device = pdf_dev)
  svg_dev <- if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else "svg"
  tryCatch(
    ggsave(paste0(out_base, ".svg"), plot, width = width_mm, height = height_mm,
           units = "mm", dpi = 300, limitsize = FALSE, device = svg_dev),
    error = function(e) warning("No se pudo crear SVG: ", conditionMessage(e))
  )
  tryCatch(
    ggsave(paste0(out_base, ".png"), plot, width = width_mm, height = height_mm,
           units = "mm", dpi = 600, limitsize = FALSE),
    error = function(e) warning("No se pudo crear PNG: ", conditionMessage(e))
  )
  if (!file.exists(paste0(out_base, ".pdf"))) {
    warning("No se pudo crear PDF: ", out_base, ".pdf")
  }
}

if (!file.exists(spec_fp)) stop("No existe spec: ", spec_fp)
spec <- read.csv(spec_fp, check.names = FALSE)
if (!("analysis_id" %in% names(spec))) stop("spec sin 'analysis_id'")
analysis_ids <- as.character(spec$analysis_id)

cat("============================================\n")
cat("miEAA GSEA bubble plots\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("spec:", spec_fp, "\n")
cat("in_root:", in_root, "\n")
cat("out_root:", out_root, "\n")
cat("run_tag:", run_tag, "\n")
cat("n_mirpathdb:", n_mirpathdb, "\n")
cat("wrap_width:", wrap_width, "\n")
cat("preset:", preset, "\n")
cat("Log:", log_fp, "\n")
cat("Log copy:", log_fp_copy, "\n")
cat("============================================\n\n")

mirpath_levels <- c(
  "GO Biological process (miRPathDB)",
  "GO Molecular function (miRPathDB)",
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

  # miRPathDB 2x2
  df_mirpath <- df[df$db %in% mirpath_levels, , drop = FALSE]
  df_mirpath$db <- factor(df_mirpath$db, levels = mirpath_levels)
  df_mirpath <- df_mirpath %>%
    group_by(db) %>%
    arrange(pval_use, .by_group = TRUE) %>%
    slice_head(n = n_mirpathdb) %>%
    ungroup()

  mirpath_counts <- table(df_mirpath$db)
  cat("miRPathDB terms:", paste(names(mirpath_counts), mirpath_counts, sep = "=", collapse = ", "), "\n")

  if (nrow(df_mirpath) > 0) {
    df_mirpath <- order_terms(df_mirpath, group_col = "db")
    size_range <- c(p$point * 0.8, p$point * 2.8)
    p_mirpath <- ggplot(df_mirpath, aes(x = logp, y = term_plot)) +
      geom_point(aes(size = observed, color = pval_use), alpha = 0.85) +
      facet_wrap(~db, ncol = 2, drop = FALSE, scales = "free_y",
                 labeller = labeller(db = label_wrap_gen(width = 16))) +
      scale_y_discrete(labels = function(x) str_wrap(gsub("\\|\\|\\|.*$", "", x), width = wrap_width)) +
      scale_color_viridis_c(option = "D", end = 0.9, name = p_sel$source) +
      scale_size(range = size_range, name = "Observed") +
      labs(x = paste0("-log10(", p_sel$source, ")"), y = NULL,
           title = paste0("miRPathDB summary: ", aid)) +
      theme_pub()

    out_base <- file.path(out_root, run_tag, paste0("Bubble_miRPathDB_2x2_", aid_fs, "_", run_tag))
    save_plot(p_mirpath, out_base, n_rows = 2, width_mult = 1.25)
  } else {
    cat("WARN", aid, "-> sin terminos miRPathDB para graficar.\n")
  }

  # QC summary
  max_label <- max(nchar(df$term), na.rm = TRUE)
  cat("QC: term_max_chars=", max_label,
      " | p_source=", p_sel$source,
      " | n_total=", nrow(df), "\n")
}
