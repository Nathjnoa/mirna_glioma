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
out_root <- get_arg("--out_root", file.path("results", "figures", "MiEAA_GSEA_categories_pvals"))
run_tag <- get_arg("--run_tag", "A_conservative")
wrap_width <- as.integer(get_arg("--wrap_width", "35"))
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
log_fp <- file.path("logs", paste0("mieaa_categories_pvals_", ts, ".txt"))
log_fp_copy <- file.path(out_root, run_tag, paste0("mieaa_categories_pvals_", ts, ".txt"))
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
      legend.title = element_text(size = p$legend),
      legend.text = element_text(size = p$legend),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(p$margin_mm, p$margin_mm, p$margin_mm, p$margin_mm, unit = "mm")
    )
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
  db
}

find_latest <- function(root, aid_fs, run_tag, kind = c("all", "top50")) {
  kind <- match.arg(kind)
  dir_path <- file.path(root, aid_fs, run_tag)
  if (!dir.exists(dir_path)) return(NA_character_)
  if (kind == "all") {
    pat <- paste0("^MiEAA_GSEA_all_", aid_fs, "_", run_tag, "_.*\\.tsv$")
  } else {
    pat <- paste0("^MiEAA_GSEA_top50_", aid_fs, "_", run_tag, "_.*\\.tsv$")
  }
  ff <- list.files(dir_path, pattern = pat, full.names = TRUE)
  if (length(ff) == 0) return(NA_character_)
  ff[order(file.info(ff)$mtime, decreasing = TRUE)][1]
}

save_plot <- function(plot, out_base) {
  width_mm <- p$width_mm
  height_mm <- p$height_mm
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
           units = "mm", dpi = 300, limitsize = FALSE, device = "png"),
    error = function(e) warning("No se pudo crear PNG: ", conditionMessage(e))
  )
}

if (!file.exists(spec_fp)) stop("No existe spec: ", spec_fp)
spec <- read.csv(spec_fp, check.names = FALSE)
if (!("analysis_id" %in% names(spec))) stop("spec sin 'analysis_id'")
analysis_ids <- as.character(spec$analysis_id)

cat("============================================\n")
cat("miEAA GSEA categories p-values\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("spec:", spec_fp, "\n")
cat("in_root:", in_root, "\n")
cat("out_root:", out_root, "\n")
cat("run_tag:", run_tag, "\n")
cat("wrap_width:", wrap_width, "\n")
cat("preset:", preset, "\n")
cat("Log:", log_fp, "\n")
cat("Log copy:", log_fp_copy, "\n")
cat("============================================\n\n")

for (aid in analysis_ids) {
  aid_fs <- safe_name(aid)
  all_fp <- find_latest(in_root, aid_fs, run_tag, kind = "all")
  use_fp <- all_fp
  used_kind <- "all"
  if (is.na(use_fp) || !file.exists(use_fp)) {
    use_fp <- find_latest(in_root, aid_fs, run_tag, kind = "top50")
    used_kind <- "top50"
  }
  if (is.na(use_fp) || !file.exists(use_fp)) {
    cat("SKIP", aid, "-> no hay archivos miEAA para run_tag=", run_tag, "\n")
    next
  }
  if (used_kind == "top50") {
    cat("WARN", aid, "-> usando top50 (distribuciones pueden estar sesgadas):", use_fp, "\n")
  } else {
    cat("Analysis:", aid, "-> usando:", use_fp, "\n")
  }

  df <- read.delim(use_fp, check.names = FALSE)
  names(df) <- trimws(names(df))
  if (!is.data.frame(df) || nrow(df) == 0) {
    cat("WARN", aid, "-> tabla vacia\n")
    next
  }

  category_col <- find_col(df, c("category", "database", "source", "collection", "set"))
  q_col <- find_col(df, c("qvalue", "qval", "fdr"))
  padj_col <- find_col(df, c("padjusted", "padj", "p\\.?adj", "adjustedp"))

  if (is.na(category_col)) {
    cat("WARN", aid, "-> no se encontro columna de categoria; usando 'Unknown'.\n")
    db_raw <- rep("Unknown", nrow(df))
  } else {
    db_raw <- as.character(df[[category_col]])
  }
  db_map <- map_db(db_raw)
  df$db <- ifelse(!is.na(db_map), db_map, db_raw)

  if (!is.na(q_col)) {
    df$qval <- as_num_safe(df[[q_col]])
  }
  if (!is.na(padj_col)) {
    df$padj <- as_num_safe(df[[padj_col]])
  }

  plot_metric <- function(metric, label, out_tag) {
    if (!metric %in% names(df)) {
      cat("WARN", aid, "-> falta columna", label, "; se omite plot.\n")
      return(invisible(NULL))
    }
    dd <- df[is.finite(df[[metric]]), , drop = FALSE]
    if (nrow(dd) == 0) {
      cat("WARN", aid, "-> sin valores validos para", label, "\n")
      return(invisible(NULL))
    }
    dd$logv <- -log10(pmax(dd[[metric]], .Machine$double.xmin))
    med <- dd %>% group_by(db) %>% summarise(med = median(logv, na.rm = TRUE), .groups = "drop")
    ord <- med %>% arrange(desc(med)) %>% pull(db)
    dd$db <- factor(dd$db, levels = ord)
    title_txt <- paste0(aid, " | ", run_tag, " | Categories p-values (metric: ", label, ")")

    p_plot <- ggplot(dd, aes(x = db, y = logv)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.25, linewidth = p$line * 0.5,
                   color = "#4D4D4D", fill = "#B3CDE3") +
      geom_jitter(width = 0.18, alpha = 0.8, size = 1.2, color = "#0072B2") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = wrap_width)) +
      labs(x = NULL, y = paste0("-log10(", label, ")"), title = title_txt) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 25, hjust = 1))

    out_base <- file.path(out_root, run_tag,
                          paste0("Categories_pvals_", out_tag, "_", aid_fs, "_", run_tag))
    save_plot(p_plot, out_base)
  }

  if (!is.na(q_col)) {
    plot_metric("qval", "Q-value", "Q")
  } else {
    cat("WARN", aid, "-> no se encontro Q-value; se omite plot Q.\n")
  }
  if (!is.na(padj_col)) {
    plot_metric("padj", "P-adjusted", "Padj")
  } else {
    cat("WARN", aid, "-> no se encontro P-adjusted; se omite plot Padj.\n")
  }
}
