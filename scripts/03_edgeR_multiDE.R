#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

# -----------------------------
# CLI parsing (sin dependencias)
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[[hit + 1]]
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

parse_kv_list <- function(s) {
  res <- list()
  if (is.null(s) || nchar(s) == 0) return(res)
  parts <- unlist(strsplit(s, ";", fixed = TRUE))
  for (p in parts) {
    kv <- unlist(strsplit(p, ":", fixed = TRUE), use.names = FALSE)
    if (length(kv) != 2) next
    key <- trimws(kv[1])
    vals <- trimws(unlist(strsplit(kv[2], ",", fixed = TRUE)))
    vals <- vals[nzchar(vals)]
    if (length(key) > 0 && length(vals) > 0) res[[key]] <- vals
  }
  res
}

parse_contrasts <- function(s) {
  res <- list()
  if (is.null(s) || nchar(s) == 0) return(res)
  parts <- unlist(strsplit(s, ";", fixed = TRUE))
  for (p in parts) {
    kv <- unlist(strsplit(p, ":", fixed = TRUE), use.names = FALSE)
    if (length(kv) != 2) next
    key <- trimws(kv[1])
    cr <- unlist(strsplit(kv[2], "vs", fixed = TRUE), use.names = FALSE)
    if (length(cr) != 2) next
    case <- trimws(cr[1]); ref <- trimws(cr[2])
    if (length(key) > 0 && nchar(case) > 0 && nchar(ref) > 0) {
      res[[key]] <- list(case = case, ref = ref)
    }
  }
  res
}

parse_range_map <- function(s) {
  if (is.null(s) || nchar(s) == 0) return(list())
  parts <- unlist(strsplit(s, ";", fixed = TRUE))
  out <- list()
  for (p in parts) {
    kv <- unlist(strsplit(p, "=", fixed = TRUE), use.names = FALSE)
    if (length(kv) != 2) next
    rng <- kv[1]
    lab <- kv[2]
    bounds <- unlist(strsplit(rng, "-", fixed = TRUE), use.names = FALSE)
    if (length(bounds) != 2) next
    mn <- suppressWarnings(as.numeric(gsub(",", ".", trimws(bounds[1]))))
    mx <- suppressWarnings(as.numeric(gsub(",", ".", trimws(bounds[2]))))
    if (is.na(mn) || is.na(mx) || nchar(lab) == 0) next
    out[[length(out) + 1]] <- list(min = mn, max = mx, label = trimws(lab))
  }
  out
}

as_numeric_robust <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
}

safe_numeric <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  x <- gsub(",", ".", as.character(x))
  suppressWarnings(as.numeric(x))
}

parse_collapse_levels <- function(s) {
  res <- list()
  if (is.null(s) || nchar(s) == 0) return(res)
  parts <- unlist(strsplit(s, ";", fixed = TRUE))
  for (p in parts) {
    kv <- unlist(strsplit(p, "=", fixed = TRUE), use.names = FALSE)
    if (length(kv) != 2) next
    lhs <- trimws(kv[1]); rhs <- trimws(kv[2])
    if (!nzchar(lhs) || !nzchar(rhs)) next
    lhs_vals <- trimws(unlist(strsplit(lhs, "\\|")))
    lhs_vals <- lhs_vals[nzchar(lhs_vals)]
    for (lv in lhs_vals) res[[lv]] <- rhs
  }
  res
}

# lectura genérica de tablas (csv/tsv autodetectado por separador de data.table)
read_table <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(path, data.table = FALSE, check.names = FALSE)
  } else {
    utils::read.csv(path, check.names = FALSE)
  }
}

counts_fp <- get_arg("--counts", file.path("data","intermediate","Gliomas_all_counts_merged.csv"))
meta_fp   <- get_arg("--meta",   file.path("data","intermediate","Metadatos_gliomas_verificados.csv"))
spec_fp   <- get_arg("--spec",   NULL)
vars_str  <- get_arg("--vars",   NULL)  # e.g. "sexo,IDH,MGMT"
drop_re   <- get_arg("--drop_regex", "^A")
min_lib   <- as.numeric(get_arg("--min_libsize", "100000"))
outdir    <- get_arg("--outdir", file.path("results","DE"))
ref_mode  <- get_arg("--ref_mode", "first")  # "first" or "last" (por si quieres referencia)
contrasts_str <- get_arg("--contrasts", "")
include_str   <- get_arg("--include_levels", "")
exclude_str   <- get_arg("--exclude_levels", "")
logcpm_fp     <- get_arg("--logcpm", "")
spec_fp       <- get_arg("--spec", "")
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

vars <- character(0)
include_map <- parse_kv_list(include_str)
exclude_map <- parse_kv_list(exclude_str)
contrasts_map <- parse_contrasts(contrasts_str)
spec_rows <- NULL
if (nchar(spec_fp) > 0) {
  if (!file.exists(spec_fp)) stop("No existe spec: ", spec_fp)
  spec_rows <- read_table(spec_fp)
  required_cols <- c("analysis_id", "mode", "var")
  missing_cols <- setdiff(required_cols, colnames(spec_rows))
  if (length(missing_cols) > 0) {
    stop("Spec faltan columnas requeridas: ", paste(missing_cols, collapse = ", "))
  }
} else {
  if (is.null(vars_str) || nchar(vars_str) == 0) {
    stop("Debes indicar --vars, por ejemplo: --vars sexo,IDH,MGMT (o usar --spec).")
  }
  vars <- trimws(strsplit(vars_str, ",")[[1]])
  vars <- vars[nzchar(vars)]
  if (length(vars) == 0) stop("No se detectaron variables en --vars.")
}

dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("data","processed"), showWarnings = FALSE, recursive = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
fig_root <- file.path("results", "figures")
tab_root <- file.path("results", "tables")
dir.create(fig_root, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_root, showWarnings = FALSE, recursive = TRUE)
logcpm_qc_dir <- file.path(fig_root, "logcpm_qc")
dir.create(logcpm_qc_dir, showWarnings = FALSE, recursive = TRUE)

log_fp <- file.path("logs", paste0("DE_multivar_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")

cat("============================================\n")
cat("edgeR multi-DE (univariada)\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("Counts:", counts_fp, "\n")
cat("Meta:", meta_fp, "\n")
cat("Vars:", paste(vars, collapse = ", "), "\n")
if (!is.null(spec_rows)) {
  cat("Spec:", spec_fp, "| filas:", nrow(spec_rows), "\n")
}
cat("drop_regex:", drop_re, "\n")
cat("min_libsize:", min_lib, "\n")
cat("outdir:", outdir, "\n")
cat("ref_mode:", ref_mode, "\n")
cat("============================================\n\n")

suppressPackageStartupMessages({
  library(edgeR)
})

plot_volcano <- function(tt, analysis_fs, tag, fig_dir, ts) {
  if (!("logFC" %in% colnames(tt))) return()
  has_fdr <- "FDR" %in% colnames(tt)
  has_p <- "PValue" %in% colnames(tt)
  if (!has_fdr && !has_p) return()

  # Load ggplot2 and ggrepel for publication-quality plots
  suppressPackageStartupMessages({
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      cat("WARN: ggplot2 not available, skipping volcano plot\n")
      return()
    }
    library(ggplot2)
  })

  # Prepare data
  df_plot <- data.frame(
    feature_id = tt$feature_id,
    logFC = as.numeric(tt$logFC),
    pval = if (has_fdr) as.numeric(tt$FDR) else as.numeric(tt$PValue),
    stringsAsFactors = FALSE
  )
  df_plot <- df_plot[is.finite(df_plot$logFC) & is.finite(df_plot$pval), ]
  if (nrow(df_plot) == 0) return()

  df_plot$neglog10p <- -log10(pmax(df_plot$pval, .Machine$double.xmin))

  # Classify genes by significance and direction
  fdr_thr <- 0.05
  lfc_thr <- 0.584  # log2(1.5) ≈ 0.584, corresponds to 1.5-fold change
  df_plot$category <- "Not Significant"
  df_plot$category[df_plot$pval < fdr_thr & df_plot$logFC > lfc_thr] <- "Up-regulated"
  df_plot$category[df_plot$pval < fdr_thr & df_plot$logFC < -lfc_thr] <- "Down-regulated"
  df_plot$category[df_plot$pval < fdr_thr & abs(df_plot$logFC) <= lfc_thr] <- "Significant (|logFC| < 0.584)"
  df_plot$category <- factor(df_plot$category,
                             levels = c("Up-regulated", "Down-regulated",
                                        "Significant (|logFC| < 0.584)", "Not Significant"))

  # Colorblind-safe palette
  colors <- c("Up-regulated" = "#D55E00",
              "Down-regulated" = "#0072B2",
              "Significant (|logFC| < 0.584)" = "#CC79A7",
              "Not Significant" = "#999999")

  # Top genes to label (top 10 by significance among significant)
  df_sig <- df_plot[df_plot$pval < fdr_thr, ]
  if (nrow(df_sig) > 0) {
    df_sig <- df_sig[order(df_sig$pval), ]
    top_genes <- head(df_sig$feature_id, 10)
    df_plot$label <- ifelse(df_plot$feature_id %in% top_genes, df_plot$feature_id, NA)
  } else {
    df_plot$label <- NA
  }

  # Count genes per category
  n_up <- sum(df_plot$category == "Up-regulated")
  n_down <- sum(df_plot$category == "Down-regulated")
  n_sig_low <- sum(df_plot$category == "Significant (|logFC| < 0.584)")

  # Build plot
  p <- ggplot(df_plot, aes(x = logFC, y = neglog10p, color = category)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = colors, name = "Category", drop = FALSE) +
    geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", color = "#333333", linewidth = 0.5) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dotted", color = "#666666", linewidth = 0.4) +
    labs(
      x = expression(log[2]~Fold~Change),
      y = if (has_fdr) expression(-log[10]~FDR) else expression(-log[10]~P-value),
      title = paste0(analysis_fs, ": ", tag),
      subtitle = sprintf("Up: %d | Down: %d | Sig (|logFC|<0.584): %d", n_up, n_down, n_sig_low)
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 12, color = "#555555", hjust = 0),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.position = "right",
      panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
      plot.margin = margin(10, 15, 10, 10, unit = "mm")
    )

  # Add labels if ggrepel available
  if (requireNamespace("ggrepel", quietly = TRUE) && any(!is.na(df_plot$label))) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = label),
      size = 3.5,
      max.overlaps = 15,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.3,
      na.rm = TRUE,
      show.legend = FALSE
    )
  }

  # Export high-resolution PNG (300 DPI)
  out_png <- file.path(fig_dir, paste0("Volcano_", analysis_fs, "_", tag, "_", ts, ".png"))
  ggsave(out_png, p, width = 10, height = 8, dpi = 300, bg = "white")
  cat("  -> Volcano:", out_png, "\n")

  # Also export PDF for publication
  out_pdf <- file.path(fig_dir, paste0("Volcano_", analysis_fs, "_", tag, "_", ts, ".pdf"))
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else pdf
  ggsave(out_pdf, p, width = 10, height = 8, device = pdf_dev)
}

# -----------------------------
# Lectura de datos
# -----------------------------
read_table <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(path, data.table = FALSE, check.names = FALSE)
  } else {
    read.csv(path, check.names = FALSE)
  }
}

if (!file.exists(counts_fp)) stop("No existe counts: ", counts_fp)
if (!file.exists(meta_fp)) stop("No existe meta: ", meta_fp)

df <- read_table(counts_fp)
meta <- read_table(meta_fp)

# limpiar filas meta con id vacío si existieran
if ("id" %in% colnames(meta)) {
  meta <- meta[!(is.na(meta$id) | trimws(meta$id) == ""), , drop = FALSE]
}

cat("Features iniciales (df):", nrow(df), "\n")
cat("Dim counts df:", nrow(df), "x", ncol(df), "\n")
cat("Dim meta:", nrow(meta), "x", ncol(meta), "\n\n")

# -----------------------------
# Construcción counts_mat base
# -----------------------------
anno_candidates <- c("id","type","entrez_id","HGNC_symbol")
anno_present <- intersect(anno_candidates, colnames(df))
count_cols <- setdiff(colnames(df), anno_present)

# drop post mortem por regex
drop_samples <- grep(drop_re, count_cols, value = TRUE)
count_cols <- setdiff(count_cols, drop_samples)

feature_id <- if ("id" %in% colnames(df)) df[["id"]] else seq_len(nrow(df))
counts_df <- df[, count_cols, drop = FALSE]
counts_df[] <- lapply(counts_df, function(x) suppressWarnings(as.numeric(x)))

counts_mat <- as.matrix(counts_df)
storage.mode(counts_mat) <- "numeric"
rownames(counts_mat) <- make.unique(as.character(feature_id))
colnames(counts_mat) <- make.unique(colnames(counts_mat))

# limpieza mínima features
keep_nz <- rowSums(counts_mat, na.rm = TRUE) > 0
counts_mat <- counts_mat[keep_nz, , drop = FALSE]
row_var <- apply(counts_mat, 1, function(x) stats::var(x, na.rm = TRUE))
keep_var <- is.finite(row_var) & row_var > 0
counts_mat <- counts_mat[keep_var, , drop = FALSE]

cat("Tras drop_regex y limpieza mínima features:\n")
cat("- counts_mat:", nrow(counts_mat), "features x", ncol(counts_mat), "samples\n\n")

# filtro global por tamaño de librería
y0 <- DGEList(counts = counts_mat)
lib0 <- y0$samples$lib.size
keep_lib <- lib0 >= min_lib
cat("Filtro global lib_size >= ", min_lib, ":\n", sep="")
cat("- samples antes:", length(lib0), "\n")
cat("- samples removidas:", sum(!keep_lib), "\n")
if (sum(!keep_lib) > 0) {
  cat("- removidas:\n")
  print(names(lib0)[!keep_lib])
}
cat("- samples después:", sum(keep_lib), "\n\n")

counts_mat <- counts_mat[, keep_lib, drop = FALSE]

# -----------------------------
# Alinear metadata al counts_mat base
# -----------------------------
if (!("id" %in% colnames(meta))) stop("El metadata debe tener columna 'id' con IDs de muestra.")

samples_final <- colnames(counts_mat)
meta_ids <- as.character(meta$id)

missing_in_meta <- setdiff(samples_final, meta_ids)
if (length(missing_in_meta) > 0) {
  cat("WARNING: muestras en counts pero no en meta$id. Se descartarán de counts:\n",
      paste(missing_in_meta, collapse = ", "), "\n")
  keep_cols <- setdiff(samples_final, missing_in_meta)
  if (length(keep_cols) == 0) {
    stop("No quedan muestras tras descartar las ausentes en metadata.")
  }
  counts_mat <- counts_mat[, keep_cols, drop = FALSE]
  samples_final <- colnames(counts_mat)
}

meta_base <- meta[match(samples_final, meta_ids), , drop = FALSE]
stopifnot(all(as.character(meta_base$id) == samples_final))

# exportar dataset base filtrado
base_counts_fp <- file.path("data","processed", paste0("counts_filtered_base_", ts, ".csv"))
base_meta_fp   <- file.path("data","processed", paste0("metadata_filtered_base_", ts, ".csv"))
write.csv(data.frame(id = rownames(counts_mat), counts_mat, check.names = FALSE),
          base_counts_fp, row.names = FALSE)
write.csv(meta_base, base_meta_fp, row.names = FALSE)

cat("Export base:\n")
cat("- ", base_counts_fp, "\n", sep="")
cat("- ", base_meta_fp, "\n\n", sep="")

# ---- Opcional: QC con matriz logCPM externa filtrada a meta ----
if (nchar(logcpm_fp) > 0) {
  if (!file.exists(logcpm_fp)) {
    cat("[logCPM QC] WARNING: no existe logCPM en ruta:", logcpm_fp, "\n\n")
  } else {
    cat("[logCPM QC] Leyendo logCPM:", logcpm_fp, "\n")
    lc_df <- read_table(logcpm_fp)
    if (ncol(lc_df) >= 2) {
      lc_features <- lc_df[[1]]
      lc_mat <- as.matrix(lc_df[, -1, drop = FALSE])
      rownames(lc_mat) <- make.unique(as.character(lc_features))
      colnames(lc_mat) <- make.unique(colnames(lc_mat))

      keep_samp <- intersect(colnames(lc_mat), samples_final)
      missing_lc <- setdiff(samples_final, colnames(lc_mat))
      if (length(missing_lc) > 0) {
        cat("[logCPM QC] WARNING: estas muestras del meta no están en logCPM:", paste(missing_lc, collapse = ", "), "\n")
      }
      lc_mat_f <- lc_mat[, keep_samp, drop = FALSE]
      cat("[logCPM QC] Dim logCPM filtrada:", nrow(lc_mat_f), "features x", ncol(lc_mat_f), "muestras\n")

      lc_filt_fp <- file.path(tab_root, paste0("logCPM_filtered_", ts, ".csv"))
      write.csv(
        cbind(feature_id = rownames(lc_mat_f), as.data.frame(lc_mat_f, check.names = FALSE)),
        lc_filt_fp,
        row.names = FALSE
      )
      cat("[logCPM QC] Exportada logCPM filtrada:", lc_filt_fp, "\n")

      if (ncol(lc_mat_f) >= 2) {
        # MDS con cmdscale
        dmat <- stats::dist(t(lc_mat_f))
        mds <- cmdscale(dmat, k = 2)
        mds_fp <- file.path(logcpm_qc_dir, paste0("MDS_logCPM_filtered_", ts, ".png"))
        png(mds_fp, width = 1400, height = 1000, res = 150)
        plot(mds[,1], mds[,2],
             xlab = "MDS1", ylab = "MDS2",
             main = "MDS logCPM (filtrado a meta)")
        text(mds[,1], mds[,2], labels = colnames(lc_mat_f), pos = 3, cex = 0.7)
        dev.off()
        cat("[logCPM QC] MDS filtrado guardado en:", mds_fp, "\n")

        # PCA
        pca <- stats::prcomp(t(lc_mat_f), center = TRUE, scale. = FALSE)
        pca_fp <- file.path(logcpm_qc_dir, paste0("PCA_logCPM_filtered_", ts, ".png"))
        png(pca_fp, width = 1400, height = 1000, res = 150)
        plot(pca$x[,1], pca$x[,2],
             xlab = paste0("PC1 (", round(100*summary(pca)$importance[2,1],1), "%)"),
             ylab = paste0("PC2 (", round(100*summary(pca)$importance[2,2],1), "%)"),
             main = "PCA logCPM (filtrado a meta)")
        text(pca$x[,1], pca$x[,2], labels = colnames(lc_mat_f), pos = 3, cex = 0.7)
        dev.off()
        cat("[logCPM QC] PCA filtrado guardado en:", pca_fp, "\n\n")
      } else {
        cat("[logCPM QC] WARNING: <2 muestras tras filtrar; se omiten MDS/PCA.\n\n")
      }
    } else {
      cat("[logCPM QC] WARNING: logCPM no tiene columnas suficientes.\n\n")
    }
  }
}

# annotation para merge en output (si existe)
anno_df <- NULL
if ("id" %in% colnames(df)) {
  keep_anno <- intersect(c("id","type","entrez_id","HGNC_symbol"), colnames(df))
  anno_df <- unique(df[, keep_anno, drop = FALSE])
}

force_numeric_df <- function(df, cols) {
  for (cc in cols) {
    if (cc %in% colnames(df)) {
      v <- df[[cc]]
      if (is.factor(v)) v <- as.character(v)
      df[[cc]] <- suppressWarnings(as.numeric(v))
    }
  }
  df
}

qc_cols <- c("analysis_id","var","mode","case","ref","cut","cut_type","collapse_map",
             "exclude_rule","event_var","n_samples","n_features_tested","groups",
             "n_sig_FDR0.05","n_sig_FDR0.1","comparison")

finalize_qc <- function(qc_row) {
  missing_cols <- setdiff(qc_cols, colnames(qc_row))
  if (length(missing_cols) > 0) {
    for (mc in missing_cols) qc_row[[mc]] <- NA
  }
  qc_row[, qc_cols, drop = FALSE]
}

decide_status_safe <- function(qlf_obj) {
  tryCatch({
    dt_fun <- get("decideTestsDGE", envir = asNamespace("edgeR"))
    dt_fun(qlf_obj, adjust.method = "BH", p.value = 0.05)
  }, error = function(e) {
    cat("WARNING: decideTestsDGE no disponible; plotMD sin status.\n")
    NULL
  })
}

run_edgeR_binary <- function(analysis_id, analysis_fs, group, counts_sub, anno_df, ts, case_label, ref_label) {
  analysis_dir <- file.path(outdir, analysis_fs)
  dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)

  y <- DGEList(counts = counts_sub, group = group)
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~ group)

  keep_feat <- filterByExpr(y, design = design)
  cat("[QC filterByExpr] keep_edgeR =", sum(keep_feat), " / ", length(keep_feat), "\n", sep = "")
  if (!any(keep_feat)) {
    cat("SKIP:", analysis_id, "-> filterByExpr eliminó todas las features.\n\n")
    return(NULL)
  }
  y <- y[keep_feat, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")  # recalcular tras filtrar features

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)

  coef_candidates <- unique(c(
    paste0("group", case_label),
    paste0("group", make.names(case_label))
  ))
  coef_idx <- which(colnames(design) %in% coef_candidates)
  if (length(coef_idx) != 1 && ncol(design) == 2) {
    # binaria típica: usa la segunda columna como contraste
    coef_idx <- 2
  }
  if (length(coef_idx) != 1) {
    cat("SKIP:", analysis_id, "-> no se encontró coeficiente para el caso. Candidatos=",
        paste(coef_candidates, collapse = ","), " | cols diseño=",
        paste(colnames(design), collapse = ","), "\n\n", sep = "")
    return(NULL)
  }

  qlf <- glmQLFTest(fit, coef = coef_idx)
  tt <- topTags(qlf, n = Inf)$table
  tt$feature_id <- rownames(tt); rownames(tt) <- NULL
  if (!is.null(anno_df)) tt <- merge(tt, anno_df, by.x = "feature_id", by.y = "id", all.x = TRUE)
  tt <- force_numeric_df(tt, c("PValue", "FDR", "logFC", "logCPM", "F"))

  n_all <- nrow(tt)
  n_fdr_na <- sum(is.na(tt$FDR))
  n_fdr05 <- sum(!is.na(tt$FDR) & tt$FDR < 0.05)
  n_fdr10 <- sum(!is.na(tt$FDR) & tt$FDR < 0.10)
  cat("   [DEBUG] ", analysis_id, " n_all=", n_all, " FDR NA=", n_fdr_na,
      " FDR<0.05=", n_fdr05, " FDR<0.1=", n_fdr10, "\n", sep = "")

  all_fp <- file.path(analysis_dir, paste0("DE_", analysis_fs, "_all_", ts, ".tsv"))
  sig05_fp <- file.path(analysis_dir, paste0("DE_", analysis_fs, "_sig_FDR0.05_", ts, ".tsv"))
  sig10_fp <- file.path(analysis_dir, paste0("DE_", analysis_fs, "_sig_FDR0.1_", ts, ".tsv"))

  write.table(tt, all_fp, sep = "\t", quote = FALSE, row.names = FALSE)
  tt_sig05 <- tt[!is.na(tt$FDR) & tt$FDR < 0.05, , drop = FALSE]
  tt_sig10 <- tt[!is.na(tt$FDR) & tt$FDR < 0.10, , drop = FALSE]
  write.table(tt_sig05, sig05_fp, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(tt_sig10, sig10_fp, sep = "\t", quote = FALSE, row.names = FALSE)

  png(file.path(analysis_dir, paste0("MDplot_", analysis_fs, "_", ts, ".png")),
      width = 1400, height = 1000, res = 150)
  plotMD(qlf, status = decide_status_safe(qlf),
         main = paste0("MD plot: ", analysis_id, " (", case_label, " vs ", ref_label, ")"))
  abline(h = c(-1, 1), lty = 2)
  dev.off()

  plot_volcano(tt, analysis_fs, paste0(case_label, "_vs_", ref_label), analysis_dir, ts)

  qc <- data.frame(
    analysis_id = analysis_id,
    n_samples = ncol(y),
    n_features_tested = nrow(y),
    groups = paste(names(table(group)), as.integer(table(group)), sep = ":", collapse = ";"),
    n_sig_FDR0.05 = n_fdr05,
    n_sig_FDR0.1 = n_fdr10,
    comparison = paste0(case_label, " vs ", ref_label),
    stringsAsFactors = FALSE
  )
  summary_fp <- file.path(analysis_dir, paste0("QC_", analysis_fs, "_summary_", ts, ".tsv"))
  write.table(qc, summary_fp, sep = "\t", quote = FALSE, row.names = FALSE)

  cat("Output (analysis_id=", analysis_id, "):\n", sep = "")
  cat("- ", all_fp, "\n", sep = "")
  cat("- ", sig05_fp, "\n", sep = "")
  cat("- ", sig10_fp, "\n\n", sep = "")

  list(
    qc = qc,
    files = list(all = all_fp, sig05 = sig05_fp, sig10 = sig10_fp)
  )
}

run_edgeR_continuous <- function(analysis_id, analysis_fs, counts_sub, x_scaled, anno_df, ts, tag) {
  analysis_dir <- file.path(outdir, analysis_fs)
  dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)

  group <- factor(rep("all", ncol(counts_sub)))
  y <- DGEList(counts = counts_sub, group = group)
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~ x_scaled)

  keep_feat <- filterByExpr(y, design = design)
  cat("[QC filterByExpr] keep_edgeR (continuous) =", sum(keep_feat), " / ", length(keep_feat), "\n", sep = "")
  if (!any(keep_feat)) {
    cat("SKIP:", analysis_id, "-> filterByExpr eliminó todas las features (continuous).\n\n")
    return(NULL)
  }
  y <- y[keep_feat, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  tt <- topTags(qlf, n = Inf)$table
  tt$feature_id <- rownames(tt); rownames(tt) <- NULL
  if (!is.null(anno_df)) tt <- merge(tt, anno_df, by.x = "feature_id", by.y = "id", all.x = TRUE)
  tt <- force_numeric_df(tt, c("PValue", "FDR", "logFC", "logCPM", "F"))

  n_all <- nrow(tt)
  n_fdr05 <- sum(!is.na(tt$FDR) & tt$FDR < 0.05)
  n_fdr10 <- sum(!is.na(tt$FDR) & tt$FDR < 0.10)
  cat("   [DEBUG] ", analysis_id, " ", tag, " n_all=", n_all,
      " FDR<0.05=", n_fdr05, " FDR<0.1=", n_fdr10, "\n", sep = "")

  all_fp <- file.path(analysis_dir, paste0("DE_", analysis_fs, "_", tag, "_all_", ts, ".tsv"))
  sig05_fp <- file.path(analysis_dir, paste0("DE_", analysis_fs, "_", tag, "_sig_FDR0.05_", ts, ".tsv"))
  sig10_fp <- file.path(analysis_dir, paste0("DE_", analysis_fs, "_", tag, "_sig_FDR0.1_", ts, ".tsv"))
  write.table(tt, all_fp, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(tt[!is.na(tt$FDR) & tt$FDR < 0.05, , drop = FALSE], sig05_fp, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(tt[!is.na(tt$FDR) & tt$FDR < 0.10, , drop = FALSE], sig10_fp, sep = "\t", quote = FALSE, row.names = FALSE)

  png(file.path(analysis_dir, paste0("MDplot_", analysis_fs, "_", tag, "_", ts, ".png")),
      width = 1400, height = 1000, res = 150)
  plotMD(qlf, status = decide_status_safe(qlf),
         main = paste0("MD plot: ", analysis_id, " (continuous ", tag, ")"))
  abline(h = c(-1, 1), lty = 2)
  dev.off()

  plot_volcano(tt, analysis_fs, tag, analysis_dir, ts)

  qc <- data.frame(
    analysis_id = analysis_id,
    n_samples = ncol(y),
    n_features_tested = nrow(y),
    groups = paste(names(table(group)), as.integer(table(group)), sep = ":", collapse = ";"),
    n_sig_FDR0.05 = n_fdr05,
    n_sig_FDR0.1 = n_fdr10,
    comparison = tag,
    stringsAsFactors = FALSE
  )
  summary_fp <- file.path(analysis_dir, paste0("QC_", analysis_fs, "_summary_", ts, ".tsv"))
  write.table(qc, summary_fp, sep = "\t", quote = FALSE, row.names = FALSE)

  cat("Output (analysis_id=", analysis_id, ", continuous):\n", sep = "")
  cat("- ", all_fp, "\n", sep = "")
  cat("- ", sig05_fp, "\n", sep = "")
  cat("- ", sig10_fp, "\n\n", sep = "")

  list(qc = qc, files = list(all = all_fp, sig05 = sig05_fp, sig10 = sig10_fp))
}

run_spec_analysis <- function(spec_row) {
  norm_string <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x)] <- ""
    x[x == "NA"] <- ""
    x
  }
  analysis_id <- norm_string(spec_row$analysis_id)
  mode <- tolower(norm_string(spec_row$mode))
  var_key <- norm_string(spec_row$var)
  case_label <- norm_string(spec_row$case)
  ref_label <- norm_string(spec_row$ref)
  cut_val <- safe_numeric(spec_row$cut)[1]
  cut_type <- tolower(norm_string(spec_row$cut_type))
  collapse_map_str <- norm_string(spec_row$collapse_map)
  exclude_rule <- norm_string(spec_row$exclude_rule)
  event_var <- norm_string(spec_row$event_var)
  exclude_samples_str <- ""
  if ("exclude_samples" %in% names(spec_row)) {
    exclude_samples_str <- norm_string(spec_row$exclude_samples)
  }

  analysis_fs <- safe_name(analysis_id)
  meta_use <- meta_base
  counts_use <- counts_mat
  if (nzchar(exclude_samples_str)) {
    exclude_ids <- trimws(unlist(strsplit(exclude_samples_str, "[,;| ]+")))
    exclude_ids <- exclude_ids[nzchar(exclude_ids)]
    if (length(exclude_ids) > 0) {
      excl_mask <- meta_use$id %in% exclude_ids
      if (any(excl_mask)) {
        cat("Excluyendo muestras (spec):", paste(meta_use$id[excl_mask], collapse = ", "), "\n")
        meta_use <- meta_use[!excl_mask, , drop = FALSE]
        counts_use <- counts_use[, !excl_mask, drop = FALSE]
      } else {
        cat("WARNING: exclude_samples no coincide con meta$id\n")
      }
    }
  }
  n_before <- nrow(meta_use)

  cat("============================================\n")
  cat("[SPEC] analysis_id=", analysis_id,
      " mode=", mode,
      " var=", var_key,
      " case=", case_label,
      " ref=", ref_label,
      " cut=", ifelse(is.na(cut_val), "NA", cut_val),
      " cut_type=", ifelse(nzchar(cut_type), cut_type, "NA"),
      " collapse_map=", ifelse(nzchar(collapse_map_str), collapse_map_str, "NA"),
      " exclude_rule=", ifelse(nzchar(exclude_rule), exclude_rule, "NA"),
      "\n", sep = "")

  if (!(var_key %in% colnames(meta_use))) {
    cat("SKIP:", analysis_id, "-> variable no existe en metadata:", var_key, "\n\n")
    return(NULL)
  }
  if (!nzchar(mode)) {
    cat("SKIP:", analysis_id, "-> modo vacío.\n\n")
    return(NULL)
  }

  group <- NULL
  keep_mask <- rep(TRUE, nrow(meta_use))

  if (identical(mode, "continuous")) {
    x_num <- safe_numeric(meta_use[[var_key]])
    keep_mask <- !is.na(x_num)
    if (sum(keep_mask) < 3) {
      cat("SKIP:", analysis_id, "-> continuous sin suficientes muestras válidas.\n\n")
      return(NULL)
    }
    counts_sub <- counts_use[, keep_mask, drop = FALSE]
    x_scaled <- x_num[keep_mask]
    scale_val <- if (is.finite(cut_val) && cut_val != 0) cut_val else 1
    x_scaled <- x_scaled / scale_val
    tag <- paste0("continuous_per", scale_val)
    res <- run_edgeR_continuous(
      analysis_id = analysis_id,
      analysis_fs = analysis_fs,
      counts_sub = counts_sub,
      x_scaled = x_scaled,
      anno_df = anno_df,
      ts = ts,
      tag = tag
    )
    if (is.null(res)) return(NULL)
    qc_row <- res$qc
    qc_row$mode <- mode
    qc_row$var <- var_key
    qc_row$case <- ""
    qc_row$ref <- ""
    qc_row$cut <- scale_val
    qc_row$cut_type <- ""
    qc_row$exclude_rule <- exclude_rule
    qc_row$event_var <- event_var
    qc_row <- finalize_qc(qc_row)
    return(qc_row)

  } else if (identical(mode, "range_vs")) {
    ranges <- parse_range_map(collapse_map_str)
    if (length(ranges) == 0) {
      cat("SKIP:", analysis_id, "-> collapse_map (rangos) vacío para range_vs.\n\n")
      return(NULL)
    }
    x_num <- as_numeric_robust(meta_use[[var_key]])
    mapped <- rep(NA_character_, length(x_num))
    for (rng in ranges) {
      hit <- !is.na(x_num) & x_num >= rng$min & x_num <= rng$max
      mapped[hit] <- rng$label
    }
    if (identical(exclude_rule, "drop_unmapped")) {
      mapped[is.na(mapped)] <- NA
    }
    group_raw <- mapped
    if (!nzchar(case_label) || !nzchar(ref_label)) {
      cat("SKIP:", analysis_id, "-> debes indicar case y ref para range_vs.\n\n")
      return(NULL)
    }
    keep_mask <- !is.na(group_raw) & group_raw %in% c(case_label, ref_label)
    group <- droplevels(factor(group_raw[keep_mask]))
    if (nlevels(group) < 2) {
      cat("SKIP:", analysis_id, "-> range_vs sin 2 niveles tras filtrar.\n\n")
      return(NULL)
    }
    group <- stats::relevel(group, ref = ref_label)
    tab <- table(group)
    if (any(tab < 2)) {
      cat("SKIP:", analysis_id, "-> range_vs algún nivel <2 muestras.\n\n")
      return(NULL)
    }
    counts_sub <- counts_use[, keep_mask, drop = FALSE]
    meta_sub <- meta_use[keep_mask, , drop = FALSE]
    names(group) <- colnames(counts_sub)

    res <- run_edgeR_binary(
      analysis_id = analysis_id,
      analysis_fs = analysis_fs,
      group = group,
      counts_sub = counts_sub,
      anno_df = anno_df,
      ts = ts,
      case_label = case_label,
      ref_label = ref_label
    )
    if (is.null(res)) return(NULL)
    qc_row <- res$qc
    qc_row$mode <- mode
    qc_row$var <- var_key
    qc_row$case <- case_label
    qc_row$ref <- ref_label
    qc_row$cut <- cut_val
    qc_row$cut_type <- cut_type
    qc_row$collapse_map <- collapse_map_str
    qc_row$exclude_rule <- exclude_rule
    qc_row$event_var <- event_var
    qc_row <- finalize_qc(qc_row)
    return(qc_row)

  } else if (identical(mode, "as_is")) {
    if (!nzchar(case_label) || !nzchar(ref_label)) {
      cat("SKIP:", analysis_id, "-> debes indicar case y ref para as_is.\n\n")
      return(NULL)
    }
    v <- meta_use[[var_key]]
    group_raw <- as.factor(v)
    keep_mask <- !is.na(group_raw) & trimws(as.character(group_raw)) != ""
    group_chr <- as.character(group_raw)
    keep_mask <- keep_mask & (group_chr %in% c(case_label, ref_label))
    group <- factor(group_chr[keep_mask])
    group <- droplevels(group)
  } else if (identical(mode, "binary_cut")) {
    if (!is.finite(cut_val)) {
      cat("SKIP:", analysis_id, "-> cut no es numérico.\n\n")
      return(NULL)
    }
    if (!cut_type %in% c("gt", "ge")) {
      cat("SKIP:", analysis_id, "-> cut_type debe ser gt o ge.\n\n")
      return(NULL)
    }
    if (!nzchar(case_label) || !nzchar(ref_label)) {
      cat("SKIP:", analysis_id, "-> debes indicar case y ref para binary_cut.\n\n")
      return(NULL)
    }
    x <- safe_numeric(meta_use[[var_key]])
    keep_mask <- !is.na(x)
    if (cut_type == "gt") {
      group_chr <- ifelse(x > cut_val, case_label, ref_label)
    } else {
      group_chr <- ifelse(x >= cut_val, case_label, ref_label)
    }
    group <- factor(group_chr[keep_mask])
    group <- droplevels(group)
  } else if (identical(mode, "collapse_levels")) {
    if (!nzchar(case_label) || !nzchar(ref_label)) {
      cat("SKIP:", analysis_id, "-> debes indicar case y ref para collapse_levels.\n\n")
      return(NULL)
    }
    cmap <- parse_collapse_levels(collapse_map_str)
    if (length(cmap) == 0) {
      cat("SKIP:", analysis_id, "-> collapse_map vacío o malformado.\n\n")
      return(NULL)
    }
    v <- as.character(meta_use[[var_key]])
    mapped <- rep(NA_character_, length(v))
    for (nm in names(cmap)) {
      mapped[!is.na(v) & v == nm] <- cmap[[nm]]
    }
    keep_mask <- !is.na(mapped)
    keep_mask <- keep_mask & mapped %in% c(case_label, ref_label)
    group <- factor(mapped[keep_mask])
    group <- droplevels(group)
  } else if (identical(mode, "survival_cut")) {
    if (!nzchar(case_label) || !nzchar(ref_label)) {
      cat("SKIP:", analysis_id, "-> debes indicar case y ref para survival_cut.\n\n")
      return(NULL)
    }
    if (!is.finite(cut_val)) {
      cat("SKIP:", analysis_id, "-> cut no es numérico para survival_cut.\n\n")
      return(NULL)
    }
    if (!(event_var %in% colnames(meta_use))) {
      cat("SKIP:", analysis_id, "-> event_var no existe en metadata: ", event_var, "\n\n", sep = "")
      return(NULL)
    }
    time <- safe_numeric(meta_use[[var_key]])
    event <- safe_numeric(meta_use[[event_var]])
    g <- rep(NA_character_, length(time))
    lt_idx <- which(!is.na(time) & !is.na(event) & time < cut_val & event == 1)
    ge_idx <- which(!is.na(time) & time >= cut_val)
    g[lt_idx] <- ref_label
    g[ge_idx] <- case_label
    if (identical(exclude_rule, "exclude_if_censored_lt_cut")) {
      cens_lt <- which(!is.na(time) & !is.na(event) & event == 0 & time < cut_val)
      if (length(cens_lt) > 0) g[cens_lt] <- NA
    }
    keep_mask <- !is.na(g)
    group <- factor(g[keep_mask])
    group <- droplevels(group)
  } else {
    cat("SKIP:", analysis_id, "-> modo no soportado:", mode, "\n\n")
    return(NULL)
  }

  n_after <- sum(keep_mask)
  cat("Samples antes:", n_before, " | después de exclusiones:", n_after, "\n")

  if (length(group) == 0 || n_after == 0) {
    cat("SKIP:", analysis_id, "-> sin muestras tras aplicar reglas.\n\n")
    return(NULL)
  }

  if (!(ref_label %in% levels(group)) || !(case_label %in% levels(group))) {
    cat("SKIP:", analysis_id, "-> niveles no encontrados tras filtrar. levels=", paste(levels(group), collapse = ","), "\n\n")
    return(NULL)
  }

  group <- stats::relevel(group, ref = ref_label)
  tab <- table(group)
  cat("Tabla grupos:\n"); print(tab); cat("\n")

  if (any(tab < 2)) {
    cat("SKIP:", analysis_id, "-> algún nivel tiene <2 muestras.\n\n")
    return(NULL)
  }

  counts_sub <- counts_use[, keep_mask, drop = FALSE]
  meta_sub <- meta_use[keep_mask, , drop = FALSE]
  if (ncol(counts_sub) != length(group)) {
    stop("Inconsistencia interna: counts_sub y group longitud distinta.")
  }
  names(group) <- colnames(counts_sub)

  res <- run_edgeR_binary(
    analysis_id = analysis_id,
    analysis_fs = analysis_fs,
    group = group,
    counts_sub = counts_sub,
    anno_df = anno_df,
    ts = ts,
    case_label = case_label,
    ref_label = ref_label
  )
  if (is.null(res)) return(NULL)

  qc_row <- res$qc
  qc_row$mode <- mode
  qc_row$var <- var_key
  qc_row$case <- case_label
  qc_row$ref <- ref_label
  qc_row$cut <- cut_val
  qc_row$cut_type <- cut_type
  qc_row$collapse_map <- collapse_map_str
  qc_row$exclude_rule <- exclude_rule
  qc_row$event_var <- event_var
  qc_row <- finalize_qc(qc_row)
  qc_row
}

# -----------------------------
# Helper: correr DE para una variable
# -----------------------------
pick_reference <- function(f) {
  f <- droplevels(f)
  if (nlevels(f) < 2) return(f)
  if (ref_mode == "last") {
    stats::relevel(f, ref = tail(levels(f), 1))
  } else {
    stats::relevel(f, ref = levels(f)[1])
  }
}

run_de_for_var <- function(varname,
                           analysis_id = NULL,
                           mode = "factor",
                           case_level = NULL,
                           ref_level = NULL,
                           collapse_map = NULL,
                           exclude_rule = NULL,
                           cut = NULL) {
  var_key <- varname
  analysis_id <- if (is.null(analysis_id) || nchar(analysis_id) == 0) var_key else analysis_id
  analysis_fs <- safe_name(analysis_id)
  cat("[INFO] var_key=", var_key, " analysis_id=", analysis_id, " analysis_fs=", analysis_fs, " mode=", mode, "\n", sep = "")

  if (!(var_key %in% colnames(meta_base))) {
    cat("SKIP: variable no existe en metadata:", var_key, "\n\n")
    return(NULL)
  }

  v <- meta_base[[var_key]]

  if (var_key %in% names(include_map)) {
    allowed <- include_map[[var_key]]
    v[!(as.character(v) %in% allowed)] <- NA
  }
  if (var_key %in% names(exclude_map)) {
    banned <- exclude_map[[var_key]]
    v[as.character(v) %in% banned] <- NA
  }

  # manejar modos especiales
  if (mode == "range_vs") {
    x_num <- as_numeric_robust(v)
    ranges <- parse_range_map(collapse_map)
    mapped <- rep(NA_character_, length(x_num))
    if (length(ranges) == 0) {
      cat("SKIP:", var_key, "-> collapse_map vacío para range_vs.\n\n")
      return(NULL)
    }
    for (rng in ranges) {
      hit <- !is.na(x_num) & x_num >= rng$min & x_num <= rng$max
      mapped[hit] <- rng$label
    }
    if (!is.null(exclude_rule) && exclude_rule == "drop_unmapped") {
      mapped[is.na(mapped)] <- NA
    }
    v <- mapped
  }

  # modo continuo
  if (mode == "continuous") {
    x_num <- as_numeric_robust(v)
    keep <- !is.na(x_num)
    if (sum(keep) < 3) {
      cat("SKIP:", var_key, "-> continuous sin suficientes muestras válidas.\n\n")
      return(NULL)
    }
    x_num <- x_num[keep]
    counts_sub <- counts_mat[, keep, drop = FALSE]
    meta_sub <- meta_base[keep, , drop = FALSE]
    group <- factor(rep("all", ncol(counts_sub)))
    y <- DGEList(counts = counts_sub, group = group)
    y <- calcNormFactors(y, method = "TMM")
    cut_val <- suppressWarnings(as.numeric(cut))
    if (is.na(cut_val) || cut_val == 0) cut_val <- 1
    x_scaled <- x_num / cut_val
    design <- model.matrix(~ x_scaled)
    keep_feat <- filterByExpr(y, design = design)
    cat("   [DEBUG] ", var_key, " continuous filterByExpr keep=", sum(keep_feat), " / ", length(keep_feat), "\n", sep = "")
    y <- y[keep_feat, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y, method = "TMM")
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = 2)
    tt <- topTags(qlf, n = Inf)$table
    tt$feature_id <- rownames(tt); rownames(tt) <- NULL
    if (!is.null(anno_df)) tt <- merge(tt, anno_df, by.x = "feature_id", by.y = "id", all.x = TRUE)
    tt <- force_numeric(tt, c("PValue", "FDR", "logFC", "logCPM", "F"))
    tag <- paste0("continuous_per", cut_val)

    out_fp <- file.path(tab_var_dir, paste0("DE_", analysis_fs, "_", tag, "_all_", ts, ".tsv"))
    write.table(tt, out_fp, sep = "\t", quote = FALSE, row.names = FALSE)

    sig_fp <- file.path(tab_var_dir, paste0("DE_", analysis_fs, "_", tag, "_sig_FDR0.05_", ts, ".tsv"))
    tt_sig <- tt[!is.na(tt$FDR) & tt$FDR < 0.05, , drop = FALSE]
    write.table(tt_sig, sig_fp, sep = "\t", quote = FALSE, row.names = FALSE)
    sig_fp10 <- file.path(tab_var_dir, paste0("DE_", analysis_fs, "_", tag, "_sig_FDR0.1_", ts, ".tsv"))
    tt_sig10 <- tt[!is.na(tt$FDR) & tt$FDR < 0.10, , drop = FALSE]
    write.table(tt_sig10, sig_fp10, sep = "\t", quote = FALSE, row.names = FALSE)

    png(file.path(fig_var_dir, paste0("MDplot_", analysis_fs, "_", tag, "_", ts, ".png")),
        width = 1400, height = 1000, res = 150)
    plotMD(qlf, status = decide_status(qlf, coef_idx = 2),
           main = paste0("MD plot: ", var_key, " ", tag))
    abline(h = c(-1, 1), lty = 2)
    dev.off()

    plot_volcano(tt, analysis_fs, tag, fig_var_dir, ts)

    qc <- data.frame(
      analysis = analysis_id,
      variable = var_key,
      mode = mode,
      n_samples = ncol(y),
      n_features_tested = nrow(y),
      levels = "",
      group_counts = "",
      comparison = tag,
      stringsAsFactors = FALSE
    )
    write.table(qc, summary_fp <- file.path(tab_var_dir, paste0("QC_", analysis_fs, "_summary_", ts, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Output tablas:", tab_var_dir, "\n")
    cat("Output figuras:", fig_var_dir, "\n\n")
    return(qc)
  }

  # tratar como factor (univariada)
  if (is.numeric(v)) {
    if (length(unique(v[!is.na(v)])) > 10 && mode == "factor") {
      cat("WARNING:", var_key, "parece continua (numérica con >10 valores). Aquí se forzará a factor.\n")
    }
  }
  group <- as.factor(v)
  # eliminar NA del grupo para este análisis
  keep <- !is.na(group) & trimws(as.character(group)) != ""
  if (sum(!keep) > 0) {
    cat("INFO:", var_key, "-> removiendo", sum(!keep), "muestras con NA/empty.\n")
  }

  group <- droplevels(group[keep])
  counts_sub <- counts_mat[, keep, drop = FALSE]
  meta_sub <- meta_base[keep, , drop = FALSE]

  # Contraste explícito (binario) si aplica
  contrast_mode <- FALSE
  tag_contrast <- NULL
  # contraste explícito por CLI spec
  if (var_key %in% names(contrasts_map) && is.null(case_level) && is.null(ref_level)) {
    case_level <- contrasts_map[[var_key]]$case
    ref_level  <- contrasts_map[[var_key]]$ref
  }
  if (!is.null(case_level) && !is.null(ref_level)) {
    contrast_mode <- TRUE
    target_levels <- c(case_level, ref_level)
    keep_contrast <- as.character(group) %in% target_levels
    dropped_contrast <- sum(!keep_contrast)
    if (dropped_contrast > 0) {
      cat("INFO:", var_key, "-> removiendo", dropped_contrast, "muestras fuera de contraste explícito.\n")
    }
    group <- droplevels(group[keep_contrast])
    counts_sub <- counts_sub[, keep_contrast, drop = FALSE]
    meta_sub <- meta_sub[keep_contrast, , drop = FALSE]
    if (nlevels(group) < 2) {
      cat("SKIP:", var_key, "-> contraste explícito sin 2 niveles tras filtrar.\n\n")
      return(NULL)
    }
    if (!(ref_level %in% levels(group)) || !(case_level %in% levels(group))) {
      cat("SKIP:", var_key, "-> niveles de contraste no presentes tras filtrado.\n\n")
      return(NULL)
    }
    group <- stats::relevel(group, ref = ref_level)
    tag_contrast <- paste0(case_level, "_vs_", ref_level)
  }

  if (nlevels(group) < 2) {
    cat("SKIP:", var_key, "-> <2 niveles tras limpiar NA.\n\n")
    return(NULL)
  }

  # replicación mínima
  tab <- table(group)
  if (any(tab < 2)) {
    cat("SKIP:", var_key, "-> algún nivel tiene <2 muestras.\n")
    print(tab); cat("\n")
    return(NULL)
  }

  if (!contrast_mode) {
    group <- pick_reference(group)
    tab <- table(group)
    ref_level <- levels(group)[1]
  } else {
    tab <- table(group)
  }

  # crear carpeta output
  tab_var_dir <- file.path(tab_root, paste0("DE_", analysis_fs))
  fig_var_dir <- file.path(fig_root, paste0("DE_", analysis_fs))
  dir.create(tab_var_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_var_dir, showWarnings = FALSE, recursive = TRUE)

  # edgeR pipeline QL
  y <- DGEList(counts = counts_sub, group = group)
  y <- calcNormFactors(y, method = "TMM")

  design <- model.matrix(~ group)
  # ---- QC: interpretar filterByExpr (umbral implícito) ----
  tab <- table(group)
  min_n <- min(tab)

  eff_lib_tmp <- y$samples$lib.size * y$samples$norm.factors
  cpm10 <- 10 / eff_lib_tmp * 1e6

  cat("[QC filterByExpr] Variable:", var_key, "\n")
  cat("[QC filterByExpr] Grupos:", paste(names(tab), tab, sep = ":", collapse = "; "), "\n")
  cat("[QC filterByExpr] min_n (grupo más pequeño):", min_n, "\n")

  qs <- stats::quantile(cpm10, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
  cat("[QC filterByExpr] cpm10 (10 reads equiv.) summary:\n")
  print(summary(cpm10))
  cat("[QC filterByExpr] cpm10 quantiles (0.05/0.25/0.5/0.75/0.95):\n")
  print(qs)
  cat("[QC filterByExpr] cpm10 range:", paste(range(cpm10, finite = TRUE), collapse = " - "), "\n")

  cpm_mat <- edgeR::cpm(y, normalized.lib.sizes = TRUE)

  keep_manual <- tryCatch({
    th <- as.numeric(stats::quantile(cpm10, 0.25, na.rm = TRUE))
    km <- rowSums(cpm_mat >= th) >= min_n
    cat("[QC filterByExpr] keep_manual (th=Q1(cpm10) =", signif(th, 4), ") -> keep=", sum(km), "\n")

    for (mc in c(5, 10, 20)) {
      cpm_mc <- mc / eff_lib_tmp * 1e6
      th_mc <- as.numeric(stats::quantile(cpm_mc, 0.25, na.rm = TRUE))
      keep_mc <- rowSums(cpm_mat >= th_mc) >= min_n
      cat("[QC filterByExpr] sensitivity min_counts=", mc,
          " th(Q1)=", signif(th_mc, 4),
          " keep=", sum(keep_mc), "\n", sep = "")
    }
    km
  }, error = function(e) {
    cat("[QC filterByExpr] WARNING: no se pudo calcular keep_manual: ", conditionMessage(e), "\n", sep = "")
    NULL
  })

  keep_feat <- filterByExpr(y, design = design)
  cat("[QC filterByExpr] keep_edgeR =", sum(keep_feat), " / ", length(keep_feat), "\n", sep = "")
  if (!is.null(keep_manual)) {
    cat("[QC filterByExpr] concordancia edgeR vs manual:\n")
    print(table(keep_edgeR = keep_feat, keep_manual = keep_manual))
  }
  cat("   [DEBUG] ", var_key, " filterByExpr keep=", sum(keep_feat), " / ", length(keep_feat), "\n", sep = "")
  y <- y[keep_feat, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")  # recalcular tras filtrar features

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)

  decide_status <- function(qlf_obj, coef_idx = NULL) {
    tryCatch({
      dt_fun <- get("decideTestsDGE", envir = asNamespace("edgeR"))
      dt_fun(qlf_obj, adjust.method = "BH", p.value = 0.05)
    }, error = function(e) {
      cat("WARNING: decideTestsDGE no disponible; plotMD sin status.\n")
      NULL
    })
  }

  force_numeric <- function(df, cols) {
    for (cc in cols) {
      if (cc %in% colnames(df)) {
        v <- df[[cc]]
        if (is.factor(v)) v <- as.character(v)
        df[[cc]] <- suppressWarnings(as.numeric(v))
      }
    }
    df
  }

  # Output summary
  summary_fp <- file.path(tab_var_dir, paste0("QC_", analysis_fs, "_summary_", ts, ".txt"))
  cat("== ", analysis_id, " (", var_key, ") ==\n", sep="")
  cat("N muestras:", ncol(y), "\n")
  cat("N features test:", nrow(y), "\n")
  cat("Grupos:\n"); print(tab)
  cat("\n")

  # Omnibus si >2 niveles
  res_list <- list()

  if (nlevels(group) > 2) {
    # test global para todos los coeficientes de grupo (excepto intercept)
    coef_idx <- 2:ncol(design)
    qlf_all <- glmQLFTest(fit, coef = coef_idx)
    tt <- topTags(qlf_all, n = Inf)$table
    tt$feature_id <- rownames(tt); rownames(tt) <- NULL
    if (!is.null(anno_df)) tt <- merge(tt, anno_df, by.x = "feature_id", by.y = "id", all.x = TRUE)
    out_fp <- file.path(tab_var_dir, paste0("DE_", analysis_fs, "_omnibus_all_", ts, ".tsv"))
    write.table(tt, out_fp, sep = "\t", quote = FALSE, row.names = FALSE)
    res_list[["omnibus"]] <- out_fp

    # pairwise vs referencia
    lev <- levels(group)
    ref <- lev[1]
    for (lv in lev[-1]) {
      # contraste: grouplv vs groupref => coef del diseño correspondiente
      # En model.matrix(~group), el coef para lv es "grouplv" (comparado con ref)
      coef_name <- paste0("group", lv)
      coef_idx <- which(colnames(design) == coef_name)
      if (length(coef_idx) != 1) next
      qlf <- glmQLFTest(fit, coef = coef_idx)
      tt <- topTags(qlf, n = Inf)$table
      tt$feature_id <- rownames(tt); rownames(tt) <- NULL
      if (!is.null(anno_df)) tt <- merge(tt, anno_df, by.x = "feature_id", by.y = "id", all.x = TRUE)
      tt <- force_numeric(tt, c("PValue", "FDR", "logFC", "logCPM", "F"))

      tag <- paste0(lv, "_vs_", ref)

      # checkpoints de conteos FDR
      n_all <- nrow(tt)
      n_fdr_na <- sum(is.na(tt$FDR))
      n_fdr05 <- sum(!is.na(tt$FDR) & tt$FDR < 0.05)
      n_fdr10 <- sum(!is.na(tt$FDR) & tt$FDR < 0.10)
      cat("   [DEBUG] ", tag, " n_all=", n_all, " FDR NA=", n_fdr_na,
          " FDR<0.05=", n_fdr05, " FDR<0.1=", n_fdr10, "\n", sep = "")
      out_fp <- file.path(tab_var_dir, paste0("DE_", analysis_fs, "_", tag, "_all_", ts, ".tsv"))
      write.table(tt, out_fp, sep = "\t", quote = FALSE, row.names = FALSE)

      sig_fp <- file.path(tab_var_dir, paste0("DE_", analysis_fs, "_", tag, "_sig_FDR0.05_", ts, ".tsv"))
      tt_sig <- tt[!is.na(tt$FDR) & tt$FDR < 0.05, , drop = FALSE]
      write.table(tt_sig, sig_fp, sep = "\t", quote = FALSE, row.names = FALSE)
      sig_fp10 <- file.path(tab_var_dir, paste0("DE_", analysis_fs, "_", tag, "_sig_FDR0.1_", ts, ".tsv"))
      tt_sig10 <- tt[!is.na(tt$FDR) & tt$FDR < 0.10, , drop = FALSE]
      write.table(tt_sig10, sig_fp10, sep = "\t", quote = FALSE, row.names = FALSE)

      # MD plot
      png(file.path(fig_var_dir, paste0("MDplot_", analysis_fs, "_", tag, "_", ts, ".png")),
          width = 1400, height = 1000, res = 150)
      plotMD(qlf, status = decide_status(qlf),
             main = paste0("MD plot: ", var_key, " ", tag))
      abline(h = c(-1, 1), lty = 2)
      dev.off()

      plot_volcano(tt, analysis_fs, tag, fig_var_dir, ts)

      res_list[[tag]] <- out_fp
    }

  } else {
    # binaria: un solo coeficiente (grupo2 vs ref)
    coef_idx <- 2
    qlf <- glmQLFTest(fit, coef = coef_idx)
    tt <- topTags(qlf, n = Inf)$table
    tt$feature_id <- rownames(tt); rownames(tt) <- NULL
    if (!is.null(anno_df)) tt <- merge(tt, anno_df, by.x = "feature_id", by.y = "id", all.x = TRUE)
    tt <- force_numeric(tt, c("PValue", "FDR", "logFC", "logCPM", "F"))

    lev <- levels(group); tag <- paste0(lev[2], "_vs_", lev[1])

    # checkpoints de conteos FDR
    n_all <- nrow(tt)
    n_fdr_na <- sum(is.na(tt$FDR))
    n_fdr05 <- sum(!is.na(tt$FDR) & tt$FDR < 0.05)
    n_fdr10 <- sum(!is.na(tt$FDR) & tt$FDR < 0.10)
    cat("   [DEBUG] ", tag, " n_all=", n_all, " FDR NA=", n_fdr_na,
        " FDR<0.05=", n_fdr05, " FDR<0.1=", n_fdr10, "\n", sep = "")

    out_fp <- file.path(tab_var_dir, paste0("DE_", analysis_fs, "_", tag, "_all_", ts, ".tsv"))
    write.table(tt, out_fp, sep = "\t", quote = FALSE, row.names = FALSE)

    sig_fp <- file.path(tab_var_dir, paste0("DE_", analysis_fs, "_", tag, "_sig_FDR0.05_", ts, ".tsv"))
    tt_sig <- tt[!is.na(tt$FDR) & tt$FDR < 0.05, , drop = FALSE]
    write.table(tt_sig, sig_fp, sep = "\t", quote = FALSE, row.names = FALSE)
    sig_fp10 <- file.path(tab_var_dir, paste0("DE_", analysis_fs, "_", tag, "_sig_FDR0.1_", ts, ".tsv"))
    tt_sig10 <- tt[!is.na(tt$FDR) & tt$FDR < 0.10, , drop = FALSE]
    write.table(tt_sig10, sig_fp10, sep = "\t", quote = FALSE, row.names = FALSE)

    # MD plot
    png(file.path(fig_var_dir, paste0("MDplot_", analysis_fs, "_", tag, "_", ts, ".png")),
        width = 1400, height = 1000, res = 150)
    plotMD(qlf, status = decide_status(qlf),
           main = paste0("MD plot: ", var_key, " ", tag))
    abline(h = c(-1, 1), lty = 2)
    dev.off()

    plot_volcano(tt, analysis_fs, tag, fig_var_dir, ts)

    res_list[[tag]] <- out_fp
  }

  # guardar resumen QC por variable
  qc <- data.frame(
    variable = var_key,
    n_samples = ncol(y),
    n_features_tested = nrow(y),
    levels = paste(names(tab), tab, sep=":", collapse=";"),
    group_counts = paste(names(tab), tab, sep=":", collapse=";"),
    comparison = if (contrast_mode && !is.null(case_level) && !is.null(ref_level)) {
      paste0(case_level, " vs ", ref_level, " (explicito)")
    } else {
      other_levels <- setdiff(levels(group), ref_level)
      paste0("ref=", ref_level, "; otros=", paste(other_levels, collapse=","), " (auto ref_mode=", ref_mode, ")")
    },
    stringsAsFactors = FALSE
  )
  write.table(qc, summary_fp, sep = "\t", quote = FALSE, row.names = FALSE)

  cat("Output tablas:", tab_var_dir, "\n")
  cat("Output figuras:", fig_var_dir, "\n\n")
  invisible(qc)
}

# -----------------------------
# Ejecutar todas las variables
# -----------------------------
spec_df <- NULL
if (!is.null(spec_fp)) {
  if (!file.exists(spec_fp)) stop("No existe spec: ", spec_fp)
  spec_df <- read.csv(spec_fp, check.names = FALSE, stringsAsFactors = FALSE)
  required_cols <- c("analysis_id", "mode", "var", "case", "ref",
                     "cut", "cut_type", "event_var", "collapse_map", "exclude_rule")
  missing_cols <- setdiff(required_cols, colnames(spec_df))
  if (length(missing_cols) > 0) {
    stop("El archivo spec carece de columnas requeridas: ", paste(missing_cols, collapse = ", "))
  }
  if (anyDuplicated(spec_df$analysis_id)) {
    dup <- spec_df$analysis_id[duplicated(spec_df$analysis_id)]
    stop("analysis_id debe ser único. Duplicados: ", paste(unique(dup), collapse = ", "))
  }
  if (length(vars) > 0) {
    cat("INFO: --spec proporcionado; se ignorará --vars.\n")
  }
  cat("Modo spec activado. Se ejecutarán ", nrow(spec_df), " análisis.\n", sep = "")
  cat("IDs:", paste(spec_df$analysis_id, collapse = ", "), "\n\n")
}

if (!is.null(spec_df)) {
  qc_list <- list()
  for (i in seq_len(nrow(spec_df))) {
    aid <- as.character(spec_df$analysis_id[i])
    qc_list[[aid]] <- tryCatch(
      run_spec_analysis(spec_df[i, , drop = FALSE]),
      error = function(e) { cat("ERROR en", aid, ":", conditionMessage(e), "\n\n"); NULL }
    )
  }
  qc_df <- do.call(rbind, qc_list)
  if (!is.null(qc_df) && nrow(qc_df) > 0) {
    summary_all_fp <- file.path(outdir, paste0("DE_summary_all_", ts, ".tsv"))
    write.table(qc_df, summary_all_fp, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Resumen global (spec):", summary_all_fp, "\n")
  }
} else {
  qc_list <- list()
  for (v in vars) {
    qc_list[[v]] <- tryCatch(run_de_for_var(v),
                             error = function(e) { cat("ERROR en", v, ":", conditionMessage(e), "\n\n"); NULL })
  }

  qc_df <- do.call(rbind, qc_list)
  if (!is.null(qc_df) && nrow(qc_df) > 0) {
    summary_all_fp <- file.path(tab_root, paste0("DE_summary_all_", ts, ".tsv"))
    write.table(qc_df, summary_all_fp, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Resumen global:", summary_all_fp, "\n")
  }
}

cat("\nFIN. Log:", log_fp, "\n")

sink(type="message"); sink(type="output"); close(zz)
