#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[[hit + 1]]
}
as_num <- function(x) {
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
latest_file <- function(dir_path, pattern) {
  ff <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(ff) == 0) return(NA_character_)
  ff[order(file.info(ff)$mtime, decreasing = TRUE)][1]
}
read_table <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(path, data.table = FALSE, check.names = FALSE)
  } else {
    read.csv(path, check.names = FALSE)
  }
}

# ---- Inputs / Params ----
counts_fp <- get_arg("--counts", file.path("data","intermediate","Gliomas_all_counts_merged.csv"))
meta_fp   <- get_arg("--meta",   file.path("data","intermediate","Metadatos_gliomas_verificados.csv"))
spec_fp   <- get_arg("--spec",   file.path("config","de_specs.csv"))
de_root   <- get_arg("--de_root", file.path("results","DE"))

drop_re   <- get_arg("--drop_regex", "^A")
min_lib   <- as_num(get_arg("--min_libsize", "100000"))
prior_cnt <- as_num(get_arg("--prior_count", "2"))

time_var  <- get_arg("--time_var", "MESES_SEGUIMIENTO_")
event_var <- get_arg("--event_var", "MUERTE")

km_cut    <- tolower(get_arg("--km_cut", "median"))  # median (por ahora)
min_group_n <- as.integer(get_arg("--min_group_n", "5"))

# Candidate selection from DE
fdr_cut <- as_num(get_arg("--fdr_cut", "0.1"))

# Outputs
fig_dir <- file.path("results","figures","KM_DE_candidates")
tab_dir <- file.path("results","tables","KM_DE_candidates")
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_fp <- file.path("logs", paste0("KM_DE_candidates_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type="output"); sink(zz, type="message")

cat("============================================\n")
cat("Kaplan–Meier (alto vs bajo) para RNAs candidatos DE\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("Counts:", counts_fp, "\n")
cat("Meta:", meta_fp, "\n")
cat("Spec:", spec_fp, "\n")
cat("DE root:", de_root, "\n")
cat("drop_regex:", drop_re, "\n")
cat("min_libsize:", min_lib, "\n")
cat("prior_count:", prior_cnt, "\n")
cat("time_var:", time_var, "\n")
cat("event_var:", event_var, "\n")
cat("km_cut:", km_cut, "\n")
cat("min_group_n:", min_group_n, "\n")
cat("FDR candidates:", fdr_cut, "\n")
cat("Fig dir:", fig_dir, "\n")
cat("Tab dir:", tab_dir, "\n")
cat("Log:", log_fp, "\n")
cat("============================================\n\n")

suppressPackageStartupMessages({
  library(edgeR)
  library(survival)
})

use_survminer <- requireNamespace("survminer", quietly = TRUE)
if (use_survminer) {
  suppressPackageStartupMessages(library(survminer))
  cat("survminer: disponible (usando ggsurvplot)\n")
} else {
  cat("survminer: NO disponible (usando plot base)\n")
}

# ---- Load counts + meta ----
if (!file.exists(counts_fp)) stop("No existe counts: ", counts_fp)
if (!file.exists(meta_fp)) stop("No existe meta: ", meta_fp)
if (!file.exists(spec_fp)) stop("No existe spec: ", spec_fp)

df <- read_table(counts_fp)
meta <- read_table(meta_fp)
meta <- meta[!(is.na(meta$id) | trimws(as.character(meta$id)) == ""), , drop = FALSE]

# columnas de anotación
anno_candidates <- c("id","type","entrez_id","HGNC_symbol")
anno_present <- intersect(anno_candidates, colnames(df))
count_cols <- setdiff(colnames(df), anno_present)

# drop samples por regex (post-mortem A*)
drop_samples <- grep(drop_re, count_cols, value = TRUE)
if (length(drop_samples) > 0) {
  cat("Drop samples (regex):", paste(drop_samples, collapse=", "), "\n")
}
count_cols <- setdiff(count_cols, drop_samples)

# counts matrix
feature_id <- if ("id" %in% colnames(df)) df[["id"]] else seq_len(nrow(df))
counts_df <- df[, count_cols, drop = FALSE]
counts_df[] <- lapply(counts_df, as_num)
counts_mat <- as.matrix(counts_df)
storage.mode(counts_mat) <- "numeric"
rownames(counts_mat) <- make.unique(as.character(feature_id))
colnames(counts_mat) <- make.unique(colnames(counts_mat))

# limpieza mínima features
counts_mat <- counts_mat[rowSums(counts_mat, na.rm = TRUE) > 0, , drop = FALSE]
rv <- apply(counts_mat, 1, function(x) stats::var(x, na.rm = TRUE))
counts_mat <- counts_mat[is.finite(rv) & rv > 0, , drop = FALSE]
cat("Tras limpieza mínima:", nrow(counts_mat), "features x", ncol(counts_mat), "samples\n")

# filtro global por libsize
y0 <- DGEList(counts_mat)
lib0 <- y0$samples$lib.size
keep_lib <- lib0 >= min_lib
cat("Filtro lib_size >=", min_lib, "-> removidas:", sum(!keep_lib), "/", length(keep_lib), "\n")
if (sum(!keep_lib) > 0) print(names(lib0)[!keep_lib])
counts_mat <- counts_mat[, keep_lib, drop = FALSE]

# alinear meta con muestras finales
samples_final <- colnames(counts_mat)
missing_in_meta <- setdiff(samples_final, as.character(meta$id))
if (length(missing_in_meta) > 0) {
  cat("WARNING: muestras en counts sin meta$id. Se descartan de counts:\n",
      paste(missing_in_meta, collapse = ", "), "\n")
  keep_cols <- setdiff(samples_final, missing_in_meta)
  if (length(keep_cols) == 0) {
    stop("No quedan muestras tras descartar las ausentes en metadata.")
  }
  counts_mat <- counts_mat[, keep_cols, drop = FALSE]
  samples_final <- colnames(counts_mat)
}
meta_base <- meta[match(samples_final, as.character(meta$id)), , drop = FALSE]
stopifnot(all(as.character(meta_base$id) == samples_final))
cat("Muestras finales:", nrow(meta_base), "\n\n")

# validar variables de sobrevida
if (!(time_var %in% colnames(meta_base))) stop("time_var no existe en meta: ", time_var)
if (!(event_var %in% colnames(meta_base))) stop("event_var no existe en meta: ", event_var)

time <- as_num(meta_base[[time_var]])
event <- as_num(meta_base[[event_var]])
ok_surv <- is.finite(time) & is.finite(event)
event <- ifelse(event > 0, 1, 0)

cat("QC survival:\n")
cat("- time finite:", sum(is.finite(time)), "/", length(time), "\n")
cat("- event finite:", sum(is.finite(as_num(meta_base[[event_var]]))), "/", length(event), "\n")
cat("- ok_surv:", sum(ok_surv), "/", length(ok_surv), "\n")
cat("- eventos (MUERTE=1) en ok_surv:", sum(event[ok_surv] == 1), "\n\n")

# ---- Build logCPM (TMM) ----
y <- DGEList(counts_mat)
y <- calcNormFactors(y, method = "TMM")
logcpm <- cpm(y, log = TRUE, prior.count = prior_cnt)
cat("logCPM matrix:", nrow(logcpm), "features x", ncol(logcpm), "samples\n\n")

# ---- Collect candidate features from DE results (FDR < fdr_cut) ----
spec <- read.csv(spec_fp, check.names = FALSE)
if (!("analysis_id" %in% colnames(spec))) stop("spec sin analysis_id")
analysis_ids <- as.character(spec$analysis_id)

candidates <- character(0)
per_analysis_counts <- list()

for (aid in analysis_ids) {
  aid_fs <- safe_name(aid)
  de_dir <- file.path(de_root, aid_fs)
  if (!dir.exists(de_dir)) {
    cat("SKIP DE:", aid, "-> no existe dir:", de_dir, "\n")
    next
  }

  # preferir archivo sig_FDR0.1 más reciente; si no, usar all más reciente y filtrar
  sig_fp <- latest_file(de_dir, pattern = "sig_FDR0\\.1_.*\\.tsv$")
  if (is.na(sig_fp) || !file.exists(sig_fp)) {
    all_fp <- latest_file(de_dir, pattern = "_all_.*\\.tsv$")
    if (is.na(all_fp) || !file.exists(all_fp)) {
      cat("SKIP DE:", aid, "-> no hay TSV\n")
      next
    }
    de <- read.delim(all_fp, check.names = FALSE)
    if (!("feature_id" %in% colnames(de))) {
      cat("SKIP DE:", aid, "-> all sin feature_id\n")
      next
    }
    if (!("FDR" %in% colnames(de))) {
      cat("SKIP DE:", aid, "-> all sin FDR (no puedo filtrar)\n")
      next
    }
    de$FDR <- as_num(de$FDR)
    de_sig <- de[is.finite(de$FDR) & de$FDR < fdr_cut, , drop = FALSE]
  } else {
    de_sig <- read.delim(sig_fp, check.names = FALSE)
  }

  if (!("feature_id" %in% colnames(de_sig))) {
    cat("SKIP DE:", aid, "-> sig sin feature_id\n")
    next
  }

  feats <- unique(as.character(de_sig$feature_id))
  feats <- feats[!is.na(feats) & nzchar(feats)]
  per_analysis_counts[[aid]] <- length(feats)
  candidates <- unique(c(candidates, feats))

  cat("Candidates from", aid, ":", length(feats), "\n")
}

cat("\nTotal candidatos (union):", length(candidates), "\n")

cand_fp <- file.path(tab_dir, paste0("DE_candidates_FDR", fdr_cut, "_", ts, ".txt"))
writeLines(candidates, cand_fp)
cat("Export candidatos:", cand_fp, "\n\n")

# mantener solo candidatos presentes en logcpm
candidates_in <- intersect(candidates, rownames(logcpm))
cat("Candidatos presentes en logCPM:", length(candidates_in), "/", length(candidates), "\n\n")
if (length(candidates_in) == 0) stop("0 candidatos presentes en la matriz de expresión.")

# ---- KM loop ----
# Publication-quality PDF settings
pdf_fp <- file.path(fig_dir, paste0("KM_curves_DE_candidates_", ts, ".pdf"))
pdf(pdf_fp, width = 10, height = 8, family = "sans")

km_rows <- list()
n_plotted <- 0
n_skipped <- 0

for (fid in candidates_in) {
  expr <- as.numeric(logcpm[fid, ])
  ok <- ok_surv & is.finite(expr)

  if (sum(ok) < (2 * min_group_n)) {
    n_skipped <- n_skipped + 1
    next
  }

  expr_ok <- expr[ok]
  t_ok <- time[ok]
  e_ok <- event[ok]

  cut_val <- NA_real_
  cut_low <- NA_real_
  cut_high <- NA_real_
  n_dropped_middle <- 0
  keep_idx <- rep(TRUE, length(expr_ok))

  if (km_cut == "median") {
    cut_val <- stats::median(expr_ok, na.rm = TRUE)
    grp <- ifelse(expr_ok > cut_val, "High", "Low")
    grp <- factor(grp, levels = c("Low","High"))
  } else if (km_cut == "terciles_extremos") {
    cut_low <- as.numeric(stats::quantile(expr_ok, 1/3, na.rm = TRUE, type = 7))
    cut_high <- as.numeric(stats::quantile(expr_ok, 2/3, na.rm = TRUE, type = 7))
    grp <- ifelse(expr_ok <= cut_low, "Low",
                  ifelse(expr_ok >= cut_high, "High", NA))
    n_dropped_middle <- sum(is.na(grp))
    keep_idx <- !is.na(grp)
    grp <- factor(grp[keep_idx], levels = c("Low","High"))
  } else if (km_cut == "cuartiles_extremos") {
    cut_low <- as.numeric(stats::quantile(expr_ok, 0.25, na.rm = TRUE, type = 7))
    cut_high <- as.numeric(stats::quantile(expr_ok, 0.75, na.rm = TRUE, type = 7))
    grp <- ifelse(expr_ok <= cut_low, "Low",
                  ifelse(expr_ok >= cut_high, "High", NA))
    n_dropped_middle <- sum(is.na(grp))
    keep_idx <- !is.na(grp)
    grp <- factor(grp[keep_idx], levels = c("Low","High"))
  } else {
    cat("SKIP ", fid, ": km_cut no soportado: ", km_cut, "\n", sep = "")
    n_skipped <- n_skipped + 1
    next
  }

  expr_use <- expr_ok[keep_idx]
  t_use <- t_ok[keep_idx]
  e_use <- e_ok[keep_idx]

  n_low <- sum(grp == "Low")
  n_high <- sum(grp == "High")
  cat("KM groups ", fid, " | ", km_cut, " | nLow=", n_low,
      " nHigh=", n_high, " dropped=", n_dropped_middle, "\n", sep = "")

  if (min(n_low, n_high) < min_group_n) {
    cat("SKIP ", fid, ": grupos demasiado pequeños tras ", km_cut,
        " (nLow=", n_low, " nHigh=", n_high, " dropped=", n_dropped_middle, ")\n", sep = "")
    n_skipped <- n_skipped + 1
    next
  }

  df_km <- data.frame(time = t_use, event = e_use, expr = expr_use, group = grp)

  sf <- survival::survfit(survival::Surv(time, event) ~ group, data = df_km)
  sd <- survival::survdiff(survival::Surv(time, event) ~ group, data = df_km)
  p_lr <- 1 - stats::pchisq(sd$chisq, df = (length(sd$n) - 1))

  # mediana por grupo (si existe)
  med_low <- NA_real_; med_high <- NA_real_
  tt <- tryCatch(summary(sf)$table, error = function(e) NULL)
  if (!is.null(tt)) {
    # cuando hay strata: filas por grupo
    if (is.matrix(tt) && nrow(tt) >= 2 && "median" %in% colnames(tt)) {
      rn <- rownames(tt)
      # rn suele ser "group=Low" / "group=High"
      if (any(grepl("Low", rn))) med_low <- tt[grep("Low", rn)[1], "median"]
      if (any(grepl("High", rn))) med_high <- tt[grep("High", rn)[1], "median"]
    } else if (is.numeric(tt) && "median" %in% names(tt)) {
      # caso raro sin strata
    }
  }

  if (km_cut == "median") {
    cut_txt <- paste0("Median = ", sprintf("%.3f", cut_val))
    lab_low <- "Low Expression"
    lab_high <- "High Expression"
  } else if (km_cut == "terciles_extremos") {
    cut_txt <- paste0("T1 = ", sprintf("%.3f", cut_low), "; T2 = ", sprintf("%.3f", cut_high))
    lab_low <- "Low (T1)"
    lab_high <- "High (T3)"
  } else if (km_cut == "cuartiles_extremos") {
    cut_txt <- paste0("Q1 = ", sprintf("%.3f", cut_low), "; Q3 = ", sprintf("%.3f", cut_high))
    lab_low <- "Low (Q1)"
    lab_high <- "High (Q4)"
  } else {
    cut_txt <- "Cut = NA"
    lab_low <- "Low"
    lab_high <- "High"
  }

  # Simplified English title for publication

main_txt <- paste0(fid, " (n = ", n_low + n_high, ")")

  # Format p-value for display
  p_txt <- if (p_lr < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", sprintf("%.3f", p_lr))
  }

  if (use_survminer) {
    gp <- suppressMessages(suppressWarnings(survminer::ggsurvplot(
      sf, data = df_km,
      risk.table = TRUE,
      risk.table.col = "strata",
      risk.table.height = 0.25,
      pval = FALSE,  # We'll add p-value manually in top-right
      conf.int = FALSE,
      palette = c("#0072B2", "#D55E00"),  # Colorblind-safe
      legend.title = "Expression",
      legend.labs = c(lab_low, lab_high),
      title = main_txt,
      xlab = "Time (months)",
      ylab = "Survival Probability",
      font.main = c(18, "bold"),
      font.x = c(16, "bold"),
      font.y = c(16, "bold"),
      font.tickslab = c(14, "plain"),
      font.legend = c(14, "plain"),
      ggtheme = ggplot2::theme_classic(base_size = 14)
    )))
    # Add p-value annotation in upper-right area (80% of x-axis)
    x_pval <- max(df_km$time, na.rm = TRUE) * 0.80
    gp$plot <- gp$plot +
      ggplot2::annotate("text", x = x_pval, y = 0.82, label = p_txt,
                        hjust = 0.5, vjust = 0.5, size = 6, fontface = "bold")
    render_plot <- function() {
      suppressMessages(suppressWarnings({
        grid::grid.newpage()
        print(gp)
      }))
    }
  } else {
    render_plot <- function() {
      par(mar = c(5, 5, 4, 2) + 0.1, cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.5)
      plot(sf, main = main_txt, xlab = "Time (months)", ylab = "Survival Probability",
           lwd = 3, col = c("#0072B2", "#D55E00"))
      legend("bottomleft", legend = c(lab_low, lab_high), lwd = 3,
             col = c("#0072B2", "#D55E00"), bty = "n", cex = 1.3)
      # p-value in top-right corner
      legend("topright", legend = p_txt, bty = "n", cex = 1.4, text.font = 2)
    }
  }
  render_plot()

  # High-resolution PNG export
  png_fp <- file.path(fig_dir, paste0("KM_", safe_name(fid), "_", ts, ".png"))
  png(png_fp, width = 3000, height = 2400, res = 300)
  render_plot()
  dev.off()

  km_rows[[fid]] <- data.frame(
    feature_id = fid,
    n_total = sum(ok),
    n_low = n_low,
    n_high = n_high,
    n_dropped_middle = n_dropped_middle,
    cut_method = km_cut,
    cut_value = cut_val,
    cut_low = cut_low,
    cut_high = cut_high,
    logrank_p = p_lr,
    median_surv_low = med_low,
    median_surv_high = med_high,
    stringsAsFactors = FALSE
  )

  n_plotted <- n_plotted + 1
}

dev.off()

km_df <- do.call(rbind, km_rows)
if (is.null(km_df) || nrow(km_df) == 0) {
  cat("WARNING: no se generaron curvas KM (todo fue SKIP).\n")
  km_df <- data.frame()
} else {
  km_df$FDR_logrank <- p.adjust(km_df$logrank_p, method = "BH")
  km_df <- km_df[order(km_df$FDR_logrank, km_df$logrank_p), ]
}

km_fp <- file.path(tab_dir, paste0("KM_summary_", ts, ".tsv"))
write.table(km_df, km_fp, sep = "\t", quote = FALSE, row.names = FALSE)

cat("============================================\n")
cat("FIN\n")
cat("KM plotted:", n_plotted, "\n")
cat("KM skipped:", n_skipped, "\n")
cat("PDF:", pdf_fp, "\n")
cat("Summary:", km_fp, "\n")
cat("Candidates list:", cand_fp, "\n")
cat("Log:", log_fp, "\n")
cat("============================================\n")

sink(type="message"); sink(type="output"); close(zz)
