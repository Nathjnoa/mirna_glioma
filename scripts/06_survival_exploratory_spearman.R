#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

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

read_table <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(path, data.table = FALSE, check.names = FALSE)
  } else {
    read.delim(path, check.names = FALSE)
  }
}

as_num_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

# ---- CLI ----
counts_fp <- get_arg("--counts", file.path("data","intermediate","Gliomas_all_counts_merged.csv"))
meta_fp   <- get_arg("--meta",   file.path("data","intermediate","Metadatos_gliomas_verificados.csv"))
spec_fp   <- get_arg("--spec",   file.path("config","de_specs.csv"))
de_dir    <- get_arg("--de_dir", file.path("results","DE"))

drop_re   <- get_arg("--drop_regex", "^A")
min_lib   <- as.numeric(get_arg("--min_libsize", "100000"))
prior_cnt <- as.numeric(get_arg("--prior_count", "2"))

fdr_thr   <- as.numeric(get_arg("--fdr", "0.1"))
max_plots <- as.integer(get_arg("--max_plots", "80"))  # limita PNGs por sanidad

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

dir.create("logs", showWarnings = FALSE, recursive = TRUE)
tab_outdir <- file.path("results","tables","survival_exploratory")
fig_outdir <- file.path("results","figures","survival_exploratory", paste0("FDR", gsub("\\.","", sprintf("%.2g", fdr_thr))))
dir.create(tab_outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_outdir, showWarnings = FALSE, recursive = TRUE)

log_fp <- file.path("logs", paste0("survival_exploratory_spearman_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type = "output"); sink(zz, type = "message")

cat("============================================\n")
cat("Exploratorio: correlación expresión vs seguimiento (Spearman)\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("Counts:", counts_fp, "\n")
cat("Meta:", meta_fp, "\n")
cat("Spec:", spec_fp, "\n")
cat("DE dir:", de_dir, "\n")
cat("drop_regex:", drop_re, "\n")
cat("min_libsize:", min_lib, "\n")
cat("prior_count:", prior_cnt, "\n")
cat("FDR threshold (union):", fdr_thr, "\n")
cat("max_plots:", max_plots, "\n")
cat("Tables:", tab_outdir, "\n")
cat("Figures:", fig_outdir, "\n")
cat("Log:", log_fp, "\n")
cat("============================================\n\n")

# ---- deps ----
suppressPackageStartupMessages({
  library(edgeR)
})

# ---- read spec ----
if (!file.exists(spec_fp)) stop("No existe spec: ", spec_fp)
spec <- read_table(spec_fp)

required_cols <- c("analysis_id","mode","var","event_var","cut","cut_type","case","ref","collapse_map","exclude_rule")
missing_cols <- setdiff(required_cols, colnames(spec))
if (length(missing_cols) > 0) stop("Faltan columnas en spec: ", paste(missing_cols, collapse=", "))

spec$analysis_id <- as.character(spec$analysis_id)
if (anyDuplicated(spec$analysis_id)) {
  dups <- unique(spec$analysis_id[duplicated(spec$analysis_id)])
  stop("analysis_id duplicados en spec: ", paste(dups, collapse=", "))
}
cat("Analyses en spec:", nrow(spec), "\n\n")

# ---- read DE tables and build union of features (FDR < fdr_thr) ----
union_feats <- character(0)
membership <- list()

for (i in seq_len(nrow(spec))) {
  aid <- as.character(spec$analysis_id[i])
  aid_fs <- safe_name(aid)
  adir <- file.path(de_dir, aid_fs)

  if (!dir.exists(adir)) {
    cat("WARN: no existe carpeta DE para", aid, "->", adir, "\n")
    next
  }

  all_files <- list.files(adir, pattern = "DE_.*_all_.*\\.tsv$", full.names = TRUE)
  if (length(all_files) == 0) {
    cat("WARN: no hay DE_*_all_*.tsv en", adir, "\n")
    next
  }
  fp <- all_files[order(file.info(all_files)$mtime, decreasing = TRUE)][1]
  de <- read.delim(fp, check.names = FALSE)

  if (!all(c("feature_id","FDR") %in% colnames(de))) {
    cat("WARN:", aid, "tabla sin feature_id/FDR:", fp, "\n")
    next
  }

  fdr <- as_num_safe(de$FDR)
  feats <- as.character(de$feature_id)[is.finite(fdr) & fdr < fdr_thr]
  feats <- unique(feats)

  membership[[aid]] <- feats
  union_feats <- unique(c(union_feats, feats))

  cat(sprintf("[DE] %-25s  sig(FDR<%.2g)=%d  file=%s\n", aid, fdr_thr, length(feats), basename(fp)))
}

cat("\nUnion de features (FDR<", fdr_thr, "): ", length(union_feats), "\n\n", sep="")

if (length(union_feats) == 0) {
  stop("No hay features con FDR < ", fdr_thr, " en ninguna comparación.")
}

# guardar membership wide (feature x analysis)
mem_df <- data.frame(feature_id = union_feats, stringsAsFactors = FALSE)
for (aid in names(membership)) {
  mem_df[[safe_name(aid)]] <- union_feats %in% membership[[aid]]
}
mem_fp <- file.path(tab_outdir, paste0("features_membership_FDR", gsub("\\.","", sprintf("%.2g", fdr_thr)), "_", ts, ".tsv"))
write.table(mem_df, mem_fp, sep="\t", quote=FALSE, row.names=FALSE)

featlist_fp <- file.path(tab_outdir, paste0("features_union_FDR", gsub("\\.","", sprintf("%.2g", fdr_thr)), "_", ts, ".txt"))
writeLines(union_feats, featlist_fp)

# ---- read counts + meta (build logCPM TMM) ----
if (!file.exists(counts_fp)) stop("No existe counts: ", counts_fp)
if (!file.exists(meta_fp)) stop("No existe meta: ", meta_fp)

df <- if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::fread(counts_fp, data.table = FALSE, check.names = FALSE)
} else read.csv(counts_fp, check.names = FALSE)

meta <- if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::fread(meta_fp, data.table = FALSE, check.names = FALSE)
} else read.csv(meta_fp, check.names = FALSE)

if (!("id" %in% colnames(meta))) stop("El metadata debe tener columna 'id'.")
meta <- meta[!(is.na(meta$id) | trimws(as.character(meta$id)) == ""), , drop = FALSE]

# columnas counts
anno_candidates <- c("id","type","entrez_id","HGNC_symbol")
anno_present <- intersect(anno_candidates, colnames(df))
count_cols <- setdiff(colnames(df), anno_present)

# drop por regex (post mortem)
drop_samples <- grep(drop_re, count_cols, value = TRUE)
count_cols <- setdiff(count_cols, drop_samples)

feature_id <- if ("id" %in% colnames(df)) df[["id"]] else seq_len(nrow(df))
counts_df <- df[, count_cols, drop = FALSE]
counts_df[] <- lapply(counts_df, function(x) suppressWarnings(as.numeric(x)))

counts_mat <- as.matrix(counts_df)
storage.mode(counts_mat) <- "numeric"
rownames(counts_mat) <- make.unique(as.character(feature_id))
colnames(counts_mat) <- make.unique(colnames(counts_mat))

# limpieza mínima
counts_mat <- counts_mat[rowSums(counts_mat, na.rm = TRUE) > 0, , drop = FALSE]
rv <- apply(counts_mat, 1, function(x) stats::var(x, na.rm = TRUE))
counts_mat <- counts_mat[is.finite(rv) & rv > 0, , drop = FALSE]

cat("counts_mat base:", nrow(counts_mat), "features x", ncol(counts_mat), "samples\n")

# filtro libsize global
y0 <- DGEList(counts_mat)
lib0 <- y0$samples$lib.size
keep_lib <- lib0 >= min_lib
counts_mat <- counts_mat[, keep_lib, drop = FALSE]
cat("Tras min_libsize:", ncol(counts_mat), "samples\n\n")

# alinear meta
samples_final <- colnames(counts_mat)
meta_ids <- as.character(meta$id)
missing_in_meta <- setdiff(samples_final, meta_ids)
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
meta2 <- meta[match(samples_final, meta_ids), , drop = FALSE]
stopifnot(all(as.character(meta2$id) == samples_final))

# columnas survival
time_var <- "MESES_SEGUIMIENTO_"
event_var <- "MUERTE"
if (!(time_var %in% colnames(meta2))) stop("Meta sin columna: ", time_var)
if (!(event_var %in% colnames(meta2))) stop("Meta sin columna: ", event_var)

time_m <- as_num_safe(meta2[[time_var]])
event_m <- as_num_safe(meta2[[event_var]])
event_m <- ifelse(is.finite(event_m) & event_m != 0, 1, ifelse(is.finite(event_m), 0, NA))

# logCPM (TMM)
y <- DGEList(counts_mat)
y <- calcNormFactors(y, method = "TMM")
logcpm <- cpm(y, log = TRUE, prior.count = prior_cnt)

# restringir features a las presentes
union_feats2 <- union_feats[union_feats %in% rownames(logcpm)]
cat("Union presentes en logcpm:", length(union_feats2), "/", length(union_feats), "\n\n")
if (length(union_feats2) == 0) stop("Ninguna feature del union está en logcpm (IDs no coinciden).")

# ---- Spearman per feature ----
spearman_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 5) return(list(rho = NA_real_, p = NA_real_, n = length(x)))
  ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = length(x))
}

res <- data.frame(
  feature_id = union_feats2,
  rho_all = NA_real_, p_all = NA_real_, n_all = NA_integer_,
  rho_event = NA_real_, p_event = NA_real_, n_event = NA_integer_,
  stringsAsFactors = FALSE
)

for (k in seq_along(union_feats2)) {
  feat <- union_feats2[k]
  expr <- as.numeric(logcpm[feat, ])
  names(expr) <- colnames(logcpm)

  # all (incluye censura) — exploratorio
  r1 <- spearman_one(expr, time_m)
  res$rho_all[k] <- r1$rho
  res$p_all[k] <- r1$p
  res$n_all[k] <- r1$n

  # solo eventos
  idx_ev <- which(is.finite(event_m) & event_m == 1)
  r2 <- spearman_one(expr[idx_ev], time_m[idx_ev])
  res$rho_event[k] <- r2$rho
  res$p_event[k] <- r2$p
  res$n_event[k] <- r2$n
}

# FDR (BH) dentro del set exploratorio
res$fdr_all <- p.adjust(res$p_all, method = "BH")
res$fdr_event <- p.adjust(res$p_event, method = "BH")

# ordenar por evidencia (all)
res <- res[order(res$fdr_all, res$p_all, na.last = TRUE), , drop = FALSE]

out_fp <- file.path(tab_outdir, paste0("spearman_expr_vs_followup_unionFDR", gsub("\\.","", sprintf("%.2g", fdr_thr)), "_", ts, ".tsv"))
write.table(res, out_fp, sep="\t", quote=FALSE, row.names=FALSE)

cat("Tabla Spearman:", out_fp, "\n")
cat("Membership:", mem_fp, "\n")
cat("Lista union:", featlist_fp, "\n\n")

# ---- Plots (limitados) ----
plot_one <- function(feat, expr, time_m, event_m, out_png, rho_all, p_all, rho_ev, p_ev) {
  ok <- is.finite(expr) & is.finite(time_m) & is.finite(event_m)
  x <- time_m[ok]
  yv <- expr[ok]
  ev <- event_m[ok]

  # estilos: evento = punto lleno; censura = círculo vacío
  pch <- ifelse(ev == 1, 16, 1)

  png(out_png, width = 1600, height = 1200, res = 160)
  par(mar = c(5,5,4,1) + 0.1)
  plot(x, yv, pch = pch, xlab = "MESES_SEGUIMIENTO_", ylab = "logCPM (TMM)",
       main = feat)
  # línea LM (solo visual)
  if (length(x) >= 5) {
    abline(stats::lm(yv ~ x), lwd = 2)
    # loess suave
    ss <- try(stats::loess(yv ~ x), silent = TRUE)
    if (!inherits(ss, "try-error")) {
      xx <- seq(min(x), max(x), length.out = 200)
      yy <- predict(ss, newdata = data.frame(x = xx))
      lines(xx, yy, lwd = 2, lty = 2)
    }
  }
  mtext(
    sprintf("Spearman(all): rho=%.3f p=%.3g | Spearman(events): rho=%.3f p=%.3g",
            rho_all, p_all, rho_ev, p_ev),
    side = 3, line = 0.2, cex = 0.95
  )
  legend("topright",
         legend = c("Evento (MUERTE=1)", "Censurado (MUERTE=0)"),
         pch = c(16,1), bty = "n")
  dev.off()
}

nplot <- min(max_plots, nrow(res))
cat("Generando plots:", nplot, "de", nrow(res), "features (limitado por --max_plots)\n\n")

for (k in seq_len(nplot)) {
  feat <- res$feature_id[k]
  expr <- as.numeric(logcpm[feat, ])
  out_png <- file.path(fig_outdir, paste0(sprintf("%03d", k), "_", safe_name(feat), "_", ts, ".png"))
  plot_one(
    feat = feat,
    expr = expr,
    time_m = time_m,
    event_m = event_m,
    out_png = out_png,
    rho_all = res$rho_all[k],
    p_all = res$p_all[k],
    rho_ev = res$rho_event[k],
    p_ev = res$p_event[k]
  )
}

cat("FIN.\n")
cat("- Log:", log_fp, "\n")
cat("- Tabla Spearman:", out_fp, "\n")
cat("- Figuras:", fig_outdir, "\n")

sink(type="message"); sink(type="output"); close(zz)
