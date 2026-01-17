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
exclude_samples_arg <- get_arg("--exclude_samples", "")

drop_re   <- get_arg("--drop_regex", "^A")
min_lib   <- as.numeric(get_arg("--min_libsize", "100000"))
prior_cnt <- as.numeric(get_arg("--prior_count", "2"))

fdr_thr   <- as.numeric(get_arg("--fdr", "0.1"))
time_var  <- get_arg("--time_var", "MESES_SEGUIMIENTO_")
event_var <- get_arg("--event_var", "MUERTE")

standardize <- tolower(get_arg("--standardize", "false")) %in% c("true","t","1","yes","y")
# standardize=TRUE -> HR por 1 SD de expresión (z-score). Si FALSE -> HR por 1 unidad logCPM.

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

tab_outdir <- file.path("results","tables","survival_cox_univariate")
dir.create(tab_outdir, showWarnings = FALSE, recursive = TRUE)

log_fp <- file.path("logs", paste0("survival_cox_univariate_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type = "output"); sink(zz, type = "message")

cat("============================================\n")
cat("Cox univariado por feature (union DE FDR<thr)\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("Counts:", counts_fp, "\n")
cat("Meta:", meta_fp, "\n")
cat("Spec:", spec_fp, "\n")
cat("DE dir:", de_dir, "\n")
cat("exclude_samples:", ifelse(nzchar(exclude_samples_arg), exclude_samples_arg, "NA"), "\n")
cat("drop_regex:", drop_re, "\n")
cat("min_libsize:", min_lib, "\n")
cat("prior_count:", prior_cnt, "\n")
cat("Union FDR<thr:", fdr_thr, "\n")
cat("time_var:", time_var, "\n")
cat("event_var:", event_var, "\n")
cat("standardize:", standardize, "\n")
cat("Output tables:", tab_outdir, "\n")
cat("Log:", log_fp, "\n")
cat("============================================\n\n")

suppressPackageStartupMessages({
  library(edgeR)
  library(survival)
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

exclude_samples_str <- exclude_samples_arg
if (!nzchar(exclude_samples_str) && "exclude_samples" %in% colnames(spec)) {
  vals <- trimws(as.character(spec$exclude_samples))
  vals <- vals[nzchar(vals)]
  if (length(vals) > 0) {
    exclude_samples_str <- paste(unique(vals), collapse = ",")
  }
}

# ---- read DE tables and union of features ----
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

# guardar membership wide
mem_df <- data.frame(feature_id = union_feats, stringsAsFactors = FALSE)
for (aid in names(membership)) {
  mem_df[[safe_name(aid)]] <- union_feats %in% membership[[aid]]
}
mem_fp <- file.path(tab_outdir, paste0("features_membership_unionFDR", gsub("\\.","", sprintf("%.2g", fdr_thr)), "_", ts, ".tsv"))
write.table(mem_df, mem_fp, sep="\t", quote=FALSE, row.names=FALSE)

featlist_fp <- file.path(tab_outdir, paste0("features_unionFDR", gsub("\\.","", sprintf("%.2g", fdr_thr)), "_", ts, ".txt"))
writeLines(union_feats, featlist_fp)

# ---- read counts + meta; build logCPM (TMM) ----
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

if (!(time_var %in% colnames(meta))) stop("Meta sin columna time_var: ", time_var)
if (!(event_var %in% colnames(meta))) stop("Meta sin columna event_var: ", event_var)

# columns counts
anno_candidates <- c("id","type","entrez_id","HGNC_symbol")
anno_present <- intersect(anno_candidates, colnames(df))
count_cols <- setdiff(colnames(df), anno_present)

# drop regex
drop_samples <- grep(drop_re, count_cols, value = TRUE)
count_cols <- setdiff(count_cols, drop_samples)

feature_id <- if ("id" %in% colnames(df)) df[["id"]] else seq_len(nrow(df))
counts_df <- df[, count_cols, drop = FALSE]
counts_df[] <- lapply(counts_df, function(x) suppressWarnings(as.numeric(x)))

counts_mat <- as.matrix(counts_df)
storage.mode(counts_mat) <- "numeric"
rownames(counts_mat) <- make.unique(as.character(feature_id))
colnames(counts_mat) <- make.unique(colnames(counts_mat))

# minimal feature cleaning
counts_mat <- counts_mat[rowSums(counts_mat, na.rm = TRUE) > 0, , drop = FALSE]
rv <- apply(counts_mat, 1, function(x) stats::var(x, na.rm = TRUE))
counts_mat <- counts_mat[is.finite(rv) & rv > 0, , drop = FALSE]
cat("counts_mat base:", nrow(counts_mat), "features x", ncol(counts_mat), "samples\n")

# libsize global filter
y0 <- DGEList(counts_mat)
lib0 <- y0$samples$lib.size
keep_lib <- lib0 >= min_lib
counts_mat <- counts_mat[, keep_lib, drop = FALSE]
cat("Tras min_libsize:", ncol(counts_mat), "samples\n\n")

# align meta
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

# excluir muestras específicas (si aplica)
if (nzchar(exclude_samples_str)) {
  exclude_ids <- trimws(unlist(strsplit(exclude_samples_str, "[,;| ]+")))
  exclude_ids <- exclude_ids[nzchar(exclude_ids)]
  if (length(exclude_ids) > 0) {
    excl_mask <- meta2$id %in% exclude_ids
    if (any(excl_mask)) {
      cat("Excluyendo muestras:", paste(meta2$id[excl_mask], collapse = ", "), "\n")
      keep <- !excl_mask
      meta2 <- meta2[keep, , drop = FALSE]
      counts_mat <- counts_mat[, keep, drop = FALSE]
    } else {
      cat("WARNING: exclude_samples no coincide con meta$id\n")
    }
  }
}

time_m <- as_num_safe(meta2[[time_var]])
event_m <- as_num_safe(meta2[[event_var]])
event_m <- ifelse(is.finite(event_m) & event_m != 0, 1, ifelse(is.finite(event_m), 0, NA))

# keep samples with valid survival info
keep_surv <- is.finite(time_m) & is.finite(event_m)
cat("Muestras con survival válido:", sum(keep_surv), "/", length(keep_surv), "\n\n")

if (sum(keep_surv) < 8) stop("Muy pocas muestras con survival válido para Cox.")

counts_mat2 <- counts_mat[, keep_surv, drop = FALSE]
meta3 <- meta2[keep_surv, , drop = FALSE]
time_m2 <- time_m[keep_surv]
event_m2 <- event_m[keep_surv]

# logCPM post-TMM
y <- DGEList(counts_mat2)
y <- calcNormFactors(y, method = "TMM")
logcpm <- cpm(y, log = TRUE, prior.count = prior_cnt)

# restrict features present
union_feats2 <- union_feats[union_feats %in% rownames(logcpm)]
cat("Union presentes en logcpm:", length(union_feats2), "/", length(union_feats), "\n\n")
if (length(union_feats2) == 0) stop("Ninguna feature del union está en logcpm (IDs no coinciden).")

# ---- Cox per feature ----
cox_one <- function(expr, time, event) {
  # expr: numeric vector
  ok <- is.finite(expr) & is.finite(time) & is.finite(event)
  expr <- expr[ok]; time <- time[ok]; event <- event[ok]
  if (length(expr) < 8) return(NULL)
  if (stats::sd(expr) == 0) return(NULL)

  dfc <- data.frame(time = time, event = event, expr = expr)
  fit <- survival::coxph(survival::Surv(time, event) ~ expr, data = dfc, ties = "efron")
  s <- summary(fit)

  beta <- s$coef[1, "coef"]
  se   <- s$coef[1, "se(coef)"]
  z    <- s$coef[1, "z"]
  p    <- s$coef[1, "Pr(>|z|)"]

  hr <- exp(beta)
  ci <- suppressWarnings(exp(confint(fit)))
  hr_lo <- ci[1,1]
  hr_hi <- ci[1,2]

  list(
    n = nrow(dfc),
    beta = beta, se = se, z = z, p = p,
    hr = hr, hr_lo = hr_lo, hr_hi = hr_hi
  )
}

res_list <- vector("list", length(union_feats2))

for (k in seq_along(union_feats2)) {
  feat <- union_feats2[k]
  expr <- as.numeric(logcpm[feat, ])
  if (standardize) {
    expr <- as.numeric(scale(expr))
  }
  out <- cox_one(expr, time_m2, event_m2)
  if (is.null(out)) {
    res_list[[k]] <- NULL
  } else {
    res_list[[k]] <- data.frame(
      feature_id = feat,
      n = out$n,
      beta = out$beta,
      se = out$se,
      z = out$z,
      p = out$p,
      HR = out$hr,
      HR_CI_low = out$hr_lo,
      HR_CI_high = out$hr_hi,
      stringsAsFactors = FALSE
    )
  }
}

res_df <- do.call(rbind, res_list)
if (is.null(res_df) || nrow(res_df) == 0) stop("No se pudo ajustar Cox para ninguna feature (revisar datos).")

res_df$FDR <- p.adjust(res_df$p, method = "BH")
res_df <- res_df[order(res_df$FDR, res_df$p, na.last = TRUE), , drop = FALSE]

# attach membership columns (si existen)
res_df2 <- merge(res_df, mem_df, by = "feature_id", all.x = TRUE)

# outputs
tag_std <- if (standardize) "HRper1SD" else "HRper1logCPM"
out_fp <- file.path(tab_outdir, paste0("cox_univariate_unionFDR", gsub("\\.","", sprintf("%.2g", fdr_thr)),
                                       "_", tag_std, "_", ts, ".tsv"))
write.table(res_df2, out_fp, sep="\t", quote=FALSE, row.names=FALSE)

cat("OK. Cox results:", out_fp, "\n")
cat("Membership:", mem_fp, "\n")
cat("Union list:", featlist_fp, "\n")
cat("\nResumen:\n")
cat("- n_features_tested:", nrow(res_df), "\n")
cat("- min FDR:", min(res_df$FDR, na.rm = TRUE), "\n")
cat("- n_FDR<0.05:", sum(res_df$FDR < 0.05, na.rm = TRUE), "\n")
cat("- n_FDR<0.1:", sum(res_df$FDR < 0.1, na.rm = TRUE), "\n")

cat("\nFIN. Log:", log_fp, "\n")

sink(type="message"); sink(type="output"); close(zz)
