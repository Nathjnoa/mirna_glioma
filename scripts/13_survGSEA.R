#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

# ============================================================
# Survival-ranked GSEA (miEAA) for smallRNA-seq miRNAs
# - Build logCPM(TMM) from counts
# - Cox univariate per miRNA: Surv(time, event) ~ expression
# - Rank by signed Wald z
# - miEAA GSEA using 4 miRPathDB expert categories
# ============================================================

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

# ---- CLI defaults (adapt to your project structure) ----
counts_fp <- get_arg("--counts", file.path("data","intermediate","Gliomas_all_counts_merged.csv"))
meta_fp   <- get_arg("--meta",   file.path("data","intermediate","Metadatos_gliomas_verificados.csv"))

out_root  <- get_arg("--out_root", "results")
run_tag   <- get_arg("--run_tag", "SurvivalRank_CoxZ_miRPathDB")

drop_regex   <- get_arg("--drop_regex", "^A")
min_libsize  <- as_num_safe(get_arg("--min_libsize", "100000"))
prior_count  <- as_num_safe(get_arg("--prior_count", "2"))

time_var  <- get_arg("--time_var",  "MESES_SEGUIMIENTO_")
event_var <- get_arg("--event_var", "MUERTE")  # expected 0/1
event_dead_value <- as_num_safe(get_arg("--event_dead_value","1"))

exclude_samples_arg <- get_arg("--exclude_samples", "") # e.g. "T11,T05"
exclude_samples <- if (nzchar(exclude_samples_arg)) trimws(strsplit(exclude_samples_arg, ",")[[1]]) else character()

scale_expr <- tolower(get_arg("--scale_expr", "true")) %in% c("true","t","1","yes","y")

# miEAA params
species <- get_arg("--species", "hsa")
mirna_type <- get_arg("--mirna_type", "mature")
p_adj_method <- get_arg("--p_adj_method", "fdr")
independent_p_adj <- tolower(get_arg("--independent_p_adj", "false")) %in% c("true","t","1","yes","y")
sig_level <- as_num_safe(get_arg("--sig_level", "1"))  # keep all by default
min_hits  <- as.numeric(get_arg("--min_hits", "10"))

# miRPathDB expert categories (mature)
categories <- c(
  "miRPathDB_GO_Biological_process_mature",
  "miRPathDB_GO_Molecular_function_mature",
  "miRPathDB_KEGG_mature",
  "miRPathDB_Reactome_mature"
)

# ---- robust paths (relative to project root) ----
get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
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
counts_fp <- resolve_path(counts_fp)
meta_fp   <- resolve_path(meta_fp)
out_root  <- resolve_path(out_root)

# ---- logging ----
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

out_tables <- file.path(out_root, "tables", "SurvivalRank_GSEA", run_tag)
dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)

log_fp <- file.path("logs", paste0("survival_rank_gsea_", run_tag, "_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type="output"); sink(zz, type="message")
on.exit({ sink(type="message"); sink(type="output"); close(zz) }, add=TRUE)

cat("============================================\n")
cat("Survival-ranked GSEA (miEAA) - miRPathDB expert\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("counts:", counts_fp, "\n")
cat("meta:", meta_fp, "\n")
cat("drop_regex:", drop_regex, "\n")
cat("min_libsize:", min_libsize, "\n")
cat("prior_count:", prior_count, "\n")
cat("exclude_samples:", if (length(exclude_samples)) paste(exclude_samples, collapse=",") else "(none)", "\n")
cat("time_var:", time_var, " | event_var:", event_var, " | event_dead_value:", event_dead_value, "\n")
cat("scale_expr:", scale_expr, "\n")
cat("miEAA species:", species, " | mirna_type:", mirna_type, "\n")
cat("p_adj_method:", p_adj_method, " | independent_p_adj:", independent_p_adj, "\n")
cat("sig_level:", sig_level, " | min_hits:", min_hits, "\n")
cat("categories:", paste(categories, collapse=", "), "\n")
cat("out_tables:", out_tables, "\n")
cat("log:", log_fp, "\n")
cat("============================================\n\n")

suppressPackageStartupMessages({
  library(edgeR)
  library(survival)
  library(rbioapi)
})

# ---- validate miEAA categories exist ----
cats_available <- tryCatch(rbioapi::rba_mieaa_cats(), error = function(e) NULL)
if (!is.null(cats_available) && is.data.frame(cats_available) && nrow(cats_available) > 0) {
  # heuristics: first column often contains category IDs
  cat_col <- intersect(c("category", "categories", "id", "ID", "name"), names(cats_available))
  if (length(cat_col) == 0) cat_col <- names(cats_available)[1]
  avail_ids <- as.character(cats_available[[cat_col[1]]])
  missing <- setdiff(categories, avail_ids)
  if (length(missing) > 0) {
    stop("miEAA categories not found: ", paste(missing, collapse=", "),
         "\nTip: inspect rba_mieaa_cats() output and update category IDs accordingly.")
  }
} else {
  cat("WARN: Could not validate categories via rba_mieaa_cats() (continuing).\n")
}

# ---- read data ----
if (!file.exists(counts_fp)) stop("Counts file not found: ", counts_fp)
if (!file.exists(meta_fp)) stop("Metadata file not found: ", meta_fp)

counts_df <- read.csv(counts_fp, check.names = FALSE)
meta <- read.csv(meta_fp, check.names = FALSE)

if (!("id" %in% names(counts_df))) stop("Counts table missing column 'id'")
if (!("type" %in% names(counts_df))) stop("Counts table missing column 'type'")

# annotation columns (common)
anno_cols <- intersect(c("id","type","entrez_id","HGNC_symbol","feature_id"), names(counts_df))
non_anno <- setdiff(names(counts_df), anno_cols)

# sample columns: usually those starting with A/T, or just all non-anno
sample_cols <- non_anno
if (any(grepl("^[AT]", sample_cols))) {
  sample_cols <- sample_cols[grepl("^[AT]", sample_cols)]
}

# drop controls (^A)
sample_cols <- sample_cols[!grepl(drop_regex, sample_cols)]
# exclude specific samples if requested
if (length(exclude_samples)) sample_cols <- setdiff(sample_cols, exclude_samples)

if (length(sample_cols) < 5) stop("Too few sample columns inferred from counts.")

# detect sample ID column in metadata by maximum overlap with sample_cols
meta_char <- meta
for (cc in names(meta_char)) meta_char[[cc]] <- as.character(meta_char[[cc]])
overlaps <- sapply(names(meta_char), function(cc) sum(meta_char[[cc]] %in% sample_cols, na.rm = TRUE))
id_col <- names(overlaps)[which.max(overlaps)]
if (max(overlaps, na.rm = TRUE) < 5) {
  stop("Could not detect sample ID column in metadata by overlap (need >=5 matches). ",
       "Please add --meta_id_col support or ensure metadata contains sample IDs.")
}
meta$.__sid__. <- as.character(meta[[id_col]])
cat("Detected metadata sample ID column:", id_col, " (matches=", max(overlaps), ")\n", sep="")

# align sample order
keep_samples <- intersect(sample_cols, meta$.__sid__.)
if (length(keep_samples) < 5) stop("After aligning metadata, too few samples remain.")
meta2 <- meta[match(keep_samples, meta$.__sid__.), , drop = FALSE]
stopifnot(all(meta2$.__sid__. == keep_samples))

# build counts matrix
counts_mat <- as.matrix(counts_df[, keep_samples, drop = FALSE])
storage.mode(counts_mat) <- "numeric"
if (any(!is.finite(counts_mat))) stop("Counts matrix has non-finite values.")
# counts should be integer-ish; warn if not
nonint <- sum(abs(counts_mat - round(counts_mat)) > 1e-6)
if (nonint > 0) cat("WARN: counts have non-integer entries (n=", nonint, "); proceeding with numeric.\n", sep="")
counts_mat <- round(counts_mat)
rownames(counts_mat) <- as.character(counts_df$id)

# ---- edgeR normalization (TMM) + logCPM ----
y <- DGEList(counts = counts_mat)

libsize <- y$samples$lib.size
keep_ls <- libsize >= min_libsize
cat("Samples before libsize filter:", ncol(y), "\n")
cat("Samples removed by libsize <", min_libsize, ":", sum(!keep_ls), "\n")

y <- y[, keep_ls, keep.lib.sizes = FALSE]
meta3 <- meta2[keep_ls, , drop = FALSE]

y <- calcNormFactors(y, method = "TMM")
logcpm <- cpm(y, log = TRUE, prior.count = prior_count)

# ---- keep miRNAs only (type == mirna) ----
type_vec <- tolower(as.character(counts_df$type))
names(type_vec) <- as.character(counts_df$id)
mir_ids <- intersect(rownames(logcpm), names(type_vec)[type_vec == "mirna"])
if (length(mir_ids) < 20) stop("Too few miRNAs after filters (n<20).")

X <- t(logcpm[mir_ids, , drop = FALSE])  # samples x miRNAs

# ---- survival variables ----
time <- as_num_safe(meta3[[time_var]])
event_raw <- meta3[[event_var]]
event <- as_num_safe(event_raw)
event <- ifelse(is.na(event), NA_real_, ifelse(event == event_dead_value, 1, 0))

ok_surv <- is.finite(time) & is.finite(event)
cat("Samples with valid time+event:", sum(ok_surv), "/", length(ok_surv), "\n")
if (sum(ok_surv) < 10) stop("Too few samples with valid time/event (n<10).")

X <- X[ok_surv, , drop = FALSE]
meta4 <- meta3[ok_surv, , drop = FALSE]
time <- time[ok_surv]
event <- event[ok_surv]

n_events <- sum(event == 1, na.rm = TRUE)
cat("Final N for Cox:", nrow(X), " | Events:", n_events, "\n")
if (n_events < 5) stop("Too few events for Cox (<5).")

# ---- Cox univariate per miRNA ----
cox_one <- function(x, time, event, scale_x = TRUE) {
  x <- as_num_safe(x)
  if (!all(is.finite(x))) return(NULL)
  if (sd(x) == 0) return(NULL)
  if (scale_x) x <- as.numeric(scale(x))

  df <- data.frame(time=time, event=event, x=x)
  fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ x, data=df, ties="efron"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  sm <- tryCatch(summary(fit), error = function(e) NULL)
  if (is.null(sm) || is.null(sm$coef)) return(NULL)

  beta <- sm$coef["x","coef"]
  se   <- sm$coef["x","se(coef)"]
  z    <- sm$coef["x","z"]
  p    <- sm$coef["x","Pr(>|z|)"]

  data.frame(beta=beta, se=se, z=z, p=p, HR=exp(beta), stringsAsFactors = FALSE)
}

res_list <- vector("list", ncol(X))
names(res_list) <- colnames(X)
warn_msgs <- character()

for (j in seq_len(ncol(X))) {
  mi <- colnames(X)[j]
  out <- withCallingHandlers(
    cox_one(X[, j], time, event, scale_x = scale_expr),
    warning = function(w) {
      warn_msgs <<- c(warn_msgs, conditionMessage(w))
    }
  )
  if (is.null(out)) next
  out$feature_id <- mi
  res_list[[j]] <- out
  if (j %% 500 == 0) cat("Cox progress:", j, "/", ncol(X), "\n")
}

cox_df <- do.call(rbind, res_list)
if (is.null(cox_df) || nrow(cox_df) == 0) stop("Cox produced no results.")

if (length(warn_msgs) > 0) {
  warn_fp <- file.path(out_tables, paste0("Cox_warnings_", ts, ".txt"))
  writeLines(warn_msgs, warn_fp)
  cat("Cox warnings (first 20):\n")
  cat(paste(head(warn_msgs, 20), collapse = "\n"), "\n")
  cat("Cox warnings file:", warn_fp, "\n")
}

cox_df$FDR <- p.adjust(cox_df$p, method = "BH")
# ranking metric: signed Wald z
cox_df$rank_stat <- cox_df$z
# interpretation helper
cox_df$direction <- ifelse(cox_df$rank_stat > 0, "higher_expr_worse_survival", "higher_expr_better_survival")

# export Cox table
cox_df <- cox_df[order(cox_df$p), ]
cox_all_fp <- file.path(out_tables, paste0("Cox_univariate_miRNA_all_", ts, ".tsv"))
write.table(cox_df, cox_all_fp, sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote Cox table:", cox_all_fp, " (n=", nrow(cox_df), ")\n", sep="")

# ranked list for miEAA (IDs only, sorted by rank_stat desc)
ranked <- cox_df[is.finite(cox_df$rank_stat) & nzchar(cox_df$feature_id), ]
ranked <- ranked[order(ranked$rank_stat, decreasing = TRUE), ]
ranked <- ranked[!duplicated(ranked$feature_id), ]
ranked_ids <- ranked$feature_id

rank_fp <- file.path(out_tables, paste0("Ranked_miRNAs_by_CoxZ_", ts, ".txt"))
writeLines(ranked_ids, rank_fp)
cat("Wrote ranked miRNAs:", rank_fp, " (n=", length(ranked_ids), ")\n", sep="")

# ---- miEAA GSEA ----
cat("\n--- miEAA GSEA (miRPathDB expert) ---\n")
res <- tryCatch(
  rbioapi::rba_mieaa_enrich(
    test_set = ranked_ids,
    mirna_type = mirna_type,
    test_type = "GSEA",
    species = species,
    categories = categories,
    p_adj_method = p_adj_method,
    independent_p_adj = independent_p_adj,
    sig_level = sig_level,
    min_hits = min_hits
  ),
  error = function(e) {
    cat("ERROR miEAA:", conditionMessage(e), "\n")
    NULL
  }
)

if (is.null(res) || !is.data.frame(res) || nrow(res) == 0) {
  stop("miEAA returned empty results.")
}

# coerce numeric columns
for (cc in c("P-value","P-adjusted","Q-value","Observed")) {
  if (cc %in% names(res)) res[[cc]] <- as_num_safe(res[[cc]])
}

# write full + top50
mieaa_all_fp <- file.path(out_tables, paste0("miEAA_GSEA_all_miRPathDB_expert_", ts, ".tsv"))
write.table(res, mieaa_all_fp, sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote miEAA all:", mieaa_all_fp, "\n")

qcol <- if ("Q-value" %in% names(res)) "Q-value" else if ("P-adjusted" %in% names(res)) "P-adjusted" else "P-value"
ord <- order(res[[qcol]], na.last = TRUE)
top <- res[ord, , drop=FALSE]
if (nrow(top) > 50) top <- top[1:50, , drop=FALSE]

mieaa_top_fp <- file.path(out_tables, paste0("miEAA_GSEA_top50_miRPathDB_expert_", ts, ".tsv"))
write.table(top, mieaa_top_fp, sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote miEAA top50:", mieaa_top_fp, "\n")

cat("\nDONE\n")
cat("Tables:", out_tables, "\n")
cat("Log:", log_fp, "\n")
