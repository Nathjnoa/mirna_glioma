#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

# =========================
# CLI parsing (simple)
# =========================
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

as_num_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

read_table <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(path, data.table = FALSE, check.names = FALSE)
  } else {
    read.csv(path, check.names = FALSE)
  }
}

parse_collapse_levels <- function(s) {
  # "0|1=01;2=2" -> map "0"->"01", "1"->"01", "2"->"2"
  res <- list()
  if (is.null(s) || !nzchar(s)) return(res)
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

parse_range_map <- function(map_str) {
  out <- list()
  parts <- strsplit(map_str, ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  for (p in parts) {
    kv <- strsplit(p, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) next
    rng <- trimws(kv[1]); lab <- trimws(kv[2])
    mm <- strsplit(rng, "-", fixed = TRUE)[[1]]
    if (length(mm) != 2) next
    mn <- as_num_safe(mm[1]); mx <- as_num_safe(mm[2])
    if (!is.finite(mn) || !is.finite(mx)) next
    out[[length(out) + 1]] <- list(min = mn, max = mx, label = lab)
  }
  out
}

# =========================
# Args
# =========================
counts_fp <- get_arg("--counts", file.path("data","intermediate","Gliomas_all_counts_merged.csv"))
meta_fp   <- get_arg("--meta",   file.path("data","intermediate","Metadatos_gliomas_verificados.csv"))
spec_fp   <- get_arg("--spec",   file.path("config","de_specs.csv"))

outdir    <- get_arg("--outdir", file.path("results","DE_miRNA_protCoding"))
drop_re   <- get_arg("--drop_regex", "^A")
min_lib   <- as.numeric(get_arg("--min_libsize", "100000"))
prior_cnt <- as.numeric(get_arg("--prior_count", "2"))

keep_types_str <- get_arg("--keep_types", "mirna,protein_coding")  # editable
keep_types <- trimws(strsplit(keep_types_str, ",")[[1]])
keep_types <- keep_types[nzchar(keep_types)]

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("data","processed"), showWarnings = FALSE, recursive = TRUE)

log_fp <- file.path("logs", paste0("DE_miRNA_protCoding_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type = "output"); sink(zz, type = "message")

cat("============================================\n")
cat("edgeR DE (miRNA + protein_coding) desde spec\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("Counts:", counts_fp, "\n")
cat("Meta:", meta_fp, "\n")
cat("Spec:", spec_fp, "\n")
cat("keep_types:", paste(keep_types, collapse = ", "), "\n")
cat("drop_regex:", drop_re, "\n")
cat("min_libsize:", min_lib, "\n")
cat("prior_count:", prior_cnt, "\n")
cat("outdir:", outdir, "\n")
cat("Log:", log_fp, "\n")
cat("============================================\n\n")

suppressPackageStartupMessages({
  library(edgeR)
})

# =========================
# Read inputs
# =========================
if (!file.exists(counts_fp)) stop("No existe counts: ", counts_fp)
if (!file.exists(meta_fp)) stop("No existe meta: ", meta_fp)
if (!file.exists(spec_fp)) stop("No existe spec: ", spec_fp)

df <- read_table(counts_fp)
meta <- read_table(meta_fp)
spec <- read_table(spec_fp)

if (!("id" %in% colnames(meta))) stop("El metadata debe tener columna 'id'.")
meta <- meta[!(is.na(meta$id) | trimws(as.character(meta$id)) == ""), , drop = FALSE]

required_cols <- c("analysis_id","mode","var","case","ref","cut","cut_type","event_var","collapse_map","exclude_rule")
miss <- setdiff(required_cols, colnames(spec))
if (length(miss) > 0) stop("Spec falta columnas: ", paste(miss, collapse=", "))
if (anyDuplicated(spec$analysis_id)) stop("analysis_id duplicados en spec.")

cat("Dim df:", nrow(df), "x", ncol(df), "\n")
cat("Dim meta:", nrow(meta), "x", ncol(meta), "\n")
cat("Spec analyses:", nrow(spec), "\n\n")

# =========================
# Build counts_mat base (samples)
# =========================
anno_candidates <- c("id","type","entrez_id","HGNC_symbol")
anno_present <- intersect(anno_candidates, colnames(df))
count_cols <- setdiff(colnames(df), anno_present)

# drop post-mortem by regex
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

cat("Tras drop_regex + limpieza mínima:", nrow(counts_mat), "features x", ncol(counts_mat), "samples\n")

# libsize filter
y0 <- DGEList(counts_mat)
lib0 <- y0$samples$lib.size
keep_lib <- lib0 >= min_lib
cat("Filtro global lib_size >= ", min_lib, ":\n", sep="")
cat("- samples antes:", length(lib0), "\n")
cat("- removidas:", sum(!keep_lib), "\n")
if (sum(!keep_lib) > 0) print(names(lib0)[!keep_lib])
cat("- samples después:", sum(keep_lib), "\n\n")
counts_mat <- counts_mat[, keep_lib, drop = FALSE]

# align meta
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

# =========================
# Filter by type (miRNA + protein_coding)
# =========================
if (!("type" %in% colnames(df))) stop("El counts CSV debe tener columna 'type' para filtrar biotipos.")
type_vec <- as.character(df$type)
type_vec[is.na(type_vec)] <- ""

# map types to rownames
# OJO: rownames(counts_mat) viene de df$id; asumimos misma fila/orden original.
# Para robustez, reconstruimos un index por la posición original.
df_ids <- make.unique(as.character(feature_id))
# counts_mat rownames también es make.unique(feature_id); deben coincidir en orden.
if (length(df_ids) != nrow(df)) stop("Inconsistencia interna: df_ids vs nrow(df).")

# Construir un data.frame auxiliar para mapear tipo por fila actual de counts_mat
aux <- data.frame(feature_id = df_ids, type = type_vec, stringsAsFactors = FALSE)
aux <- aux[match(rownames(counts_mat), aux$feature_id), , drop = FALSE]
stopifnot(all(aux$feature_id == rownames(counts_mat)))

keep_type <- aux$type %in% keep_types
cat("Filtro por type:", paste(keep_types, collapse=", "), "\n")
cat("- features antes:", nrow(counts_mat), "\n")
cat("- keep_type TRUE:", sum(keep_type), "\n")

if (sum(keep_type) == 0) {
  cat("WARNING: 0 features con type en keep_types. Tipos disponibles (top 20):\n")
  print(head(sort(table(aux$type), decreasing = TRUE), 20))
  stop("No quedan features tras filtrar por type.")
}

# Warn si protein_coding no existe
if (!("protein_coding" %in% unique(aux$type)) && ("protein_coding" %in% keep_types)) {
  cat("WARNING: 'protein_coding' no aparece en columna type. Se continuará con los tipos presentes.\n")
}

counts_mat <- counts_mat[keep_type, , drop = FALSE]
cat("- features después:", nrow(counts_mat), "\n\n")

# Export base (para reproducibilidad de este análisis)
base_counts_fp <- file.path("data","processed", paste0("counts_filtered_base_miRNA_protCoding_", ts, ".csv"))
base_meta_fp   <- file.path("data","processed", paste0("metadata_filtered_base_miRNA_protCoding_", ts, ".csv"))
write.csv(data.frame(feature_id = rownames(counts_mat), counts_mat, check.names = FALSE),
          base_counts_fp, row.names = FALSE)
write.csv(meta_base, base_meta_fp, row.names = FALSE)

cat("Export base (miRNA/protCoding):\n")
cat("- ", base_counts_fp, "\n", sep="")
cat("- ", base_meta_fp, "\n\n", sep="")

# annotation merge (si existe)
anno_df <- NULL
if ("id" %in% colnames(df)) {
  keep_anno <- intersect(c("id","type","entrez_id","HGNC_symbol"), colnames(df))
  anno_df <- unique(df[, keep_anno, drop = FALSE])
}

force_numeric_df <- function(d, cols) {
  for (cc in cols) if (cc %in% colnames(d)) d[[cc]] <- as_num_safe(d[[cc]])
  d
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

plot_volcano <- function(tt, analysis_fs, tag, fig_dir, ts) {
  if (!("logFC" %in% colnames(tt))) return()
  has_fdr <- "FDR" %in% colnames(tt)
  has_p <- "PValue" %in% colnames(tt)
  if (!has_fdr && !has_p) return()
  yvals <- if (has_fdr) tt$FDR else tt$PValue
  yvals <- pmax(yvals, .Machine$double.xmin)
  neglog <- -log10(yvals)
  sig <- if (has_fdr) tt$FDR < 0.05 else rep(FALSE, nrow(tt))
  png(file.path(fig_dir, paste0("Volcano_", analysis_fs, "_", tag, "_", ts, ".png")),
      width = 1400, height = 1000, res = 150)
  plot(tt$logFC, neglog,
       xlab = "logFC",
       ylab = if (has_fdr) "-log10(FDR)" else "-log10(PValue)",
       main = paste0("Volcano: ", analysis_fs, " ", tag),
       pch = ifelse(sig, 16, 1))
  abline(h = -log10(0.10), lty = 2)
  abline(v = c(-1, 1), lty = 3)
  dev.off()
}

# =========================
# Run one analysis from spec
# =========================
run_one <- function(row) {
  norm_string <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x)] <- ""
    x[x == "NA"] <- ""
    x
  }
  analysis_id <- norm_string(row$analysis_id)
  mode <- tolower(norm_string(row$mode))
  var_key <- norm_string(row$var)
  case_label <- norm_string(row$case)
  ref_label <- norm_string(row$ref)
  cut_val <- as_num_safe(row$cut)[1]
  cut_type <- tolower(norm_string(row$cut_type))
  event_var <- norm_string(row$event_var)
  collapse_map_str <- norm_string(row$collapse_map)
  exclude_rule <- norm_string(row$exclude_rule)

  if (!(var_key %in% colnames(meta_base))) {
    cat("SKIP:", analysis_id, "-> var no existe:", var_key, "\n\n"); return(NULL)
  }

  analysis_fs <- safe_name(analysis_id)
  analysis_dir <- file.path(outdir, analysis_fs)
  dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)

  keep_mask <- rep(TRUE, nrow(meta_base))
  group <- NULL

  if (mode == "continuous") {
    x <- as_num_safe(meta_base[[var_key]])
    keep_mask <- is.finite(x)
    if (sum(keep_mask) < 8) {
      cat("SKIP:", analysis_id, "-> continuous con pocas muestras.\n\n"); return(NULL)
    }
    scale_val <- if (is.finite(cut_val) && cut_val != 0) cut_val else 1
    x_scaled <- x[keep_mask] / scale_val
    counts_sub <- counts_mat[, keep_mask, drop = FALSE]

    y <- DGEList(counts_sub)
    y <- calcNormFactors(y, method = "TMM")
    design <- model.matrix(~ x_scaled)

    keep_feat <- filterByExpr(y, design = design)
    if (!any(keep_feat)) { cat("SKIP:", analysis_id, "-> filterByExpr dejó 0.\n\n"); return(NULL) }
    y <- y[keep_feat, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y, method = "TMM")
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = 2)

    tt <- topTags(qlf, n = Inf)$table
    tt$feature_id <- rownames(tt); rownames(tt) <- NULL
    if (!is.null(anno_df)) tt <- merge(tt, anno_df, by.x="feature_id", by.y="id", all.x=TRUE)
    tt <- force_numeric_df(tt, c("logFC","logCPM","F","PValue","FDR"))

    tag <- paste0("continuous_per", scale_val)

    all_fp <- file.path(analysis_dir, paste0("DE_", analysis_fs, "_", tag, "_all_", ts, ".tsv"))
    sig05_fp <- file.path(analysis_dir, paste0("DE_", analysis_fs, "_", tag, "_sig_FDR0.05_", ts, ".tsv"))
    sig10_fp <- file.path(analysis_dir, paste0("DE_", analysis_fs, "_", tag, "_sig_FDR0.1_", ts, ".tsv"))
    write.table(tt, all_fp, sep="\t", quote=FALSE, row.names=FALSE)
    write.table(tt[!is.na(tt$FDR) & tt$FDR < 0.05, , drop=FALSE], sig05_fp, sep="\t", quote=FALSE, row.names=FALSE)
    write.table(tt[!is.na(tt$FDR) & tt$FDR < 0.10, , drop=FALSE], sig10_fp, sep="\t", quote=FALSE, row.names=FALSE)

    png(file.path(analysis_dir, paste0("MDplot_", analysis_fs, "_", tag, "_", ts, ".png")),
        width=1400, height=1000, res=150)
    plotMD(qlf, status = decide_status_safe(qlf), main = paste0("MD: ", analysis_id, " ", tag))
    abline(h=c(-1,1), lty=2); dev.off()

    plot_volcano(tt, analysis_fs, tag, analysis_dir, ts)

    qc <- data.frame(
      analysis_id = analysis_id,
      mode = mode,
      var = var_key,
      n_samples = ncol(y),
      n_features_tested = nrow(y),
      groups = NA_character_,
      comparison = tag,
      n_sig_FDR0.05 = sum(!is.na(tt$FDR) & tt$FDR < 0.05),
      n_sig_FDR0.1  = sum(!is.na(tt$FDR) & tt$FDR < 0.10),
      stringsAsFactors = FALSE
    )
    return(qc)
  }

  # ---- non-continuous: build group ----
  v <- meta_base[[var_key]]

  if (mode == "as_is") {
    g <- as.character(v); g[trimws(g)==""] <- NA
    keep_mask <- !is.na(g) & (g %in% c(case_label, ref_label))
    group <- factor(g[keep_mask])

  } else if (mode == "binary_cut") {
    x <- as_num_safe(v)
    keep_mask <- is.finite(x)
    if (!is.finite(cut_val)) { cat("SKIP:", analysis_id, "-> cut inválido.\n\n"); return(NULL) }
    if (!cut_type %in% c("gt","ge")) { cat("SKIP:", analysis_id, "-> cut_type inválido.\n\n"); return(NULL) }
    g <- if (cut_type=="gt") ifelse(x > cut_val, case_label, ref_label) else ifelse(x >= cut_val, case_label, ref_label)
    keep_mask <- keep_mask & !is.na(g)
    group <- factor(g[keep_mask])

  } else if (mode == "collapse_levels") {
    cmap <- parse_collapse_levels(collapse_map_str)
    if (length(cmap) == 0) { cat("SKIP:", analysis_id, "-> collapse_map vacío.\n\n"); return(NULL) }
    raw <- as.character(v); raw[trimws(raw)==""] <- NA
    mapped <- rep(NA_character_, length(raw))
    for (nm in names(cmap)) mapped[!is.na(raw) & raw==nm] <- cmap[[nm]]
    keep_mask <- !is.na(mapped) & (mapped %in% c(case_label, ref_label))
    group <- factor(mapped[keep_mask])

  } else if (mode == "range_vs") {
    if (!nzchar(collapse_map_str)) { cat("SKIP:", analysis_id, "-> collapse_map vacío en range_vs.\n\n"); return(NULL) }
    ranges <- parse_range_map(collapse_map_str)
    if (length(ranges) == 0) { cat("SKIP:", analysis_id, "-> collapse_map inválido en range_vs.\n\n"); return(NULL) }
    x <- as_num_safe(v)
    g <- rep(NA_character_, length(x))
    for (rr in ranges) {
      idx <- is.finite(x) & x >= rr$min & x <= rr$max
      g[idx] <- rr$label
    }
    keep_mask <- !is.na(g)
    if (nzchar(case_label) && nzchar(ref_label)) {
      keep_mask <- keep_mask & g %in% c(case_label, ref_label)
    }
    g <- g[keep_mask]
    if (length(g) == 0) { cat("SKIP:", analysis_id, "-> range_vs dejó 0 muestras.\n\n"); return(NULL) }
    group <- factor(g)

  } else if (mode == "survival_cut") {
    if (!nzchar(event_var) || !(event_var %in% colnames(meta_base))) {
      cat("SKIP:", analysis_id, "-> event_var inválido.\n\n"); return(NULL)
    }
    time <- as_num_safe(v)
    ev <- as_num_safe(meta_base[[event_var]])
    if (!is.finite(cut_val)) { cat("SKIP:", analysis_id, "-> cut inválido.\n\n"); return(NULL) }
    g <- rep(NA_character_, length(time))
    g[is.finite(time) & is.finite(ev) & ev==1 & time < cut_val] <- ref_label
    g[is.finite(time) & time >= cut_val] <- case_label
    if (exclude_rule == "exclude_if_censored_lt_cut") {
      cens_lt <- is.finite(time) & is.finite(ev) & ev==0 & time < cut_val
      g[cens_lt] <- NA
    }
    keep_mask <- !is.na(g)
    group <- factor(g[keep_mask])

  } else {
    cat("SKIP:", analysis_id, "-> mode no soportado:", mode, "\n\n"); return(NULL)
  }

  group <- droplevels(group)
  if (nlevels(group) < 2) { cat("SKIP:", analysis_id, "-> <2 niveles.\n\n"); return(NULL) }
  if (!(ref_label %in% levels(group)) || !(case_label %in% levels(group))) {
    cat("SKIP:", analysis_id, "-> no están case/ref tras filtrar.\n\n"); return(NULL)
  }
  group <- relevel(group, ref = ref_label)
  tab <- table(group)
  if (any(tab < 2)) { cat("SKIP:", analysis_id, "-> algún grupo <2.\n\n"); return(NULL) }

  counts_sub <- counts_mat[, keep_mask, drop = FALSE]

  # ---- edgeR QL ----
  y <- DGEList(counts_sub, group = group)
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~ group)

  keep_feat <- filterByExpr(y, design = design)
  if (!any(keep_feat)) { cat("SKIP:", analysis_id, "-> filterByExpr dejó 0.\n\n"); return(NULL) }
  y <- y[keep_feat, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)

  tt <- topTags(qlf, n = Inf)$table
  tt$feature_id <- rownames(tt); rownames(tt) <- NULL
  if (!is.null(anno_df)) tt <- merge(tt, anno_df, by.x="feature_id", by.y="id", all.x=TRUE)
  tt <- force_numeric_df(tt, c("logFC","logCPM","F","PValue","FDR"))

  tag <- paste0(case_label, "_vs_", ref_label)

  all_fp <- file.path(analysis_dir, paste0("DE_", analysis_fs, "_", tag, "_all_", ts, ".tsv"))
  sig05_fp <- file.path(analysis_dir, paste0("DE_", analysis_fs, "_", tag, "_sig_FDR0.05_", ts, ".tsv"))
  sig10_fp <- file.path(analysis_dir, paste0("DE_", analysis_fs, "_", tag, "_sig_FDR0.1_", ts, ".tsv"))
  write.table(tt, all_fp, sep="\t", quote=FALSE, row.names=FALSE)
  write.table(tt[!is.na(tt$FDR) & tt$FDR < 0.05, , drop=FALSE], sig05_fp, sep="\t", quote=FALSE, row.names=FALSE)
  write.table(tt[!is.na(tt$FDR) & tt$FDR < 0.10, , drop=FALSE], sig10_fp, sep="\t", quote=FALSE, row.names=FALSE)

  png(file.path(analysis_dir, paste0("MDplot_", analysis_fs, "_", tag, "_", ts, ".png")),
      width=1400, height=1000, res=150)
  plotMD(qlf, status = decide_status_safe(qlf), main = paste0("MD: ", analysis_id, " (", tag, ")"))
  abline(h=c(-1,1), lty=2); dev.off()

  plot_volcano(tt, analysis_fs, tag, analysis_dir, ts)

  qc <- data.frame(
    analysis_id = analysis_id,
    mode = mode,
    var = var_key,
    n_samples = ncol(y),
    n_features_tested = nrow(y),
    groups = paste(names(tab), as.integer(tab), sep=":", collapse=";"),
    comparison = tag,
    n_sig_FDR0.05 = sum(!is.na(tt$FDR) & tt$FDR < 0.05),
    n_sig_FDR0.1  = sum(!is.na(tt$FDR) & tt$FDR < 0.10),
    stringsAsFactors = FALSE
  )
  qc
}

# =========================
# Run all analyses
# =========================
qc_list <- list()
for (i in seq_len(nrow(spec))) {
  aid <- as.character(spec$analysis_id[i])
  cat("---- Running:", aid, "----\n")
  qc_list[[aid]] <- tryCatch(run_one(spec[i, , drop = FALSE]),
                             error = function(e) { cat("ERROR en", aid, ":", conditionMessage(e), "\n\n"); NULL })
}

qc_df <- do.call(rbind, Filter(Negate(is.null), qc_list))
if (!is.null(qc_df) && nrow(qc_df) > 0) {
  summary_fp <- file.path(outdir, paste0("DE_summary_all_miRNA_protCoding_", ts, ".tsv"))
  write.table(qc_df, summary_fp, sep="\t", quote=FALSE, row.names=FALSE)
  cat("\nResumen global:", summary_fp, "\n")
}

cat("\nFIN. Log:", log_fp, "\n")
sink(type="message"); sink(type="output"); close(zz)
