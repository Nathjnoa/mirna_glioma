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
as_num_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

# ---- inputs ----
spec_fp   <- get_arg("--spec", file.path("config","de_specs.csv"))
de_root   <- get_arg("--de_root", file.path("results","DE"))
out_root  <- get_arg("--out_root", file.path("results","tables","MiEAA_GSEA"))
species   <- get_arg("--species", "hsa")
mirna_type <- get_arg("--mirna_type", "mature")

min_hits_cli  <- as.numeric(get_arg("--min_hits", "5"))
sig_level <- 1  # importante: 1 para no filtrar (se fuerza)
p_adj_method <- get_arg("--p_adj_method", "fdr")
independent_p_adj_cli <- tolower(get_arg("--independent_p_adj", "true")) %in% c("true","t","1","yes","y")

# Si quieres cambiar categorías desde CLI: --categories "A,B,C"
cats_arg <- get_arg("--categories", "")
default_cats <- c(
  "miRPathDB_GO_Biological_process_mature",
  "miRPathDB_GO_Molecular_function_mature",
  "miRPathDB_KEGG_mature",
  "miRPathDB_Reactome_mature"
)
categories <- if (nzchar(cats_arg)) trimws(strsplit(cats_arg, ",")[[1]]) else default_cats
categories <- categories[nzchar(categories)]

# ranking score choice
# Nota: para glmQLFTest, recomendado usar signed_sqrtF; signed_logp puede sobreponderar p-valores extremos en GSEA ponderado.
rank_mode <- get_arg("--rank_mode", "signed_sqrtF") # signed_sqrtF | signed_logp | signed_F | signed_logFC

# corrida definitiva (A_conservative)
run_tag <- "A_conservative"
run_independent_p_adj <- FALSE
run_min_hits <- 5

# ---- path resolution (robusto a CWD) ----
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
de_root <- resolve_path(de_root)
out_root <- resolve_path(out_root)

is_abs_path <- function(p) {
  grepl("^(/|[A-Za-z]:)", p)
}
if (!is_abs_path(spec_fp)) spec_fp <- file.path(proj_root, spec_fp)
if (!is_abs_path(de_root)) de_root <- file.path(proj_root, de_root)
if (!is_abs_path(out_root)) out_root <- file.path(proj_root, out_root)

has_de_files <- function(root) {
  if (!dir.exists(root)) return(FALSE)
  length(list.files(root, pattern = "^DE_.*_all_.*\\.tsv$", recursive = TRUE)) > 0
}
if (!has_de_files(de_root)) {
  alt_root <- file.path(proj_root, de_root)
  if (alt_root != de_root && has_de_files(alt_root)) {
    de_root <- alt_root
  }
}

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

log_fp <- file.path("logs", paste0("mieaa_gsea_all_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type="output"); sink(zz, type="message")

cat("============================================\n")
cat("miEAA GSEA - todas las comparaciones (desde spec)\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("spec:", spec_fp, "\n")
cat("de_root:", de_root, "\n")
cat("out_root:", out_root, "\n")
cat("species:", species, " | mirna_type:", mirna_type, "\n")
cat("categories:", paste(categories, collapse=", "), "\n")
cat("min_hits CLI (ignorado):", min_hits_cli, "\n")
cat("sig_level:", sig_level, "\n")
cat("rank_mode:", rank_mode, "\n")
cat("independent_p_adj CLI (ignorado):", independent_p_adj_cli, "\n")
cat("run_tag:", run_tag, "\n")
cat("run_independent_p_adj:", run_independent_p_adj, "\n")
cat("run_min_hits:", run_min_hits, "\n")
cat("Log:", log_fp, "\n")
cat("============================================\n\n")

suppressPackageStartupMessages({
  library(rbioapi)
})

mieaa_cats <- rbioapi::rba_mieaa_cats(species = species, mirna_type = mirna_type)
mieaa_ids <- if (is.data.frame(mieaa_cats)) {
  id_col <- intersect(c("category_id", "category", "id", "name"), names(mieaa_cats))
  if (length(id_col) == 0) as.character(mieaa_cats[[1]]) else as.character(mieaa_cats[[id_col[1]]])
} else {
  as.character(mieaa_cats)
}
mieaa_ids <- mieaa_ids[!is.na(mieaa_ids) & nzchar(mieaa_ids)]
missing_cats <- setdiff(categories, mieaa_ids)
if (length(missing_cats) > 0) {
  stop(
    "Categorías miEAA no encontradas: ", paste(missing_cats, collapse = ", "),
    "\nSugerencia: imprime rba_mieaa_cats() y filtra por 'miRPathDB' y 'GO'."
  )
}

# --- helpers ---
latest_de_file <- function(dir_path) {
  ff <- list.files(dir_path, pattern = "^DE_.*_all_.*\\.tsv$", full.names = TRUE)
  if (length(ff) == 0) return(NA_character_)
  ff[order(file.info(ff)$mtime, decreasing = TRUE)][1]
}

find_de_file_fallback <- function(de_root, aid_fs) {
  ff <- list.files(de_root, pattern = "^DE_.*_all_.*\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(ff) == 0) return(NA_character_)
  sep <- .Platform$file.sep
  by_dir <- ff[grepl(paste0(sep, aid_fs, sep), ff, fixed = TRUE)]
  by_name <- ff[grepl(paste0("DE_", aid_fs), basename(ff), fixed = TRUE)]
  cand <- unique(c(by_dir, by_name))
  if (length(cand) == 0) return(NA_character_)
  cand[order(file.info(cand)$mtime, decreasing = TRUE)][1]
}

build_ranked_list <- function(de, rank_mode = "signed_sqrtF") {
  # require columns: feature_id, type, logFC, PValue (and optionally F for signed_F/signed_sqrtF)
  if (!("feature_id" %in% names(de))) stop("DE table sin 'feature_id'")
  if (!("type" %in% names(de))) stop("DE table sin 'type'")
  de_mi <- de[tolower(de$type) == "mirna", , drop = FALSE]

  # columnas necesarias
  if (!("logFC" %in% names(de_mi))) stop("DE miRNA sin 'logFC'")
  if (rank_mode == "signed_logp" && !("PValue" %in% names(de_mi))) stop("DE miRNA sin 'PValue'")
  if (rank_mode %in% c("signed_F", "signed_sqrtF") && !("F" %in% names(de_mi))) stop("DE miRNA sin 'F'")

  # limpiar
  de_mi$logFC <- as_num_safe(de_mi$logFC)
  if ("PValue" %in% names(de_mi)) de_mi$PValue <- as_num_safe(de_mi$PValue)
  if ("F" %in% names(de_mi)) de_mi$F <- as_num_safe(de_mi$F)

  ok <- is.finite(de_mi$logFC)
  if (rank_mode == "signed_logp") ok <- ok & is.finite(de_mi$PValue)
  if (rank_mode %in% c("signed_F", "signed_sqrtF")) ok <- ok & is.finite(de_mi$F)

  de_mi <- de_mi[ok, , drop = FALSE]
  if (nrow(de_mi) < 20) stop("Muy pocos miRNAs con stats válidas (n<20).")

  # score
  if (rank_mode == "signed_logp") {
    de_mi$score <- sign(de_mi$logFC) * (-log10(pmax(de_mi$PValue, .Machine$double.xmin)))
  } else if (rank_mode == "signed_sqrtF") {
    de_mi$score <- sign(de_mi$logFC) * sqrt(pmax(de_mi$F, 0))
  } else if (rank_mode == "signed_F") {
    de_mi$score <- sign(de_mi$logFC) * de_mi$F
  } else if (rank_mode == "signed_logFC") {
    de_mi$score <- de_mi$logFC
  } else {
    stop("rank_mode inválido: ", rank_mode, " (esperado: signed_sqrtF | signed_logp | signed_F | signed_logFC)")
  }

  # deduplicar por feature_id (si existiera)
  de_mi <- de_mi[order(abs(de_mi$score), decreasing = TRUE), ]
  de_mi <- de_mi[!duplicated(de_mi$feature_id), ]

  ranked <- de_mi$feature_id[order(de_mi$score, decreasing = TRUE)]
  ranked <- ranked[!is.na(ranked) & nzchar(ranked)]
  ranked
}

coerce_results <- function(res) {
  if (is.null(res) || !is.data.frame(res) || nrow(res) == 0) return(res)
  # normalizar nombres por si vienen distintos
  # En tu salida: P-value, P-adjusted, Q-value, Observed (strings)
  for (cc in c("P-value","P-adjusted","Q-value","Observed")) {
    if (cc %in% names(res)) res[[cc]] <- as_num_safe(res[[cc]])
  }
  if ("Observed" %in% names(res)) res[["Observed"]] <- as.integer(res[["Observed"]])
  res
}

# --- read spec ---
if (!file.exists(spec_fp)) stop("No existe spec: ", spec_fp)
spec <- read.csv(spec_fp, check.names = FALSE)
if (!("analysis_id" %in% names(spec))) stop("spec sin 'analysis_id'")
analysis_ids <- as.character(spec$analysis_id)

# outputs summary
summary_rows <- list()

for (aid in analysis_ids) {
  aid_fs <- safe_name(aid)
  de_dir <- file.path(de_root, aid_fs)
  if (!dir.exists(de_dir)) {
    cat("WARN", aid, "-> no existe dir:", de_dir, "\n")
  }
  de_fp <- if (dir.exists(de_dir)) latest_de_file(de_dir) else NA_character_
  if (is.na(de_fp) || !nzchar(de_fp) || !file.exists(de_fp)) {
    de_fp <- find_de_file_fallback(de_root, aid_fs)
  }
  if (is.na(de_fp) || !nzchar(de_fp) || !file.exists(de_fp)) {
    cat("SKIP", aid, "-> no hay DE_*_all_*.tsv\n")
    next
  }

  cat("\n--------------------------------------------\n")
  cat("Analysis:", aid, "\n")
  cat("DE file:", de_fp, "\n")

  de <- read.delim(de_fp, check.names = FALSE)

  # construir ranking
  ranked_mirnas <- tryCatch(build_ranked_list(de, rank_mode = rank_mode),
                            error = function(e) { cat("ERROR ranking:", conditionMessage(e), "\n"); NULL })

  if (is.null(ranked_mirnas)) next

  cat("Ranked miRNAs:", length(ranked_mirnas), "\n")

  run_dir <- file.path(out_root, aid_fs, run_tag)
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

  cat("Run:", run_tag,
      "| independent_p_adj:", run_independent_p_adj,
      "| min_hits:", run_min_hits,
      "| categories:", paste(categories, collapse=", "), "\n")

  # Export ranking (mature)
  rank_mature_fp <- file.path(run_dir, paste0("MiEAA_ranked_mature_", aid_fs, "_", ts, ".txt"))
  writeLines(ranked_mirnas, rank_mature_fp)
  cat("Export ranking mature:", rank_mature_fp, " n=", length(ranked_mirnas), "\n")

  # run GSEA (no filtrar)
  res_all <- tryCatch({
    rbioapi::rba_mieaa_enrich(
      test_set = ranked_mirnas,
      mirna_type = mirna_type,
      test_type = "GSEA",
      species = species,
      categories = categories,
      p_adj_method = p_adj_method,
      independent_p_adj = run_independent_p_adj,
      sig_level = sig_level,
      min_hits = run_min_hits
    )
  }, error = function(e) {
    cat("ERROR miEAA:", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(res_all) || !is.data.frame(res_all) || nrow(res_all) == 0) {
    cat("WARN: resultados vacíos (res_all)\n")
    summary_rows[[aid]] <- data.frame(
      analysis_id = aid,
      run_tag = run_tag,
      de_file = basename(de_fp),
      n_ranked_miRNAs = length(ranked_mirnas),
      n_terms = 0,
      n_q_lt_0.05 = 0,
      min_q = NA_real_,
      min_p_adjusted = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }

  # export full unfiltered results
  all_fp <- file.path(run_dir, paste0("MiEAA_GSEA_all_", aid_fs, "_", run_tag, "_", ts, ".tsv"))
  write.table(res_all, all_fp, sep = "\t", quote = FALSE, row.names = FALSE)

  res_all <- coerce_results(res_all)
  n_terms <- nrow(res_all)

  if ("Q-value" %in% names(res_all)) {
    qvals <- res_all$`Q-value`
    q_source <- "Q-value"
  } else if ("P-adjusted" %in% names(res_all)) {
    qvals <- res_all$`P-adjusted`
    q_source <- "P-adjusted"
    cat("WARN: Q-value no disponible; usando P-adjusted para filtros.\n")
  } else {
    qvals <- rep(NA_real_, n_terms)
    q_source <- NA_character_
    cat("WARN: Q-value/P-adjusted no disponibles; filtros pueden quedar vacíos.\n")
  }
  padj <- if ("P-adjusted" %in% names(res_all)) res_all$`P-adjusted` else rep(NA_real_, n_terms)

  res_q05 <- res_all[is.finite(qvals) & qvals < 0.05, , drop = FALSE]

  # top50 por Q-value (fallback a P-adjusted si Q no existe/finito)
  if (any(is.finite(qvals))) {
    ord <- order(qvals, padj, na.last = TRUE)
  } else if (any(is.finite(padj))) {
    ord <- order(padj, na.last = TRUE)
    cat("WARN: Q-value no disponible; top50 ordenado por P-adjusted.\n")
  } else {
    ord <- seq_len(n_terms)
    cat("WARN: Q-value/P-adjusted no disponibles; top50 mantiene orden original.\n")
  }
  res_top50 <- res_all[ord, , drop = FALSE]
  if (nrow(res_top50) > 50) res_top50 <- res_top50[seq_len(50), , drop = FALSE]

  q05_fp <- file.path(run_dir, paste0("MiEAA_GSEA_Qlt0.05_", aid_fs, "_", run_tag, "_", ts, ".tsv"))
  top50_fp <- file.path(run_dir, paste0("MiEAA_GSEA_top50_", aid_fs, "_", run_tag, "_", ts, ".tsv"))
  write.table(res_q05, q05_fp, sep="\t", quote=FALSE, row.names=FALSE)
  write.table(res_top50, top50_fp, sep="\t", quote=FALSE, row.names=FALSE)

  n_q05 <- nrow(res_q05)
  min_q <- if (any(is.finite(qvals))) min(qvals[is.finite(qvals)]) else NA_real_
  min_padj <- if (any(is.finite(padj))) min(padj[is.finite(padj)]) else NA_real_

  cat("Terms:", n_terms, " | Q<0.05 (", q_source, "):", n_q05,
      " | min_q:", min_q, " | min_padj:", min_padj, "\n", sep="")

  summary_rows[[aid]] <- data.frame(
    analysis_id = aid,
    run_tag = run_tag,
    de_file = basename(de_fp),
    n_ranked_miRNAs = length(ranked_mirnas),
    n_terms = n_terms,
    n_q_lt_0.05 = n_q05,
    min_q = min_q,
    min_p_adjusted = min_padj,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
summary_fp <- file.path(out_root, paste0("MiEAA_GSEA_summary_", ts, ".tsv"))
write.table(summary_df, summary_fp, sep="\t", quote=FALSE, row.names=FALSE)

cat("\n============================================\n")
cat("FIN\n")
cat("Summary:", summary_fp, "\n")
cat("Log:", log_fp, "\n")
cat("============================================\n")

sink(type="message"); sink(type="output"); close(zz)
