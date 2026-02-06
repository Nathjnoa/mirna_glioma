#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# 17_generate_pathway_mirna_reports.R
# Generate reports of top pathways with miRNAs annotated by rank and statistic
#
# Output format for miRNAs column:
#   hsa-miR-451a (#3, logFC=-0.85); hsa-miR-148a-3p (#6, logFC=-0.72); ...
#
# Two report types:
#   1. DE-based GSEA: ranks miRNAs as used for miEAA GSEA (from script 09),
#      and reports logFC for interpretability
#   2. Survival-based GSEA: ranks miRNAs as used for miEAA GSEA (from script 13),
#      using Cox z-score
# -----------------------------------------------------------------------------

options(stringsAsFactors = FALSE)

# Mirror safe_name() used in other pipeline scripts (e.g., 09_miEAA_GSEA_all_comparisons.R)
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

# --- CLI parsing ---
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[[hit + 1]]
}

# Parameters
spec_fp      <- get_arg("--spec", "config/de_specs.csv")
de_root      <- get_arg("--de_root", "results/DE")
gsea_root    <- get_arg("--gsea_root", "results/tables/MiEAA_GSEA")
surv_root    <- get_arg("--surv_root", "results/tables/SurvivalRank_GSEA")
out_dir      <- get_arg("--out_dir", "results/tables/reports")
run_tag      <- get_arg("--run_tag", "A_conservative")
top_n        <- as.integer(get_arg("--top_n", "12"))
surv_run_tag <- get_arg("--surv_run_tag", "SurvivalRank_CoxZ_miRPathDB")

# Create output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Timestamp
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Log setup
log_dir <- "logs"
dir.create(log_dir, showWarnings = FALSE)
log_fp <- file.path(log_dir, paste0("pathway_mirna_reports_", ts, ".txt"))
log_con <- file(log_fp, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("=== Pathway-miRNA Reports Generator ===\n")
cat("Timestamp:", ts, "\n")
cat("Parameters:\n")
cat("  spec_fp:", spec_fp, "\n")
cat("  de_root:", de_root, "\n")
cat("  gsea_root:", gsea_root, "\n")
cat("  surv_root:", surv_root, "\n")
cat("  out_dir:", out_dir, "\n")
cat("  run_tag:", run_tag, "\n")
cat("  top_n:", top_n, "\n")
cat("  surv_run_tag:", surv_run_tag, "\n\n")

# --- Helper functions ---

read_table <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(path, data.table = FALSE, check.names = FALSE)
  } else {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

# Find most recent file matching pattern
find_latest_file <- function(dir_path, pattern) {
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  # Sort by modification time, get most recent
  info <- file.info(files)
  files[which.max(info$mtime)]
}

# Extract timestamp (YYYYmmdd_HHMMSS) from filenames
extract_ts <- function(filename) {
  m <- regmatches(filename, regexpr("[0-9]{8}_[0-9]{6}", filename))
  if (length(m) == 0 || is.na(m) || !nzchar(m)) return(NA_character_)
  m
}

# Format miRNA with rank and statistic
format_mirna_with_stats <- function(mirnas_str, rank_df, stat_col, stat_name) {
  if (is.na(mirnas_str) || mirnas_str == "") return("")

  mirnas <- trimws(unlist(strsplit(mirnas_str, ";", fixed = TRUE)))
  mirnas <- mirnas[nzchar(mirnas)]

  formatted <- sapply(mirnas, function(m) {
    idx <- which(rank_df$feature_id == m)
    if (length(idx) == 0) {
      # Try matching without hsa- prefix or with variations
      idx <- which(grepl(m, rank_df$feature_id, fixed = TRUE))
    }
    if (length(idx) == 0) {
      return(m)  # Return as-is if not found
    }
    idx <- idx[1]  # Take first match
    rank_val <- rank_df$rank[idx]
    stat_val <- rank_df[[stat_col]][idx]
    sprintf("%s (#%d, %s=%.3f)", m, rank_val, stat_name, stat_val)
  })

  paste(formatted, collapse = "; ")
}

# -----------------------------------------------------------------------------
# PART 1: DE-based GSEA Reports
# -----------------------------------------------------------------------------

cat("\n========== DE-based GSEA Reports ==========\n\n")

# Read spec file
if (!file.exists(spec_fp)) {
  cat("WARNING: Spec file not found:", spec_fp, "\n")
  spec_rows <- data.frame()
} else {
  spec_rows <- read_table(spec_fp)
  cat("Loaded", nrow(spec_rows), "comparisons from spec\n")
}

# Process each comparison
de_results_list <- list()

for (i in seq_len(nrow(spec_rows))) {
  analysis_id <- spec_rows$analysis_id[i]
  analysis_fs <- safe_name(analysis_id)
  cat("\n--- Processing:", analysis_id, "---\n")

  # Find DE results directory
  de_dir <- file.path(de_root, analysis_id)
  if (!dir.exists(de_dir)) de_dir <- file.path(de_root, analysis_fs)
  if (!dir.exists(de_dir)) {
    cat("  DE directory not found (tried:", analysis_id, "and", analysis_fs, "), skipping\n")
    next
  }

  # Find GSEA results
  gsea_dir <- file.path(gsea_root, analysis_fs, run_tag)
  if (!dir.exists(gsea_dir)) gsea_dir <- file.path(gsea_root, analysis_id, run_tag)

  gsea_file <- find_latest_file(gsea_dir, "^MiEAA_GSEA_Qlt.*\\.tsv$")
  if (is.null(gsea_file)) {
    # Try all results if filtered not available
    gsea_file <- find_latest_file(gsea_dir, "^MiEAA_GSEA_all.*\\.tsv$")
  }

  if (is.null(gsea_file)) {
    cat("  GSEA file not found in", gsea_dir, ", skipping\n")
    next
  }
  cat("  GSEA file:", basename(gsea_file), "\n")

  # Pick DE file: most recent at/just before this GSEA run (fallback: latest overall)
  de_files <- list.files(de_dir, pattern = "^DE_.*_all_.*\\.tsv$", full.names = TRUE)
  if (length(de_files) == 0) {
    cat("  DE file not found, skipping\n")
    next
  }
  de_info <- file.info(de_files)
  gsea_mtime <- file.info(gsea_file)$mtime
  idx_before <- which(de_info$mtime <= gsea_mtime)
  if (length(idx_before) > 0) {
    de_file <- de_files[idx_before][which.max(de_info$mtime[idx_before])]
  } else {
    de_file <- de_files[which.max(de_info$mtime)]
  }
  cat("  DE file:", basename(de_file), "\n")

  # Read DE results (for logFC lookup)
  de_df <- read_table(de_file)
  if (is.null(de_df) || nrow(de_df) == 0) {
    cat("  Empty DE results, skipping\n")
    next
  }
  if (!("feature_id" %in% names(de_df))) {
    cat("  WARNING: DE table missing 'feature_id'; cannot annotate miRNAs, skipping\n")
    next
  }
  if (!("logFC" %in% names(de_df))) {
    cat("  WARNING: DE table missing 'logFC'; cannot annotate miRNAs, skipping\n")
    next
  }

  de_mi <- if ("type" %in% names(de_df)) {
    de_df[tolower(as.character(de_df$type)) == "mirna", , drop = FALSE]
  } else {
    cat("  WARNING: DE table missing 'type' column; using full table for logFC lookup\n")
    de_df
  }
  de_mi$logFC <- as_num_safe(de_mi$logFC)
  de_mi <- de_mi[is.finite(de_mi$logFC) & nzchar(de_mi$feature_id), , drop = FALSE]
  if (nrow(de_mi) == 0) {
    cat("  No valid miRNA logFC values found, skipping\n")
    next
  }
  logfc_map <- setNames(de_mi$logFC, de_mi$feature_id)

  # Determine miRNA rank list used for this GSEA run (prefer the exported ranked list)
  gsea_ts <- extract_ts(basename(gsea_file))
  rank_fp <- if (!is.na(gsea_ts)) {
    file.path(gsea_dir, paste0("MiEAA_ranked_mature_", analysis_fs, "_", gsea_ts, ".txt"))
  } else {
    NA_character_
  }

  ranked_ids <- NULL
  if (!is.na(rank_fp) && file.exists(rank_fp)) {
    ranked_ids <- trimws(readLines(rank_fp, warn = FALSE))
    ranked_ids <- ranked_ids[nzchar(ranked_ids)]
  } else {
    # Fallback: compute rank list like script 09 (signed_sqrtF)
    if (!("F" %in% names(de_mi))) {
      cat("  WARNING: Missing 'F' column; cannot reconstruct signed_sqrtF ranking. Skipping\n")
      next
    }
    de_mi$F <- as_num_safe(de_mi$F)
    de_mi <- de_mi[is.finite(de_mi$F), , drop = FALSE]
    if (nrow(de_mi) < 20) {
      cat("  Too few miRNAs with valid F/logFC to rank (n<20), skipping\n")
      next
    }
    de_mi$score <- sign(de_mi$logFC) * sqrt(pmax(de_mi$F, 0))
    de_mi <- de_mi[order(abs(de_mi$score), decreasing = TRUE), , drop = FALSE]
    de_mi <- de_mi[!duplicated(de_mi$feature_id), , drop = FALSE]
    ranked_ids <- de_mi$feature_id[order(de_mi$score, decreasing = TRUE)]
  }

  rank_df <- data.frame(
    feature_id = ranked_ids,
    rank = seq_along(ranked_ids),
    stringsAsFactors = FALSE
  )
  rank_df$logFC <- logfc_map[rank_df$feature_id]
  cat("  miRNAs ranked for GSEA:", nrow(rank_df), "\n")

  # Read GSEA results
  gsea_df <- read_table(gsea_file)
  if (is.null(gsea_df) || nrow(gsea_df) == 0) {
    cat("  Empty GSEA results, skipping\n")
    next
  }

  # Exclude GO Molecular Function (not biologically relevant for miRNA analysis)
  gsea_df <- gsea_df[!grepl("Molecular function", gsea_df$Category, ignore.case = TRUE), ]

  # Get top N per category and enrichment direction
  categories <- unique(gsea_df$Category)

  for (cat_name in categories) {
    for (enrich_dir in c("enriched", "depleted")) {
      subset_df <- gsea_df[gsea_df$Category == cat_name &
                           gsea_df$Enrichment == enrich_dir, ]

      if (nrow(subset_df) == 0) next

      # Sort by Q-value and take top N
      subset_df <- subset_df[order(subset_df$`Q-value`), ]
      subset_df <- head(subset_df, top_n)

      # Format miRNAs with rank and logFC
      mirna_col <- grep("miRNA|precursor", names(subset_df),
                        ignore.case = TRUE, value = TRUE)[1]
      if (is.na(mirna_col)) mirna_col <- "miRNAs/precursors"

      subset_df$miRNAs_annotated <- sapply(subset_df[[mirna_col]], function(m) {
        format_mirna_with_stats(m, rank_df, "logFC", "logFC")
      })

      # Build output row
      for (j in seq_len(nrow(subset_df))) {
        de_results_list[[length(de_results_list) + 1]] <- data.frame(
          comparison = analysis_id,
          Category = cat_name,
          Subcategory = subset_df$Subcategory[j],
          Enrichment = enrich_dir,
          Q_value = subset_df$`Q-value`[j],
          Observed = subset_df$Observed[j],
          miRNAs_with_rank = subset_df$miRNAs_annotated[j],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  cat("  Added pathways to report\n")
}

# Combine and save DE-based report
if (length(de_results_list) > 0) {
  de_report <- do.call(rbind, de_results_list)

  out_file <- file.path(out_dir,
                        paste0("GSEA_DE_top", top_n, "_annotated_", ts, ".csv"))
  write.csv(de_report, out_file, row.names = FALSE)
  cat("\n\nDE-based report saved:", out_file, "\n")
  cat("Total rows:", nrow(de_report), "\n")
} else {
  cat("\nNo DE-based results to report\n")
}

# -----------------------------------------------------------------------------
# PART 2: Survival-based GSEA Reports
# -----------------------------------------------------------------------------

cat("\n\n========== Survival-based GSEA Reports ==========\n\n")

# Find Cox results
surv_dir <- file.path(surv_root, surv_run_tag)
if (!dir.exists(surv_dir)) {
  cat("Survival directory not found:", surv_dir, "\n")
} else {
  cox_file <- find_latest_file(surv_dir, "^Cox_univariate_miRNA_all.*\\.tsv$")

  if (is.null(cox_file)) {
    cat("Cox results file not found\n")
  } else {
    cat("Cox file:", basename(cox_file), "\n")

    # Read Cox results
    cox_df <- read_table(cox_file)

    if (!is.null(cox_df) && nrow(cox_df) > 0) {
      # Rank by |z| descending
      cox_df$z <- as_num_safe(cox_df$z)
      cox_df <- cox_df[is.finite(cox_df$z) & nzchar(cox_df$feature_id), , drop = FALSE]
      if (nrow(cox_df) == 0) {
        cat("Empty Cox results after filtering non-finite z, skipping\n")
      } else {
        # Rank list used for GSEA: signed z, descending (as in script 13)
        cox_ts <- extract_ts(basename(cox_file))
        rank_fp <- if (!is.na(cox_ts)) {
          file.path(surv_dir, paste0("Ranked_miRNAs_by_CoxZ_", cox_ts, ".txt"))
        } else {
          NA_character_
        }

        ranked_ids <- NULL
        if (!is.na(rank_fp) && file.exists(rank_fp)) {
          ranked_ids <- trimws(readLines(rank_fp, warn = FALSE))
          ranked_ids <- ranked_ids[nzchar(ranked_ids)]
        } else {
          ranked_tmp <- cox_df[order(cox_df$z, decreasing = TRUE), , drop = FALSE]
          ranked_tmp <- ranked_tmp[!duplicated(ranked_tmp$feature_id), , drop = FALSE]
          ranked_ids <- ranked_tmp$feature_id
        }

        z_map <- setNames(cox_df$z, cox_df$feature_id)
        cox_rank <- data.frame(
          feature_id = ranked_ids,
          rank = seq_along(ranked_ids),
          z = z_map[ranked_ids],
          stringsAsFactors = FALSE
        )

        cat("Loaded", nrow(cox_df), "features from Cox analysis\n")
        cat("Top z-score:", round(max(cox_rank$z, na.rm = TRUE), 3), "\n")

        # Find reduced GSEA results (after rrvgo)
        reduced_dir <- file.path(surv_dir, "reduced_rrvgo")

        # Try to find GSEA results
        gsea_file <- NULL
        if (dir.exists(reduced_dir)) {
          gsea_file <- find_latest_file(reduced_dir, "^miEAA_GSEA.*reduced.*\\.tsv$")
        }
        if (is.null(gsea_file)) {
          gsea_file <- find_latest_file(surv_dir, "^miEAA_GSEA.*\\.tsv$")
        }

	        if (!is.null(gsea_file)) {
	          cat("GSEA file:", basename(gsea_file), "\n")

          gsea_df <- read_table(gsea_file)

	          if (!is.null(gsea_df) && nrow(gsea_df) > 0) {
	            surv_results_list <- list()

            # Exclude GO Molecular Function
            gsea_df <- gsea_df[!grepl("Molecular function", gsea_df$Category, ignore.case = TRUE), ]

            categories <- unique(gsea_df$Category)

            for (cat_name in categories) {
              for (enrich_dir in c("enriched", "depleted")) {
                subset_df <- gsea_df[gsea_df$Category == cat_name &
                                     gsea_df$Enrichment == enrich_dir, ]

                if (nrow(subset_df) == 0) next

                # Sort by Q-value and take top N
                q_col <- grep("Q-value|q-value|Q_value", names(subset_df),
                             ignore.case = TRUE, value = TRUE)[1]
                if (is.na(q_col)) q_col <- "Q-value"

                subset_df <- subset_df[order(subset_df[[q_col]]), ]
                subset_df <- head(subset_df, top_n)

                # Format miRNAs with rank and z-score (rank from GSEA list)
                mirna_col <- grep("miRNA|precursor", names(subset_df),
                                 ignore.case = TRUE, value = TRUE)[1]
                if (is.na(mirna_col)) mirna_col <- "miRNAs/precursors"

                subset_df$miRNAs_annotated <- sapply(subset_df[[mirna_col]], function(m) {
                  format_mirna_with_stats(m, cox_rank, "z", "z")
                })

                # Build output row
                for (j in seq_len(nrow(subset_df))) {
                  surv_results_list[[length(surv_results_list) + 1]] <- data.frame(
                    analysis = "SurvivalRank_CoxZ",
                    Category = cat_name,
                    Subcategory = subset_df$Subcategory[j],
                    Enrichment = enrich_dir,
                    Q_value = subset_df[[q_col]][j],
                    Observed = subset_df$Observed[j],
                    miRNAs_with_rank = subset_df$miRNAs_annotated[j],
                    stringsAsFactors = FALSE
                  )
                }
              }
            }

            # Combine and save survival-based report
	            if (length(surv_results_list) > 0) {
	              surv_report <- do.call(rbind, surv_results_list)

              out_file <- file.path(out_dir,
                                    paste0("GSEA_survGSEA_top", top_n, "_annotated_", ts, ".csv"))
              write.csv(surv_report, out_file, row.names = FALSE)
              cat("\nSurvival-based report saved:", out_file, "\n")
              cat("Total rows:", nrow(surv_report), "\n")
	            }
	          }
	        } else {
	          cat("GSEA results file not found for survival analysis\n")
	        }
    }
  }
}
}

# Close log
cat("\n\nDone. Log saved:", log_fp, "\n")
sink(type = "message")
sink()
close(log_con)
